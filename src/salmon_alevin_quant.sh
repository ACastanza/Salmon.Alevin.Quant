#! /bin/bash
# Using getopt
set -e

trap abort ERR PROF
abort()
{

[[ -f "transcriptome" ]] && rm -rf transcriptome
[[ -f "GTF_rrnaGenes.txt" ]] && rm GTF_rrnaGenes.txt
[[ -f "GTF_mtGenes.txt" ]] && rm GTF_mtGenes.txt
[[ -f "GTF_txp2gene.tsv" ]] && rm GTF_txp2gene.tsv

[[ -d "salmon_index" ]] && rm -rf salmon_index

    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

while getopts ":b:q:i:c:l:t:z:u:w:m:n:r:s:x:e:f:d:a:j:" opt; do
    case $opt in
        b)
            inbarcode=`realpath $OPTARG`
            ;;
        q)
            infile=`realpath $OPTARG`
            ;;
        i)
            index=`realpath $OPTARG`
            ;;
        c)
            method="$OPTARG"
            ;;
        l)
            lib="$OPTARG"
            ;;
        t)
          if [[ $OPTARG =~ ^[^-]+$ ]];then
            tgMap=`realpath $OPTARG`
            echo "--tgMap = $tgMap"
          elif [[ $OPTARG =~ ^-. ]];then
            tgMap=""
            let OPTIND=$OPTIND-1
          else
            tgMap=`realpath $OPTARG`
            echo "--tgMap = $tgMap"
          fi
            ;;
        z)
            basename="$OPTARG"
            ;;
        u)
            idtype="$OPTARG"
            ;;
        w)
          if [[ $OPTARG =~ ^[^-]+$ ]];then
            whitelist=`realpath $OPTARG`
            echo "--whitelist = $whitelist"
          elif [[ $OPTARG =~ ^-. ]];then
            whitelist=""
            let OPTIND=$OPTIND-1
          else
            whitelist=`realpath $OPTARG`
            echo "--whitelist = $whitelist"
          fi
            ;;
        m)
          if [[ $OPTARG =~ ^[^-]+$ ]];then
            mtgenes=`realpath $OPTARG`
            echo "--mtgenes = $mtgenes"
          elif [[ $OPTARG =~ ^-. ]];then
            mtgenes=""
            let OPTIND=$OPTIND-1
          else
            mtgenes=`realpath $OPTARG`
            echo "--mrna = $mtgenes"
          fi
            ;;
        n)
            mtbuild="$OPTARG"
            ;;
        r)
          if [[ $OPTARG =~ ^[^-]+$ ]];then
            rrna=`realpath $OPTARG`
            echo "--rrna = $rrna"
          elif [[ $OPTARG =~ ^-. ]];then
            rrna=""
            let OPTIND=$OPTIND-1
          else
            rrna=`realpath $OPTARG`
            echo "--rrna = $rrna"
          fi
            ;;
        s)
            rbuild="$OPTARG"
            ;;
        x)
            mtx="$OPTARG"
            ;;
        e)
            expect="$OPTARG"
            ;;
        f)
            force="$OPTARG"
            ;;
        d)
            features="$OPTARG"
            ;;
        a)
            numboot="$OPTARG"
            ;;
        j)
            threads="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            abort
            ;;
    esac
done

 mkdir -p salmon_index
 tar -zxvf $index -C salmon_index

 barcodes=$(cat $inbarcode | sort)
 reads=$(cat $infile | sort)

tgcol=$(awk -F'\t' '{print NF; exit}' "$tgMap")
if [[ tgcol -eq 2 ]]; then
	txp2gene=$tgMap
elif [[ tgcol -ne 2 ]]; then
	apt-get install less file -qq > /dev/null
 	## Less and file are needed to parse the GTF into the correct format this will install them (silently?)
	if (file $tgMap | grep -q compressed); then
	gunzip -ck $tgMap > transcriptome
	tgMap=transcriptome
	fi

	test=$(zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | head -n 1| cut -f9 | tr -s ";" " ")
	txid=$(($(echo $test | awk -v b="transcript_id" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}') + 1))

	if [[ $txid == "4" && $idtype == "gene_id" ]]; then
		zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$2}' | sort | uniq | sed 's/\"//g' > "GTF_txp2gene.tsv"
		txp2gene="GTF_txp2gene.tsv"

			if [[ $mtbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript" && ($1=="M" || $1=="chrM" || $1=="MT")' | cut -f9 | tr -s ";" " " | awk '{print$2}' | sort | uniq | sed 's/\"//g' > "GTF_mtGenes.txt"
			mtgenes=GTF_mtGenes.txt
			fi

			if [[ $rbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '$6=="\"rRNA\""' | awk '{print$2}' | sort | uniq | sed 's/\"//g' > "GTF_rrnaGenes.txt"
			rrna=GTF_rrnaGenes.txt
			fi

		elif [[ $txid == "6"  && $idtype == "gene_id" ]]; then
		echo "Separate version field (ensembl, non-gencode transcriptome, eg. rat, etc)"
		zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$6 "."  $8"\t"$2 "." $4}' | sort | uniq | sed 's/\"//g' > "GTF_txp2gene.tsv"
		txp2gene="GTF_txp2gene.tsv"

			if [[ $mtbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript" && ($1=="M" || $1=="chrM" || $1=="MT")' | cut -f9 | tr -s ";" " " | awk '{print$2 "."  $4}' | sort | uniq | sed 's/\"//g' > "GTF_mtGenes.txt"
			mtgenes=GTF_mtGenes.txt
			fi

			if [[ $rbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '$14=="\"rRNA\""' | awk '{print$2 "."  $4}' | sort | uniq | sed 's/\"//g' > "GTF_rrnaGenes.txt"
			rrna=GTF_rrnaGenes.txt
			fi

		elif [[ $txid == "4" && $idtype == "gene_id_noversion" ]]; then
		echo "Omitting gene id versions..."
		zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$2}' | awk 'BEGIN { OFS=FS="\t" } { sub("\\..*", "", $2); print }' | sort | uniq | sed 's/\"//g' > "GTF_txp2gene.tsv"
		txp2gene="GTF_txp2gene.tsv"

			if [[ $mtbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript" && ($1=="M" || $1=="chrM" || $1=="MT")' | cut -f9 | tr -s ";" " " | awk '{print$2}' | awk 'BEGIN { OFS=FS="\t" } { sub("\\..*", "", $1); print }' | sort | uniq | sed 's/\"//g' > "GTF_mtGenes.txt"
			mtgenes=GTF_mtGenes.txt
			fi

			if [[ $rbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '$6=="\"rRNA\""' | awk '{print$2}' | awk 'BEGIN { OFS=FS="\t" } { sub("\\..*", "", $1); print }' | sort | uniq | sed 's/\"//g' > "GTF_rrnaGenes.txt"
			rrna=GTF_rrnaGenes.txt
			fi

		elif [[ $txid == "6"  && $idtype == "gene_id_noversion" ]]; then
		echo "Separate version field detected (ensembl, non-gencode transcriptome, eg. rat, etc). Omitting it..."
		zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$6 "."  $8"\t"$2}' | sort | uniq | sed 's/\"//g' > "GTF_txp2gene.tsv"
		txp2gene="GTF_txp2gene.tsv"

			if [[ $mtbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript" && ($1=="M" || $1=="chrM" || $1=="MT")' | cut -f9 | tr -s ";" " " | awk '{print$2}' | sort | uniq | sed 's/\"//g' > "GTF_mtGenes.txt"
			mtgenes=GTF_mtGenes.txt
			fi

			if [[ $rbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '$14=="\"rRNA\""' | awk '{print$2}' | sort | uniq | sed 's/\"//g' > "GTF_rrnaGenes.txt"
			rrna=GTF_rrnaGenes.txt
			fi

		elif [[ $txid == "4" && $idtype == "gene_symbol" ]]; then
		zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$8}' | sort | uniq | sed 's/\"//g' > "GTF_txp2gene.tsv"
		txp2gene="GTF_txp2gene.tsv"

			if [[ $mtbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript" && ($1=="M" || $1=="chrM" || $1=="MT")' | cut -f9 | tr -s ";" " " | awk '{print$8}' | sort | uniq | sed 's/\"//g' > "GTF_mtGenes.txt"
			mtgenes=GTF_mtGenes.txt
			fi

			if [[ $rbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '$6=="\"rRNA\""' | awk '{print$8}' | sort | uniq | sed 's/\"//g' > "GTF_rrnaGenes.txt"
			rrna=GTF_rrnaGenes.txt
			fi

		elif [[ $txid == "6"  && $idtype == "gene_symbol" ]]; then
		echo "Separate version field (ensembl, non-gencode transcriptome, eg. rat, etc)"
		zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$6 "."  $8"\t"$10}' | sort | uniq | sed 's/\"//g' > "GTF_txp2gene.tsv"
		txp2gene="GTF_txp2gene.tsv"

			if [[ $mtbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript" && ($1=="M" || $1=="chrM" || $1=="MT")' | cut -f9 | tr -s ";" " " | awk '{print$10}' | sort | uniq | sed 's/\"//g' > "GTF_mtGenes.txt"
			mtgenes=GTF_mtGenes.txt
			fi

			if [[ $rbuild == "TRUE" ]]; then
			zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '$14=="\"rRNA\""' | awk '{print$10}' | sort | uniq | sed 's/\"//g' > "GTF_rrnaGenes.txt"
			rrna=GTF_rrnaGenes.txt
			fi

		else
		echo "Error Parsing GTF"
	fi

elif [[ tgcol -eq 2 ]]; then
	txp2gene=$tgMap
else
	echo "Error Parsing tgMap"
fi

outdir=$basename
mkdir -p "$outdir" ;


params=()
[[ $method == "dropseq" ]] && params+=(--dropseq)
[[ $method == "chromium" ]] && params+=(--chromium)
[[ $method == "chromiumV3" ]] && params+=(--chromiumV3)
[[ $method == "citeseq" ]] && params+=(--citeseq)
[[ $method == "celseq" ]] && params+=(--celseq)
[[ $method == "celseq2" ]] && params+=(--celseq2)
[[ $method == "quartzseq2" ]] && params+=(--quartzseq2)
[[ $method == "sciseq3" ]] && params+=(--sciseq3)
[[ $method == "indropV2" ]] && params+=(--indropV2)
[[ $method == "splitSeqV1" ]] && params+=(--splitSeqV1)
[[ $method == "splitSeqV2" ]] && params+=(--splitSeqV2)

[[ $method == "rhapsody" ]] && params+=(--umi-geometry \'1[53-60]\')
[[ $method == "rhapsody" ]] && params+=(--bc-geometry \'1[1-9,22-30,44-52]\')
[[ $method == "rhapsody" ]] && params+=(--read-geometry \'2[1-end]\')


[[ -e "$whitelist" ]] && params+=(--whitelist $whitelist)
[[ -e "$mtgenes" ]] && params+=(--mrna $mtgenes)
[[ -e "$rrna" ]] && params+=(--rrna $rrna)
[[ $mtx == "TRUE" ]] && params+=(--dumpMtx)
[[ $expect -gt "0" ]] && params+=(--expectCells $expect)
[[ $force -gt "0" ]] && params+=(--forceCells $force)
[[ $features == "TRUE" ]] && params+=(--dumpFeatures)
[[ $numboot -gt "0" ]] && params+=(--numCellBootstraps $numboot)

# [[ $CONDITION == true ]] && params+=(--param)

salmon alevin \
      --no-version-check \
      -i "salmon_index" \
      -p $threads \
      -l $lib \
      -1 $barcodes \
      -2 $reads \
      "${params[@]}" \
      --tgMap $txp2gene \
      -o $outdir ;

[[ -e GTF_txp2gene.tsv ]] && cp GTF_txp2gene.tsv $outdir
[[ -e GTF_mtGenes.txt ]] && cp GTF_mtGenes.txt $outdir
[[ -e GTF_rrnaGenes.txt ]] && cp GTF_rrnaGenes.txt $outdir

echo "Compressing Alevin output"
tar -czvf $outdir.tar.gz -C $outdir .

echo "Cleaning up"

[[ -f "transcriptome" ]] && rm -rf transcriptome
[[ -f "GTF_rrnaGenes.txt" ]] && rm GTF_rrnaGenes.txt
[[ -f "GTF_mtGenes.txt" ]] && rm GTF_mtGenes.txt
[[ -f "GTF_txp2gene.tsv" ]] && rm GTF_txp2gene.tsv

[[ -d "salmon_index" ]] && rm -rf salmon_index
[[ -d "$outdir" ]] && rm -rf $outdir

trap : 0
echo >&2 '
***********************
*** DONE Processing ***
***********************
'
