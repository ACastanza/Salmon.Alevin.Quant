#! /bin/bash
# Using getopt
set -e

abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

while getopts ":b:q:i:c:t:l:w:m:r:x:j:" opt; do
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
            method=="$OPTARG"
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
        l)
            lib="$OPTARG"
            ;;
        w)
          if [[ $OPTARG =~ ^[^-]+$ ]];then
            whitelist=`realpath $OPTARG`
            echo "-w <whitelist> = $whitelist"
          elif [[ $OPTARG =~ ^-. ]];then
            whitelist=""
            let OPTIND=$OPTIND-1
          else
            whitelist=`realpath $OPTARG`
            echo "-w <whitelist> = $whitelist"
          fi          
            ;;
        m)
          if [[ $OPTARG =~ ^[^-]+$ ]];then
            mrna=`realpath $OPTARG`
            echo "-m <mrna> = $mrna"
          elif [[ $OPTARG =~ ^-. ]];then
            mrna=""
            let OPTIND=$OPTIND-1
          else
            mrna=`realpath $OPTARG`
            echo "-m <mrna> = $mrna"
          fi          
            ;;
        r)
          if [[ $OPTARG =~ ^[^-]+$ ]];then
            rrna=`realpath $OPTARG`
            echo "-r <rrna> = $rrna"
          elif [[ $OPTARG =~ ^-. ]];then
            rrna=""
            let OPTIND=$OPTIND-1
          else
            rrna=`realpath $OPTARG`
            echo "-r <rrna> = $rrna"
          fi          
            ;;
        x)
            mtx="$OPTARG"
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
	apt install less -y#Less is needed to parse the GTF into the correct format
	apt install file -y
	if (file $tgMap | grep -q compressed); then
	gunzip -ck $tgMap > transcriptome
	tgMap=transcriptome
	fi
	test=$(zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | head -n 1| cut -f9 | tr -s ";" " " | awk '{print$3}' | sort | uniq |  sed 's/\"//g')
	if [[ $test == "transcript_id" ]]; then
		zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$2}' | sort | uniq | sed 's/\"//g' > "txp2gene.tsv"
		txp2gene="txp2gene.tsv"
		elif [[ $test == "gene_version" ]]; then
		echo "Separate version field (ensembl, non-gencode transcriptome, eg. rat, etc)"
		zless -S $tgMap | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$6 "."  $8"\t"$2 "." $4}' | sort | uniq | sed 's/\"//g' > "txp2gene.tsv"
		txp2gene="txp2gene.tsv"
		else
		echo "Error Parsing GTF"
	fi

elif [[ tgcol -eq 2 ]]; then
	txp2gene=$tgMap
else
	echo "Error Parsing tgMap"
fi

outdir="alevin"
mkdir -p "$outdir" ;


params=()
[[ $method == "dropseq" ]] && params+=(--dropseq)
[[ $method == "chromium" ]] && params+=(--chromium)
[[ $method == "chromiumV3" ]] && params+=(--chromiumV3)

[[ -e "$whitelist" ]] && params+=(--whitelist $whitelist)
[[ -e "$mrna" ]] && params+=(--mrna $mrna)
[[ -e "$rrna" ]] && params+=(--mrna $rrna)
[[ $mtx == "TRUE" ]] && params+=(--dumpMtx)


# [[ $CONDITION == true ]] && params+=(--param)

salmon alevin \
      --no-version-check \
      -i "salmon_index" \
      -p $threads \
      -l $lib \
      -1 $barcodes \
      -2 $reads \
      "${params[@]}" \
      -o $outdir ;

tar -czvf $outdir.alevin.tar.gz -C $outdir .
rm -rf $outdir

    echo "--Done." ;

rm -rf salmon_index txp2gene.tsv transcriptome
