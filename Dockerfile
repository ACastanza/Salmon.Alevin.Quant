### copyright 2017-2021 Regents of the University of California and the Broad Institute. All rights reserved.
FROM combinelab/salmon\:1.5.2

MAINTAINER Barbara Hill <bhill@broadinstitute.org>


# install conda
RUN mkdir /conda && \
    cd /conda && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh && \
    bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -b -p /opt/conda
ENV PATH="/opt/conda/bin:${PATH}"

RUN conda config --add channels defaults && \
    conda config --add channels bioconda  && \ 
    conda config --add channels conda-forge

# install alevin-fry from conda
RUN conda install --yes -c bioconda alevin-fry

# clean up conda caches
RUN conda clean -afy

# docker build --rm https://github.com/genepattern/Salmon.Alevin.Quant.git#develop -f Dockerfile -t genepattern/salmon-alevin-quant:<tag>
# make sure this repo and tag match the manifest & don't forget to docker push!
# docker push genepattern/salmon-alevin-quant:<tag>

# you can use this command to run Docker and iterate locally (update for your paths and module name, of course)
# docker run --rm -it -v /c/Users/MyUSER/PathTo/Salmon.Alevin.Quant:/mnt/mydata:rw genepattern/Salmon.Alevin.Quant:<tag> bash