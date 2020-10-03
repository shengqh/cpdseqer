FROM shengqh/bioinfo:base

ENV BOWTIE2_VERSION="2.5.1"
ENV PATH=/opt/bowtie2-${BOWTIE2_VERSION}-linux-x86_64:$PATH
RUN cd /opt; \
    wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip; \
    unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip ; \
    rm -rf bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip; \
    bowtie2 --version


RUN apt-get update && apt-get install -y tabix

RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir --upgrade \
    pysam \
    biopython \
    pytabix

RUN R -e "BiocManager::install('knitr');library('knitr')"
RUN R -e "BiocManager::install('reshape2');library('reshape2')"
RUN R -e "BiocManager::install('gplots');library('gplots')"
RUN R -e "BiocManager::install('psych');library('psych')"
RUN R -e "BiocManager::install('data.table');library('data.table')"
RUN R -e "BiocManager::install('edgeR');library('edgeR')"

ENV CPDSEQER_VERSION="0.3.4"
RUN pip3 install git+git://github.com/shengqh/cpdseqer.git
RUN cpdseqer -h
