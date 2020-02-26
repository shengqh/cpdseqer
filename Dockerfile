FROM shengqh/bioinfo:base

ENV BOWTIE2_VERSION="2.4.0"
ENV PATH=/opt/bowtie2-${BOWTIE2_VERSION}-linux-x86_64:$PATH
RUN cd /opt; \
    wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}_new/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip; \
    unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip ; \
    rm -rf bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip; \
    bowtie2 --version

RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir --upgrade \
    pysam \
    biopython \
    pytabix

RUN apt-get install -y tabix
    
ENV CPDSEQER_VERSION="0.1.2"
RUN pip3 install git+git://github.com/shengqh/cpdseqer.git
RUN cpdseqer -h
