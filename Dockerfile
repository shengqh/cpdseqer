FROM shengqh/bioinfo:base

RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir --upgrade \
    pysam \
    biopython \
    pytabix

RUN apt-get install -y tabix

ENV BOWTIE2_VERSION="2.3.5.1"
ENV PATH=/opt/bowtie2-${BOWTIE2_VERSION}-linux-x86_64:$PATH
RUN cd /opt; \
    wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip; \
    unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip ; \
    rm -rf bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip; \
    bowtie2 --version
    
ENV CPDSEQER_VERSION="0.0.8"
RUN pip3 install git+git://github.com/shengqh/cpdseqer.git
RUN cpdseqer -h
