FROM shengqh/bioinfo:java1.8.0_perl5.22.1_r3.6.1_python3.7.3

RUN pip install git+git://github.com/shengqh/annotateGenome.git

RUN annotateGenome -h
