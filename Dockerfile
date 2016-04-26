
FROM pditommaso/dkrbase
MAINTAINER Maria Chatzou <mxatzou@gmail.com>


RUN apt-get update -y --fix-missing && apt-get install -y \
    git \
    cmake \
    libargtable2-dev


#
# Get perl scripts needed 
#
RUN mkdir Shootstrap
RUN cd Shootstrap
RUN git clone git@github.com:cbcrg/shootstrap.git .
RUN cp bin/* /usr/local/bin/
RUN cd ..

#
# Get and compile FastTree
#
RUN mkdir FastTree
RUN mv FastTree
RUN wget http://meta.microbesonline.org/fasttree/FastTree.c
RUN gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
RUN mv FastTree /usr/local/bin/FastTree
RUN cd ..

ENTRYPOINT ["/usr/local/bin/FastTree"]

#
# Get and compile SeqBoot
#
RUN mkdir SeqBoot
RUN mv SeqBoot
RUN wget http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz
RUN tar -xvzf phylip-3.696.tar.gz
RUN cd phylip-3.696/src
RUN cp Makefile.unx Makefile
RUN make install
RUN mv ../exe/seqboot /usr/local/bin/

ENTRYPOINT ["/usr/local/bin/seqboot"]

