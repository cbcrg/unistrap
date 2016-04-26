
FROM pditommaso/dkrbase
MAINTAINER Maria Chatzou <mxatzou@gmail.com>


RUN apt-get update -y --fix-missing && apt-get install -y \
    git \
    cmake \
    libargtable2-dev


#RUN git clone https://github.com/cbcrg/shootstrap.git
#RUN cd shootstrap \
#    cp bin/* ../ \
#    cd ..


#
# Get and compile T-Coffee
#
RUN wget -q http://www.tcoffee.org/Packages/Stable/Version_11.00.8cbe486/linux/T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz && \
  tar xf T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz && \
  mv T-COFFEE_installer_Version_11.00.8cbe486_linux_x64 /opt/tcoffee && \
  rm -rf T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz
    
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/blast/bin:/opt/sate/sate-core/bin:/opt/tcoffee/bin:/opt/tcoffee/plugins/linux/
ENV TEMP /tmp
ENV DIR_4_TCOFFEE /opt/tcoffee
ENV EMAIL_4_TCOFFEE tcoffee.msa@gmail.com
ENV CACHE_4_TCOFFEE /tmp/cache/
ENV LOCKDIR_4_TCOFFEE /tmp/lck/
ENV TMP_4_TCOFFEE /tmp/tmp/  


#
# Get and compile ClustalO
#
RUN wget http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz
RUN tar -xvzf clustal-omega-1.2.1.tar.gz
RUN cd clustal-omega-1.2.1 &&\
    ./configure &&\
    make &&\
    cp src/clustalo /usr/local/bin/    
    #mv ../exe/seqboot /usr/local/bin/



#
# Get and compile FastTree
#
RUN mkdir FastTree_dir
RUN cd FastTree_dir
RUN wget http://meta.microbesonline.org/fasttree/FastTree.c
RUN gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
RUN mv FastTree /usr/local/bin/
RUN cd ..

#ENTRYPOINT ["/usr/local/bin/FastTree"]

#
# Get and compile SeqBoot
#
RUN wget http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz
RUN tar -xvzf phylip-3.696.tar.gz
RUN cd phylip-3.696/src \
    Makefile.unx Makefile \
    make install \
    mv ../exe/seqboot /usr/local/bin/

#ENTRYPOINT ["/usr/local/bin/seqboot"]

