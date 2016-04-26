FROM pditommaso/dkrbase:1.2
MAINTAINER Maria Chatzou <mxatzou@gmail.com>


RUN apt-get update -y --fix-missing && apt-get install -y \
    git \
    cmake \
    libargtable2-dev


#
# Download T-Coffee pre-built package
#
RUN wget -q http://www.tcoffee.org/Packages/Stable/Version_11.00.8cbe486/linux/T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz && \
  tar xf T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz && \
  mv T-COFFEE_installer_Version_11.00.8cbe486_linux_x64 /opt/tcoffee && \
  rm -rf T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz
    
ENV PATH=$PATH:/opt/tcoffee/bin:/opt/tcoffee/plugins/linux/ TEMP=/tmp DIR_4_TCOFFEE=/opt/tcoffee EMAIL_4_TCOFFEE=tcoffee.msa@gmail.com CACHE_4_TCOFFEE=/tmp/cache/ LOCKDIR_4_TCOFFEE=/tmp/lck/ TMP_4_TCOFFEE=/tmp/tmp/  


#
# Get and compile ClustalO
#
RUN wget http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz -q -O- | tar -xz &&\
  cd clustal-omega-1.2.1 &&\
  ./configure &&\
  make &&\
  cp src/clustalo /usr/local/bin/ &&\
  rm -rf /clustal-omega-1.2.1   


#
# Get and compile FastTree
#
RUN mkdir FastTree_dir &&\
 cd FastTree_dir &&\
 wget http://meta.microbesonline.org/fasttree/FastTree.c -q &&\
 gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm &&\
 mv FastTree /usr/local/bin/ &&\
 rm -rf /FastTree_dir

#
# Get and compile SeqBoot
#
RUN wget http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz -q -O- | tar -xz &&\
 cd phylip-3.696/src &&\
 make -f Makefile.unx install &&\
 mv ../exe/seqboot /usr/local/bin/ &&\
 rm -rf /phylip-3.696

