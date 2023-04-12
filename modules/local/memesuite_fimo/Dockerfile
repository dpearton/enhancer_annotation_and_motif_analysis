FROM  nfcore/base:2.1

LABEL authors="alexthiery1@gmail.com" \
      description="Docker image containing memesuite"

RUN apt-get update && \
    apt-get install -y build-essential && \
    apt-get install -y make

# Install cellranger
RUN cd /opt/ && \
    wget -O meme-5.5.2.tar.gz "https://meme-suite.org/meme/meme-software/5.5.2/meme-5.5.2.tar.gz" && \	
    tar zxf meme-5.5.2.tar.gz && \	
    cd meme-5.5.2 && \	
    ./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt && \	
    make && \	
    make test && \	
    make install && \
    export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.2:$PATH

