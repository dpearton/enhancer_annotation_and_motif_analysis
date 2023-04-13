# Motif Enhancer Screening
Mini nextflow pipeline for filtering enhancers associated to differentially expressed genes, and then screening for motifs of interest




docker run -it -d --name memesuite -v ${PWD}:/home/data memesuite/memesuite /bin/bash
docker exec memesuite jaspar2meme -bundle -pfm test_data/test.pfm



docker run -it -v $PWD:/temp quay.io/biocontainers/gtfparse:1.2.1--pyh864c0ab_0 /bin/bash


