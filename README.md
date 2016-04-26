Shootstrap-NF
===================

A workflow for simultenious MSA and tree generation. 
Shootstrap-NF generates alignment models, that allow it to estimate the phylogenetic tree having the most supported topology. 
In addition it "dresses" that tree with shootstrap support values, which takes into account boostrapping and alignmnet uncertainty effects.

How to use
-----------
    
You can run RAxML using the following command: 

    nextflow run cbcrg/shootstrap [raxml command line options]

For example: 

     nextflow run cbcrg/shootstrap -with-docker
    
    
Bonus
------

You can access the container by running this command: 

	docker run -ti --entrypoint /bin/bash cbcrg/shootstrap-nf

