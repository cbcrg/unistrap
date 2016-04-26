Shootstrap-NF
===================

A workflow for simultaneous MSA and tree generation, and tree reliability estimation. 
Shootstrap-NF generates alignment models, that allow it to estimate the phylogenetic tree having the most supported topology. 
In addition it "dresses" that tree with shootstrap support values, which takes into account boostrapping and alignmnet uncertainty effects.

How to use
-----------
    
You can run Shootstrap-NF using the following command: 

    nextflow run cbcrg/shootstrap [shootstrap-nf command line options]

For example: 

     nextflow run cbcrg/shootstrap -with-docker

Command line options:
---------------------

	--in_dir	Directory containing one or more datasets in fasta format to be analysed

	--out_dir	Directory to output the results

	--in_tree	The tree for which the shootstrap support values will be estimated 
				(default: Shootstrap-NF estimates the tree having the most supported tree topology, across the 
				different trees relsulting from the different MSA replicates)

	--rep_num	The number of MSA replicates to be produced (default: 100)

	--seed		Random number to be used to reproduce the shuffling procedure in the exact way (default: 10)


    
    
Bonus
------

You can access the container by running this command: 

	docker run -ti --entrypoint /bin/bash cbcrg/shootstrap-nf

