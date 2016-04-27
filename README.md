Shootstrap-NF
===================

Given a protein sequence dataset, Shootstrap-NF estimates a reference phylogenetic tree and the shootstrap support of every branch in that tree. The shootsrap support is a generalized version of Felsenstein branch bootstrap support measure that simultaneously takes into account the sequence input-order effect (shuffling) and the column sampling effect (Felsenstein regular Bootstrap) to provide a support value for every branch in a reference tree. 

In practice, the procedure involves generating N shuffled input sequence datasets, from which as many replicate MSAs are generated and used to generate as many trees. These trees are compared using RF and the one with the highest average similarity is selected as a refence. The replicate shuffled MSAs are then used to draw bootstrap replicates. These replicates are used to estimate shootstrap  support values for every branch of the reference tree. 

Installation
-----------



Usage
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

