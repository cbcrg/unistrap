Unistrap-NF
===================

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.18.0-brightgreen.svg)](http://nextflow.io)
![CircleCI status](https://circleci.com/gh/cbcrg/unistrap.png?style=shield)


Given a protein sequence dataset, Unistrap-NF estimates a reference phylogenetic tree and the unistrap support of every branch in that tree. The unistrap support is a generalized version of Felsenstein branch bootstrap support measure that simultaneously takes into account the alignment uncerainty effect and the column sampling effect (Felsenstein regular Bootstrap) to provide a support value for every branch in a reference tree. 

In practice, the procedure involves generating N shuffled input sequence datasets, from which as many replicate MSAs are generated and used to generate as many trees. These trees are compared using Robinson and Foulds (RF) and the one with the highest average similarity (most suported tree topology) is selected as a refence. The replicate shuffled MSAs are then used to draw bootstrap replicates. These replicates are used to estimate unistrap support values for every branch of the reference tree. 

Installation
-----------

To run unistrap-NF you need to install nextflow, by simply checking if you have Java 7+ and if yes, by then typing the following command:

	curl -fsSL get.nextflow.io | bash

If you want to replicate exactly the pipeline and/or not install all the dependencies unistrap-NF has, then you also need to install Docker and run unistrap-NF with the '-with-docker' flag, as demonstrated below. Otherwise all the dependencies of unistrap-NF have to be installed and put in the PATH.


Usage
-----------
    
You can run Unistrap-NF using the following command: 

    nextflow run cbcrg/unistrap [unistrap-nf command line options]

For example: 

     nextflow run cbcrg/unistrap -with-docker

This command will run (by automatically cloning the respository in your workstation, due to the [Nextflow integration with github](http://www.nextflow.io/docs/latest/sharing.html)) the analysis for the "hip.fa" example dataset contained in the data/dataset directory, and will generate a reference tree, with estimated unistrap support values for every branch in it.

     nextflow run unistrap.nf -with-docker --in_dir data/dataset/hip.fa 
     nextflow run unistrap.nf -with-docker --in_dir data/dataset/hip.fa --in_tree data/tree/hip.tree 

The 1st command will run the analysis for the "hip.fa" example dataset contained in the data/dataset directory, and will generate a reference tree, with estimated unistrap support values for every branch in it.
The 2nd command will estimate the unistrap support values for every branch in the given "hip.tree" example tree, contained in the data/tree directory. Note that in order for these commands to run the repository has to first be cloned and then to launch them from inside it. 

Command line options:
---------------------

	--in_dir	Directory containing one or more protein sequence datasets in fasta format to be analysed

	--out_dir	Output directory

	--in_tree	The tree in newick format for which the unistrap support values will be estimated 
				(default: Unistrap-NF estimates the tree having the most supported tree topology, across the 
				different trees relsulting from the different MSA replicates)

	--rep_num	The number of MSA replicates to be produced (default: 100)

	--seed		Seed for random number generation, (default: time )

Pipeline flow
--------------

![alt unistrap](unistrap_schematic_workflow.png)

Dependencies:
-------------

Unistrap-NF uses the following external programs (which have to exist in the PATH if the pipeline is executed without docker):

	Nextflow	http://www.nextflow.io/
	clustalo	http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz
	t-coffee	http://www.tcoffee.org/Packages/Stable/Version_11.00.8cbe486/linux/T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz
	seqboot		http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz
	FastTree	http://meta.microbesonline.org/fasttree/FastTree.c
        
Bonus
------

You can access the container by running this command: 

	docker run -ti cbcrg/unistrap-nf bash

