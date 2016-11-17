#!/usr/bin/env nextflow
/*
 * Copyright (c) 2016, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'Shootstrap-NF'.
 *
 *   Shootstrap-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Shootstrap-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Shootstrap-NF.  If not, see <http://www.gnu.org/licenses/>.
 */


/* 
 * Main Shootstrap-NF script
 *
 * @authors
 * Maria Chatzou <mxatzou@gmail.com>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

params.in_dir="$baseDir/data/dataset/*"
params.rep_num=100
params.seed=10
params.aligner="clustalo"
params.in_tree="" 
params.stats=false 


Channel
	.fromPath(params.in_dir)
	.ifEmpty { error "Cannot find any data -- Check the path specified: `${params.in_dir}`" }
    .set { file_names }
    
    
in_tree_file = params.in_tree ? file(params.in_tree) : null
if( in_tree_file ) assert in_tree_file.exists(), "The tree file does not exist: $in_tree_file !!!" 



def get_prefix(name) { 
  def m = name =~ /(.*)_\d+\.fa\.aln.*/
  m.matches() ? m.group(1) : 'unknown_tree_name' 
}

def get_tree_prefix(name) { 
  def m = name =~ /(.*)_\d+\..*\.tree/
  m.matches() ? m.group(1) : 'unknown_tree_name' 
}


def combine_trees( allFiles ) {

  def prefix = get_prefix(allFiles[0].name)
  def bigTree = java.nio.file.Files.createTempFile(prefix,".tree")
 
  def result = []
  allFiles.each {  
    bigTree << it.text  
    result << [ it, bigTree ]
  }
  result
}



/*
 * Step 1. Shuffle replicates generation. Creates X shuffle replicates for each input-dataset.
 */
process get_shuffle_replicates{

  input:
      file(seq_file) from file_names
  output:
      file "*.fa" into shuffle_replicates mode flatten
  
  shell:
  '''
      tmp_name=`basename !{seq_file} | awk -F. '{print $1}'`
      seq_shuffle.pl !{seq_file} ${tmp_name} !{params.seed} !{params.rep_num}
  '''

}



/*
 * Step 2. Replicate MSAs generation. Aligns the shuffle replicates, using one the template aligners (default : Clustalo).
 */
process get_msa_replicates{

  input:
      file(seq_file) from shuffle_replicates
  output:
      file "${seq_file}.aln" into msa_replicates, msa_replicates2, msa_replicates3, msa_replicates4
  
  script:
      template "${params.aligner}_msa_command"
}



/*
 * Step 3. Replicate trees generation. Estimates a phylogenetic tree from each one of the replicate MSAs, using one the template tree estimator programs (default : FastTree).
 */
process get_msa_trees{

  input:
      file(seq_file) from msa_replicates
  output:
      file "${seq_file}.msa_tree" into msa_trees, msa_trees2
  
  script:
      output_tree="${seq_file}.msa_tree"
      template "tree_commands"
}



/*
 * Step 4. Shuffle trees generation. For each one of the replicate trees the shuffle tree is estimated, (i.e. the replicate tree dressed with the shuffle support values 
 *         - that correspond to the percentage of how many times each branch in the tree was seen in the collection of replicate trees).
 *         This is done using the "Fast Tree-Comparison Tools".
 */
msa_trees.map { file -> tuple(get_prefix(file.name), file) }
         .groupTuple()
         .flatMap { prefix, files -> combine_trees(files) }
         .set{ msa_all_trees }

process get_shuffle_trees{

  input:
      set file(in_tree_file), file(all_tree_file) from msa_all_trees
      
  output:
      set file("*.stable.tree") into stable_trees
      set stdout, file(in_tree_file) into most_stable_tree
  
  shell:
  '''
      tmp_name=`basename !{in_tree_file} | awk -F. '{print $1}'`
      perl -I !{baseDir}/bin/ !{baseDir}/bin/CompareToBootstrap.pl -tree !{in_tree_file} -boot !{all_tree_file} > $tmp_name.stable.tree
      cat  $tmp_name.stable.tree | sed 's/)/\\n/g' | sed 's/;//g' | awk -F: '{print $1}' | grep -v "(" | grep -v "^$" | awk '{ sum=sum+$1 ; sumX1+=(($2)^2)} END { avg=sum/NR; printf "%f\\n", avg}' 
  '''
}



/*
 * Step 5. Bootstrap replicates generation. Create 1 bootstrap replicate from each MSA replicate, using "seqboot" from PHYLIP package.
 */
process get_seqboot_replicates{

  input:
      file(seq_file) from msa_replicates2
  output:
      file "${seq_file}.replicate" into boot_replicates
  
  """
      echo -e "$seq_file\nR\n1\nY\n31\n"|seqboot
      mv outfile ${seq_file}.replicate
  """

}



/*
 * Step 6. Bootstrap trees generation. Create 1 bootstrap tree from each bootstrap replicate, using one the template tree estimator programs (default : FastTree).
 */
process get_bootstrap_trees{

  input:
      file(seq_file) from boot_replicates
  output:
      file "${seq_file}.tree" into trees
  
  script:
      output_tree="${seq_file}.tree"
      template "tree_commands"
} 



/*
 * Step 7. Shootstrap tree generation. Create 1 bootstrap tree from each bootstrap replicate, using "seqboot" from PHYLIP package.
 */
trees.map { file -> tuple(get_prefix(file.name), file) }
         .groupTuple()
         .flatMap { prefix, files -> combine_trees(files) }
         .set{ all_trees }

process get_shootstrap_tree{

  input:
      set file(in_tree_file), file(all_tree_file) from all_trees
      
  output:
      file "*.shootstrap.tree" into shootstrap_tree
  
  shell:
  '''
      tmp_name=`basename !{in_tree_file} | awk -F. '{print $1}'`
      perl -I !{baseDir}/bin/ !{baseDir}/bin/CompareToBootstrap.pl -tree !{in_tree_file} -boot !{all_tree_file} > $tmp_name.shootstrap.tree
  '''
} 



/*
 * Step 8. MSA stability estimation. Estimates the MSA stability of each dataset. 
 */
process get_msa_stability_stats{

  when:
  params.stats 

  input:
      set prefix,file(msa) from msa_replicates3.map{ file -> tuple(get_prefix(file.name), file) }.groupTuple()

  output:
      file "${prefix}_insta.txt" into insta

  shell:
  '''
      echo -n !{prefix}"\t" >> !{prefix}_insta.txt 
      for j in {1..!{params.rep_num}};{ for r in {1..!{params.rep_num}};{ [ $j != $r ] && msaa -a !{prefix}_$j.fa.aln -r !{prefix}_$r.fa.aln; }; } |  awk '{ sum += $2 } END { if(NR>0) print sum/NR }' >> !{prefix}_insta.txt  
  '''
} 

insta.collectFile(name:'MSA_INSTABILITY.stats', seed:"Name\tmsaInstability\n", storeDir:params.out_dir) 



/*
 * Step 9. Tree stability estimation. Estimates the Tree stability of each the replicate trees. 
 */
process get_tree_stability_stats{

  when:
  params.stats 

  input:
      set prefix,file(msa) from stable_trees.map{ file -> tuple(get_tree_prefix(file.name), file) }.groupTuple()

  output:
      file "${prefix}_insta.txt" into tree_insta

  shell:
  '''
      echo -n !{prefix}"\t" >> !{prefix}_insta.txt 
      for r in {1..!{params.rep_num}}; { cat !{prefix}_$r.stable.tree | sed 's/);//g' |sed 's/):0.0//g' | sed 's/)/\\n/g' | awk -F: '{print "|"$1}' | grep -v "(" | awk -F"|" '{ sum+=$2} END { print sum/NR }';  } | awk '{ sum+=$1} END { printf sum/NR"\\n"}' >> !{prefix}_insta.txt  
  '''
} 

tree_insta.collectFile(name:'TREE_INSTABILITY.stats', seed:"Name\ttreeInstability\n", storeDir:params.out_dir) 



/*
 * Step 10. MSA identity and number of sequences estimation. Estimates the MSA stability of each dataset and reports how many sequences each dataset contains.
 */
process get_general_stats{

  when:
  params.stats 

  input:
      set prefix,file(msa) from msa_replicates4.map{ file -> tuple(get_prefix(file.name), file) }.groupTuple().map{ prefix,files -> tuple(prefix,files[0]) }

  output:
      file "${prefix}_insta.txt" into genStats

  shell:
  '''
      echo -n !{prefix}"\t" >> !{prefix}_insta.txt  
      grep -P '^\\s+\\d+\\s+\\d+' !{prefix}_*.fa.aln | awk '{printf $2"\t"}'  >> !{prefix}_insta.txt 
      msaa -a !{prefix}_*.fa.aln -i | grep aln | awk '{print $2}' >> !{prefix}_insta.txt 
  '''
} 

genStats.collectFile(name:'GENERAL.stats', seed:"Name\tseqNum\tmsaIdentity\n", storeDir:params.out_dir) 

