#!/usr/bin/env nextflow
/*
 * Copyright (c) 2016, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'Unistrap-NF'.
 *
 *   Unistrap-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Unistrap-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Unistrap-NF.  If not, see <http://www.gnu.org/licenses/>.
 */


/* 
 * Main Unistrap-NF script
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
params.boot=false 
params.boot_datasets="all"
params.boot_num=100

Channel
	.fromPath(params.in_dir)
	.ifEmpty { error "Cannot find any data -- Check the path specified: `${params.in_dir}`" }
    .set { file_names }
    
    
in_tree_file = params.in_tree ? file(params.in_tree) : null
assert !in_tree_file || in_tree_file.exists(), "The tree file does not exist: $in_tree_file !!!" 


def get_prefix(name) { 
  def m = name =~ /(.*)_\d+\.fa\.aln.*/
  m.matches() ? m.group(1) : 'unknown_tree_name' 
}

def get_tree_prefix(name) { 
  def m = name =~ /(.*)_\d+\..*\.tree/
  m.matches() ? m.group(1) : 'unknown_tree_name' 
}

def get_bootTree_prefix(name) { 
  def m = name =~ /(.*_\d+)\..*\.tree/
  m.matches() ? m.group(1) : 'unknown_tree_name' 
}


def combine_msa_trees( allFiles ) {

  final prefix = get_prefix(allFiles[0].name)
  final bigTree = cacheableFile(allFiles, "${prefix}.all_msa_trees")
  final isCached = bigTree.exists()
  final result = []
  allFiles.each {  
    if(!isCached) bigTree << it.text  
    result << [ it, bigTree ]
  }
  result
}

def combine_boot_trees( allFiles ) {

  final prefix = get_prefix(allFiles[0].name)
  final bigTree = cacheableFile(allFiles, "${prefix}.all_boot_trees") 
  final isCached = bigTree.exists()
  final result = []
  allFiles.each {  
    if(!isCached) bigTree << it.text  
    result << [ it, bigTree ]
  }
  result
}

def combine_rep_trees_by_prefix( allFiles ) {

  final prefix = get_prefix(allFiles[0].name)
  final bigTree = cacheableFile(allFiles, "${prefix}.all_replicate_trees")
  final isCached = bigTree.exists()
  final result = []
  allFiles.each {  
    if(!isCached) bigTree << it.text  
    result << [ prefix, bigTree ]
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
      file "${seq_file}.aln" into msa_replicates, msa_replicates2, msa_replicates3, msa_replicates4, msa_replicates5
  
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
      file "${seq_file}.msa_tree" into msa_trees, msa_trees3
  
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
         .flatMap { prefix, files -> combine_msa_trees(files) }
         .set{ msa_all_trees }

process get_shuffle_trees{

  input:
      set file(in_tree_file), file(all_tree_file) from msa_all_trees
      
  output:
      file("*.stable.tree") into stable_trees
      set stdout, file(in_tree_file) into most_stable_tree
      file "${all_tree_file}"
  
  shell:
  '''
      tmp_name=`basename !{in_tree_file} | awk -F. '{print $1}'`
      perl -I !{baseDir}/bin/ !{baseDir}/bin/CompareToBootstrap.pl -tree !{in_tree_file} -boot !{all_tree_file} > $tmp_name.stable.tree
      cat  $tmp_name.stable.tree | sed 's/)/\\n/g' | sed 's/;//g' | awk -F: '{print $1}' | grep -v "(" | grep -v "^$" | awk '{ sum=sum+$1 ; } END { avg=-1; if( NR>0 ) avg=sum/NR; printf "%f\\n", avg; }' 
  '''
}


/*
 * Step 5. Bootstrap replicates MSAs generation. Create 1 bootstrap replicate from each MSA replicate, using "seqboot" from PHYLIP package.
 */
process get_1_seqboot_replicates{

  when:         
  params.boot

  input:
      file(seq_file) from msa_replicates2
  output:
      file "${seq_file}.replicate" into boot_replicates   
  

  script:
  """
      shuf_num=`ls $seq_file| awk -F"." '{print \$1}' | awk -F"_" '{print \$NF}'`
      [ \$((shuf_num%2)) -eq 0 ] && shuf_num=\$((\$shuf_num + 10001))
      echo -e "$seq_file\nR\n1\nY\n\$shuf_num\n"|seqboot
      mv outfile ${seq_file}.replicate
  """
}


/*
 * Step 6. Bootstrap replicate trees generation. Create 1 bootstrap tree from each bootstrap replicate, using one the template tree estimator programs (default : FastTree).
 */
process get_1_bootstrap_replicate_trees{

  input:
      file(seq_file) from boot_replicates
  output:
      file "${seq_file}.tree" into boot_trees
  
  script:
      output_tree="${seq_file}.tree"
      template "tree_commands"
} 


/*
 * Step 7. Unistrap values estimation. Calculates the unistrap values for each replicate tree, using "Fast Tree-Comparison Tools" from FastTree package.
 */
boot_trees.map { file -> tuple(get_prefix(file.name), file) }
          .groupTuple()
          .flatMap { prefix, files -> combine_rep_trees_by_prefix(files) }
          .set{ all_boot_trees }


msa_trees3.map { file -> tuple(get_prefix(file.name), file) }
          .set{ all_msa_trees }


all_boot_trees
      .phase(all_msa_trees)
      .map { foo, bar -> tuple( foo[0], foo[1], bar[1] ) }
      .set{ all_msa_boot_trees }


process get_unistrap_trees{

  input:
      set prefix, file(all_tree_file), file(in_tree_file) from all_msa_boot_trees
      
  output:
      file "*.unistrap.tree" into unistrap_trees
      file "${all_tree_file}"
  
  shell:
  '''
      tmp_name=`basename !{in_tree_file} | awk -F. '{print $1}'`
      perl -I !{baseDir}/bin/ !{baseDir}/bin/CompareToBootstrap.pl -tree !{in_tree_file} -boot !{all_tree_file} > $tmp_name.unistrap.tree
  '''
} 




/**************************************************************************
 *                                                                        *
 *        Processes 8, 9 & 10 estimate the normal bootstrap values.       *
 *           They only run when the option "boot" is specified.           *
 *                                                                        *
 **************************************************************************/


/*
 * Step 8. Bootstrap replicates MSAs generation. Create X (default: 100) bootstrap replicates from each MSA replicate, using "seqboot" from PHYLIP package. 
 */

process get_100_seqboot_replicates{

  when:         
  params.boot

  input:
      file(seq_file) from msa_replicates5
  output:
      file "${seq_file}.boot_replicate" into norm_boot_replicates
  
  script:
  """
      if [ ${params.boot_datasets} == "all" ] 
      then
          echo -e "$seq_file\nR\n$params.boot_num\nY\n31\n"|seqboot
          mv outfile ${seq_file}.boot_replicate
      else
          for seq_set in $params.boot_datasets
          {
                echo -e "$seq_file\nR\n$params.boot_num\nY\n31\n"|seqboot
                mv outfile ${seq_file}.boot_replicate
          }    
      fi    
  """
}


/*
 * Step 9. Normal bootstrap trees generation. Create 100 bootstrap trees from the 100 bootstrap replicates, using one the template tree estimator programs (default : FastTree).
 */
process get_100_bootstrap_rep_trees{

  when:
  params.boot
  
  input:
      file(seq_file) from norm_boot_replicates
  output:
      file "${seq_file}.bootstrap.tree" into norm_boot_rep_trees
  
  script:
      output_tree="${seq_file}.bootstrap.tree"
      boot_num="${params.boot_num}"
      template "tree_bootstrap_commands"
} 


/*
 * Step 10. Bootstrap values estimation. Calculates the bootstrap values for each replicate tree, using "Fast Tree-Comparison Tools" from FastTree package.
 */
norm_boot_rep_trees.map { file -> tuple(get_bootTree_prefix(file.name), file) }
         .groupTuple()
         .flatMap { prefix, files -> combine_boot_trees(files) }
         .set{ all_norm_boot_rep_trees }

process get_bootstrap_trees{

  when:
  params.boot

  input:
      set file(in_tree_file), file(all_tree_file) from all_norm_boot_rep_trees
      
  output:
      file "*.bootstrap.tree" into bootstrap_trees
      file "${all_tree_file}"
  
  shell:
  '''
      tmp_name=`basename !{in_tree_file} | awk -F. '{print $1}'`
      perl -I !{baseDir}/bin/ !{baseDir}/bin/CompareToBootstrap.pl -tree !{in_tree_file} -boot !{all_tree_file} > $tmp_name.bootstrap.tree
  '''
} 




/**************************************************************************
 *                                                                        *
 *            The next processes estimate different statistics.           *
 *           They only run when the option "stats" is specified.          *
 *                                                                        *
 **************************************************************************/


/*
 * Step 11. MSA stability estimation. Estimates the MSA stability of each dataset. 
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
 * Step 12. Tree stability estimation. Estimates the Tree stability of each the replicate trees. 
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
 * Step 13. Average bootstrap value estimation. The average unistrap value from all the unistrap trees is estimated. 
 */
process get_tree_bootstrap_stats{

  when:
  params.stats && params.boot

  input:
      set prefix,file(msa) from bootstrap_trees.map{ file -> tuple(get_tree_prefix(file.name), file) }.groupTuple()

  output:
      file "${prefix}_insta.txt" into tree_avg_bootstrap

  shell:
  '''
      echo -n !{prefix}"\t" >> !{prefix}_insta.txt 
      for r in {1..!{params.rep_num}}; { cat !{prefix}_$r.bootstrap.tree | sed 's/);//g' |sed 's/):0.0//g' | sed 's/)/\\n/g' | awk -F: '{print "|"$1}' | grep -v "(" | awk -F"|" '{ sum+=$2} END { print sum/NR }';  } | awk '{ sum+=$1} END { printf sum/NR"\\n"}' >> !{prefix}_insta.txt  
  '''
} 

tree_avg_bootstrap.collectFile(name:'TREE_AVG_BOOTSTRAP.stats', seed:"Name\ttreeAvgBootstrap\n", storeDir:params.out_dir) 



/*
 * Step 14. Average unistrap value estimation. The average unistrap value from all the unistrap trees is estimated. 
 */
process get_tree_unistrap_stats{

  when:
  params.stats 

  input:
      set prefix,file(msa) from unistrap_trees.map{ file -> tuple(get_tree_prefix(file.name), file) }.groupTuple()

  output:
      file "${prefix}_insta.txt" into tree_avg_unistrap

  shell:
  '''
      echo -n !{prefix}"\t" >> !{prefix}_insta.txt 
      for r in {1..!{params.rep_num}}; { cat !{prefix}_$r.unistrap.tree | sed 's/);//g' |sed 's/):0.0//g' | sed 's/)/\\n/g' | awk -F: '{print "|"$1}' | grep -v "(" | awk -F"|" '{ sum+=$2} END { print sum/NR }';  } | awk '{ sum+=$1} END { printf sum/NR"\\n"}' >> !{prefix}_insta.txt  
  '''
} 

tree_avg_unistrap.collectFile(name:'TREE_AVG_UNISTRAP.stats', seed:"Name\ttreeAvgUnistrap\n", storeDir:params.out_dir) 



/*
 * Step 15. MSA identity and number of sequences estimation. Estimates the MSA stability of each dataset and reports how many sequences each dataset contains.
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
      grep -P '^\\s+\\d+\\s+\\d+' !{prefix}_*.fa.aln | awk '{printf $1"\t"}'  >> !{prefix}_insta.txt 
      msaa -a !{prefix}_*.fa.aln -i | grep aln | awk '{print $2}' >> !{prefix}_insta.txt 
  '''
} 

genStats.collectFile(name:'GENERAL.stats', seed:"Name\tseqNum\tmsaIdentity\n", storeDir:params.out_dir) 


workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
