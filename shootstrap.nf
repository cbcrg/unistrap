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
params.rep_num=2
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



def get_aln_prefix(name) { 
  def m = name =~ /(.*)_\d+\.fa\.aln/
  m.matches() ? m.group(1) : 'unknown' 
}

def get_tree_prefix(name) { 
  def m = name =~ /(.*)_\d+\.fa\.aln.*/
  m.matches() ? m.group(1) : 'unknown_tree_name' 
}



def combine_trees( allFiles ) {
  def prefix = get_tree_prefix(allFiles[0].name)
  println(prefix)
  def bigTree = java.nio.file.Files.createTempFile(prefix,".tree")
 
  def result = []
  allFiles.each {  
    bigTree << it.text  
    result << [ it, bigTree ]
  }
  result
}





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

process get_msa_replicates{

  input:
      file(seq_file) from shuffle_replicates
  output:
      file "${seq_file}.aln" into msa_replicates, msa_replicates2, msa_replicates3
  
  script:
      template "${params.aligner}_msa_command"
}


process get_msa_trees{

  input:
      file(seq_file) from msa_replicates
  output:
      file "${seq_file}.msa_tree" into msa_trees, msa_trees2
  
  script:
      output_tree="${seq_file}.msa_tree"
      template "tree_commands"
}


msa_trees.map { file -> tuple(file.name[0], file) }
         .groupTuple()
         .flatMap { prefix, files -> combine_trees(files) }
         .view()
         .set{ msa_all_trees }


process get_stable_msa_trees{

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


// Create X replicates for each MSA replicate
process get_seqboot_replicates{

  input:
      file(seq_file) from msa_replicates2
  output:
      file "${seq_file}.replicate" into replicates
  
  """
      echo -e "$seq_file\nR\n1\nY\n31\n"|seqboot
      mv outfile ${seq_file}.replicate
  """

}


process get_replicate_trees{

  input:
      file(seq_file) from replicates
  output:
      file "${seq_file}.tree" into trees
  
  script:
      output_tree="${seq_file}.tree"
      template "tree_commands"
} 



trees.map { file -> tuple(get_tree_prefix(file.name), file) }
         .groupTuple()
         .flatMap { prefix, files -> combine_trees(files) }
         .view()
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


process get_stability_stats{

  when:
  params.stats 

  input:
      set prefix,file(msa) from msa_replicates3.map{ file -> tuple(get_aln_prefix(file.name), file) }.groupTuple()

  output:
      file "${prefix}_insta.txt" 
  

  shell:
  '''
      echo -n !{prefix}" " >> !{prefix}_insta.txt 
      for j in {1..!{params.rep_num}};{ for r in {1..!{params.rep_num}};{ [ $j != $r ] && msaa -a !{prefix}_$j.fa.aln -r !{prefix}_$r.fa.aln; }; } |  awk '{ sum += $2 } END { if(NR>0) print sum/NR }' >> !{prefix}_insta.txt  
  '''
} 
