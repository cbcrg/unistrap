#!/usr/bin/env nextflow

//nextflow run shootstrap.nf --in_tree hip_clustalo_100.phy.replicate.tree

params.in_dir="$baseDir/data/dataset/*"
params.out_dir="Shootstrap_Analysis_Results"
params.rep_num=2
params.seed=10
params.aligner="clustalo"
params.in_tree="" 

file_names=Channel.fromPath(params.in_dir)
in_tree_file = params.in_tree ? file(params.in_tree) : null
if( in_tree_file ) assert in_tree_file.exists(), "The tree file does not exist: $in_tree_file !!!" 

process get_shuffle_replicates{
  publishDir params.out_dir, mode: 'copy'

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
  publishDir params.out_dir, mode: 'copy'

  input:
      file(seq_file) from shuffle_replicates
  output:
      file "${seq_file}.aln" into msa_replicates, msa_replicates2
  
  script:
      template "${params.aligner}_msa_command"
}

process get_msa_trees{
  publishDir params.out_dir, mode: 'copy'

  input:
      file(seq_file) from msa_replicates
  output:
      file "${seq_file}.msa_tree" into msa_trees, msa_trees2
  
  script:
      output_tree="${seq_file}.msa_tree"
      template "tree_commands"
}


process get_stable_msa_trees{
  publishDir params.out_dir, mode: 'copy'

  input:
      file(in_tree_file) from msa_trees
      file(all_tree_file) from msa_trees2.collectFile(name: 'big.tree').first()
      
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
  publishDir params.out_dir, mode: 'copy'

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
  publishDir params.out_dir, mode: 'copy'

  input:
      file(seq_file) from replicates
  output:
      file "${seq_file}.tree" into trees
  
  script:
      output_tree="${seq_file}.tree"
      template "tree_commands"
} 


process get_shootstrap_tree{
  publishDir params.out_dir, mode: 'copy'

  input:
      file(all_tree_file) from trees.collectFile(name: 'big.tree')
      set val(score), file(tree) from in_tree_file ? Channel.value([100, in_tree_file]) : most_stable_tree.max { it[0].trim().toFloat() }
  output:
      file "*.shootstrap.tree" into shootstrap_tree
  
  shell:
  '''
      tmp_name=`basename !{tree} | awk -F. '{print $1}'`
      perl -I !{baseDir}/bin/ !{baseDir}/bin/CompareToBootstrap.pl -tree !{tree} -boot !{all_tree_file} > $tmp_name.shootstrap.tree
  '''
} 
