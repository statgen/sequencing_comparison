pairs = Channel.from(file(params.pairs_list_path).readLines()).map { line -> fields = line.split();  [ fields[0], fields[1] ] }
files1 = Channel.from(file(params.study1_files_list_path).readLines()).map { line -> fields = line.split(); [ fields[0], fields[1], fields[2] ] }
files2 = Channel.from(file(params.study2_files_list_path).readLines()).map { line -> fields = line.split(); [ fields[0], fields[1], fields[2] ] }


process compare {
   errorStrategy "ignore"

   input:
   set val(sample1), val(sample2), val(dp1), val(gt1), val(vep1), val(dp2), val(gt2), val(vep2) from pairs.combine(files1.merge(files2))

   output:
   set file("${sample1}.summary.gz"), file("${sample1}.summary.log") into compared

   publishDir "results", pattern: "*.summary.*"

   """
   ${params.compare} -s1 ${sample1} -s2 ${sample2} -d1 ${dp1} -g1 ${gt1} -a1 ${vep1} -d2 ${dp2} -g2 ${gt2} -a2 ${vep2} -o ${sample1}.summary.gz > ${sample1}.summary.log 
   """
}
