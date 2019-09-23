vcf_with_index = Channel.fromPath(params.coverage_files_path).map { f -> [ f, f + params.coverage_files_index_suffix ]}

process histograms {
   errorStrategy "ignore"

   input:
   set file(vcf), file(index) from vcf_with_index
   each file(gtf) from Channel.fromPath(params.gencode_gtf_path)

   output:
   file "${vcf.baseName}.histograms.json.gz"

   publishDir "result", pattern: "*.histograms.json.gz"

   """
   ${params.histograms} -g ${gtf} -d ${vcf} -o ${vcf.baseName}.histograms.json.gz
   """
}
