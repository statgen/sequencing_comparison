
process subset {
   errorStrategy "ignore"

   input:
   file vcf from Channel.fromPath(params.vcfs)
   each file(bams_list) from Channel.fromPath(params.bams_list_path)
   each file(reference) from Channel.fromPath(params.reference_path)
   each file(reference_index) from Channel.fromPath(params.reference_path + ".fai")

   output:
   set file("${vcf.baseName}.subset.bcf"), file("${vcf.baseName}.subset.bcf.csi") into subsetted

   publishDir "results/", pattern: "*.subset.bcf*"

   """
   cut -f1 -d" " ${bams_list} > samples.list
   ${params.bcftools} view -S samples.list -c 1 ${vcf} | ${params.bcftools} norm -m -any | ${params.bcftools} norm -f ${reference} | ${params.bcftools} view -v snps,indels -c 1 -Ob -o ${vcf.baseName}.subset.bcf
   ${params.tabix} ${vcf.baseName}.subset.bcf 
   """
}
