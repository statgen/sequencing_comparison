process vep {
   errorStrategy "ignore"

   input:
   file vcf from Channel.fromPath(params.vcfs)

   output:
   file "${vcf}.vep.vcf.gz"
   file "${vcf}.vep.vcf.gz.tbi"

   publishDir "result/vep", pattern: "*.vep.vcf.gz*"

   """
   ${params.bcftools} view -G ${vcf} | ${params.bcftools} annotate -x INFO | ${params.vep} -o ${vcf}.vep.vcf.gz --format vcf --cache --offline --vcf --compress_output bgzip --no_stats ${params.vep_flags}
   ${params.tabix} ${vcf}.vep.vcf.gz
   """
}
