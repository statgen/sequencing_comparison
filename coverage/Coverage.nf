bams = Channel.from(file(params.bams_list_path).readLines()).map { line -> fields = line.split();  [ file(fields[1]), file(fields[1] + (params.cram ? ".crai" : ".bai"))] }
intervals = Channel.fromPath(params.intervals).map { file -> tokens = file.baseName.split('-'); [ tokens[0], tokens[1], tokens[2], file ] }

grouped = intervals.combine(bams).groupTuple(by: [0, 1, 2, 3])

process mpileup {
   errorStrategy "ignore"

   input:
   set val(chrom), val(start), val(stop), file(intervals), file(bams), file(bai) from grouped
   each file(bams_list) from Channel.fromPath(params.bams_list_path)
   each file(reference) from Channel.fromPath(params.reference_path)
   each file(reference_index) from Channel.fromPath(params.reference_path + ".fai")

   output:
   set file("${chrom}-${start}-${stop}.DP.bcf"), file("${chrom}-${start}-${stop}.DP.bcf.csi") into interval_bcfs

   publishDir "results/pileup", pattern: "*.DP.bcf*"

   """
   samples=`cut -f1 -d" " ${bams_list} | tr "\n" " "`
   bam_names=`cut -f2 -d" " ${bams_list} | xargs -L1 basename | tr "\n" " "`
   ${params.samtools} mpileup -f ${reference} -E -q 20 -Q 20 --ff 0x0704 -d ${params.max_depth} -r ${chrom}:${start}-${stop} -l ${intervals} \${bam_names} | ${params.mpileup_subset} -p - -s \${samples} -o ${chrom}-${start}-${stop}.DP.bcf
   ${params.tabix} ${chrom}-${start}-${stop}.DP.bcf
   """
}

process merge {
   errorStrategy "ignore"

   input:
   set val(chrom), file(bcfs), file(csi) from interval_bcfs.map { item -> [ item[0].baseName.split('-')[0], item[0], item[1] ] }.groupTuple(by: 0)

   output:
   set file("${chrom}.merged.DP.bcf"), file("${chrom}.merged.DP.bcf.csi") into merged_bcfs

   publishDir "results/merged", pattern: "*.merged.DP.bcf*"

   """
   find . -name "${chrom}-*.DP.bcf" | sort -V > files.list
   ${params.bcftools} concat -f files.list -Ob -o ${chrom}.merged.DP.bcf
   ${params.tabix} ${chrom}.merged.DP.bcf
   """
}
