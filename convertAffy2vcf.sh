#!/bin/bash

ref="/gpfs/group/ibb3/default/genome_resources/coral/Adig_NCBI/NCBI_assembly/GCF_000222465.1_Adig_1.1_genomic.fna"
dir="/gpfs/group/ibb3/default/shallow_genome_SNP_data/chipRawData_database/PRO100175_PSU175_SAX_b01/cluster_coral/"
platenum=PRO100175_PSU175_SAX_b01
annot_file=/gpfs/group/ibb3/default/shallow_genome_SNP_data/chipRawData_database/PRO100175_PSU175_SAX_b01/AxiomReference/Axiom_AcropSNP_coral_Annotation.r1.csv
out="/gpfs/group/ibb3/default/shallow_genome_SNP_data/chipRawData_database/PRO100175_PSU175_SAX_b01/"

/gpfs/group/ibb3/default/tools/tmp/bcftools/bcftools +/gpfs/group/ibb3/default/tools/tmp/bcftools/plugins/affy2vcf.so --no-version -Ov --fasta-ref $ref \
  --annot $annot_file \
  --threads 6 \
  --snp-posteriors $dir/AxiomGT1.snp-posteriors.txt \
  --summary $dir/AxiomGT1.summary.txt \
  --report $dir/AxiomGT1.report.txt \
  --calls $dir/AxiomGT1.calls.txt \
  --confidences $dir/AxiomGT1.confidences.txt \
  --output $out/$platenum.vcf
