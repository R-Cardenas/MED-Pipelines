for f in 'GY1204-S12' 'NN5301-S53' 'S0102P-S01' 'S0108-S01' 'S0501-S05' 'S0502-S05' 'S0505-S05' 'S0801-S08' 'S0803-S08' 'S0901-S09' 'S0904-S09' 'S1001-S10' 'S1003-S10' 'S1101-S11' 'S1103-S11' 'S1201-S12' 'S5303-S53' 'WX0701-WX' 'WX0703-WX'; do bcftools concat -a -O z -o $f-merge.vcf.gz $f-caller.indels.slivar.vcf.gz $f-callermerged.snps.slivar.vcf.gz; done
bcftools merge -O z -o choleB2-all-merge.vcf.gz *merge.vcf.gz

for f in 'S204-S2' 'S205-S2'; do bcftools concat -a -O z -o $f-merge.vcf.gz $f-caller.indels.slivar.vcf.gz $f-callermerged.snps.slivar.vcf.gz; done
bcftools merge -O z -o choleB1-all-merge.vcf.gz *merge.vcf.gz
