#to get header
bcftools head eg/ERR031940.bcf
bcftools head eg/ERR031940.bcf > header.vcf
#list of samples
 bcftools query -l ERR031940.bcf
 