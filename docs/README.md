## Summary of the genotyping workflow

Graphical overview of the genotyping.nf pipeline:

<img src="genotyping.png" title="plot of chunk plot1" alt="plot of chunk plot1" width="95%" />

---

Short summary of the SNP filtering prcess:

| Stage	| SNPs (n) | Comment |
|--------:| --------:| --------:|
| intermediate.vcf.gz | 599636340 | raw GATK output (all Variants) |
| raw_var_sites.vcf.gz | 48638754 | raw GATK output (only SNPs) |
| intermediate.filterd.vcf.gz | 31046775 | filtered using GATK stats |
| filterd_bi-allelic.vcf.gz | 28912732 | filtered by VCFtools (max-missing: 17, max-alleles: 2) |
| phased.vcf.gz | 27537293 | raw phased |
| phased_mac2.vcf.gz | 23344192 | phased and filterd (mac >= 2) |

---

<center><img src="../logo.svg" alt="logo" width="150"/></center>