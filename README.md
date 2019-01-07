# The **chapter2** repository

This project contains the analysis of phylogenetic relationships between hamlet species (*Hypoplectrus* spp) from Belize, Honduras and Panama.

Steps involved include:

- Genotyping from raw sequencing reads (`*.fq.gz` ==> `phased_mac2.vcf.gz`) using [**GATK4**](https://software.broadinstitute.org/gatk/) and phasing using [**shapeit2**](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
- Analysis of divergence (*F<sub>ST</sub>*, *D<sub>xy</sub>*, G x P) using [**VCFtools**](https://vcftools.github.io/examples.html), **popgenWindows.py** (by [Simon Martin](https://github.com/simonhmartin/genomics_general)) and [**GEMMA**](https://github.com/genetics-statistics/GEMMA)
- cluster analysis (*PCA* and *phylogenetic trees*) using the [**R**](https://cran.r-project.org/) package [**SNPrelate**](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html), [**PhyML**](http://www.atgc-montpellier.fr/index.php?type=bn)/[**Fasttree2**](http://www.microbesonline.org/fasttree/) in combination with [**twisst**](https://github.com/simonhmartin/twisst), and [**structure**](https://web.stanford.edu/group/pritchardlab/structure.html)/[**admixture**](https://www.genetics.ucla.edu/software/admixture/index.html) (broken link)

The whole process is monitored using [**Nextflow**](https://www.nextflow.io/index.html) and distributed over several pipelines:

- genotyping (including varian calling and phasing) is done in `./genotyping.nf`
- population genetic analysis happens in `./analysis.nf`

Additional information about the input data is located in the `./metadata` folder,
a summary of the nextflow processes can be found in the `./docs` folder.

---

<center><img src="logo.svg" alt="logo" width="150"/></center>
