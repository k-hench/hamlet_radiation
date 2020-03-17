# The resources folder

This folder contains mostly *external data* which are not novel but inherited from previous work.

The hamlet genome files are originally from the Hench *et al* 2019 ([doi: 10.1038/s41559-019-0814-5](https://doi.org/10.1038/s41559-019-0814-5)),
while the mappability masks are from Moran *et al* 2019 ([doi: 10.1111/mec.15110](https://doi.org/10.1111/mec.15110)).

The subfolder `other_studies` contains data from studies by other authors:
The file `han_etal_2017.tsv` contains *F<sub>ST</sub>* data of Darwin’s finches published by Han *et al.*, while the file `stankowski_etal_2019_spec_short.tsv` contains the populations assignment of the samples from Stankowski *et al* (2019) study.
This folder should also contain the genotype file `all_9_taxa_G1_vars.vcf.gz` form that study which is not tracked here but can be downloaded from the authors dryad repository ([doi: 10.5061/dryad.18j0j5q](https://doi.org/10.5061/dryad.18j0j5q)).
Prior to usage I un- and rezipped this file using `bgzip` to create `all_9_taxa_G1_vars.bgziped.vcf.gz` as input for the pipeline `nf/analysis_stankowski_etal_2019/analysis_stankowski_etal_2019.nf`.

Apart from this, the `ressources` folder also contains some intermediate results from this study that were reused a later stages (this was mainly an awkward fix to solve some cluster incompatibilities with some software used in some pipelines or a manual resurrection of intermediate steps after a pipeline had finished).
This refers specifically to the sub-folders `coverage_masks`, `genotypes`, `indel_masks` and `plugin`.


These files are not tracked to minimize the repository size but the following structure is assumed:

```
.
|-- HP_genome_unmasked_01.dict
|-- HP_genome_unmasked_01.fa
|-- HP_genome_unmasked_01.fa.amb
|-- HP_genome_unmasked_01.fa.ann
|-- HP_genome_unmasked_01.fa.bwt
|-- HP_genome_unmasked_01.fa.fai
|-- HP_genome_unmasked_01.fa.gz
|-- HP_genome_unmasked_01.fa.pac
|-- HP_genome_unmasked_01.fa.sa
|-- README.md
|-- bams
|-- coverage_masks
|   |-- 16_21-30nigpan.LG01.coverage_mask.bed.gz
|   |-- ...
|   `-- PL17_95maybel.LG24.coverage_mask.bed.gz
|-- genotypes
|   |-- gvcfs
|   |   |-- 16_21-30nigpan.g.vcf.gz
|   |   |-- ...
|   |   `-- s_tort_3torpan.g.vcf.gz
|   |-- raw_var_sites.vcf.gz
|   `-- raw_var_sites.vcf.gz.tbi
|-- img
|   |-- bars_r.c.svg
|   |-- finch.c.svg
|   |-- flycatcher.c.svg
|   |-- generic_hamlet.c.svg
|   |-- heliconius.c.svg
|   |-- monkeyflower.c.svg
|   |-- peduncle_r.c.svg
|   |-- snout_r.c.svg
|   `-- sunflower.c.svg
|-- indel_masks
|   |-- indel_mask.LG01.bed.gz
|   |-- ...
|   `-- indel_mask.LG24.bed.gz
|-- mappability_masks
|   |-- LG01.mapmask.bed.txt.gz
|   |-- ...
|   |-- LG24.mapmask.bed.txt.gz
|   `-- LG_M.mapmask.bed.txt.gz
|-- other_studies
|   |-- all_9_taxa_G1_vars.bgziped.vcf.gz
|   |-- all_9_taxa_G1_vars.bgziped.vcf.gz.tbi
|   |-- han_etal_2017.tsv
|   |-- stankowski_etal_2019.tsv
|   `-- stankowski_etal_2019_spec_short.tsv
`-- plugin
    |-- indel_masks
    |   |-- indel_mask.LG01.bed.gz
    |   |-- ...
    |   `-- indel_mask.LG24.bed.gz
    `-- trees
        |-- bel
        |   |-- bel.LG01.w200.phyml_bionj.data.tsv
        |   |-- bel.LG01.w200.phyml_bionj.trees.gz
        |   |-- ...
        |   |-- bel.LG24.w50.phyml_bionj.data.tsv
        |   |-- bel.LG24.w50.phyml_bionj.trees.gz
        |   |-- bel.phyml_bionj.data.tsv
        |   `-- bel.phyml_bionj.trees.gz
        `-- hon
            |-- hon.LG01.w200.phyml_bionj.data.tsv
            |-- hon.LG01.w200.phyml_bionj.trees.gz
            |-- ...
            |-- hon.LG24.w50.phyml_bionj.data.tsv
            |-- hon.LG24.w50.phyml_bionj.trees.gz
            |-- hon.phyml_bionj.data.tsv
            `-- hon.phyml_bionj.trees.gz

13 directories, 4448 files
```

References:

Han, F. *et al.* Gene flow, ancient polymorphism, and ecological adaptation shape the genomic landscape of divergence among darwin’s finches. Genome research 27, 1004–1015 (2017).

Stankowski, S. *et al.* Widespread selection and gene flow shape the genomic landscape during a radiation of monkeyflowers. PLoS biology 17, e3000391 (2019).

---

<p align="center"><img src="../logo.svg" alt="logo" width="150"/></p>