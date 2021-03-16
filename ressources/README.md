# The resources folder

This folder contains mostly *external data* which are not novel but inherited from previous work.

The hamlet genome files are originally from the Hench *et al* 2019 ([doi: 10.1038/s41559-019-0814-5](https://doi.org/10.1038/s41559-019-0814-5)),
while the mappability masks are from Moran *et al* 2019 ([doi: 10.1111/mec.15110](https://doi.org/10.1111/mec.15110)).

Furthermore, panel **a** of **figure 1** is based on data and scripts by Rabosky *et al.* (2018). For this, the corresponding Dryad repository ([doi: 10.5061/dryad.fc71cp4](doi.org/10.5061/dryad.fc71cp4)) was downloaded and decompressed into the sub-folder `ressources/Rabosky_etal_2018`.


Apart from this, the `ressources` folder also contains some intermediate results from this study that were reused a later stages (this was mainly an awkward fix to solve some cluster incompatibilities with some software used in some pipelines or a manual resurrection of intermediate steps after a pipeline had finished).
This refers specifically to the sub-folders `coverage_masks`, `genotypes`, `indel_masks` and `plugin`.


These files are not tracked to minimize the repository size but the following structure is assumed:

```
.
├── HP_genome_unmasked_01.dict
├── HP_genome_unmasked_01.fa
├── HP_genome_unmasked_01.fa.amb
├── HP_genome_unmasked_01.fa.ann
├── HP_genome_unmasked_01.fa.bwt
├── HP_genome_unmasked_01.fa.fai
├── HP_genome_unmasked_01.fa.gz
├── HP_genome_unmasked_01.fa.pac
├── HP_genome_unmasked_01.fa.sa
├── README.md
├── bams
├── coverage_masks
|   ├── 16_21-30nigpan.LG01.coverage_mask.bed.gz
|   ├── ...
|   └── PL17_95maybel.LG24.coverage_mask.bed.gz
├── genotypes
|   ├── gvcfs
|   |   ├── 16_21-30nigpan.g.vcf.gz
|   |   ├── ...
|   |   └── s_tort_3torpan.g.vcf.gz
|   ├── raw_var_sites.vcf.gz
|   └── raw_var_sites.vcf.gz.tbi
├── indel_masks
|   ├── indel_mask.LG01.bed.gz
|   ├── ...
|   └── indel_mask.LG24.bed.gz
├── mappability_masks
|   ├── LG01.mapmask.bed.txt.gz
|   ├── ...
|   ├── LG24.mapmask.bed.txt.gz
|   └── LG_M.mapmask.bed.txt.gz
├── plugin
│   ├── indel_masks
│   |   ├── indel_mask.LG01.bed.gz
│   |   ├── ...
│   |   └── indel_mask.LG24.bed.gz
│   └── trees
│       ├── bel
│       |   ├── bel.LG01.w200.phyml_bionj.data.tsv
│       |   ├── bel.LG01.w200.phyml_bionj.trees.gz
│       |   ├── ...
│       |   ├── bel.LG24.w50.phyml_bionj.data.tsv
│       |   ├── bel.LG24.w50.phyml_bionj.trees.gz
│       |   ├── bel.phyml_bionj.data.tsv
│       |   └── bel.phyml_bionj.trees.gz
│       └── hon
│           ├── hon.LG01.w200.phyml_bionj.data.tsv
│           ├── hon.LG01.w200.phyml_bionj.trees.gz
│           ├── ...
│           ├── hon.LG24.w50.phyml_bionj.data.tsv
│           ├── hon.LG24.w50.phyml_bionj.trees.gz
│           ├── hon.phyml_bionj.data.tsv
│           └── hon.phyml_bionj.trees.gz
├── Rabosky_etal_2018
│   ├── dataFiles
│   │   └── ...
│   ├── __MACOSX
│   │   └── ...
│   ├── readme.txt
│   └── scripts
│       └── ...
├── README.md
├── vcf2gp.spid
├── vcf2nh.spid
└── windows_1kb.bed.gz
```

---

**References:**

Hench, K., Vargas, M., Höppner, M. P., McMillan, W. O., & Puebla, O. (2019). Inter-chromosomal coupling between vision and pigmentation genes during genomic divergence. **Nature Ecology & Evolution**, 3 (4), 657–667. https://doi.org/10.1038/s41559-019-0814-5

Hench, K., Vargas, M., Höppner, M. P., McMillan, W. O., & Puebla, O. (2019). Data from: Inter-chromosomal coupling between vision and pigmentation genes during genomic divergence. **Dryad**. https://doi.org/10.5061/DRYAD.PG8Q56G

Moran, B. M., Hench, K., Waples, R. S., Höppner, M. P., Baldwin, C. C., McMillan, W. O., & Puebla, O. (2019). The evolution of microendemism in a reef fish (Hypoplectrus maya). **Molecular Ecology**, 28 (11), 2872–2885. https://doi.org/10.1111/mec.15110

Moran, B. M., Hench, K., Waples, R. S., Höppner, M. P., Baldwin, C. C., McMillan, W. O., & Puebla, O. (2019). Data from: The evolution of microendemism in a reef fish (Hypoplectrus maya). **Dryad**. https://doi.org/10.5061/DRYAD.HP388DM

Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kaschner, K., Garilao, C., Near, T. J., Coll, M., & Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. **Nature**, 559 (7714), 392–395. https://doi.org/10.1038/s41586-018-0273-1

Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kaschner, K., Garilao, C., Near, T. J., Coll, M., & Alfaro, M. E. (2019). Data from: An inverse latitudinal gradient in speciation rate for marine fishes. **Dryad**. https://doi.org/10.5061/DRYAD.FC71CP4

---

<p align="center"><img src="../logo.svg" alt="logo" width="150"/></p>
