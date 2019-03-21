# The **Nextflow** folder

This folder contains the pipelines that actually bring everything else together.
Each Nextflow script is placed within its own sub directory to allow praralell execution without the processes blocking the respective usage of *logistic* files/folders (`.nextflow`,`.nextflow.log` & `trace.txt`)

Of course, the **genotyping\*** pipelines need to be run *before* the **analysis\*** pipelines.

---

<p align="center"><img src="../logo.svg" alt="logo" width="150"/></p>