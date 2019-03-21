# The **python** folder

This folder contains analysis run in **python**.
The scripts are run from within the **nextflow** pipelines in the root folder.

## Something is missing here:

For `bamHamletCaller.py` to run, the script `utils.py` from the repo [**msmc-tools**](https://github.com/stschiff/msmc-tools) by Stephan Schiffels needs to be located here.
Since I do not want to redistribute this myself, you will need to download it on your own and crate a link at this location.

This could for example look like this:

```sh
# assuming $SOFTWARE_DIR is the location you would usually install software:
cd $SOFTWARE_DIR
git clone https://github.com/stschiff/msmc-tools.git

cd -
ln -s $SOFTWARE_DIR/msmc-tools/utils.py ./
```
---

<p align="center"><img src="../logo.svg" alt="logo" width="150"/></p>