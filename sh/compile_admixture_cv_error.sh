cd 2_analysis/admixture
echo -e "gid\tk\tcv_error" > cv_errors.tsv
grep CV log*.out | sed "s/.[0-9]*.out:CV error (K=/\t/; s/log.//; s/): /\t/" >> cv_errors.tsv
