#!/usr/bin/python
import sys, getopt
import numpy as np
import scipy
import pandas as pd
import allel

# define error message
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# the main function
def main(argv):
	vcffile = ''
	popfile = ''
	try:
		opts, args = getopt.getopt(argv,"hv:p:",["vcffile=","popfile="])
	except getopt.GetoptError:
		print( 'multipop_fst_per_snp.py -v <vcf-file> -p <pop-file> ')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print( 'multipop_fst_per_snp.py -v <vcf-file> -p <pop-file>');
			sys.exit()
		elif opt in ("-v", "--vcf"):
			vcffile = arg
		elif opt in ("-p", "--pop"):
			popfile = arg

	# the "welcome screen"
	eprint("\n==== \033[1;31mmultipop_fst\033[1;m ====")
	eprint("======= (SNPs) =======")
	eprint("== \033[1;31mmulti_fst\033[1;m (SNPs) ==")
	eprint("Genotypes: ----------> \033[1;34m",vcffile,"\033[1;m", sep ='')
	eprint("Pop-File: -----------> \033[1;34m",popfile,"\033[1;m", sep ='')
	eprint("======================")
	# open the vcf file
	callset = allel.read_vcf(vcffile)
	# open the pop file
	pops = pd.read_csv(popfile, sep='\t', names=['ind','pop'])
	# determine all pops present in popfile
	pops_list = pd.Series(pops['pop'], dtype="category")
	pops_val = pops_list.cat.categories
	# init the sample indices for each pop
	subpops =  list()
	# fill the sample indices for each pop
	for pop in pops_val:
		eprint("pop: \033[1;34m",pop,"\033[1;m", sep ='')
		current_pop = [i for i, e in enumerate(pops_list) if e == pop]
		subpops.append( current_pop );

	eprint("======================")
	# store linkage groups, positions and genotypes
	chrom = callset['variants/CHROM']
	pos = callset['variants/POS']
	g = callset['calldata/GT']
	# compute f stats per SNP
	a, b, c = allel.weir_cockerham_fst(g, subpops)
	# summarise global fst
	global_fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))
	# compute per SNP fst
	fst = (np.sum(a, axis=1) /(np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))
	# define the output column order
	sorter = ['chrom', 'pos', 'fst']
	# compose output table
	df = pd.DataFrame({'chrom':chrom, 'pos':pos, 'fst':[round(item,5) for item in fst]})[sorter]
	# print ouput
	df.to_csv(sys.stdout, index=False, sep='\t')
	# report global fst to stderr
	eprint("Global Fst:    ",round(global_fst,5),sep='')
	# goodbye message
	eprint("======== \033[1;31mDone\033[1;m! =======")

if __name__ == "__main__":
	main(sys.argv[1:])



