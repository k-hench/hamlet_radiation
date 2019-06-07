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
	win = [];
	step = [];
	try:
		opts, args = getopt.getopt(argv,"hv:p:w:s:",["vcffile=","popfile="])
	except getopt.GetoptError:
		print( 'multipop_fst_windows.py -v <vcf-file> -p <pop-file> -w <window-size> -s <step/increment-size>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print( 'multipop_fst_windows.py -v <vcf-file> -p <pop-file> -w <window-size> -s <step/increment-size>');
			sys.exit()
		elif opt in ("-v", "--vcf"):
			vcffile = arg
		elif opt in ("-p", "--pop"):
			popfile = arg
		elif opt in ("-w", "--window"):
			win = int(arg)
		elif opt in ("-s", "--step"):
			step = int(arg)

	# the "welcome screen"
	eprint("\n==== \033[1;31mmultipop_fst\033[1;m ====")
	eprint("====== (Windows) =====")
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
	# determine all present linkage groups
	chrom_list = pd.Series(callset['variants/CHROM'], dtype="category")
	chrom_val = chrom_list.cat.categories

	# define the output column order
	sorter = ['chrom', 'start', 'end', 'n_snps', 'fst']
	# init the output data frame
	df = pd.DataFrame({'chrom':[], 'start':[], 'end':[], 'n_snps':[], 'fst':[]})[sorter]
	# init the output data frame column precursor
	chroms =  list()
	starts = list();
	ends = list();
	n_snpss = list();
	fsts = list();
	# compute the windowed fst for each linkage group
	for chrom in chrom_val:
		# report current lg
		eprint("Working on: \033[1;34m",chrom,"\033[1;m", sep='')
		# determine relevant indices
		current_rows = [callset['variants/CHROM'] == chrom];
		# select relevant positions
		pos = callset['variants/POS'][tuple(current_rows)]
		# select relevant genotypes
		g = callset['calldata/GT'][tuple(current_rows)]
		# compute the window statistics
		aw, bw, cw = allel.windowed_weir_cockerham_fst(pos, g, subpops, size=win, start=0, step=step)
		# derive the number of windows
		n_windows = len(aw)
		# append the lg specific results
		chroms.append( [chrom] * n_windows );
		starts.append(bw[:, 0]+1);
		ends.append(bw[:, 1]+1);
		n_snpss.append(cw);
		fsts.append(aw);

	# "flatten" the individual columns and compose output table
	df = pd.DataFrame({'chrom':[item for sublist in chroms for item in sublist],
	'start':[item for sublist in starts for item in sublist],
	'end':[item for sublist in ends for item in sublist],
	'n_snps':[item for sublist in n_snpss for item in sublist],
	'fst':[round(item,5) for sublist in fsts for item in sublist]})[sorter]
	# print ouput
	df[df.n_snps != 0].to_csv(sys.stdout, index=False, sep='\t')
	# goodbye message
	eprint("======== \033[1;31mDone\033[1;m! =======")

if __name__ == "__main__":
	main(sys.argv[1:])