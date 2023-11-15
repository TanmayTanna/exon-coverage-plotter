### Tanmay Tanna #

import HTSeq
import numpy
from matplotlib import pyplot
import argparse
import csv
import pandas as pd
import math
import os

parser = argparse.ArgumentParser()

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')  # change grouping of argument to required or optional in help through this chunk of code

### inputs ###

# required 
required.add_argument('-i','--inbamlist',  nargs='+', help='array of bam files with path', required=True)
required.add_argument('-g', '--genome', help='genome annotation file path (gff/gtf)', required=True)


# optional
optional.add_argument('-o', '--outPath', help='path to output directory.', default='.')
optional.add_argument('-w', '--winwidth', help='distance upstream and downstream of exon to be displayed Default=300 bp', type=int, dest='winwidth', default=300)
optional.add_argument('-e', '--exonsize', help='standardized exon length - to be used for condensing coverage for all exons to a single length for output plots Default=1000', type=int, dest='exonsize', default=1000)
optional.add_argument('-f', '--feature', help='gff3 feature for plotting coverage', default='gene', dest='feature')
#initialize values

parser._action_groups.append(optional) 
args = parser.parse_args()

bamlist=args.inbamlist
gtffile = HTSeq.GFF_Reader(str(args.genome))
outdir=str(args.outPath)
feature=str(args.feature)
winwidth = int(args.winwidth)
exonsize=int(args.exonsize)
if not os.path.exists(outdir):
	os.mkdir(outdir)

x_coord=numpy.arange(-winwidth, exonsize+winwidth) # user defined or default size of coverage plot
# x_coord_smooth = numpy.linspace( -winwidth, exonsize+winwidth - 1,300) # x coordinates used for smoothing
# adding exon start and end to smoothing coordinates
# if 0 not in x_coord_smooth:
#     x_coord_smooth=numpy.append(x_coord_smooth, [0], axis=0)
# if exonsize not in x_coord_smooth:
#     x_coord_smooth=numpy.append(x_coord_smooth, [exonsize], axis=0)
# x_coord_smooth = numpy.sort(x_coord_smooth)

cvg_list_global = []
filenumber=0 # file counter

### looping through for each bam file ###
for file in bamlist:
	print("Processing: " + file)
	filenumber+=1
	bamfile = HTSeq.BAM_Reader(str(file))
	coverage = HTSeq.GenomicArray("auto", stranded=False, typecode="i") #initializing coverage array
	outfile= file.split("/")[-1]

	readnumber=0 # number of reads
	cvg_list = []

	for almnt in bamfile:
		readnumber+=1
		if almnt.aligned:
			coverage[almnt.iv] += 1  # add 1 to coverage for each aligned "fragment"

	exonpos = {}
	gene_names = []
	for feat in gtffile:
		if feat.type == feature:
			exonpos[feat.attr['gene']] = [feat.iv.start_d_as_pos, feat.iv.end_d_as_pos] # read gene/exon coordinates
	for gene in exonpos:
		p = exonpos[gene]
		if p[0].pos < p[1].pos: # positive strand
			exon = HTSeq.GenomicInterval(p[0].chrom, p[0].pos, p[1].pos, ".")
			upstream = HTSeq.GenomicInterval(p[0].chrom, p[0].pos - winwidth, p[0].pos, ".")
			downstream = HTSeq.GenomicInterval(p[0].chrom, p[1].pos, p[1].pos + winwidth, ".")
			k = (p[1].pos - p[0].pos) # exon length
		else: # negative strand
			exon = HTSeq.GenomicInterval(p[0].chrom, p[1].pos+1, p[0].pos+1, ".")
			upstream = HTSeq.GenomicInterval(p[0].chrom, p[1].pos - winwidth+1, p[1].pos+1, ".")
			downstream = HTSeq.GenomicInterval(p[0].chrom, p[0].pos+1, p[0].pos + winwidth +1, ".")
			k = (p[0].pos - p[1].pos)
		if k < 100:
			continue
		exoncvg = numpy.fromiter(coverage[exon], dtype='f', count=k) # coverage over exon 
		exoncvg_scaled = numpy.zeros(exonsize, dtype='f') # scaling coverage over each exon to user defined or default exonsize for comparability
		k = (k-1)/exonsize # increment to be used for loop
		i = 0 # loop iterators
		j = 0
		while j < exonsize:
			if round(i) != round(i+k):
				exoncvg_scaled[j] = numpy.mean(exoncvg[int(round(i)):int(round(i+k))])
			else:
				exoncvg_scaled[j] = exoncvg[int(round(i))]
			i += k
			j += 1
		# concatenating coverage over upstream and downstream window with scaled exon coverage
		wincvg = numpy.concatenate((numpy.fromiter(coverage[upstream], dtype='f', count=winwidth), exoncvg_scaled), axis=0)  
		wincvg = numpy.concatenate((wincvg, numpy.fromiter(coverage[downstream], dtype='f', count=winwidth)), axis=0)

		# add coverage array to local coverage profile
		if p[1].strand == "+":
			cvg_list.append(wincvg)
			cvg_list_global.append(wincvg)
		else:
			cvg_list.append(wincvg[::-1])
			cvg_list_global.append(wincvg[::-1])
		gene_names.append(gene)

	cvg_df = pd.concat([pd.Series(x) for x in cvg_list], axis=1)
	cvg_df_norm = cvg_df/cvg_df.max()
	profile = cvg_df_norm.mean(axis=1)
	profile_se = cvg_df_norm.sem(axis=1)
	cvg_df_norm.columns = gene_names
	cvg_df_norm.insert(0, "coordinate", x_coord)
	cvg_df_norm.to_csv(outdir+"/"+outfile+"_exoncoverage.csv", index=False)

	### Normalization ###

	# MinMax normalization
	# p_max=numpy.amax(profile_smooth)
	# p_min=numpy.amin(profile_smooth)
	# profile_smooth = (profile_smooth-p_min)/(p_max-p_min)

	# Normalization by max
	# p_max=numpy.amax(profile_smooth)
	# profile_smooth = profile_smooth/p_max

	# Normalization by mean
	# profile *= len(profile)/sum(profile)  

	###This part is for plotting profiles from individual bams to separate files ### Comment out if these aren't needed ###

	pyplot.style.use('seaborn-white')
	pyplot.plot(x_coord, profile, color="blue")

	readnumber=0 # number of reads
	cvg_list = []

	for almnt in bamfile:
		readnumber+=1
		if almnt.aligned:
			coverage[ almnt.iv ] += 1  # add 1 to coverage for each aligned "fragment"

	exonpos = {}
	gene_names = []
	for feat in gtffile:
		if feat.type == feature:
			exonpos[feat.attr['gene']] = [feat.iv.start_d_as_pos, feat.iv.end_d_as_pos] # read gene/exon coordinates
	for gene in exonpos:
		p = exonpos[gene]
		if p[0].pos < p[1].pos: # positive strand
			exon = HTSeq.GenomicInterval(p[0].chrom, p[0].pos, p[1].pos, ".")
			upstream = HTSeq.GenomicInterval(p[0].chrom, p[0].pos - winwidth, p[0].pos, ".")
			downstream = HTSeq.GenomicInterval(p[0].chrom, p[1].pos, p[1].pos + winwidth, ".")
			k = (p[1].pos - p[0].pos) # exon length
		else: # negative strand
			exon = HTSeq.GenomicInterval(p[0].chrom, p[1].pos+1, p[0].pos+1, ".")
			upstream = HTSeq.GenomicInterval(p[0].chrom, p[1].pos - winwidth+1, p[1].pos+1, ".")
			downstream = HTSeq.GenomicInterval(p[0].chrom, p[0].pos+1, p[0].pos + winwidth+1, ".")
			k = (p[0].pos - p[1].pos)
		if k < 100:
			continue
		exoncvg = numpy.fromiter(coverage[exon], dtype='f', count=k) # coverage over exon 
		exoncvg_scaled = numpy.zeros(exonsize, dtype='f') # scaling coverage over each exon to user defined or default exonsize for comparability
		k = (k-1)/exonsize # increment to be used for loop
		i = 0 # loop iterators
		j = 0
		while j < exonsize:
			if round(i) != round(i+k):
				exoncvg_scaled[j] = numpy.mean(exoncvg[int(round(i)):int(round(i+k))])
			else:
				exoncvg_scaled[j] = exoncvg[int(round(i))]
			i += k
			j += 1
		# concatenating coverage over upstream and downstream window with scaled exon coverage
		wincvg = numpy.concatenate((numpy.fromiter(coverage[upstream], dtype='f', count=winwidth), exoncvg_scaled), axis=0)  
		wincvg = numpy.concatenate((wincvg, numpy.fromiter(coverage[downstream], dtype='f', count=winwidth)), axis=0)

		# add coverage array to local coverage profile
		if p[1].strand == "+":
			cvg_list.append(wincvg)
			cvg_list_global.append(wincvg)
		else:
			cvg_list.append(wincvg[::-1])
			cvg_list_global.append(wincvg[::-1])
		gene_names.append(gene)

	cvg_df = pd.concat([pd.Series(x) for x in cvg_list], axis=1)
	cvg_df_norm = cvg_df/cvg_df.max()
	profile = cvg_df_norm.mean(axis=1)
	profile_se = cvg_df_norm.sem(axis=1)
	cvg_df_norm.columns = gene_names
	cvg_df_norm.insert(0, "coordinate", x_coord)
	cvg_df_norm.to_csv(outdir+"/"+outfile+"_exoncoverage.csv", index=False)

	# profile = profile/readnumber # normalize by readnumber
	# interpolate = interp1d(x_coord, profile)
	# profile_smooth = interpolate(x_coord_smooth)

	### Normalization ###

	# MinMax normalization
	# p_max=numpy.amax(profile_smooth)
	# p_min=numpy.amin(profile_smooth)
	# profile_smooth = (profile_smooth-p_min)/(p_max-p_min)

	# Normalization by max
	# p_max=numpy.amax(profile_smooth)
	# profile_smooth = profile_smooth/p_max

	# Normalization by mean
	# profile *= len(profile)/sum(profile)  


	###This part is for plotting profiles from individual bams to separate files ### Comment out if these aren't needed ###

	pyplot.style.use('seaborn-white')
	pyplot.plot(x_coord, profile, color="blue")
	profile_se = cvg_df_norm.sem(axis=1)
	cvg_df_norm.columns = gene_names
	cvg_df_norm.insert(0,"coordinate", x_coord)
	cvg_df_norm.to_csv(outdir+"/"+outfile+"_exoncoverage.csv", index=False)

	# profile = profile/readnumber # normalize by readnumber
	# interpolate = interp1d(x_coord, profile)
	# profile_smooth = interpolate(x_coord_smooth)

    ### Normalization ###

	# MinMax normalization
	# p_max=numpy.amax(profile_smooth)
	# p_min=numpy.amin(profile_smooth)
	# profile_smooth = (profile_smooth-p_min)/(p_max-p_min)

	# Normalization by max
	# p_max=numpy.amax(profile_smooth)
	# profile_smooth = profile_smooth/p_max

	# Normalization by mean
	# profile *= len(profile)/sum(profile)  

	
	###This part is for plotting profiles from individual bams to separate files ### Comment out if these aren't needed ###

	pyplot.style.use('seaborn-white')
	pyplot.plot( x_coord, profile, color="blue")
	pyplot.fill_between(x_coord, profile-profile_se, profile+profile_se,
	alpha=0.5,  facecolor="blue")
	pyplot.title(outfile)
	pyplot.axvline(x=0, ls="-.", lw="2", color="black")
	pyplot.axvline(x=exonsize, ls="-.", lw="2", color="black")
	pyplot.xlabel('position')
	pyplot.ylabel('normalized mean density')
	pyplot.ylim( ymax=max(profile)+0.01)
	pyplot.xlim( xmax=1350, xmin=-350)
	pyplot.xticks([-winwidth, -winwidth/2, 0, exonsize/4, exonsize/2, 3*exonsize/4, exonsize, exonsize + (winwidth/2), exonsize +winwidth], [-winwidth, -round(winwidth/2),'ESS', '25%', '50%', '75%', 'EES', round(winwidth/2), winwidth])
	pyplot.savefig(outdir+"/"+outfile+"_exoncoverage.png", dpi=300)
	pyplot.close()
	
	### This part is to output coverage data as csvs for individual bam files ### Comment out if these aren't needed ###

	print(file + ' completed.')  #to check progress when dealing with a lot of files



### Normalization ###
#MinMax Normalization
# p_max=numpy.amax(profile_global)
# p_min=numpy.amin(profile_global)
# profile_global = (profile_global-p_min)/(p_max-p_min)
# p_max=numpy.amax(profile_global_max)
# p_min=numpy.amin(profile_global_max)
# profile_global_max = (profile_global_max-p_min)/(p_max-p_min)
# p_max=numpy.amax(profile_global_min)
# p_min=numpy.amin(profile_global_min)
# profile_global_min = (profile_global_min-p_min)/(p_max-p_min)
# Normalize by mean
# profile_global *= len(profile_global)/sum(profile_global)  # Normalize average coverage profile by mean coverage per base over profile
# profile_global_min *= len(profile_global_min)/sum(profile_global_min)  # Normalize average coverage profile by mean coverage per base over profile
# profile_global_max *= len(profile_global_max)/sum(profile_global_max)  # Normalize average coverage profile by mean coverage per base over profile

### This part is for plotting global profile ### Comment out if this isn't needed ###
cvg_df_global = pd.concat([pd.Series(x) for x in cvg_list_global], axis=1)
cvg_df_norm_global = cvg_df_global/cvg_df_global.max()
profile_global = cvg_df_norm_global.mean(axis=1)
profile_se = cvg_df_norm_global.sem(axis=1)

pyplot.style.use('seaborn-white')
pyplot.plot( x_coord, profile, color="blue")
pyplot.fill_between(x_coord, profile-profile_se, profile+profile_se,
alpha=0.5,  facecolor="blue")
pyplot.title("global mean")
pyplot.axvline(x=0, ls="-.", lw="2")
pyplot.axvline(x=exonsize, ls="-.", lw="2")
pyplot.xlabel('position')
pyplot.ylabel('normalized mean density')
pyplot.ylim( ymax=max(profile_global)+0.01)
pyplot.xlim( xmax=1350, xmin=-350)
pyplot.xticks([-winwidth, -winwidth/2, 0, exonsize/4, exonsize/2, 3*exonsize/4, exonsize, exonsize + (winwidth/2), exonsize +winwidth], [-winwidth, -round(winwidth/2),'ESS', '25%', '50%', '75%', 'EES', round(winwidth/2), winwidth])
pyplot.savefig(outdir+"/global_exoncoverage.png", dpi=300)
pyplot.close()

### This part is to output global coverage data as csv ### Comment out if this isn't needed ###

# coverage_list=pd.DataFrame({"coordinate" : x_coord, "normalized coverage" : profile_global, })
# coverage_list.to_csv(outdir+"/global_exoncoverage.csv", index=False)

# coverage_list_max=pd.DataFrame({"coordinate" : x_coord, "normalized coverage" : profile_global_max, })
# coverage_list_max.to_csv(outdir+"/global_exoncoverage_maximum.csv", index=False)

# coverage_list_min=pd.DataFrame({"coordinate" : x_coord, "normalized coverage" : profile_global_min, })
# coverage_list_min.to_csv(outdir+"/global_exoncoverage_minimum.csv", index=False)
