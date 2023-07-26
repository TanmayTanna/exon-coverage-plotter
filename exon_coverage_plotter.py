### Tanmay Tanna ### 14-06-2018 ###

from __future__ import division
import HTSeq
import numpy
from matplotlib import pyplot
import argparse
from scipy.interpolate import interp1d
import csv
import pandas as pd

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

#initialize values

parser._action_groups.append(optional) 
args = parser.parse_args()

bamlist=args.inbamlist
gtffile = HTSeq.GFF_Reader(str(args.genome))
outdir=str(args.outPath)
winwidth = int(args.winwidth)
exonsize=int(args.exonsize)
fragmentsize = 200 # This sets size of alignment to 200 as 200 bp is more representative of standard fragment length than the observed readlength 



x_coord=numpy.arange( -winwidth, exonsize+winwidth) # user defined or default size of coverage plot
x_coord_smooth = numpy.linspace( -winwidth, exonsize+winwidth - 1,300) # x coordinates used for smoothing
# adding exon start and end to smoothing coordinates
if 0 not in x_coord_smooth:
	x_coord_smooth=numpy.append(x_coord_smooth, [0], axis=0)
if exonsize not in x_coord_smooth:
	x_coord_smooth=numpy.append(x_coord_smooth, [exonsize], axis=0)
x_coord_smooth = numpy.sort(x_coord_smooth)

profile_global = numpy.zeros( numpy.size(x_coord_smooth) , dtype='f' )  # average global (for all bams) coverage profile
profile_global_min = numpy.zeros( numpy.size(x_coord_smooth) , dtype='f' ) # minimum coverage values at each coordinate
profile_global_max = numpy.zeros( numpy.size(x_coord_smooth) , dtype='f' ) # maximum coverage values at each coordinate
profile_global_min[0:len(profile_global_min)]=100000 # initializing minimum to high value to allow initializing with first value in loop

filenumber=0 # file counter

### looping through for each bam file ###
for file in bamlist:
	filenumber+=1
	bamfile = HTSeq.BAM_Reader(str(file))
	coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" ) #initializing coverage array
	outfile= file.split("/")[-1]

	readnumber=0 # number of reads

	for almnt in bamfile:
	    readnumber+=1
	    if almnt.aligned:
	    	if almnt.iv.end > fragmentsize:
	    		almnt.iv.length = fragmentsize
	    	else:
	    		almnt.iv.length = almnt.iv.end
	    	coverage[ almnt.iv ] += 1  # add 1 to coverage for each aligned "fragment"

	exonpos = list()
	for feature in gtffile:
	    if feature.type == "gene":
	       exonpos.append( [feature.iv.start_d_as_pos, feature.iv.end_d_as_pos] ) # read gene/exon coordinates

	profile = numpy.zeros( 2*winwidth + exonsize , dtype='f' ) # average local (for current bam) coverage profile

	for p in exonpos:
	    if p[0].pos<p[1].pos: # positive strand
	        exon = HTSeq.GenomicInterval( p[0].chrom, p[0].pos , p[1].pos , "." )
	        upstream = HTSeq.GenomicInterval( p[0].chrom, p[0].pos - winwidth, p[0].pos, "." )
	        downstream = HTSeq.GenomicInterval( p[0].chrom, p[1].pos, p[1].pos + winwidth, "." )
	        k = (p[1].pos - p[0].pos) # exon length
	    else: # negative strand
	        exon = HTSeq.GenomicInterval( p[0].chrom, p[1].pos+1 , p[0].pos+1 , "." )
	        upstream = HTSeq.GenomicInterval( p[0].chrom, p[1].pos - winwidth+1, p[1].pos+1, "." )
	        downstream = HTSeq.GenomicInterval( p[0].chrom, p[0].pos+1, p[0].pos + winwidth +1,  "." )
	        k = (p[0].pos - p[1].pos)
	    if k<100:
	    	continue
	    exoncvg = numpy.fromiter( coverage[exon], dtype='f', count=k) # coverage over exon 
	    exoncvg_scaled = numpy.zeros( exonsize, dtype='f' ) # scaling coverage over each exon to user defined or default exonsize for comparability
	    k=(k-1)/exonsize # increment to be used for loop
	    i=0 # loop iterators
	    j=0
	    while j<exonsize:
	        if round(i)!=round(i+k):
		        exoncvg_scaled[j]=numpy.mean(exoncvg[int(round(i)):int(round(i+k))])
	        else:
		        exoncvg_scaled[j]=exoncvg[int(round(i))]
	        i+=k
	        j+=1
	    # concatenating coverage over upstream and downstream window with scaled exon coverage
	    wincvg = numpy.concatenate((numpy.fromiter(coverage[upstream], dtype='f', count=winwidth), exoncvg_scaled), axis = 0)  
	    wincvg = numpy.concatenate((wincvg, numpy.fromiter(coverage[downstream], dtype='f', count=winwidth)), axis = 0)
	  
	    # add coverage array to local coverage profile
	    if p[1].strand == "+":
	        profile += wincvg
	    else:
	        profile += wincvg[::-1] # reverse the coverage array if transcription direction is opposite

	profile = profile/readnumber # normalize by readnumber
	interpolate = interp1d(x_coord, profile)
	profile_smooth = interpolate(x_coord_smooth)

	profile_global = profile_global + profile_smooth # Add profile to global profile

	# Check if any point in profile has the maximum or minimum value compared to other profiles
	for i in numpy.arange(0,len(profile_smooth)):
		if profile_smooth[i] < profile_global_min[i]:
			profile_global_min[i]=profile_smooth[i]
		if profile_smooth[i] > profile_global_max[i]:
			profile_global_max[i]=profile_smooth[i]
    
    ### Normalization ###

	# MinMax normalization
	p_max=numpy.amax(profile_smooth)
	p_min=numpy.amin(profile_smooth)
	profile_smooth = (profile_smooth-p_min)/(p_max-p_min)

	# Normalization by max
	# p_max=numpy.amax(profile_smooth)
	# profile_smooth = profile_smooth/p_max

	# Normalization by mean
	# profile *= len(profile)/sum(profile)  

	
	###This part is for plotting profiles from individual bams to separate files ### Comment out if these aren't needed ###

	pyplot.style.use('ggplot')
	pyplot.plot( x_coord_smooth, profile_smooth, color="blue")
	pyplot.title(outfile)
	pyplot.axvline(x=0, ls="-.", lw="2")
	pyplot.axvline(x=exonsize, ls="-.", lw="2")
	pyplot.xlabel('ESS = exon start site     Position     EES = exon end site')
	pyplot.ylabel('Normalized Density')
	pyplot.ylim( ymax=max(profile)+1)
	pyplot.xticks([-winwidth, -winwidth/2, 0, exonsize/4, exonsize/2, 3*exonsize/4, exonsize, exonsize + (winwidth/2), exonsize +winwidth], [-winwidth, -round(winwidth/2),'ESS', '25%', '50%', '75%', 'EES', round(winwidth/2), winwidth])
	pyplot.savefig(outdir+"/"+outfile+"_exoncoverage.png")
	pyplot.close()
	
	### This part is to output coverage data as csvs for individual bam files ### Comment out if these aren't needed ###

	coverage_list=pd.DataFrame({"coordinate" : x_coord_smooth, "normalized coverage" : profile_smooth})
	coverage_list.to_csv(outdir+"/"+outfile+"_exoncoverage.csv", index=False)
	print('1 (more) instance complete')  #to check progress when dealing with a lot of files


profile_global=profile_global/filenumber # divide sum of all profiles by number of profiles

### Normalization ###
#MinMax Normalization
p_max=numpy.amax(profile_global)
p_min=numpy.amin(profile_global)
profile_global = (profile_global-p_min)/(p_max-p_min)

p_max=numpy.amax(profile_global_max)
p_min=numpy.amin(profile_global_max)
profile_global_max = (profile_global_max-p_min)/(p_max-p_min)
p_max=numpy.amax(profile_global_min)
p_min=numpy.amin(profile_global_min)
profile_global_min = (profile_global_min-p_min)/(p_max-p_min)

# Normalize by mean
# profile_global *= len(profile_global)/sum(profile_global)  # Normalize average coverage profile by mean coverage per base over profile
# profile_global_min *= len(profile_global_min)/sum(profile_global_min)  # Normalize average coverage profile by mean coverage per base over profile
# profile_global_max *= len(profile_global_max)/sum(profile_global_max)  # Normalize average coverage profile by mean coverage per base over profile

### This part is for plotting global profile ### Comment out if this isn't needed ###

pyplot.style.use('ggplot')
pyplot.plot( x_coord_smooth, profile_global, color="blue")
pyplot.title("global average")
pyplot.axvline(x=0, ls="-.", lw="2")
pyplot.axvline(x=exonsize, ls="-.", lw="2")
pyplot.xlabel('ESS = exon start site     Position     EES = exon end site')
pyplot.ylabel('Normalized Density')
pyplot.ylim( ymax=max(profile_global)+1)
pyplot.xticks([-winwidth, -winwidth/2, 0, exonsize/4, exonsize/2, 3*exonsize/4, exonsize, exonsize + (winwidth/2), exonsize +winwidth], [-winwidth, -round(winwidth/2),'ESS', '25%', '50%', '75%', 'EES', round(winwidth/2), winwidth])
pyplot.savefig(outdir+"/global_exoncoverage.png")
pyplot.close()

### This part is to output global coverage data as csv ### Comment out if this isn't needed ###

coverage_list=pd.DataFrame({"coordinate" : x_coord_smooth, "normalized coverage" : profile_global, })
coverage_list.to_csv(outdir+"/global_exoncoverage.csv", index=False)

coverage_list_max=pd.DataFrame({"coordinate" : x_coord_smooth, "normalized coverage" : profile_global_max, })
coverage_list_max.to_csv(outdir+"/global_exoncoverage_maximum.csv", index=False)

coverage_list_min=pd.DataFrame({"coordinate" : x_coord_smooth, "normalized coverage" : profile_global_min, })
coverage_list_min.to_csv(outdir+"/global_exoncoverage_minimum.csv", index=False)
