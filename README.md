### Sequence Sampler

This python program takes single/multiple bam file(s) as input along with a reference feature list file (gff/gtf), and generates exon coverage plots (individual and mean) using these bam files. The [HTSeq package](https://htseq.readthedocs.io/en/release_0.11.1/) is required. 

The program standardizes each exon in the feature list to a user-defined size, and includes a window (upstream and downstream) to the exon in the coverage plots. It generates plots for each individual input bam file, a 'global' coverage plot that represents the mean coverage in all bam files. In addition to this, it also generates CSV files with coverage data as output. 

Usage example:\
`python exon_coverage_plotter.py --inbamlist a.bam b.bam c.bam --genome genome.gff3 --outPath output_path --winwidth 300 --exonsize 1000` 

The inputs are as follows:

```
Required inputs:
--inbamlist	input bam file(s) separated by spaces
--genome	input reference feature list file (gff/gtf)

Optional inputs:
--outPath 	output path
--winwidth	window (upstream and downstream) to the exon to be included in the coverage plots
--exonsize	size (in bp) to standardize all exon to
```

### Reference:
This script is adapted to a significant extent from https://htseq.readthedocs.io/en/release_0.11.1/tss.html


