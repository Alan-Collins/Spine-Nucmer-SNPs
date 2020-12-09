# Spine-Nucmer-SNPs
 A collection of scripts to process Spine and Nucmer outputs to analyse SNPs from core genome alignments.

## Workflow overview

### TL;DR

1. Spine used to find core genome of each assembly.
2. Nucmer aligns cores and outputs SNPs.
3. snps2fasta.py processes SNPs and output aligned fasta.
4. fasta2diffmat.py does pairwise comparisons of all sequences to create pairwise SNP distances and optionally plot them.
5. get_snps_support_MP.py uses SAM files to identify positions in reads that map to SNP locations to find and summarise the support for SNP calls in reads.

### Workflow upstream of these scripts

These scripts were written to accept the output from [Spine](https://github.com/egonozer/Spine) and [Nucmer](http://mummer.sourceforge.net/) as their input. However, they could be adapted to be used on nucmer produced SNPs of other alignments too.

Spine is used to find the core genome of input assemblies. These core genomes are then aligned against a reference using Nucmer (Spine produces a file called output.backbone.fasta which can be used as this reference and will be in the below example workflow). The resulting alignments are processed using Nucmer programs to produce a list of the SNPs identified in each core genome relative to the reference genome.

### Use of the scripts in this repository

#### snps2fasta.py

The first of the scripts provided here (snps2fasta.py) incorporates all of the SNP lists for each core genome and builds a list of all the positions that are variant in at least one core genome. It records the SNP that is present in every core genome that is variant at that position and adds the base present in the reference genome to any core genomes that are not variant at that site. This SNP information is processed in such a way that all the sequences are aligned with one another and can be output as an aligned fasta of just the variant positions (i.e. variant in at least one genome in the dataset) or of the entire reference core genome with invariant positions filled in using the reference genome.

#### fasta2diffmat.py

The second script provided here (fasta2diffmat.py) performs all pairwise comparisons of the fasta sequences output by snps2fasta.py to create a python dict object of the format {("genome1", "genome2") : #snps} where #snps is an integer count of the positions that differ between genome1 and genome2. This dict is pickled and optionally compressed using gzip. In addition to outputing a pairwise snp distance dict, this script can plot the distribution of snp distances as a histogram of either the entire dataset, or just those pairs of genomes who differ by fewer than a user-provided SNP threshold. These histograms are output if the user provides a filename destination for them. However, other parameters are hardcoded.

#### get_snp_support_MP.py

The third script provided here (get_snp_support_MP.py) queries the reads that were used to create your assemblies to ascertain the support for SNP calls in the final assembly. In short, it works by working through the .snps file output by Nucmer and using the .coords file (also output by Nucmer) to identify the position in the original assembly that is being considered to contain a core genome SNP. Once this position is identified, the .sam file associated with the assembly is used to find all reads that map to that position. The position in each of these reads that corresponds to the SNP is determined and adjusted using the CIGAR string in the sam file. Once all the reads mapping to agiven SNP have been processed in this way, the support information is reported by listing the base calls in each read at that position as well as their quality score. This information is summarized in a single "consensus score" by dividing the sum of the quality scores that agree with the base call by the sum of all quality scores at that position. For example, for a position called as A with 5 reads mapping to it, with 3 A and 2 T calls where all have a quality score of 10, consensus score = (3 * 10)/((3 * 10) + (2 * 10)) = 0.6

The consensus score is reported in an output file along with other information to help with interpretation of the consensus score (e.g. coverage and the base calls and quality scores that are used to calculate the consensus score) as well as position information in the core and original assembly, and a list of the reads that were identified as mapping to the SNP location. A header line is provided to indicate what each column in the file represents.


## Example workflow

### Intro and setup

In this example I will describe a simple use of Spine and Nucmer as well as the included scripts on dummy data. This assumes that you have assemblies and associated .sam files before  beginning. In all cases that commands are executed, the path will be indicated relative to the following dummy directory structure (for clarity directories are named in upper-case letters while files are lower-case):

```
~/
|-- ASSEMBLIES/
|   |-- genome1.fna
|   |-- genome2.fna
|   `-- ...
|-- SAMS/
|   |-- genome1.sam
|   |-- genome2.sam
|   `-- ...
|-- SPINE/
|-- NUCMER/
|-- SNPS/
`-- SNP_SUPPORT/
```

### Spine

Spine takes as input a text file with 3 columns: path/to/file	ID	fasta
The config file to run spine can be created with the following command (assuming only assembly files are present in your ASSEMBLIES/ dir):

```bash
~/ASSEMBLIES$ ls | awk 'BEGIN { FS="\t"; OFS="\t" } { print "../ASSEMBLIES/"$1, $1, "fasta" }' > ../SPINE/config.txt
```

Spine can then be run with as many cores as you like using something like the following command:

```bash
~/SPINE$ spine.pl -f config.txt -t 40
```