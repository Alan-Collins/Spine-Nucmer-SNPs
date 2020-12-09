# Spine-Nucmer-SNPs
 A collection of scripts to process Spine and Nucmer outputs to analyse SNPs from core genome alignments.

## Workflow overview

### TL;DR

1. Spine used to find core genome of each assembly.
2. Nucmer aligns cores and outputs SNPs.
3. snps2fasta.py processes SNPs and output aligned fasta.
4. fasta2diffmat.py does pairwise comparisons of all sequences to create pairwise SNP distances and optionally plot them.
5. get_snps_support_MP.py uses SAM files to identify positions in reads that map to SNP locations to find and summarise the support for SNP calls in reads.

N.B. All of the scripts in this repository have a hopefully clear description of their functionality and arguments they take if you call them with -h (e.g. snps2fasta.py -h)

### Workflow upstream of these scripts

These scripts were written to accept the output from [Spine](https://github.com/egonozer/Spine) and [Nucmer](http://mummer.sourceforge.net/) as their input. However, they could be adapted to be used on nucmer produced SNPs of other alignments too.

Spine is used to find the core genome of input assemblies. These core genomes are then aligned against a reference using Nucmer (Spine produces a file called output.backbone.fasta which can be used as this reference and will be in the below example workflow). The resulting alignments are processed using Nucmer programs to produce a list of the SNPs identified in each core genome relative to the reference genome.

### Use of the scripts in this repository

#### snps2fasta.py

The first of the scripts provided here (snps2fasta.py) incorporates all of the SNP lists for each core genome and builds a list of all the positions that are variant in at least one core genome. It records the SNP that is present in every core genome that is variant at that position and adds the base present in the reference genome to any core genomes that are not variant at that site. This SNP information is processed in such a way that all the sequences are aligned with one another and can be output as an aligned fasta of just the variant positions (i.e. variant in at least one genome in the dataset) or of the entire reference core genome with invariant positions filled in using the reference genome.

#### fasta2diffmat.py

The second script provided here (fasta2diffmat.py) performs all pairwise comparisons of the fasta sequences output by snps2fasta.py to create a python dict object of the format {("genome1", "genome2") : #snps} where #snps is an integer count of the positions that differ between genome1 and genome2. This dict is pickled and optionally compressed using gzip. In addition to outputing a pairwise snp distance dict, this script can plot the distribution of snp distances as a histogram of either the entire dataset, or just those pairs of genomes who differ by fewer than a user-provided SNP threshold. These histograms are output if the user provides a filename destination for them. However, other parameters are hardcoded.

#### get_snp_support_MP.py

The third script provided here (get_snp_support_MP.py) queries the reads that were used to create your assemblies to ascertain the support for SNP calls in the final assembly. In short, it works by working through the .snps file output by Nucmer and using the .core_coords file (output by Spine) to identify the position in the original assembly that is being considered to contain a core genome SNP. Once this position is identified, the .sam file associated with the assembly is used to find all reads that map to that position. The position in each of these reads that corresponds to the SNP is determined and adjusted using the CIGAR string in the sam file. Once all the reads mapping to agiven SNP have been processed in this way, the support information is reported by listing the base calls in each read at that position as well as their quality score. This information is summarized in a single "consensus score" by dividing the sum of the quality scores that agree with the base call by the sum of all quality scores at that position. For example, for a position called as A with 5 reads mapping to it, with 3 A and 2 T calls where all have a quality score of 10, consensus score = (3 * 10)/((3 * 10) + (2 * 10)) = 0.6

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

Spine can then be run with as many cores as you like (-t option) using something like the following command:

```bash
~/SPINE$ spine.pl -f config.txt -t 40
```
The SPINE/ directory should look something like this once it has finished running:

```
~/SPINE/
 |-- config.txt
 |-- output.backbone_coords.txt
 |-- output.backbone.fasta
 |-- output.coords.txt
 |-- output.delta
 |-- output.genome1.fna.accessory_coords.txt
 |-- output.genome1.fna.accessory.fasta
 |-- output.genome1.fna.core_coords.txt
 |-- output.genome1.fna.core.fasta
 |-- output.genome2.fna.accessory_coords.txt
 |-- output.genome2.fna.accessory.fasta
 |-- output.genome2.fna.core_coords.txt
 |-- output.genome2.fna.core.fasta
 `-- ...
```

Spine will also output summary statistics of the overall core to your terminal which can also be found in the output.statistics.txt file.

### Nucmer

We will now use Nucmer to align each of the .core.fasta files with the output.backbone.fasta and use some other Nucmer programs so that we have SNP information for each core genome relative to the same reference.

All of the Nucmer steps can be run with the following line of code:

```bash
~/NUCMER$ ls ../SPINE/*.core.fasta | while read i; do acc=${i%.core*}; acc=${acc#../SPINE/output.}; nucmer --prefix=${acc}_core ../SPINE/output.backbone.fasta $i; delta-filter -r -q ${acc}_core.delta > ${acc}_core.filter; show-snps -Clr ${acc}_core.filter > ${acc}_core.snps; done
```

Once that has finished your NUCMER/ directory will look something like this:

```
~/NUCMER/
 |-- genome1.fna_core.delta
 |-- genome1.fna_core.filter
 |-- genome1.fna_core.snps
 |-- genome2.fna_core.delta
 |-- genome2.fna_core.filter
 |-- genome2.fna_core.snps
 `-- ...
```

### snps2fasta.py

Now that we have all of our SNPs listed in something resembling variant call format, we need to process that into a fasta sequence for downstream comparisons and tree making etc. This can be done with something like the following command (N.B. if you want the entire core genome rather than just the variant positions, use -whole):

```bash
~/SNPS$ python3 /PATH/TO/snps2fasta.py -r ../SPINE/output.backbone.fasta -f variant_core.fasta -p '(.*)_core\.snps' ../NUCMER/*.snps
```

### fasta2diffmat.py

The aligned fasta sequence ouput by snps2fasta.py can then be used to count the SNPs between each assembly in your dataset. fasta2diffmat.py can be run with something like the following command:

```bash
~/SNPS$ python3 /PATH/TO/fasta2diffmat.py -f variant_core.fasta -d diff_dict.pkl -z -t 20 -g SNP_dist_hist.png -c under_2500_SNP_dist_hist.png -ct 2500
```

The above command will output a compressed SNP dict as well as two histograms: one showing the SNP distance between all pairs of core genomes, and one showing just those with fewer than 2500 SNPs between them.

### get_snps_support_MP.py

Finally, given that all of the above information is ultimately based on the SNP calls made during the assembly of your genomes, it is worth checking how confident you can be in those SNPs. get_snps_support_MP.py checks the base call and confidence scores in the reads that contributed to SNP calls (as described in a bit more detail above). Before you run it, you need to make a config file that describes where all the corresponding files are that this script needs. That config needs to have the following columns:

```
PATH/TO/snps_files	PATH/TO/coords_files	PATH/TO/assemblies	PATH/TO/sam_files	ID
```

You can make the config file in our example with the following command:

```bash
~/SNP_SUPPORT$ paste <(ls ../NUCMER/*.snps) <(ls ../SPINE/*.core_coords.txt) <(ls ../ASSEMBLIES/*.fna) <(ls ../SAMS/*.sam) <(ls ../SAMS/) | awk '{gsub("../SAMS/","",$5)}1 {gsub(".sam","",$5)}1' | sed 's/ /\t/g' > config.txt
```

get_snps_support_MP.py can then be run with something like the following command:

```bash
~/SNP_SUPPORT$ get_snps_support.py -config config.txt -outdir ./ -threads 20 -chunksize 10 -pickle
```
 N.B. This will run quickest if your chunksize ~ number of genomes / number of threads.

 This script outputs one file per genome where each line is the support information for a single SNP. You can use the consensus score column (column 9) to quickly identify positions with low support.