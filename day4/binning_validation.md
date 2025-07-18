# Binning & Validation Tutorial

For this tutorial we are going to be working in the `checkm2` environment in conda. That environment has both `metabat2` and `checkm2` pre-installed:

```
conda activate checkm2
```

We also need to install MetaBAT2, although it's possible that it has already been installed, so we will check first. To see if the command `metabat2` is installed in our environment, we can use the `which` command:

```
which metabat2
# Should print --> /opt/miniconda3/envs/checkm2/bin/metabat2
```

If that `which` command doesn't print a path to a `metabat2` file, then we will need to install it. If it does, you can skip the following command:

```
conda install -c bioconda metabat2
```

If that succeeded, you can make sure by testing with the `which` command again:
```
which metabat2
```

If that printed a path to a `metabat2` file (e.g. `/opt/miniconda3/envs/checkm2/bin/metabat2`), then skip the next step and go right to "Collecting the Data from Last Time"

**NOTE**: It is possible that that instalation may run into an error (related to libboost, but that's for another day). If it does, you'll have to create a second environment for metabat2 and switch between them as needed. As of the 7/25 morning session the command above appears to be working, but in case it doesn't the following will create a separate conda environment for MetaBAT2 that you can switch in and out of below. (See the lines with `# NOTE:` at the start of them.)

```
conda deactivate
conda create -n metabat2 python=3.9 metabat2
conda activate metabat2
```

## Collecting the Data from Last Time

Let's start by collecting the data that we worked with last time, during the [QC & Assembly](qc_assembly.md) tutorial. If you recall, during that tutorial we went through the following steps:

1. Downloading & Sanity-checking the Raw Reads
2. Running FastQC on the reads
3. Trimming adapters off the reads using BBduk
4. Assembled the reads using MEGAHIT
	* This step created an output file called `final.contigs.fa` which is hopefully still in your folders. 
5. Binned the contigs by:
	* Mapping reads back to the contigs using a combination of `bowtie2` and `samtools`. This step(s) created an output file called `tara.bam` (Unless you named it something else). That file provided the coverage information needed for binning.
	* Running MetaBAT using the script `runMetaBat.sh` on the contigs together with the coverage (`.bam`) file, which created an output folder containing (hopefully) 10 bins, each in a file with the extension `.fa`. If you ran all the commands in the tutorial exactly, those should be in a subfolder called `metabat`.

<details>
<summary><b>If you have not run the QC & Assembly tutorial... </b></summary>

below is a set of commands you should be able to run all at once to do it:

```
# 1) Make assembly conda environment:
conda create -y -n assembly -c conda-forge -c bioconda -c defaults fastqc bbmap megahit metabat2 bowtie2 samtools
conda activate assembly

# 2) Get the Data:
mkdir -p ~/data && cd ~/data
cp /opt/shared/assembly-data/tara_reads_R* .
chmod -w tara_reads_R*

# 3) Adapter trimming
mkdir -p ~/results && cd ~/results
wget https://sourceforge.net/projects/bbmap/files/BBMap_39.08.tar.gz
tar -xzf BBMap_39.08.tar.gz
bbmappath=$(pwd)/bbmap
$bbmappath/bbduk.sh in1=~/data/tara_reads_R1.fastq.gz in2=~/data/tara_reads_R2.fastq.gz out1=~/data/tara_reads_R1_trimmed.fastq.gz out2=~/data/tara_reads_R2_trimmed.fastq.gz ref=$bbmappath/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

# 4) Assembly, Contig-Coverage & Binning:
megahit -1 ~/data/tara_reads_R1_trimmed.fastq.gz -2 ~/data/tara_reads_R2_trimmed.fastq.gz -o tara_assembly
ln -s tara_assembly/final.contigs.fa .
# ...coverage:
bowtie2-build final.contigs.fa final.contigs
bowtie2 -x final.contigs -1 ~/data/tara_reads_R1_trimmed.fastq.gz -2 ~/data/tara_reads_R2_trimmed.fastq.gz | samtools view -bS -o tara_to_sort.bam
samtools sort tara_to_sort.bam -o tara.bam
samtools index tara.bam
# ...binning:
runMetaBat.sh -m 1500 final.contigs.fa tara.bam
mv final.contigs.fa.metabat-bins1500* metabat

# *** DONE ***
```
</details>

Some enterprising students went beyond this and followed the link to the original source material, which involved running CheckM (version 1) on the bins, although this ran into some errors because CheckM-1 is not supported currently.

#### Making Variables to Point to the Outputs

We're going to use a Linux shell trick here to make the commands much neater. Specifically, we're going to set some environment variables that contain *absolute* paths to the outputs from the last tutorial. By "absolute" paths, we mean paths that start with either `~/` or with a `/`, so they represent the same folder no matter what directory you might be in when you call them. Here are the commands to set these variables using the absolute paths based on the *default* locations from the last tutorial:

```
contigs=~/results/final.contigs.fa
bam_unsorted=~/results/tara_to_sort.bam
bam_sorted=~/results/tara.bam
metabat_output=~/results/metabat/
```

#### If you dont have them at the locations Titus added them to a shared directory
```
mkdir -p ~/results
cd ~/results
cp /opt/shared/assembly-data/metabat-data/final.contigs.fa .
cp /opt/shared/assembly-data/metabat-data/tara.bam .
cp -r /opt/shared/assembly-data/metabat-data/metabat metabat
```
If you have saved these in some other folder than `~/results/...`, just change everything after the `=` there to be the absolute path to the location where you saved them. (We don't actually need the path to the unsorted `.bam` file so if you want to omit that one, go ahead.) Note that the last of these is the path to the *folder* containing all of the bins output by MetaBAT, which we renamed in the last tutorial to simply `metabat`. But that folder should have 10 files in it that all end in `.#.fa` where `#` goes from 1 to 10. 

These are now variables that are stored in the shell environment, and you can now refer to the stuff on the right using the names on the left (but always preceded by a `$` character. Let's see how that works by testing them with the `echo` command, which just prints the first argument right back to the console. Run the following:

```
echo $contigs
echo $bam_unsorted
echo $bam_sorted
echo $metabat_output
```

After each one of those you should see the path you gave earlier printed right after the command:

```
(checkm2) stamps@149.165.174.116:~/binning$ echo $contigs
echo $bam_unsorted
echo $bam_sorted
echo $metabat_output
/home/stamps/results/final.contigs.fa
/home/stamps/results/tara_to_sort.bam
/home/stamps/results/tara.bam
/home/stamps/results/metabat/
```

Note a couple things here. First, the console has replaced the `~` character with `/home/stamps`, which is normal and because `~` is itself just a special environment variable representing that home folder. Let's see if we can list the MetaBAT bin files using this environment variable. Try running the following command and if it doesn't list the files, make sure you've put the correct absolute path in the variable assignment above:

```
ls $metabat_output -l
```

In the remainder of this tutorial, we will refer to these paths using the environment variables we've set here, and we may set and use a few more as well :-)

## Running CheckM2 on the Previous Bins

First we are going to `cd` into the results folder (or whatever your equivalent was), then we are going to create a folder to hold the CheckM2 output:

```
cd ~/results
mkdir metabat_checkm2_output
```

Note that the path in that `mkdir` command was a *relative* path rather than an absolute one. So the absolute path to the folder we just created would be `~/results/metabat_checkm2_output`, and we could just as easily have done `mkdir ~/results/metabat_checkm2_output` and gotten the same result. In general absolute paths are nice because you don't have to worry about where you are when you execute the commands, but they are a bit more work. In this case it is for illustration.

OK, Now we can run CheckM2. (Note that for the steps below where we run MetaBAT2, we will name equivalent folders with `metabat2` instead of just `metabat`.) This command will run CheckM2, which should take a minute or less (roughly):

```
# NOTE: if you had to make a separate environment for MetaBAT2, run the following command (without the # character):
#conda activate checkm2	
checkm2 predict --threads 16 --input $metabat_output -x fa --output-directory metabat_checkm2_output/
```

After that command is finished running, there should be a couple of files and also a couple of subfolders in the folder `metabat_checkm2_output`. Look at the file called `quality_report.tsv`. You can print it to the console using `cat ~/results/metabat_checkm2_output/quality_report.tsv`, though you might want to pipe that to the `head` command so it doesn't print the entire thing. (Do you remember how to pipe commands together?) Or you can open it in Rstudio or download it to your laptop and open it in Excel, whatever floats your boat... 

When you open the report, what do you think about these bins? Are they good? Bad? Roughly equal size or highly variable?

## Running MetaBAT2 to Create New Bins

MetaBAT2 requires that the coverage information be converted to a very specific format first. Fortunately, it comes with a script to do that conversion automatically from a `.bam` file. In the following command, we're first going to set a variable to the path (again, *absolute* path) where want the output file to be located. (Note that this is a path to a file that doesn't yet exist, but that's ok. Feel free to point this to whatever folder and filename you would like.) Then in the second command we are going to refer to that output as well as the variable `$bam_sorted` from earlier:

```
# NOTE: if you had to make a separate environment for MetaBAT2, run the following command (without the # character):
#conda activate metabat2
depths_output=~/results/tara_depths.txt
jgi_summarize_bam_contig_depths --outputDepth $depths_output $bam_sorted
```
If the above code gave you some issues try running a similar analysis using this:
```
# Make coverm conda env
conda create -y -n coverm -c bioconda coverm

# Activate it
conda activate coverm

# Run coverm on the files
metabat_depths=~/results/tara__metabat_depths.txt
coverm genome -t 15 -x .fa --bam-files $bam_sorted --genome-fasta-directory $metabat_output/ -o $metabat_depths
```

If the jgi_summarize_bam_contig_depths step didn't work then neither will the metabat2 without more debugging than I'm able to right now so lets try this other bin tool instead!
```
# Make semibin2 conda env
conda create -y -n semibin2 -c conda-forge -c bioconda semibin

# Activate it
conda activate semibin2

# Run semibin2 on the files
SemiBin2 single_easy_bin -t 15 -i $contigs -b $bam_sorted --environment ocean -o ~/results/semibin2_output
```
#### Challenge: Try to get these new bins analyzed with checkm2 and coverm by adapting the above scripts
This will allow us to compare the binning tools

Now we can run MetaBAT2. First we'll make a specific directory for the output, though we're going to use a variable containing an absolute path this time to do it:

```
metabat2_output=~/results/metabat2
mkdir $metabat2_output
metabat2 -i $contigs -o $metabat2_output -a $depths_output
```

This is a good example of why making variables is helpful: if we want to re-run this, for example on another sample, we only have to change the folder location in one place and it propogates through the code every time the variable is used later. Let's do the same thing to make an output folder for CheckM2 when we run it on these bins, then let's run it:

```
metabat2_checkm2=~/results/metabat2_checkm2_output
mkdir $metabat2_checkm2
# NOTE: if you had to make a separate environment for MetaBAT2, run the following command (without the # character):
#conda activate checkm2	
checkm2 predict --threads 16 --input $metabat2_output -x fa --output-directory $metabat2_checkm2
```

#### Compare the bins from MetaBAT with MetaBAT2

Remember that MetaBAT2 is ostensibly a major overhaul of the algorithm for doing contig binning. This Tara dataset though is relatively simple, so it's an interesting question to compare results of the old and the new methods on these data 

* Is the number of bins the same? (Recall that the data was supposedly sampled from 20 organisms)
* Can you tell which bins in the old version are equivalent to which in the new version? Which ones are bigger?
* Which version seems to score better on quality metrics like completeness or N50?
	
## (Optional) Run CheckM2 on Contigs from a Human Gut Sample

I have run MetaBAT2 on one of the assemblies that we talked about in Monday's lectures (and which Titus also talked about by way of comparing MEGAHIT and metaSPAdes assemblies). The folder containing the output of MetaBAT2 has been zipped and uploaded to the server, so first make a copy of it into your local `~/results` folder in a sensibly named subfolder. Note that the first line in the commands below sets another environment variable to the location of the zip file. (That line has been updated to have to correct location to the zip file, and should work.)

```
gut_metabat_zip_location=/opt/shared/assembly-data/PRJNA1049470_SRR27117388_megahit_metabat_out.zip
human_gut=~/results/SRR27117388
mkdir $human_gut
cd $human_gut
cp $gut_metabat_zip_location ./
unzip PRJNA1049470_SRR27117388_megahit_metabat_out.zip
```

Now the bins should be in a subfolder called `metabat_out`, so let's make a variable containing an *absolute* path to the bins folder for this sample. Note that we can actually use previous variables to define new ones, so we'll do that here for illustration:

```
human_gut_metabat2=$human_gut/metabat_out
echo $human_gut_metabat2
# Should be: /home/stamps/results/SRR27117388/metabat_out
```

Now let's make a subfolder for the CheckM2 output and run it on these bins:

```
# NOTE: if you had to make a separate environment for MetaBAT2, run the following command (without the # character):
#conda activate checkm2	
human_gut_checkm2=$human_gut/checkm2_out
checkm2 predict --threads 16 --input $human_gut_metabat2 -x fa --output-directory $human_gut_checkm2
```

## (Double Optional) Run Sourmash on Human Gut Bins

This is an exercise for you to figure out how to run Sourmash on these bins. For each bin, you'll want to search the kmer content in that fasta file against GTDB; see [the tutorial section here](https://hackmd.io/9ORFRJGaTOiOdEAY-Aih2A?view#Generating-taxonomic-classifications-for-metagenomes-with-sourmash)). To assign a taxon to that bin, you can use `sourmash tax genome` instead of `tax metagenome`. Here is an incomplete bash script that you could possibly use to automate this:

```
for fi in $human_gut_metabat2/*.fa
do
	printf "$fi: "
	<sourmash command> | sed -n '#p'
done
```

This is a `for` loop in `bash`, where we are iterating over all of the files in our bins folder (i.e. `$human_gut_metabat2/*.fa`, which will yield *absolute* paths). At each step of the loop, the path is stored in a temporary variable called `$fi`, which you can feed into the appropriate sourmash command. Since that command outputs to the console (IIRC), you can pipe the output to `sed` to pick only a specified line number. So for example if the sourmash command outputs the top genome containing the bin on the 3rd line, you would change `#` to `3`. See if you can figure out how to do that :-).
