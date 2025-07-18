# Strain-level analyses with ParSNP (STAMPS 2024)

This tutorial is to go over how to use ParSNP for strain-level analyses of genomes. The first dataset is a MERS coronavirus outbreak dataset involving 49 isolates. The second dataset is a selected set of 31 Streptococcus pneumoniae genomes. For reference, both of these datasets should run on modestly equipped laptops in a few minutes or less.	

Before starting, you will need to activate the `parsnp` conda environment and double check that `parsnp` is installed:

```
conda activate parsnp
which parsnp
```

## Example 1: 49 MERS Coronavirus genomes

#### Download genomes: 

First we will pull these genomes from the parsnp Git repository. These genomes are included with the documentation for ParSNP.

```
mkdir parsnp_demo1
cd parsnp_demo1
wget https://github.com/marbl/harvest/raw/master/docs/content/parsnp/mers49.tar.gz
tar -xvf mers49.tar.gz
```

#### Run parsnp with default parameters 

```
parsnp -r ./mers49/England1.fna -d ./mers49 -c
```
 
* Command-line output:

![merscmd](https://github.com/marbl/harvest/raw/master/docs/content/parsnp/run_mers.cmd1.png?raw=true)

#### Visualize with Gingr 

##### MacOS Install
+ Get your [Gingr download file](https://harvest.readthedocs.io/en/latest/content/gingr.html)
+ Go to your downloads folder and unzip the gingr-OSX64-v1.3.app.zip file
+ Drag the app file to your applications folder
+ In your applications folder ctrl+click then click open (if this works the bioinfo gods have blessed you today)
	+ If not please let one of us know

![mers1](https://github.com/marbl/harvest/raw/master/docs/content/parsnp/run_mers.gingr1.png?raw=true)

#### Configure parameters
+ 95% of the reference is covered by the alignment. This is <100% mainly due to a 1kbp unaligned region from 26kbp to 27kbp.
+ To force alignment across large collinear regions, use the `-C` maximum distance between two collinear MUMs::

```
parsnp -r ./mers49/England1.fna -d ./mers49 -C 2000
```
	
####  Visualize again with Gingr
+ By adjusting the `-C` parameter, this region is no longer unaligned, boosting the reference coverage to 97%.

![mers2](https://github.com/marbl/harvest/raw/master/docs/content/parsnp/run_mers.gingr2.png?raw=true)

####  Zoom in with Gingr for nucleotide view of region:
+ On closer inspection, a large stretch of N's in Jeddah isolate C7569 was the culprit
 
![mers3](https://github.com/marbl/harvest/raw/master/docs/content/parsnp/run_mers.gingr3.png?raw=true)
 
#### Inspect Output:
+ Multiple alignment: `XMFA parsnp.xmfa` 
+ SNPs: `VCF parsnp.vcf` -- this can be created by adding the --vcf command when running parsnp (not necessary today)
+ Phylogeny: `Newick parsnp.tree`

## Example 2: 31 Streptococcus pneumoniae genomes

#### 1) Download genomes:
```
cd $HOME
mkdir parsnp_demo2
cd parsnp_demo2
wget https://github.com/marbl/harvest/raw/master/docs/content/parsnp/strep31.tar.gz
tar -xvf strep31.tar.gz
```

#### 2) Run parsnp:
  
```
parsnp -r ./strep31/NC_011900.fna -d ./strep31 -p 15
```

#### 3) Force inclusion of all genomes (-c):
  
```
parsnp -r ./strep31/NC_011900.fna -d ./strep31 -p 15 -c
```

#### 5) Inspect Output:
  
 * Multiple alignment: `parsnp.xmfa`
 * Phylogeny: `parsnp.tree`

This last step requires you to download software and is to highlight the ability to inspect strain-level differences within genomes assembled from metagenomic samples.

## Use AliView (fuhgeddaboudit)

#### Download AliView:

Aliview can be downloaded at: [https://ormbunkar.se/aliview/downloads/](https://ormbunkar.se/aliview/downloads/)

 * Download MFA file:
```
wget https://obj.umiacs.umd.edu/stamps2019/aliview.input.mfa
```

 * Open AliView
      
 * Load MFA file:
	+ File->Open File

