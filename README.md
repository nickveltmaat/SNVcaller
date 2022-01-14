# SNVcaller
#### **Pipeline to call SNV's with 4 tools ([`VarDict`](https://pubmed.ncbi.nlm.nih.gov/27060149/), [`LoFreq`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526318/), [`Mutect2`](https://www.biorxiv.org/content/10.1101/861054v1) & [`SiNVICT`](https://pubmed.ncbi.nlm.nih.gov/27531099/))**

DISCLAIMER: This project is under construction and valid results are therefore not guaranteed

### Prerequisites:
 * A process-ready BAM file (or folder containing them), e.g. pre-processed following [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)
 * [GATK](https://gatk.broadinstitute.org/hc/en-us) > [4.1.4.1](https://github.com/broadinstitute/gatk/releases/tag/4.1.4.1)
 * [Vardict Java](https://bioconda.github.io/recipes/vardict-java/README.html) > 1.8.2
 * [SAMtools](http://www.htslib.org/) > [1.9](http://www.htslib.org/download/)
 * BCFtools > [1.11](http://www.htslib.org/download/)
 * HTSlib > [1.11](http://www.htslib.org/download/)
 * Python > [3.9](https://www.python.org/downloads/release/python-390/)
 * R > [4.0](https://cran.r-project.org/bin/windows/base/)
 * *LoFreq, [vt](https://github.com/atks/vt) (for post-processing), and SiNVICT with its prerequisites are pre-built, more info on this in* [__'Installation'__](https://github.com/nickveltmaat/SNVcaller/blob/main/README.md#installation)

## General Description
This is a pipeline made to reliably generate calls for somatic mutations in Low Variant Allele Frequencies (VAF) samples in specific regions, such as NGS data from cfDNA. This is done by analyzing `.BAM` files with 4 different tools (`VarDict`, `LoFreq`, `Mutect2` & `SiNVICT`). The pipeline will output variants that are called with at least an `x` amount of tools (this can be set from 1-4). Of course, the higher the number, the lower False Positive call rate, the higher the reliability of the call, but also the higher the chance you'll miss relevant somatic variants. 

The general workflow in the pipeline is as follows: 

`.BAM` and `.bed` files are copied to a temporary folder, where the processing happens. Ather that, the 4 tools will run in parralel, generating preliminary results which are also stored in the temporary folder. Since `VarDict` and `SiNVICT` don't output regular `.vcf` files, this data first needs to be processed in order to compare the overlapping variants in the `.vcf` files. This processing consists of sorting variants and generating `.vcf` files, which is done with custom Python and R scripts. Then, all `.vcf` files are [decomposed](https://genome.sph.umich.edu/wiki/Vt#Decompose), [normalized](https://genome.sph.umich.edu/wiki/Vt#Normalization), gunzipped and indexed. Finally, with all `.vcf` files processed, the variants can be compared on overlapping variants. All variants called with `x` or more tools will be saved. Also a venn diagram of mutation calls per tool is generated, together with histograms of amount of mutations with a certain VAF & Read Depth. VAF & Read Depth are calculated with the data from `VarDict`, `LoFreq`, `Mutect2`, since `SiNVICT` doensn't output this data. A folder contatining normal samples can be provided to generate a Panal of Normals (PoN) a.k.a. a blacklist. This can be used to filter out SNP's and/or technical artifacts (when the same library prep and sequencing methods are performed as in tumor samples). Finally, the remaining variants will be annotated using [openCRAVAT](https://opencravat.org/). Finally, all mutations (if provided blacklisted- and non blacklisted) are annotated using [OpenCravat](https://open-cravat.readthedocs.io/en/latest/quickstart.html). This is a wrapper around multiple well-known annotating tools, such as [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [dbSNP](https://www.ncbi.nlm.nih.gov/snp/), [COSMIC](https://cancer.sanger.ac.uk/cosmic), [gnomAD](https://gnomad.broadinstitute.org/) and many more. All annotated mutations are saved in an excel file. 

*n.b. : A single `.bam` file or a directory containing `.bam` files can be given as arguments. When a directory is given, the process above will loop over all files, generating output folders for each `.bam` file* 



## Installation
**1. Clone the repo**

`git clone https://github.com/nickveltmaat/SNVcaller`

**2. Set working directory to the repo**

`cd /path/to/SNVCaller`

**3. Create python virtual environment (env)**

`python3 -m venv ./env`

**4. Install needed packages in env with pip3**

```
source ./env/bin/activate
pip3 install numpy
pip3 install pandas
pip3 install venn
pip3 install matplotlib
pip3 install pandas_bokeh
pip3 install glob
pip install open-cravat
oc module install-base
oc module ls -a -t annotator  (this generartes a list of available annotators that can be downloaded)
oc module install clinvar cosmic dbsnp ...  (see https://open-cravat.readthedocs.io/en/latest/1.-Installation-Instructions.html for more detailed instructions)
deactivate
```
**5. [Download](https://drive.google.com/drive/folders/1QBt0NdPqjQU_y-A7omxoyiPfl1DL65Xn?usp=sharing) and copy the pre-built tools to `/path/to/SNVCaller/` and unzip**
 
 `unzip ./tools.zip`

## Usage
Once all tools and pre-requisites are installed correctly, the pipeline can be called with: 

`bash ./SNVcaller.sh ARGUMENTS`

**Required** arguments:
- `-I` Input:              **String**   --> example: `/path/to/input.bam` Either one-file or directory
- `-R` Reference:          **String**    --> example: `/path/to/reference.fa`
- `-L` Regions List:       **String**    --> example: `/path/to/panel.bed`
- `-D` minimum Read Depth:  **Int**       --> example: `100`
- `-V` minimum VAF:         __float *[0-1]*__  --> example: `0.002`
- `-C` minimum Calls:       __Int *[1-4]*__ --> example: `2`
- `-P` Panel of Normal:     **String** --> example: `/path/to/PoN/directory/` Optional


Output will be generated in `/path/to/SNVcaller/output/name_of_.bam_file/`
