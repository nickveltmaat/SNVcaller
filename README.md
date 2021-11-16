# SNVcaller
#### **Pipeline to call SNV's with 4 tools (VarDict, LoFreq, Mutect2 &amp; SiNVICT)**

### Prerequisites:
 * A process-ready BAM file (or folder containing them), e.g. pre-processed following [GATK best practices workflows; Data pre-processing for variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)
 * [GATK](https://gatk.broadinstitute.org/hc/en-us) > [4.1.4.1](https://github.com/broadinstitute/gatk/releases/tag/4.1.4.1)
 * [Vardict Java](https://bioconda.github.io/recipes/vardict-java/README.html) > 1.8.2
 * [SAMtools](http://www.htslib.org/) > [1.9](http://www.htslib.org/download/)
 * BCFtools > [1.11](http://www.htslib.org/download/)
 * HTSlib > [1.11](http://www.htslib.org/download/)
 * Python > [3.9](https://www.python.org/downloads/release/python-390/)
 * R > [4.0](https://cran.r-project.org/bin/windows/base/)
 * *LoFreq, vt (for post-processing) and SiNVICT with its prerequisites are pre-installed, more info on this in* [__'INSTALLATION'__](https://github.com/nickveltmaat/SNVcaller/blob/main/README.md#installation)

## General Description
This is a pipeline made to reliably generate calls for somatic mutations in Low-VAF samples, such as NGS data from cfDNA. This is done by analyzing .BAM files with 4 different tools (VarDict, LoFreq, Mutect2 & SiNVICT). The pipeline will output variants that are called with at least an `x` amount of tools (this can be set from 1-4). Of course, the higher the number, the lower False Positive calls, the higher the reliability of the call, but also the higher the chance you'll miss relevant somatic variants. The general workflow is as follows: 


## Installation
