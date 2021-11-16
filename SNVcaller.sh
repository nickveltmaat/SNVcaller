# -R = reference = hg19.fa
# -L = list of regions = UCSC_Combined.bed
# -I = Inputfile.bam
# -O = output dir (location)

while getopts "R:L:I:O:" arg; do
  case $arg in
    R) Reference=$OPTARG;;
    L) Listregions=$OPTARG;;
    I) Inputbam=$OPTARG;;
    O) Outputdir=$OPTARG;;
  esac
done

# Creating temp dirs 
mkdir SNVcallingPipeline/temp
mkdir SNVcallingPipeline/temp/LF && mkdir SNVcallingPipeline/temp/M2 && mkdir SNVcallingPipeline/temp/SV && mkdir SNVcallingPipeline/temp/VD
mkdir SNVcallingPipeline/temp/output-readcount && mkdir SNVcallingPipeline/temp/output-sinvict

### SINVICT
echo -e '\nRunning SiNVICT: Creating bam-readcount files\n\n'
sinvict/bam-readcount/build/bin/bam-readcount -l $Listregions -w 1 -f $Reference $Inputbam > SNVcallingPipeline/temp/output-readcount/output.readcount
echo -e '\nReadcount file generated; running SiNVICT\n\n'
sinvict/sinvict -t SNVcallingPipeline/temp/output-readcount/ -o SNVcallingPipeline/temp/output-sinvict/

Rscript SNVcallingPipeline/sort_sinvict.R
python SNVcallingPipeline/sinvict_to_vcf.py SNVcallingPipeline/temp/output-sinvict/calls_level1_sorted.sinvict SNVcallingPipeline/temp/SV/sinvict-vcf-summary.vcf
bgzip SNVcallingPipeline/temp/SV/sinvict-vcf-summary.vcf 
bcftools index SNVcallingPipeline/temp/SV/sinvict-vcf-summary.vcf.gz

### Mutect 2
echo -e '\nRunning Mutect2: \n\n'
gatk/gatk Mutect2 -R $Reference -I $Inputbam -O SNVcallingPipeline/temp/M2/M2Unfiltered.vcf --native-pair-hmm-threads 40 -L $Listregions 
gatk/gatk FilterMutectCalls -R $Reference -V SNVcallingPipeline/temp/M2/M2Unfiltered.vcf -O SNVcallingPipeline/temp/M2/M2filtered.vcf
bgzip SNVcallingPipeline/temp/M2/M2filtered.vcf
bcftools index SNVcallingPipeline/temp/M2/M2filtered.vcf.gz

### LoFreq
echo -e '\nRunning LoFreq: \n\n'
lofreq_star-2.1.5_linux-x86-64/bin/lofreq call-parallel --pp-threads 40 --call-indels -f $Reference -o SNVcallingPipeline/temp/LF/output_lofreq_bed.vcf -l $Listregions $Inputbam
bgzip SNVcallingPipeline/temp/LF/output_lofreq_bed.vcf
bcftools index SNVcallingPipeline/temp/LF/output_lofreq_bed.vcf.gz

### Vardict 
echo -e '\nRunning Vardict: \n\n'
vardict -G $Reference -b $Inputbam -c 1 -S 2 -E 3 $Listregions > SNVcallingPipeline/temp/VD/vardict_raw.vcf
Rscript SNVcallingPipeline/sort_vardict.R
python SNVcallingPipeline/vardict_vcf.py
bgzip SNVcallingPipeline/temp/VD/vardict_out.vcf
bcftools index SNVcallingPipeline/temp/VD/vardict_out.vcf.gz

### Comparing Mutation call overlaps
echo -e '\nComparing SNV Tools: \n'
bcftools isec -p SNVcallingPipeline/output/ -O v -n +2 SNVcallingPipeline/temp/SV/sinvict-vcf-summary.vcf.gz SNVcallingPipeline/temp/M2/M2filtered.vcf.gz SNVcallingPipeline/temp/LF/output_lofreq_bed.vcf.gz SNVcallingPipeline/temp/VD/vardict_out.vcf.gz

echo -e '\nDONE, GREAT SUCCES\n\n'


rm -rf SNVcallingPipeline/temp
