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
mkdir /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp
mkdir /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp/LF && mkdir /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp/M2 && mkdir SNVcallingPipeline/temp/SV && mkdir /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp/VD
mkdir /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp/output-readcount && mkdir /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp/output-sinvict

#TODO
#############
## Preprocessing .bam and .bed (removing 'chr') 
#samtools view -H $Inputbam | sed 's/chr//g' > /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp/header.sam
#samtools reheader /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp/header.sam $Inputbam > /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp/newbam.bam
#samtools index /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp/newbam.bam
#sed 's/chr//g' $Listregions > /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/temp/newbed.bed
#############



### Running 4 tools simultaneously
VarDictJava/build/libs/VarDict-1.8.2.jar -G $Reference -b $Inputbam -c 1 -S 2 -E 3 $Listregions > SNVcallingPipeline/temp/VD/vardict_raw.vcf &
sinvict/bam-readcount/build/bin/bam-readcount -l $Listregions -w 1 -f $Reference $Inputbam > SNVcallingPipeline/temp/output-readcount/output.readcount &
gatk/gatk Mutect2 -R $Reference -I $Inputbam -O SNVcallingPipeline/temp/M2/M2Unfiltered.vcf --native-pair-hmm-threads 8 -L $Listregions &
lofreq_star-2.1.5_linux-x86-64/bin/lofreq call-parallel --pp-threads 8 --call-indels -f $Reference -o SNVcallingPipeline/temp/LF/output_lofreq_bed.vcf -l $Listregions $Inputbam &
wait


### SINVICT
echo -e '\nRunning SiNVICT: \n\n'
sinvict/sinvict -t SNVcallingPipeline/temp/output-readcount/ -o SNVcallingPipeline/temp/output-sinvict/
Rscript SNVcallingPipeline/sort_sinvict.R
python SNVcallingPipeline/sinvict_to_vcf.py SNVcallingPipeline/temp/output-sinvict/calls_level1_sorted.sinvict SNVcallingPipeline/temp/SV/sinvict-vcf-summary.vcf
bgzip SNVcallingPipeline/temp/SV/sinvict-vcf-summary.vcf 
bcftools index SNVcallingPipeline/temp/SV/sinvict-vcf-summary.vcf.gz

### Processing Mutect 2
echo -e '\nRunning Mutect2: \n\n'
gatk/gatk FilterMutectCalls -R $Reference -V SNVcallingPipeline/temp/M2/M2Unfiltered.vcf -O SNVcallingPipeline/temp/M2/M2filtered.vcf
bgzip SNVcallingPipeline/temp/M2/M2filtered.vcf
bcftools index SNVcallingPipeline/temp/M2/M2filtered.vcf.gz

### Processing LoFreq
echo -e '\nRunning LoFreq: \n\n'
bgzip SNVcallingPipeline/temp/LF/output_lofreq_bed.vcf
bcftools index SNVcallingPipeline/temp/LF/output_lofreq_bed.vcf.gz

### Processing Vardict 
echo -e '\nProcessing Vardict data: \n\n'
Rscript SNVcallingPipeline/sort_vardict.R
python SNVcallingPipeline/vardict_vcf.py
bgzip SNVcallingPipeline/temp/VD/vardict_out.vcf
bcftools index SNVcallingPipeline/temp/VD/vardict_out.vcf.gz

### Comparing Mutation call overlaps
echo -e '\nComparing SNV Tools: \n'
bcftools isec -p SNVcallingPipeline/output/ -O v -n +1 SNVcallingPipeline/temp/SV/sinvict-vcf-summary.vcf.gz SNVcallingPipeline/temp/M2/M2filtered.vcf.gz SNVcallingPipeline/temp/LF/output_lofreq_bed.vcf.gz SNVcallingPipeline/temp/VD/vardict_out.vcf.gz

echo -e '\nDONE, GREAT SUCCES\n\n'


# rm -rf SNVcallingPipeline/temp
