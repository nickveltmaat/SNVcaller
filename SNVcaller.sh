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
mkdir ./temp && mkdir ./temp/M2 && mkdir ./temp/VD && mkdir ./temp/LF && mkdir ./temp/SV
mkdir ./temp/output-readcount && mkdir ./temp/output-sinvict

echo -e '\nPre-process bam and bed file: \n\n'
## Preprocessing .bam and .bed (removing 'chr')
samtools view -H $Inputbam | sed 's/xxxxx//g' > ./temp/header.sam
samtools reheader ./temp/header.sam $Inputbam > ./temp/newbam.bam
samtools index ./temp/newbam.bam
sed 's/xxxxx//g' $Listregions > ./temp/newbed.bed

############# TOOLS: 
## Running all 4 tools simultaneously
#Lofreq
./tools/lofreq/src/lofreq/lofreq call-parallel --pp-threads 8 --call-indels --no-default-filter -f $Reference -o ./temp/LF/output_lofreq_bed.vcf -l ./temp/newbed.bed ./temp/newbam.bam && \
./tools/lofreq/src/lofreq/lofreq filter -i ./temp/LF/output_lofreq_bed.vcf -o ./temp/LF/output_lofreq.vcf -v 100 -a 0.01 &
#Mutect2
gatk Mutect2 --verbosity ERROR --callable-depth 100 --minimum-allele-fraction 0.01 -R $Reference -I ./temp/newbam.bam -O ./temp/M2/M2Unfiltered.vcf --native-pair-hmm-threads 8 -L ./temp/newbed.bed && \
gatk FilterMutectCalls --verbosity ERROR -R $Reference -V ./temp/M2/M2Unfiltered.vcf -O ./temp/M2/M2filtered.vcf &
#VarDict
~/.conda/pkgs/vardict-java-1.8.2-hdfd78af_3/bin/vardict-java -th 8 -G $Reference -b ./temp/newbam.bam -c 1 -S 2 -E 3 ./temp/newbed.bed > ./temp/VD/vardict_raw.vcf &
#SiNVICT
./tools/sinvict/bam-readcount/build/bin/bam-readcount -l ./temp/newbed.bed -w 1 -f $Reference ./temp/newbam.bam > ./temp/output-readcount/output.readcount && ./tools/sinvict/sinvict -m 100 -f 0.01 -t ./temp/output-readcount/ -o ./temp/output-sinvict/ &
wait

## Processing LoFreq:
echo -e '\nProcessing LoFreq data: \n\n'
./tools/vt/vt decompose -s -o ./temp/LF/output_lofreq-decomposed.vcf ./temp/LF/output_lofreq.vcf
./tools/vt/vt normalize -q -n -m -o ./temp/LF/output_lofreq-decomposed-normalized.vcf -r ../hg19upper.fna ./temp/LF/output_lofreq-decomposed.vcf
bgzip ./temp/LF/output_lofreq-decomposed-normalized.vcf
bcftools index ./temp/LF/output_lofreq-decomposed-normalized.vcf.gz


## Processing Mutect2:
echo -e '\nProcessing Mutect2: \n\n'
./tools/vt/vt decompose -s -o ./temp/M2/M2filtered-decomposed.vcf ./temp/M2/M2filtered.vcf
./tools/vt/vt normalize -q -n -m -o ./temp/M2/M2filtered-decomposed-normalized.vcf -r ../hg19upper.fna ./temp/M2/M2filtered-decomposed.vcf
bgzip ./temp/M2/M2filtered-decomposed-normalized.vcf
bcftools index ./temp/M2/M2filtered-decomposed-normalized.vcf.gz

## Processing VarDict
echo -e '\nProcessing Vardict data: \n\n'
Rscript ./sort_vardict.R
source ./env/bin/activate
python3 ./vardict_vcf.py
deactivate
./tools/vt/vt decompose -s -o ./temp/VD/vardict_out-decomposed.vcf ./temp/VD/vardict_out.vcf
./tools/vt/vt normalize -q -n -m -o ./temp/VD/vardict_out-decomposed-normalized.vcf -r ../hg19upper.fna ./temp/VD/vardict_out-decomposed.vcf
bgzip ./temp/VD/vardict_out-decomposed-normalized.vcf
bcftools index ./temp/VD/vardict_out-decomposed-normalized.vcf.gz


## Processing SiNVICT
echo -e '\nProcessing SiNVICT data: \n\n'
Rscript ./sort_sinvict.R
source ./env/bin/activate
python ./sinvict_to_vcf2.py ./temp/output-sinvict/calls_level1_sorted.sinvict ./temp/SV/sinvict-vcf-summary.vcf
deactivate
./tools/vt/vt decompose -s -o ./temp/SV/sinvict-vcf-summary-decomposed.vcf  ./temp/SV/sinvict-vcf-summary.vcf 
./tools/vt/vt normalize -q -n -m -o ./temp/SV/sinvict-vcf-summary-decomposed-normalized.vcf  -r ../hg19upper.fna ./temp/SV/sinvict-vcf-summary-decomposed.vcf 
bgzip ./temp/SV/sinvict-vcf-summary-decomposed-normalized.vcf
bcftools index ./temp/SV/sinvict-vcf-summary-decomposed-normalized.vcf.gz


### Comparing Mutation call overlaps
wait

echo -e '\nData processed with all 4 tools\n\n'
