#TODO:
# Stringdoc alle scriptjes

while getopts "R:L:I:O:V:D:C:" arg; do
  case $arg in
    R) Reference=$OPTARG;;
    L) Listregions=$OPTARG;;
    I) Inputbam=$OPTARG;;
    O) Outputdir=$OPTARG;;
    V) VAF=$OPTARG;;
    D) RDP=$OPTARG;;
    C) Calls=$OPTARG;;
  esac
done

cd /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/

echo -e '\nLoading modules: \n'
module load GATK/4.1.4.0-Java-8-LTS
module load SAMtools/1.9-GCCcore-7.3.0
module load BCFtools/1.11-GCCcore-7.3.0
module load HTSlib/1.11-GCCcore-7.3.0
module load Python/3.9.1-GCCcore-7.3.0-bare
module load R/4.0.3-foss-2018b-bare
module list

#rm -rf ./temp
mkdir ./output


postprocess_vcf() {
  TOOL=$1
  ./tools/vt/vt decompose -s -o ./temp/${TOOL}/${TOOL}-decomposed.vcf ./temp/${TOOL}/${TOOL}.vcf
  ./tools/vt/vt normalize -q -n -m -o ./temp/${TOOL}/${TOOL}-decomposed-normalized.vcf -r ../hg19upper.fna ./temp/${TOOL}/${TOOL}-decomposed.vcf
  bgzip ./temp/${TOOL}/${TOOL}-decomposed-normalized.vcf
  bcftools index ./temp/${TOOL}/${TOOL}-decomposed-normalized.vcf.gz  
}

pythonscript() {
  source ./env/bin/activate
  python3 $1 $2 $3
  deactivate
}

run_tools() {
  # Creating temp dirs
  mkdir ./temp && mkdir ./temp/M2 && mkdir ./temp/VD && mkdir ./temp/LF && mkdir ./temp/SV
  mkdir ./temp/output-readcount && mkdir ./temp/output-sinvict  
  echo -e '\nPre-process bam and bed file: \n\n'
  ## Preprocessing .bam and .bed (removing 'chr')
  samtools view -H $1 | sed 's/xxxxx//g' > ./temp/header.sam
  samtools reheader ./temp/header.sam $1 > ./temp/newbam.bam
  samtools index ./temp/newbam.bam
  sed 's/xxxxx//g' $Listregions > ./temp/newbed.bed
  
  ############# TOOLS: 
  ## Running all 4 tools simultaneously
  #Lofreq
  ./tools/lofreq/src/lofreq/lofreq call-parallel \
     --pp-threads 8 \
     --call-indels \
     --no-default-filter \
     -f $Reference \
     -o ./temp/LF/output_lofreq_bed.vcf \
     -l ./temp/newbed.bed \
     ./temp/newbam.bam && \
  #Filter Lofreq     
  ./tools/lofreq/src/lofreq/lofreq filter \
    -i ./temp/LF/output_lofreq_bed.vcf \
      -o ./temp/LF/LF.vcf \
      -v $RDP \
      -a $VAF & \
      
  #Mutect2
  gatk Mutect2 \
    --verbosity ERROR \
    --callable-depth $RDP \
    --minimum-allele-fraction $VAF \
    -R $Reference \
    -I ./temp/newbam.bam \
    -O ./temp/M2/M2Unfiltered.vcf \
    --native-pair-hmm-threads 8 \
    -L ./temp/newbed.bed && \
  #Filter Mutect2    
  gatk FilterMutectCalls \
    --verbosity ERROR \
    --min-allele-fraction $VAF \
    -R $Reference \
    -V ./temp/M2/M2Unfiltered.vcf \
    -O ./temp/M2/M2.vcf & \
    
  #VarDict
  ~/.conda/pkgs/vardict-java-1.8.2-hdfd78af_3/bin/vardict-java \
    -f $VAF \
    -th 8 \
    -G $Reference \
    -b ./temp/newbam.bam \
    -c 1 -S 2 -E 3 \
    ./temp/newbed.bed > ./temp/VD/vardict_raw.vcf & \
    
  #SiNVICT
  ./tools/sinvict/bam-readcount/build/bin/bam-readcount \
    -l ./temp/newbed.bed \
    -w 1 \
    -f $Reference \
    ./temp/newbam.bam > ./temp/output-readcount/output.readcount && \
    
    ./tools/sinvict/sinvict \
    -m $RDP \
    -f $VAF \
    -t ./temp/output-readcount/ \
    -o ./temp/output-sinvict/ & \
  wait
  
  ## Processing LoFreq:
  echo -e '\n\n\nProcessing LoFreq data: \n\n'
  postprocess_vcf "LF"
  
  ## Processing Mutect2:
  echo -e '\n\n\nProcessing Mutect2: \n\n'
  postprocess_vcf "M2"
  
  ## Processing VarDict
  echo -e '\n\n\nProcessing Vardict data: \n\n'
  Rscript ./sort_vardict.R
  pythonscript ./vardict_vcf.py
  postprocess_vcf "VD"
  
  ## Processing SiNVICT
  echo -e '\n\n\nProcessing SiNVICT data: \n\n'
  Rscript ./sort_sinvict.R
  pythonscript ./sinvict_to_vcf.py ./temp/output-sinvict/calls_level1_sorted.sinvict ./temp/SV/SV.vcf
  postprocess_vcf "SV"
  
  ### Comparing Mutation call overlaps
  wait
  echo -e '\nData processed with all 4 tools\n\n'
}

process_bam() {
  a=$1
  xbase=${a##*/}
  xpref=${xbase%.*}
  echo sample= ${xpref}
  mkdir ./output/${xpref}  
  run_tools $a
  echo -e '\nComparing SNV Tools: \n'
  bcftools isec \
    -p ./output/${xpref} \
    -O v \
    -n +$Calls \
    ./temp/SV/SV-decomposed-normalized.vcf.gz \
    ./temp/M2/M2-decomposed-normalized.vcf.gz \
    ./temp/LF/LF-decomposed-normalized.vcf.gz \
    ./temp/VD/VD-decomposed-normalized.vcf.gz  
  pythonscript ./create_plots.py ${xpref}
  rm -rf ./temp  
}

#Directory:
if [[ -d $Inputbam ]]; then
    for entry in "$Inputbam"/*.bam
    do
      process_bam $entry
    done
    
#One-File:     
elif [[ -f $Inputbam ]]; then
    process_bam $Inputbam

#Else (error)
else
    echo "$Inputbam is not valid"
    exit 1
fi

echo -e '\nFinished! \n'

