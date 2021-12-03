#TODO:
# R sort scripts samenvoegen met argumenten
# Sinvict to vcf. py Write loopjes in functie
# Default values for parameters (if argument is given > overwrite..

while getopts "R:L:I:O:V:D:C:P:" arg; do
  case $arg in
    R) Reference=$OPTARG;;      # "/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta"
    L) Listregions=$OPTARG;;    # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/2020_covered_probe_notranslocation.bed"
    I) Inputbam=$OPTARG;;       # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/103584/"
    O) Outputdir=$OPTARG;;
    V) VAF=$OPTARG;;            # 0.004
    D) RDP=$OPTARG;;            # 100
    C) Calls=$OPTARG;;          # 1
    P) PoN=$OPTARG;;            # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/104103_normals/"
  esac
done

cd /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/

echo -e "\nSettings: \n\nInput bam(s): $Inputbam \nReference: $Reference \nPanel: $Listregions\nminimum VAF: $VAF\nminimum Read Depth: $RDP\nminimum overlapping calls: $Calls\nPoN:$PoN"
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
  ./tools/vt/vt normalize -q -n -m -o ./temp/${TOOL}/${TOOL}-decomposed-normalized.vcf -r $Reference ./temp/${TOOL}/${TOOL}-decomposed.vcf
  bgzip ./temp/${TOOL}/${TOOL}-decomposed-normalized.vcf
  bcftools index ./temp/${TOOL}/${TOOL}-decomposed-normalized.vcf.gz  
}

pythonscript() {
  source ./env/bin/activate
  python3 $1 $2 $3 $4
  deactivate
}

create_pon(){
  #Step 1; Tumor only mode on all normal bam files
  for entry in "$PoN"/*.bam
  do
    echo -e '\n\nGenerating VCF for normal control: '
    echo normalcontrol= $(basename $entry .bam)
    gatk Mutect2 --verbosity ERROR -R $Reference -I $entry --max-mnp-distance 0 -O ./PoN/normal_$(basename $entry .bam).vcf.gz --native-pair-hmm-threads 9 -L ./PoN/newbed.bed
  done
  
  #Step 2 merge pon data
  echo -e '\nMerging normal vcfs into GenomicsDB \n\n'
  ls ./PoN/*vcf.gz > ./PoN/normals.dat
  sed -i -e 's/^/ -V /' ./PoN/normals.dat
  xargs -a ./PoN/normals.dat gatk --java-options "-Xmx8g" GenomicsDBImport --verbosity ERROR --genomicsdb-workspace-path ./PoN/controls_pon_db_chr -R $Reference -L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y

  #Step 3 create PoN VCF file 
  echo -e '\nGenerating Panel Of Normals VCF file: \n\n'
  gatk --java-options "-Xmx8g" CreateSomaticPanelOfNormals \
    --verbosity ERROR -R $Reference \
    -V gendb://$PWD/PoN/controls_pon_db_chr \
    -O ./PoN/merged_PoN.vcf
  
  echo -e '\ndecomposing and normalizing PoN VCF file: \n\n' 
  ./tools/vt/vt decompose -s -o ./PoN/merged_PoN-decomposed.vcf ./PoN/merged_PoN.vcf
  ./tools/vt/vt normalize -q -n -m -o ./PoN/merged_PoN-decomposed-normalized.vcf -r $Reference ./PoN/merged_PoN-decomposed.vcf
}


run_tools() {
  # Creating temp dirs
  mkdir ./temp && mkdir ./temp/M2 && mkdir ./temp/VD && mkdir ./temp/LF && mkdir ./temp/SV
  mkdir ./temp/output-readcount && mkdir ./temp/output-sinvict  
  echo -e '\nPre-processing bam file... \n\n'
  ## Preprocessing .bam and (removing 'chr')
  sed 's/chr//g' $Listregions > ./temp/newbed.bed
  samtools view -H $1 | sed 's/chr//g' > ./temp/header.sam
  samtools reheader ./temp/header.sam $1 > ./temp/newbam.bam
  samtools index ./temp/newbam.bam

  ############# TOOLS: 
  echo -e '\nRunning all 4 tools simultaneously...\n'
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
  pythonscript ./sinvict_to_vcf.py ./temp/output-sinvict/calls_level1_sorted.sinvict ./temp/SV/SV.vcf $Reference
  postprocess_vcf "SV"
  
  ### Comparing Mutation call overlaps
  wait
  echo -e '\nData processed with all 4 tools\n\n'
}

process_bam() {
  a=$1
  xbase=${a##*/}
  xpref=${xbase%.*}
  echo -e "\nProcessing sample $xpref\n"  
  run_tools $a
  mkdir ./output/${xpref} 
  echo -e '\nComparing SNV Tools: \n'
  bcftools isec \
    -p ./output/${xpref} \
    -O v \
    -n +$Calls \
    ./temp/SV/SV-decomposed-normalized.vcf.gz \
    ./temp/M2/M2-decomposed-normalized.vcf.gz \
    ./temp/LF/LF-decomposed-normalized.vcf.gz \
    ./temp/VD/VD-decomposed-normalized.vcf.gz  
  
  if [ -z "$PoN" ]
  then
    echo ""
  else
    echo -e "source for PoN .bam files = $PoN \n\nBlacklisting...\n"
    pythonscript ./pon_blacklist.py ${xpref}
    pythonscript ./create_plots.py ${xpref} "sites_PoN.txt"
  fi    
  pythonscript ./create_plots.py ${xpref} "sites.txt"
  
  ##Test ANNOTATION (make it pon or not)
  echo 'annotating variants...'
  source ./env/bin/activate
  oc run ./output/${xpref}/sites.txt \
    -l hg19 \
    -n annotated_SNVs --silent \
    -a clinvar civic_gene cgc cgl \
    -t excel
  deactivate 
  ## EndTest
  rm -rf ./temp  
  echo -e "\nAnalysis of $xpref is complete!\n\n\n"
}

# Create PoN if argument is given:
if [ -z "$PoN" ]
then
  echo -e "No PoN mode\n\n"
else
  echo -e "\nPoN included: source for PoN .bam files = $PoN"
  echo -e 'Generating PoN from normal samples: \n\n'
  mkdir ./PoN/
  sed 's/chr//g' $Listregions > ./PoN/newbed.bed
  create_pon
fi

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

echo -e '\nFinished run! \n'
