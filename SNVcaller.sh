#!/bin/bash
#SBATCH --job-name=SNVcaller_1_1
#SBATCH --output="SNVtest.out"
#SBATCH --error="SNVtest.err"
#SBATCH --time=16:00:00
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --get-user-env=L60
#SBATCH --export=NONE

#TODO:
# R sort scripts samenvoegen met argumenten
# Sinvict to vcf. py Write loopjes in functie
# Default values for parameters (if argument is given > overwrite..
# Filter Lofreq functie
# Add number of CPU argument (X ofzo)
# StrandBias toevoegen

# PoN when newly generated: Filter P/LP & on genomic location (intronic / synonymous / etc.) annotated blacklist can be used for this
# Save list from P/LP filterd-out variants from PoN


while getopts "R:L:I:O:V:D:C:P:Q:B:M:" arg; do 
  case $arg in
    R) Reference=$OPTARG;;      # "/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta" OR "/groups/umcg-pmb/tmp01/apps/data/reference_sequences/Homo_sapiens_assembly19_1000genomes_decoy.fasta"
    L) Listregions=$OPTARG;;    # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/2020_covered_probe_notranslocation.bed"
    I) Inputbam=$OPTARG;;       # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/103584/"
    O) Outputdir=$OPTARG;;
    V) VAF=$OPTARG;;            # 0.004
    D) RDP=$OPTARG;;            # 100
    C) Calls=$OPTARG;;          # 1
    P) PoN=$OPTARG;;            # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/104103_normals/"
    Q) Qual=$OPTARG;;           # BaseQuality (18)
    B) Bias=$OPTARG;;           # Strand Bias (5-95)
    M) MRD=$OPTARG;;            # Mimimal Mutant Read depth: 3
  esac
done

cd /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/

echo -e "\nSettings: \n\nInput bam(s): $Inputbam \nReference: $Reference \nPanel: $Listregions\nminimum VAF: $VAF\nminimum Read Depth: $RDP\nminimum overlapping calls: $Calls\nPoN:$PoN"
echo -e '\nLoading modules: \n'
module load GATK/4.1.4.1-Java-8-LTS
module load libjpeg-turbo/2.0.2-GCCcore-7.3.0
module load picard
module load SAMtools/1.9-GCCcore-7.3.0
module load BCFtools/1.11-GCCcore-7.3.0
module load HTSlib/1.11-GCCcore-7.3.0
#module load Python/3.9.1-GCCcore-7.3.0-bare
module load R/4.0.3-foss-2018b-bare

module list

mkdir ./output

postprocess_vcf() {
  Input=$1
  ./tools/vt/vt decompose -s -o $(dirname $Input)/$(basename $Input .vcf)-decomposed.vcf ${Input}
  ./tools/vt/vt normalize -q -n -m -o $(dirname $Input)/$(basename $Input .vcf)-decomposed-normalized.vcf -r $Reference $(dirname $Input)/$(basename $Input .vcf)-decomposed.vcf
  bgzip $(dirname $Input)/$(basename $Input .vcf)-decomposed-normalized.vcf
  bcftools index $(dirname $Input)/$(basename $Input .vcf)-decomposed-normalized.vcf.gz
}

pythonscript() {
  source ./env/bin/activate
  python3 $1 $2 $3 $4 $5 $6 $7
  deactivate
}

annotate() {
  source ./env/bin/activate
  oc run $1 \
    -l hg19 \
    -n $2 --silent \
    -a clinvar cgc cadd chasmplus clinpred cosmic dbsnp mutation_assessor\
        thousandgenomes thousandgenomes_european vest cadd_exome gnomad3  \
    -t excel
  deactivate 
}

mutect() { # $1 = sample, $2 = output, $3 = bed, $4 --max-mnp-distance 0 (PoN mode)
  gatk Mutect2 \
    --verbosity ERROR \
    --callable-depth $RDP \
    --minimum-allele-fraction $VAF \
    --min-base-quality-score $Qual \
    -R $Reference \
    -I $1 \
    -O $2 \
    --native-pair-hmm-threads 9 \
    -L $3 --QUIET $4
}

filter_mutect() { # $1 = sample, $2 = output
gatk FilterMutectCalls \
  --verbosity ERROR \
  --min-allele-fraction $VAF \
  -R $Reference \
  -V $1 \
  -O $2 --QUIET
}

vardict() { # $1 = sample, $2 = output, $3 = bed
  ~/.conda/pkgs/vardict-java-1.8.2-hdfd78af_3/bin/vardict-java \
    -f $VAF \
    -th 9 \
    -G $Reference \
    -b $1 \
    -q $Qual \
    -c 1 -S 2 -E 3 \
    $3 > $2
}

lofreq() { # $1 = sample, $2 = output, $3 = bed
  ./tools/lofreq/src/lofreq/lofreq call-parallel \
     --pp-threads 9 \
     --call-indels \
     --no-default-filter \
     -f $Reference \
     -o $2 \
     -l $3 -q $Qual -Q $Qual \
     $1
}

filter_lofreq() { # $1 = sample, $2 = output
  ./tools/lofreq/src/lofreq/lofreq filter \
    -i $1 \
    -o $2 \
    -v $RDP \
    -a $VAF
}

sinvict() { # $1 = sample, $2 = output, $3 = bed, $4 temp_readcount
./tools/sinvict/bam-readcount/build/bin/bam-readcount \
  -l $3 \
  -w 1 \
  -b $Qual \
  -f $Reference \
  $1 > $4/output.readcount && \
    
./tools/sinvict/sinvict \
  -m $RDP \
  -f $VAF \
  -t $4 \
  -o $2
}


create_pon() {
  echo -e "\nPoN included: source for PoN .bam files = $PoN"
  #Step 1; Tumor only mode on all normal bam files
  for entry in "$PoN"/*.bam
  do
    echo -e '\n\nGenerating VCF for normal control: '
    echo -e normalcontrol= $(basename $entry .bam) '\n'
    #Mutect2
    mutect $entry ./PoN/normal_$(basename $entry .bam)_Mutect_uf.vcf ./PoN/newbed.bed "--max-mnp-distance 0" && \
    filter_mutect ./PoN/normal_$(basename $entry .bam)_Mutect_uf.vcf ./PoN/normal_$(basename $entry .bam)_Mutect.vcf && \
    pythonscript ./hardfilter.py ${MRD} 'PoN' 'Mutect' ./PoN/normal_$(basename $entry .bam)_Mutect_F.vcf $RDP $(basename $entry .bam) && \
    gatk IndexFeatureFile -I ./PoN/normal_$(basename $entry .bam)_Mutect_F.vcf --QUIET --verbosity ERROR & \
    ########
     
    #Vardict
    vardict $entry ./PoN/normal_$(basename $entry .bam)_Vardict1.vcf ./PoN/newbed.bed && \
    Rscript ./sort_vardict.R ./PoN/normal_$(basename $entry .bam)_Vardict1.vcf ./PoN/normal_$(basename $entry .bam)_Vardict2.vcf && \
    pythonscript ./vardict_vcf.py ./PoN/normal_$(basename $entry .bam)_Vardict2.vcf ./PoN/normal_$(basename $entry .bam)_Vardict.vcf && \
    pythonscript ./hardfilter.py ${MRD} 'PoN' 'Vardict' ./PoN/normal_$(basename $entry .bam)_Vardict_F.vcf $RDP $(basename $entry .bam) && \
    bgzip ./PoN/normal_$(basename $entry .bam)_Vardict_F.vcf && \
    bcftools index ./PoN/normal_$(basename $entry .bam)_Vardict_F.vcf.gz & \
    
    #Lofreq
    lofreq $entry ./PoN/normal_$(basename $entry .bam)_LofreqUnfiltered.vcf ./PoN/newbed.bed && \
    ./tools/lofreq/src/lofreq/lofreq filter \
      -i ./PoN/normal_$(basename $entry .bam)_LofreqUnfiltered.vcf \
      -o ./PoN/normal_$(basename $entry .bam)_Lofreq.vcf  \
      -v $RDP \
      -a $VAF && \
    pythonscript ./hardfilter.py ${MRD} 'PoN' 'Lofreq' ./PoN/normal_$(basename $entry .bam)_Lofreq_F.vcf $RDP $(basename $entry .bam) && \
    bgzip ./PoN/normal_$(basename $entry .bam)_Lofreq_F.vcf && \
    bcftools index ./PoN/normal_$(basename $entry .bam)_Lofreq_F.vcf.gz & \
    
    #SiNVICT
    mkdir ./PoN/$(basename $entry .bam)output-readcount/ && mkdir ./PoN/$(basename $entry .bam)output-sinvict/ && \
    sinvict $entry ./PoN/$(basename $entry .bam)output-sinvict/ ./PoN/newbed.bed  ./PoN/$(basename $entry .bam)output-readcount/ && \
    rm -rf ./PoN/$(basename $entry .bam)output-readcount/ && \
    Rscript ./sort_sinvict.R ./PoN/$(basename $entry .bam)output-sinvict/calls_level1.sinvict ./PoN/$(basename $entry .bam)output-sinvict/calls_level1_sorted.sinvict && \
    pythonscript ./create_vcf.py ./PoN/$(basename $entry .bam)output-sinvict/calls_level1_sorted.sinvict ./PoN/normal_$(basename $entry .bam)_Sinvict.vcf $Reference 'sinvict' && \
    pythonscript ./hardfilter.py ${MRD} 'PoN' 'Sinvict' ./PoN/normal_$(basename $entry .bam)_Sinvict_F.vcf $RDP $(basename $entry .bam) && \
    bgzip ./PoN/normal_$(basename $entry .bam)_Sinvict_F.vcf && \
    bcftools index ./PoN/normal_$(basename $entry .bam)_Sinvict_F.vcf.gz && \
    wait
  done
  
  #Step 2 merge pon data per tool
  echo -e '\nMerging normal vcfs into GenomicsDB \n\n'
  ls ./PoN/*Mutect_F.vcf > ./PoN/M2normals.dat
  ls ./PoN/*Vardict_F.vcf.gz > ./PoN/VDnormals.dat
  ls ./PoN/*Lofreq_F.vcf.gz > ./PoN/LFnormals.dat
  ls ./PoN/*Sinvict_F.vcf.gz > ./PoN/SVnormals.dat
  sed -i -e 's/^/ -V /' ./PoN/M2normals.dat
  sed -i -e 's/^/ /' ./PoN/VDnormals.dat
  sed -i -e 's/^/ /' ./PoN/LFnormals.dat
  sed -i -e 's/^/ /' ./PoN/SVnormals.dat
  xargs -a ./PoN/M2normals.dat gatk --java-options "-Xmx8g" GenomicsDBImport --verbosity ERROR --genomicsdb-workspace-path ./PoN/M2controls_pon_db_chr -R $Reference -L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y --QUIET

  #Step 3 create merged PoN VCF file 
  echo -e '\nGenerating Panel Of Normals VCF file: \n'
  gatk --java-options "-Xmx8g" CreateSomaticPanelOfNormals \
    --verbosity ERROR -R $Reference \
    --QUIET \
    -V gendb://$PWD/PoN/M2controls_pon_db_chr \
    -O ./PoN/merged_PoN_Mutect2.vcf --min-sample-count 1 ###CHANGED (n samples)
  
  xargs -a ./PoN/VDnormals.dat bcftools isec -o ./PoN/Xmerged_PoN_Vardict.vcf -O v -n +1 #CHANGED(n samples)
  xargs -a ./PoN/LFnormals.dat bcftools isec -o ./PoN/Xmerged_PoN_Lofreq.vcf -O v -n +1  #CHANGED(n samples)
  xargs -a ./PoN/SVnormals.dat bcftools isec -o ./PoN/Xmerged_PoN_Sinvict.vcf -O v -n +1 #CHANGED(n samples)
  echo -e '\ndecomposing and normalizing PoN VCF file: \n\n' 

  pythonscript ./create_vcf.py ./PoN/Xmerged_PoN_Vardict.vcf ./PoN/merged_PoN_Vardict.vcf $Reference 'bcftools_isec'
  pythonscript ./create_vcf.py ./PoN/Xmerged_PoN_Lofreq.vcf ./PoN/merged_PoN_Lofreq.vcf $Reference 'bcftools_isec'
  pythonscript ./create_vcf.py ./PoN/Xmerged_PoN_Sinvict.vcf ./PoN/merged_PoN_Sinvict.vcf $Reference 'bcftools_isec'
  postprocess_vcf "./PoN/merged_PoN_Vardict.vcf"
  postprocess_vcf "./PoN/merged_PoN_Lofreq.vcf"
  postprocess_vcf "./PoN/merged_PoN_Sinvict.vcf"
  postprocess_vcf "./PoN/merged_PoN_Mutect2.vcf"
  
  bcftools isec -o ./PoN/BLACKLIST_uf.txt -O v -n +1 ./PoN/merged_PoN_Sinvict-decomposed-normalized.vcf.gz ./PoN/merged_PoN_Mutect2-decomposed-normalized.vcf.gz ./PoN/merged_PoN_Lofreq-decomposed-normalized.vcf.gz ./PoN/merged_PoN_Vardict-decomposed-normalized.vcf.gz  #CHANGED(n tools)

  #annotating & Filtering the PoN blacklist  
  echo -e '\n Annotating unfiltered PoN Blacklist...\n'
  annotate "./PoN/BLACKLIST_uf.txt" "annotated_blacklist_uf"
  echo -e '\nCleaning PoN: Removing P/LP variants & Synonymous, Intronic, UTR and Up-/downstream variants from PoN...\n'
  pythonscript ./clean_pon.py "./PoN/annotated_blacklist_uf.xlsx" "./PoN/BLACKLIST.txt"
  echo -e '\n Annotating PoN Blacklist...\n'
  annotate "./PoN/BLACKLIST.txt" "annotated_blacklist"
  echo -e '\nDone! PoN is ready for use. Now on to analyzing samples...\n\n'
}


run_tools() {
  # Creating temp dirs
  mkdir ./temp && mkdir ./temp/M2 && mkdir ./temp/VD && mkdir ./temp/LF && mkdir ./temp/SV
  mkdir ./temp/output-readcount && mkdir ./temp/output-sinvict  
  echo -e '\nPre-processing bam file... \n\n'
  ## Preprocessing .bam and (removing 'chr')
  sed 's/chr//g' $Listregions > ./temp/newbed.bed
  
  #samtools view -H $1 | sed 's/chr//g' > ./temp/header.sam
  #samtools reheader ./temp/header.sam $1 > ./temp/newbam.bam
  #cp $1 ./temp/newbam.bam
  #samtools index ./temp/newbam.bam
  samtools index $1
  echo -e 'Collecting Alignment Metrics... \n'
  java -jar ${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics QUIET=true VERBOSITY=ERROR \
    R=$Reference \
    I=./temp/newbam.bam \
    O=./output/${xpref}/alignmentsummarymetrics.txt


  ############# TOOLS: 
  echo -e '\n\n\nRunning all 4 tools simultaneously...\n'
  
  #Lofreq
  lofreq $1 ./temp/LF/output_lofreq_bed.vcf ./temp/newbed.bed && \
  #Filter Lofreq
  filter_lofreq ./temp/LF/output_lofreq_bed.vcf ./temp/LF/LF_1.vcf & \
      
  #Mutect2
  mutect $1 ./temp/M2/M2Unfiltered.vcf ./temp/newbed.bed && \
  #Filter Mutect2
  filter_mutect ./temp/M2/M2Unfiltered.vcf ./temp/M2/M2_1.vcf & \
    
  #VarDict
  vardict $1 ./temp/VD/vardict_raw.vcf ./temp/newbed.bed & \
    
  #SiNVICT
  sinvict $1 ./temp/output-sinvict/ ./temp/newbed.bed ./temp/output-readcount/ & \
  wait
  
  ## Processing LoFreq:
  echo -e '\nProcessing LoFreq data: \n'
  pythonscript ./hardfilter.py ${MRD} 'sample' 'Lofreq' './temp/LF/LF.vcf' $RDP
  postprocess_vcf "./temp/LF/LF.vcf"
  
  ## Processing Mutect2:
  echo -e '\nProcessing Mutect2: \n'
  pythonscript ./hardfilter.py ${MRD} 'sample' 'Mutect' './temp/M2/M2.vcf' $RDP
  postprocess_vcf "./temp/M2/M2.vcf"
  
  ## Processing VarDict
  echo -e '\nProcessing Vardict data: \n'
  Rscript ./sort_vardict.R ./temp/VD/vardict_raw.vcf ./temp/VD/vardict_output_sorted.vcf
  pythonscript ./vardict_vcf.py ./temp/VD/vardict_output_sorted.vcf ./temp/VD/VD_1.vcf
  pythonscript ./hardfilter.py ${MRD} 'sample' 'Vardict' './temp/VD/VD.vcf' $RDP
  postprocess_vcf "./temp/VD/VD.vcf"
  
  ## Processing SiNVICT
  echo -e '\nProcessing SiNVICT data: \n'
  Rscript ./sort_sinvict.R ./temp/output-sinvict/calls_level1.sinvict ./temp/output-sinvict/calls_level1_sorted.sinvict
  pythonscript ./create_vcf.py ./temp/output-sinvict/calls_level1_sorted.sinvict ./temp/SV/SV_1.vcf $Reference 'sinvict'
  pythonscript ./hardfilter.py ${MRD} 'sample' 'Sinvict' './temp/SV/SV.vcf' $RDP
  postprocess_vcf "./temp/SV/SV.vcf"

  wait 
  echo -e '\nData processed with all 4 tools\n'
}

##RUN SAMPLE
process_bam() {
  a=$1
  xbase=${a##*/}
  xpref=${xbase%.*}
  echo -e "Processing sample $xpref\n" 
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

  #Merging and Annotating SNV's
  echo 'Merging variants stats & M2 & VD filtering steps...'   
  pythonscript ./merge_variants.py ${xpref} "sites.txt" $RDP
  echo 'annotating variants...'
  annotate ./output/${xpref}/vcfs_merged.txt "annotated_SNVs"
  
  #Plotting and Annotating blacklisted SNV's
  if [ -z "$PoN" ]
  then
    echo ""
  else
    echo -e "source for PoN .bam files = $PoN \nBlacklisting..."   
    #Filtering, Plotting and Annotating SNV's 
    echo 'post-filtering variants filterd by pon...' 
    pythonscript ./post_filtering.py ${xpref} "PoN"
    echo 'plotting variants filterd by pon...'
    pythonscript ./plotting.py ${xpref} "filtered_variants_PoN.xlsx"
  fi 
  
  echo 'post-filtering variants...'
  pythonscript ./post_filtering.py ${xpref} "sample"
  echo 'plotting variants...'
  pythonscript ./plotting.py ${xpref} "filtered_variants.xlsx"
  
  #Cleaning
  mkdir ./output/${xpref}/annotations
  mv ./output/${xpref}/anno*.* ./output/${xpref}/annotations
  mkdir ./output/${xpref}/plots
  mv ./output/${xpref}/*.html ./output/${xpref}/plots
  mv ./output/${xpref}/*.png ./output/${xpref}/plots
  mkdir ./output/${xpref}/filtering_info
  mv ./output/${xpref}/filtering_info*.* ./output/${xpref}/filtering_info
  mv ./output/${xpref}/removed_variants* ./output/${xpref}/filtering_info
  mkdir ./output/${xpref}/intermediate_files
  mv ./output/${xpref}/*.txt ./output/${xpref}/intermediate_files
    
  rm -rf ./temp  
  echo -e "\nAnalysis of $xpref is complete!\n\n\n"
}


## PON INITIATION:
# Create PoN if argument is given:
if [ -z "$PoN" ]
then
  #No PoN Argument
  echo -e "\nNo PoN argument given\n"
else
  # If PoN Argument is given:
  echo -e "\nPoN argument is given!!\n"
  mkdir ./PoN/
  
  #Directory PoN:
  if [[ -d $PoN ]]; then
    echo 'PoN is directory'
    sed 's/chr//g' $Listregions > ./PoN/newbed.bed
    create_pon
    
  #One-File: a.k.a pre-made blacklist:      
  elif [[ -f $PoN ]]; then
    echo -e "PoN is given as pre-made blacklist: \n$PoN \n"
    cp $PoN ./PoN/BLACKLIST.txt
  
  #Else (error)
  else
    echo "Invalid PoN input. Input folder with .bam files or PoN.txt "
    exit 1
  fi
fi


## RUNNING SAMPLES FROM DIRECTORY OR SINGLE MODE
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

echo -e 'Finished run! \n\n\n'

