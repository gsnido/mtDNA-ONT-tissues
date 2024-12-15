#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Specify sample id (e.g. 6226-CT-0010)"
    exit 2
fi

sample_id=$1

fastqdir="/data1/romain/ont/fastqfiles"
outputdir="/data1/romain/ont/results"
bindir="/data1/romain/ont/bin"

fastqfile="${fastqdir}/${sample_id}.fastq.gz"
if [ ! -s ${fastqfile} ]; then
    echo "Fastq file \"${fastqfile}\" not found"
fi

echo "Input FASTQ file: ${fastqfile}"
echo "Output folder:: ${outputdir}"

# 1) TRIM
echo "1) TRIM"
trimmedfile="${outputdir}/${sample_id}_trimmed.fastq.gz"
## CMD="gunzip -c ${fastqfile} | ${bindir}/chopper -q 10 --headcrop 75 --tailcrop 20 --minlength 1000 --threads 16 | gzip -nc > ${trimmedfile}"
## echo "   ${CMD}"
echo -ne "   Trimming... "
if [ -s ${trimmedfile} ]; then
    echo "[ SKIP ] - File \"${trimmedfile}\" exists"
else
    gunzip -c ${fastqfile} | ${bindir}/chopper -q 10 --headcrop 75 --tailcrop 20 --minlength 1000 --threads 16 | gzip -nc > ${trimmedfile}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi

# 2) ALIGN

## ./minimap2 -x map-ont -d Homo_sapiens_assembly38.chrM.mmi Homo_sapiens_assembly38.chrM.fasta
## ./minimap2 -x map-ont -d Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.mmi Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
reference1="/ref/MT/Homo_sapiens_assembly38.chrM.mmi"
fasta1="/ref/MT/Homo_sapiens_assembly38.chrM.fasta"
reference2="/ref/MT/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.mmi"
fasta2="/ref/MT/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
reference5="/ref/MT/chr22.mmi"
fasta5="/ref/MT/chr22.fasta"

bamfile1="${outputdir}/${sample_id}_chrM.bam"
bamfile3="${outputdir}/${sample_id}_chrM_shifted.bam"
bamfile9="${outputdir}/${sample_id}_chr22.bam"
bamfile9small="${outputdir}/${sample_id}_chr22_small.bam"

stderr="${outputdir}/${sample_id}_aln.stderr"
cat /dev/null > ${stderr}

echo "2) ALIGN AGAINST chrM"

echo -ne "   Aligning with minimap2... "
if [ -s ${bamfile1} ]; then
echo "[ SKIP ] - File \"${bamfile1}\" exists"
else
    minimap2 -t 8 --secondary=no --MD --eqx -a -R "@RG\tID:RG1\tSM:${sample_id}\tPL:Illumina\tLB:Library" --sam-hit-only ${reference1} ${trimmedfile} 2>> ${stderr} | \
        awk '$1 ~ /^@/ || $2 == "0" || $2 == "16"' | \
        samtools view -@ 8 -bS - | \
        samtools sort - > ${bamfile1}
    samtools index ${bamfile1}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi



echo "3) Aligning against SHIFTED chrM"

echo -ne "   Aligning with minimap2... "
if [ -s ${bamfile3} ]; then
    echo "[ SKIP ] - File \"${bamfile3}\" exists"
else
    minimap2 -t 8 --secondary=no --MD --eqx -a -R "@RG\tID:RG1\tSM:${sample_id}\tPL:Illumina\tLB:Library" --sam-hit-only ${reference2} ${trimmedfile} 2>> ${stderr} | \
        awk '$1 ~ /^@/ || $2 == "0" || $2 == "16"' | \
        samtools view -@ 8 -bS - | \
        samtools sort - > ${bamfile3}
    samtools index ${bamfile3}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi



echo "4) Aligning against chr22"
echo -ne "   Aligning with minimap2... "
if [ -s ${bamfile9} ]; then
    echo "[ SKIP ] - File \"${bamfile9}\" exists"
else
    minimap2 -t 8 --secondary=no --MD --eqx -a \
        -R "@RG\tID:RG1\tSM:${sample_id}\tPL:Illumina\tLB:Library" \
        --sam-hit-only ${reference5} ${trimmedfile} 2>> ${stderr} | \
        awk '$1 ~ /^@/ || $2 == "0" || $2 == "16"' | \
        samtools view -@ 8 -bS - | \
        samtools sort - > ${bamfile9}
    samtools index ${bamfile9}
    ${bindir}/samtools-1.9/samtools view -bh -@ 8 ${bamfile9} chr22:20000000-20100000 > ${bamfile9small}
    ${bindir}/samtools-1.9/samtools index ${bamfile9small}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi




### Depth ###
mkdir -p ${outputdir}/depth/

echo "5) Depth of coverage... "

depth1="${outputdir}/depth/${sample_id}_chrM.txt"
depth3="${outputdir}/depth/${sample_id}_chrM_shifted.txt"
depth8="${outputdir}/depth/${sample_id}_chr22.txt"


echo -ne "   chrM... "
if [ -s ${depth1} ]; then
    echo "[ SKIP ] - File \"${depth1}\" exists"
else
    ${bindir}/samtools-1.9/samtools depth -a -r chrM:1-16569 -d 0 ${bamfile1} > ${depth1}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi
echo -ne "   shifted chrM... "
if [ -s ${depth3} ]; then
    echo "[ SKIP ] - File \"${depth3}\" exists"
else
    ${bindir}/samtools-1.9/samtools depth -a -r chrM:1-16569 -d 0 ${bamfile3} > ${depth3}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi

echo -ne "   chr22... "
if [ -s ${depth8} ]; then
    echo "[ SKIP ] - File \"${depth8}\" exists"
else
    ${bindir}/samtools-1.9/samtools depth -a -r chr22:20000000-20100000 -d 0 ${bamfile9} > ${depth8}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi


### Depth Per Nucleotide###

echo "6) Depth of coverage per position and per nucleotide... "

depth9="${outputdir}/depth/${sample_id}_per_nuc_depth_chrM.txt"
depth11="${outputdir}/depth/${sample_id}_per_nuc_depth_chr22.txt"


echo -ne "   chrM and chrM_shifted... "
if [ -s ${depth9} ]; then
    echo "[ SKIP ] - File \"${depth9}\" exists"
else
    bash ${bindir}/depth_per_nucleotide_MT.sh ${sample_id} 2> /dev/null >> ${depth9}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi

echo -ne "   chr22... "
if [ -s ${depth11} ]; then
    echo "[ SKIP ] - File \"${depth11}\" exists"
else
    bash ${bindir}/depth_per_nucleotide.sh ${sample_id} 2> /dev/null >> ${depth11}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi




### Matches / mismatches ###

echo "7) Matches / mismatches per read... "

mm1="${outputdir}/snvs_per_read/${sample_id}.tsv"
mm3="${outputdir}/snvs_per_read/${sample_id}_chr22.tsv"


echo -ne "   chrM, chrM_shifted and chr22... "
if [ -s ${mm1} ] && [ -s ${mm3} ]; then
    echo "[ SKIP ] - Output files exist"
else
    bash ${bindir}/count_matches_mismatches_per_read.sh ${sample_id}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi





### GC content ###

echo "8) GC content per read... "

gc1="${outputdir}/gc_content/${sample_id}_chrM.tsv"
gc3="${outputdir}/gc_content/${sample_id}_chr22.tsv"


echo -ne "   chrM, chrM_shifted and chr22... "
if [ -s ${gc1} ] && [ -s ${gc3} ]; then
    echo "[ SKIP ] - Output files exist"
else
    bash ${bindir}/calc_GC.sh ${sample_id}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi





### Get deletions from CIGAR ###

echo "9) Deletion breakpoints from CIGAR... "

dels1="${outputdir}/deletions/${sample_id}_chrM.tsv"
dels2="${outputdir}/deletions/${sample_id}_chrM_shifted.tsv"
dels3="${outputdir}/deletions/${sample_id}_chr22.tsv"


echo -ne "   chrM, chrM_shifted and chr22... "
if [ -s ${dels1} ] && [ -s ${dels2} ] && [ -s ${dels3} ]; then
    echo "[ SKIP ] - Output files exist"
else
    bash ${bindir}/get_indels.sh ${sample_id}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi


exit 1






####################
#  Variant calling #
####################

echo "10) Variant calling (MUTSERVE)"

output1="${outputdir}/${sample_id}_chrM_variants.txt"
output2="${outputdir}/${sample_id}_chrM_shifted_variants.txt"
# output5="${outputdir}/${sample_id}_chr22_variants.txt"

echo -ne "   Call variants for chrM... "
if [ -s ${output1} ]; then
    echo "[ SKIP ] - File \"${output1}\" exists"
else
    ${bindir}/mutserve call --reference=${fasta1} \
        --output=${output1} --threads=8 --level=0.01 ${bamfile1} >> ${stderr} 2>&1
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi


echo -ne "   Call variants for chrM shifted... "
if [ -s ${output2} ]; then
    echo "[ SKIP ] - File \"${output2}\" exists"
else
    ${bindir}/mutserve call --reference=${fasta2} \
        --output=${output2} --threads=8 --level=0.01 ${bamfile3} >> ${stderr} 2>&1
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi


# echo -ne "   Call variants for chr22... "
# if [ -s ${output5} ]; then
#     echo "[ SKIP ] - File \"${output5}\" exists"
# else
#     ${bindir}/mutserve call --reference=${fasta5} \
#         --output=${output5} --threads=8 --level=0.01 ${bamfile9} >> ${stderr} 2>&1
#     if [[ $? -eq 0 ]]; then
#         echo "[ OK ]"
#     else
#         echo "[ ERROR! ]"
#         #exit 2
#     fi
# fi




echo "11) Variant calling (MUTECT2)"
vcf1="${outputdir}/${sample_id}_chrM.vcf"
vcf2="${outputdir}/${sample_id}_chrM_shifted.vcf"
vcf5="${outputdir}/${sample_id}_chr22.vcf"

echo -ne "   Call variants for chrM... "
if [ -s ${vcf1} ]; then
    echo "[ SKIP ] - File \"${vcf1}\" exists"
else
    ${bindir}/gatk-4.6.0.0/gatk Mutect2 -R ${fasta1} \
        --mitochondria-mode \
        -I ${bamfile1} \
        -L chrM:4000-11999 \
        -O ${vcf1} >> ${stderr} 2>&1
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi


echo -ne "   Call variants for chrM shifted... "
if [ -s ${vcf2} ]; then
    echo "[ SKIP ] - File \"${vcf2}\" exists"
else
    ${bindir}/gatk-4.6.0.0/gatk Mutect2 -R ${fasta2} \
        --mitochondria-mode \
        -I ${bamfile3} \
        -L chrM:4000-12568 \
        -O ${vcf2} >> ${stderr} 2>&1
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi



echo -ne "   Call variants for chr22... "
if [ -s ${vcf5} ]; then
    echo "[ SKIP ] - File \"${vcf5}\" exists"
else
    ${bindir}/gatk-4.6.0.0/gatk Mutect2 -R ${fasta5} -L chr22:20000000-20100000 \
        --mitochondria-mode \
        -I ${bamfile9} \
        -O ${vcf5} >> ${stderr} 2>&1
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi


## LIFTOVER AND COMBINE VCFs
vcf="${outputdir}/${sample_id}_chrM_combined.vcf"
combined_stats="${outputdir}/${sample_id}_combined.stats"
shift_back_chain="/ref/MT/ShiftBack.chain"
echo "9) Liftover and combine VCFs... "
if [ -s ${vcf} ]; then
    echo "[ SKIP ] - File \"${vcf}\" exists"
else
    ${bindir}/gatk-4.6.0.0/gatk LiftoverVcf \
      -I ${vcf2} \
      -O ${vcf2}.shifted_back.vcf \
      -R ${fasta1} \
      -C ${shift_back_chain} \
      --REJECT ${vcf2}.rejected.vcf >> ${stderr} 2>&1

    ${bindir}/gatk-4.6.0.0/gatk MergeVcfs \
      -I ${vcf2}.shifted_back.vcf \
      -I ${vcf1} \
      -O ${vcf} >> ${stderr} 2>&1

    ${bindir}/gatk-4.6.0.0/gatk MergeMutectStats \
        --stats ${vcf1}.stats --stats ${vcf2}.stats \
        -O ${combined_stats} >> ${stderr} 2>&1

    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
        rm ${vcf2}.shifted_back.vcf
        rm ${vcf2}.shifted_back.vcf.idx
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi


# FILTERING
vcf_tmp="${outputdir}/${sample_id}_chrM_combined_tmp.vcf"
vcf_filt="${outputdir}/${sample_id}_chrM_combined_filt.vcf"
vcf_split_tmp="${outputdir}/${sample_id}_chrM_combined_filt_splt_tmp.vcf"
vcf_split="${outputdir}/${sample_id}_chrM_combined_filt_splt.vcf"
blacklisted_sites="/ref/MT/blacklist_sites.hg38.chrM.bed"
echo -ne "   Filtering of variants... "
if [ -s ${vcf_split} ]; then
    echo "[ SKIP ] - File \"${vcf_split}\" exists"
else
    ${bindir}/gatk-4.6.0.0/gatk FilterMutectCalls \
        -V ${vcf} \
        -R ${fasta1} \
        -O ${vcf_tmp} \
        --stats ${combined_stats} \
        --max-alt-allele-count 4 \
        --mitochondria-mode \
        --min-allele-fraction 0 >> ${stderr} 2>&1

    ${bindir}/gatk-4.6.0.0/gatk VariantFiltration \
        -V ${vcf_tmp} \
        -O ${vcf_filt} \
        --apply-allele-specific-filters \
        --mask ${blacklisted_sites} \
        --mask-name "blacklisted_site" >> ${stderr} 2>&1

    ${bindir}/gatk-4.6.0.0/gatk LeftAlignAndTrimVariants \
      -R ${fasta1} \
      -V ${vcf_filt} \
      -O ${vcf_split_tmp} \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac >> ${stderr} 2>&1

    ${bindir}/gatk-4.6.0.0/gatk SelectVariants \
        -V ${vcf_split_tmp} \
        -O ${vcf_split} \
        --exclude-filtered >> ${stderr} 2>&1

    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
        rm ${vcf_tmp}
        rm ${vcf_tmp}.idx
        rm ${vcf_split_tmp}
        rm ${vcf_split_tmp}.idx
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi


# CONTAMINATION
contamination="${outputdir}/${sample_id}_contamination.tsv"

echo -ne "   Check contamination... "
if [ -s ${contamination} ]; then
    echo "[ SKIP ] - File \"${contamination}\" exists"
else
    ${bindir}/haplocheck \
        --out ${contamination} \
        ${vcf_split} >> ${stderr} 2>&1
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        exit 2
    fi
fi








### Deletions from CIGAR strings ###
mkdir -p ${outputdir}/deletions/

echo "10) Deletions in single reads... "

deletions1="${outputdir}/deletions/${sample_id}_chrM.txt"
deletions3="${outputdir}/deletions/${sample_id}_chrM_shifted.txt"
deletions8="${outputdir}/deletions/${sample_id}_chr22.txt"

bam22_small="${outputdir}/${sample_id}_chr22_small.bam"


echo -ne "   chrM... "
if [ -s ${deletions1} ]; then
    echo "[ SKIP ] - File \"${deletions1}\" exists"
else
    python ${bindir}/get_indels_CIGAR.py ${bamfile1} > ${deletions1}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi
echo -ne "   shifted chrM... "
if [ -s ${deletions3} ]; then
    echo "[ SKIP ] - File \"${deletions3}\" exists"
else
    python ${bindir}/get_indels_CIGAR.py ${bamfile3} > ${deletions3}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi

echo -ne "   chr22... "
if [ ! -s ${bam22_small} ]; then
    ${bindir}/samtools-1.9/samtools view -bh -@ 8 ${bamfile9} chr22:20000000-20100000 > ${bam22_small}
    ${bindir}/samtools-1.9/samtools index ${bam22_small}
fi

if [ -s ${deletions8} ]; then
    echo "[ SKIP ] - File \"${deletions8}\" exists"
else
    python ${bindir}/get_indels_CIGAR.py ${bam22_small} > ${deletions8}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi




### Read lengths ###
mkdir -p ${outputdir}/readlengths/

echo "11) Read lengths... "

lengths1="${outputdir}/readlengths/${sample_id}_chrM.txt"
lengths3="${outputdir}/readlengths/${sample_id}_chrM_shifted.txt"
lengths8="${outputdir}/readlengths/${sample_id}_chr22.txt"


echo -ne "   chrM... "
if [ -s ${lengths1} ]; then
    echo "[ SKIP ] - File \"${lengths1}\" exists"
else
    python ${bindir}/get_read_lengths.py ${bamfile1} > ${lengths1}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi

echo -ne "   shifted chrM... "
if [ -s ${lengths3} ]; then
    echo "[ SKIP ] - File \"${lengths3}\" exists"
else
    python ${bindir}/get_read_lengths.py ${bamfile3} > ${lengths3}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi


echo -ne "   chr22... "
if [ -s ${lengths8} ]; then
    echo "[ SKIP ] - File \"${lengths8}\" exists"
else
    python ${bindir}/get_read_lengths.py ${bam22_small} > ${lengths8}
    if [[ $? -eq 0 ]]; then
        echo "[ OK ]"
    else
        echo "[ ERROR! ]"
        #exit 2
    fi
fi




