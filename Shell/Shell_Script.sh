### HuxinLiang

### 2023-08-23

############################## 文件下载及预处理 ########################################
###使用prefetch命令从SRA数据库下载文件
prefetch -v -p --option-file Aim_file_ID.txt output-directory -O /Directory/to/save/

###拆分SRA文件
for i in *.sra; do fastq-dump --gzip --split-3 -O /home/hxl/sequence/sample -A $i; done

###检查拆分后的fastq文件质量
time fastqc -q -t 4 -o /home/hxl/result/qc/ *.fq.gz

###去除低质量的碱基和接头序列
for f in $(ls *.fastq.gz | sed -e 's/_1.fastq.gz//' -e 's/_2.fastq.gz//' | sort -u)
do
        time fastp -i ${f}_1.fastq.gz -I ${f}_2.fastq.gz -o /home/hxl/sequence/clean/${f}_clean_1.fasq.gz -O /home/hxl/sequence/clean/${f}_clean_2.fasq.gz
done


############################### 使用metaphlan3进行物种注释 ###############################
for file in $(ls *.fastq.gz | sed -e 's/_1.fastq.gz//' -e 's/_2.fastq.gz//' | sort -u)
do
        time metaphlan ${file}_1.fastq.gz,${file}_2.fastq.gz --bowtie2out /home/hxl/result/tmp/${file}.bt2out --input_type fastq
        metaphlan /home/hxl/result/tmp/${file}.bt2out --input_type bowtie2out > /home/hxl/result/tmp/${file}_profile.txt
done

#提取物种水平的丰度信息
grep -E "s__|clade" merged_abundance_table.txt | sed 's/^.*s__//g' | cut -f1,3-12 | sed -e 's/clade_name/sample/g' > merged_abundance_table_species.txt


###################################### 序列组装 ###########################################
###使用megahit进行序列组装
for file in $(ls *.fastq.gz | sed -e 's/_clean_1.fastq.gz//' -e 's/_clean_2.fastq.gz//' | sort -u)
do
        time megahit -t 14 -m 0.9 --k-min 29 --min-contig-len 1000 -1 ${file}_clean_1.fastq.gz -2 ${file}_clean_2.fastq.gz --out-dir ~/result/Spanish/Spanish_contig/${file} --out-prefix ${file}
done

###Prodigal预测ORF
for file in $(ls *.fa | sed -e 's/_clean.contigs.fa//' | sort -u); do time prodigal -i ${file}_clean.contigs.fa -d ~/work_space2/prodigal_test/${file}.fa -o ~/work_space2/prodigal_test/${file}.gff -p meta -f gff > ~/work_space2/prodigal_test/${file}.log; done

###CD-hit对基因进行聚类
for file in $(ls *.fa | sed -e 's/.fa//' | sort -u); do time cd-hit-est -i ${file}.fa -o ~/work_space2/${file}_protein.fa -aS 0.9 -c 0.95 -G 0 -g 0 -T 0 -M 0; done


###########################################################################################
##########                      genne quantify                              ###############
###########################################################################################

### bwa
for file in $(ls *.fa | sed -e 's/.fa//' | sort -u)
do
        time bwa index -p ${file} ${file}.fa
        time bwa mem -t 16 ${file} /media/hxl/My_Book/China_clean/${file}_clean_1.fastq.gz /media/hxl/My_Book/China_clean/${file}_clean_2.fastq.gz | samtools sort --threads 15 -O bam -o ~/work_space6/bwa_result/${file}.bam
        rm ${file}.amb
        rm ${file}.ann
        samtools index ~/work_space6/bwa_result/${file}.bam
        samtools idxstats ~/work_space6/bwa_result/${file}.bam -@ 16 > ~/work_space6/gene_count/${file}.txt
        rm ${file}.bw
        rm ${file}.pac
        rm ${file}.sa
        rm ~/work_space6/bwa_result/${file}.bam
done

#### extract mapped rate from ${file}_stat.txt
for file in $(ls *.txt | sed -e 's/_stat.txt//' | sort -u)
do
        sed -n '/+ 0 mapped/ p' ${file}_stat.txt | sed 's/ /\t/' | cut -f2 | sed 's/+ 0 mapped (//' | sed 's/:/\t/' | cut -f1 > ${file}_mapped.tmp
        sed -i '1 i '${file}' ' ${file}_mapped.tmp
        cat ${file}_mapped.tmp | tr "\n" " " | sed 's/$/\n/' > ${file}.mapped
        rm ${file}_mapped.tmp
done

###########################################################################################
##########                      EGGNOG-mapper 功能注释                              #######
###########################################################################################

###将聚类后的序列翻译为氨基酸序列
for file in $(ls *.fa | sed -e 's/.fa//' | sort -u); do transeq -sequence ${file}.fa -outseq ~/work_space2/protein/${file}_protein.fa -trim Y; done
#调整蛋白质序列ID
for file in $(ls *.fa | sed -e 's/.fa//' | sort -u); do sed -i 's/_1 / /' ${file}.fa; done

for file in $(ls *.fa | sed -e 's/_protein.fa//' | sort -u); do time emapper.py --no_annot --no_file_comments --override --data_dir ~/database/eggnog_db/ -i ${file}_protein.fa --cpu 14 -m diamond -o ~/work_space2/tmp/${file}; done

for file in $(ls *.seed_orthologs | sed -e 's/.emapper.seed_orthologs//' | sort -u); do time emapper.py --annotate_hits_table ${file}.emapper.seed_orthologs --data_dir ~/database/eggnog_db/ --cpu 14 --no_file_comments --override -o ~/work_space2/emapper_result/${file}; done


########################### KEGG quantify
#####extract KO massage
for file in $(ls *.annotations | sed -e 's/.emapper.annotations//' | sort -u)
do
        cat ${file}.emapper.annotations | cut -f 1,12 > ~/work_space6/kos/${file}.tmp
        sed 's/,/\t/g' ~/work_space6/kos/${file}.tmp | sed '/-/ d' | sed 's/_1\t/\t/g' | cut -f 1,2 > ~/work_space6/kos/${file}_kos.txt
        rm ~/work_space6/kos/${file}.tmp
done

#####merge KO table and gene quantification table.
for file in $(ls *.txt | sed -e 's/.txt//' | sort -u)
do
        awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next}NR>FNR{if($1 in a)print a[$1],a[$2],$3}' ~/work_space6/${file}.annotations.txt ${file}.txt | sed '/-/ d' > ${file}.merge.txt
done
#${file}.annotations.txt ，kegg注释文件
#${file}.txt，非冗余基因丰度表

####对重复的行值求和
for file in $(ls *.txt | sed -e 's/.txt//' | sort -u)
do
        cat ${file}.txt | cut -f 2,4 > ${file}.tmp
        awk 'BEGIN{OFS="\t"} { seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' ${file}.tmp > ${file}.counts
        rm ${file}.tmp
done

#################################### Abricate 功能注释 ####################################
#vfdb
for file in $(ls *.fa | sed 's/.fa//' | sort -u); do time abricate --db vfdb --quiet --threads 14 ${file}.fa > ${file}.txt; done
#CARD
for file in $(ls *.fa | sed 's/.fa//' | sort -u); do time abricate --db card --quiet --threads 14 ${file}.fa > ${file}.txt; done


###########################################################################################
##########                               FishTaco                                   #######
###########################################################################################

### FishTaco with no inference
run_fishtaco.py -ta ./examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab -fu ./examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab -l ./examples/SAMPLE_vs_CLASS.tab -gc ./examples/METAPHLAN_taxa_vs_KO_only_K00001.tab -op fishtaco_out_no_inf -map_function_level none -functional_profile_already_corrected_with_musicc -assessment single_taxa -log

### FishTAco with prior-based inference
run_fishtaco.py -ta ./examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab -fu ./examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab -l ./examples/SAMPLE_vs_CLASS.tab -gc ./examples/METAPHLAN_taxa_vs_KO_only_K00001.tab -op fishtaco_out_prior_based_inf -map_function_level none -functional_profile_already_corrected_with_musicc -inf -assessment single_taxa -log

### FishTaco with de novo inference
#multi_taxa
time run_fishtaco.py -ta ./data/Diff_bacteria_abundance.txt -fu ./data/Diff_kos_abundance.txt -l ./data/Sample_Class.txt -op fishtaco_out_de_novo_inf -map_function_level none -functional_profile_already_corrected_with_musicc -inf -assessment multi_taxa -log


###########################################################################################
##########                             Virome Analysis                            #########
###########################################################################################

#Run Virsorter
for file in $(ls *.fa | sed -e 's/.fa//' | sort -u); do time virsorter run -w ~/result/data4/virsort_result/${file} -i ${file}.fa --include-groups "dsDNAphage,ssDNA" --min-length 5000 -j 8 all; done

#CheckV
#cluster to vOTU
for file in $(ls *.fa | sed -e 's/.fa//' | sort -u); do time cd-hit-est -i ${file}.fa -o ~/work_space5/${file}_cluster.fa -aS 0.8 -c 0.95 -T 0 -M 0; done

################ reads mapping ############
#build bowtie2 reference database
time bowtie2-build all-7_vir_cluster.fa meta-vir-db
#mapping
for file in $(ls *.fastq.gz | sed -e 's/_1.fastq.gz//' -e 's/_2.fastq.gz//' | sort -u); do time bowtie2 -p 16 -N 0 --no-unal -x ~/database/vir_bow_db/meta-vir-db -1 ${file}_1.fastq.gz -2 ${file}_2.fastq.gz -S /home/hxl/work_space5/tmp3/${file}.sam; done

#output bam files derectly
for file in $(ls *.fastq.gz | sed -e 's/_1.fastq.gz//' -e 's/_2.fastq.gz//' | sort -u); do time bowtie2 -p 15 -N 0 --no-unal -x ~/database/vir_bow_db/Virome_50_db -1 ${file}_1.fastq.gz -2 ${file}_2.fastq.gz | samtools sort --threads 15 -O bam -o ~/result/China/China_sort_bam/${file}.sort.bam; done

###count reads
#1.change format
for file in $(ls *.sam | sed -e 's/_clean.sam//' | sort -u); do time samtools view --threads 16 -b -S ${file}_clean.sam > ~/result/China/China_bam/${file}.bam; done
#2.sort
samtools sort --threads 12 ERR260132.sra_clean.bam > ERR260132.sra_clean.sort.bam
#3.index
samtools index file.sort.bam
for file in $(ls *.bam | sort -u); do samtools index ${file}; done

#4.merge (not necessary)
samtools merge mergedfile.bam *.bam

#5.stats, mapped reads
for file in $(ls *.bam | sed -e 's/.bam//' | sort -u); do samtools idxstats -@ 16 ${file}.bam -@ 16 > ./${file}.txt; done


























