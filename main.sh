path=/path/to/data    #change to the real path of raw data.
host_genome_path=/path/to/host_genome/      #change to the real path of host genome (.fa or .fasta).
emapper_path=/path/to/emapper_db   #change to the real path of eggnog-mapper database. 
mmseq_gtdb_db=path/to/mmseqs/gtdb/database  #change to the real path of mmseq gtdb database. 


##Trim data and remove adaptor using fastp
conda activate fastp
mkdir -p $path/trimmed
for file in *1.fastq.gz
do
name2=${file%1.fastq.gz}2.fastq.gz
name=$(echo $file | cut -d'_' -f2)
fastp \
-i $path/$file \
-I $path/$name2 \
-o $path/${name}_trimmed_R1.fastq.gz \
-O $path/${name}_trimmed_R2.fastq.gz \
--thread 10 \
-h $path/trimmed/${name}_report.html \
-j $path/trimmed/${name}_report.json
done

## Host removal using bowtie2
conda activate bowtie2
# Build host index
bowtie2-build $host_genome_path/hg38-mm10.fa $host_genome_path/hg38_mm10_index
# align
mkdir -p $path/host_removed 
cd $path/trimmed    #use trimmed data as input file
for file in *R1.fastq.gz
do
name2=${file%R1.fastq.gz}R2.fastq.gz
name=${file%_trimmed_R1.fastq.gz}
bowtie2 -x $host_genome_path/hg38_mm10_index  -1 $file -2 $name2 -S $path/host_removed/${name}.sam --un-conc $path/host_removed/${name}_host_removed -p 50
done


##Rename host removed fastq files


##Assemble using megahit
conda activate megtahit
mkdir -p $path/megahit_output
cd $path/host_removed 
for file in *1.fastq
do
megahit -1 $file -2 ${file%1.fastq}2.fastq -o $path/megahit_output/${file%_host_removed_1.fastq}_megahit_output -t 50
done

##QC of contigs using quast
conda activate quast
cd $path/megahit_output
for dir in *megahit_output
do
cd $path/megahit_output/$dir
quast -o ${dir%_megahit_output}_quast_output final.contigs.fa
done

##gene prediction using prodigal
conda activate prodigal
mkdir -p $path/prodigal_output
cd $path/megahit_output
for dir in *megahit_output
do
mkdir -p $path/prodigal_output/${dir%_megahit_output}
prodigal -i $path/megahit_output/$dir/final.contigs.fa -o $path/prodigal_output/${dir%_megahit_output}.gff -f gff \
-a $path/prodigal_output/${dir%_megahit_output}.faa -d $path/prodigal_output/${dir%_megahit_output}.fna -p meta
done



##gene counts using feature counts 
conda activate subread
mkdir -p $path/gene_counts
cd $path/prodigal_output/
for dir in *
mkdir -p $path/gene_counts/$dir 
do
bowtie2-build $path/megahit_output/${dir}_megahit_output/final.contigs.fa  $path/gene_counts/$dir/${dir}.fa_index
bowtie2 -x $path/gene_counts/$dir/${dir}.fa_index -1 $path/host_removed/${dir}_host_removed_1.fastq -2 $path/host_removed/${dir}_host_removed_2.fastq -S $path/gene_counts/$dir/output.sam -p 50
samtools view -bS $path/gene_counts/$dir/output.sam > $path/gene_counts/$dir/output.bam
samtools sort $path/gene_counts/$dir/output.bam -o $path/gene_counts/$dir/sorted_output.bam
samtools index $path/gene_counts/$dir/sorted_output.bam
featureCounts -a $path/prodigal_output/${dir}.gff -o $path/gene_counts/$dir/counts.txt -t CDS -g ID -p -B $path/gene_counts/$dir/sorted_output.bam
done


##gene functional annotation using eggnog-mapper
conda activate eggnog-mapper
mkdir -p $path/emapper_output
cd $path/prodigal_output
for dir in *
do 
emapper.py -i  $path/prodigal_output/${dir}.fna  -o $path/emapper_output/$dir --cpu 50   --data_dir $emapper_path
done


##gene taxonomy classification using mmseqs2+gtdb database
conda activate mmseqs
mkdir -p $path/taxonomy_contig_level
cd $path/prodigal_output
for dir in *
do
mkdir -p $path/taxonomy_contig_level/$dir 
mmseqs easy-taxonomy $path/prodigal_output/$dir/${dir}.faa  $mmseq_gtdb_db  $path/taxonomy_contig_level/$dir/${dir}_taxonomy $path/taxonomy_contig_level/$dir/tmp --tax-lineage 1  --threads 50
done


## binning using metawrap
conda activate metawrap 
mkdir $path/metawrap_binning
cd $path/megahit_output
for dir in *_megahit_output
do
metawrap binning -o $path/metawrap_binning/${dir%_megahit_output} -t 50 -a $path/$dir/final.contigs.fa \
--metabat2 --maxbin2 --concoct \
 $path/host_removed/${dir%_megahit_output}_host_removed_1.fastq  $path/host_removed/${dir%_megahit_output}_host_removed_2.fastq
done 


## bin refinement using metawrap
cd $path/metawrap_binning
for dir in *
do cd $path/metawrap_binning/$dir
metawrap bin_refinement -o $path/metawrap_binning/$dir/refinement_output -t 50 \
-A $path/metawrap_binning/$dir/metabat2_bins/ \
-B $path/metawrap_binning/$dir/maxbin2_bins/ \
-C $path/metawrap_binning/$dir/concoct_bins/ \
-c 70 -x 10
done


## quant bins using metawrap 
cd $path/metawrap_binning
for dir in *
do
cd $path/metawrap_binning/$dir/
mkdir -p quant_bins
metaWRAP quant_bins -b $path/metawrap_binning/$dir/refinement_output/metawrap_70_10_bins -o $path/metawrap_binning/$dir/quant_bins \
-a $path/megahit_output/${dir}_megahit_output/final.contigs.fa \
$path/host_removed/${dir%}_host_removed_1.fastq  $path/host_removed/${dir%}_host_removed_2.fastq
done


## bin taxonomy classification using gtdbtk
conda activate gtdbtk
mkdir -p $path/metawrap_binning/gtdbtk
cd $path/metawrap_binning
for dir in *
do
mkdir -p $path/metawrap_binning/gtdbtk/$dir
gtdbtk classify_wf --genome_dir $path/metawrap_binning/$dir/refinement_output/metawrap_70_10_bins  \
--out_dir $path/metawrap_binning/gtdbtk/$dir   --extension fa   --skip_ani_screen  -cpus 50
done



## bin functional annotation using prokka
conda activate prokka
mkdir -p $path/metawrap_binning/prokka
cd $path/metawrap_binning
for dir in *
do
mkdir -p $path/metawrap_binning/prokka/$dir
cd $path/metawrap_binning/$dir/refinement_output/metawrap_70_10_bins
  for file in *.fa
  do prokka --outdir $path/metawrap_binning/prokka/$dir $file --cpus 50
  done
done


conda deactivate
