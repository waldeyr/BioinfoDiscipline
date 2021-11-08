# :school: Práticas da Disciplina de Bioinformática (Bioinformatics Discipline Practices)

This is a repository for the Bioinformatics Discipline practices.

## :point_right: Setup the environment

#### How to install Conda?

* Download Anaconda software from  https://www.anaconda.com/products/individual
* Install it (https://docs.anaconda.com/anaconda/install)

#### How to set up the environment using conda?

##### Option 01:

Download the file [environment.yml](https://raw.githubusercontent.com/waldeyr/BioinfoDiscipline/main/environment.yml)

`conda env create -f environment.yml`

##### Option 02:

`conda create -n pipelines python=3.6 r=3.6 -y`

#### How to enter in the conda environment?

`conda activate pipelines`

#### How to setup the channels (repositories) with the needed tools?

```
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### How to install the needed tools into the DisciplinaBioinfo environment?

`conda install pandas numpy jupyterlab jupyter nano readline=6.2 sra-tools entrez-direct bwa fastqc fastp spades quast star htseq seqtk samtools bcftools r-xml freebayes bedtools vcflib rtg-tools matplotlib -y`

* R packages (run it from the R prompt):

`install.packages('IRkernel')`

`install.packages("BiocManager")`
    
`BiocManager::install(c("limma","edgeR","Glimma","data.table","org.Mm.eg.db", "statmod"))`


## :notebook_with_decorative_cover: Practice 01 - De novo assembly of a Brazillian isolate of Sars-Cov-2 Genome

### Enter in the conda environmet

`conda activate pipelines`


### Obtaining the raw material

* SARS-CoV-2 genome sequencing Rio Grande do Sul / Brazil, Dec 2020; Total RNA from SARS-CoV-2 positive samples was converted to cDNA. Viral whole-genome amplification was performed according to the Artic Network. Available in [SRR13510367](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13510367).

`fastq-dump --accession SRR13510367 --split-files --outdir rawdata -v`

### Generating a quality report for the sequenced reads

`fastqc rawdata/SRR13510367_1.fastq rawdata/SRR13510367_2.fastq`

### Filtering data with fastp

`mkdir filtered_data`

`fastp --thread 4 -p -q 30 -i rawdata/SRR13510367_1.fastq -I rawdata/SRR13510367_2.fastq -o filtered_data/SRR13510367_1_FILTERED.fastq -O filtered_data/SRR13510367_2_FILTERED.fastq &`

`mv fastp.* filtered_data/`

### Running the asembly with Spades

`spades.py -t 4 -o assembly --careful -k 21,33,55 -1 filtered_data/SRR13510367_1_FILTERED.fastq -2 filtered_data/SRR13510367_2_FILTERED.fastq &`

### generate a report about the assembly with Quast

* Download the Sars-Cov-2 reference genome

`esearch -db nucleotide -query "NC_045512.2" | efetch -format fasta > NC_045512.2.fasta`

* Performe the report

`quast assembly/scaffolds.fasta -R NC_045512.2.fasta`

## :notebook_with_decorative_cover: Practice 02 - [Transcriptome of human T cell stimulated with anti-CD3 antibody (Sousa, 2019)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112899)

Briefing: human peripheral blood mononuclear cells were purified from healthy volunteers blood and were cultured in the presence of the monoclonal antibody OKT3 or a recombinant fragment of humanized anti-CD3 (FvFcR) or recombinant fragment chimeric anti-CD3 (FvFcM).

* For the tutorial, we only will use the chromossome 22 and a monoclonal antibody OKT3 sample to make it feasible in a personal computer.


### Enter in the conda environmet

`conda activate pipelines`

### Obtaining the raw material

* [Homo sapiens reference genome](http://www.ensembl.org/info/data/ftp/index.html)
* [Chromossome 22 DNA sequence](http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz)
* [Chromossome 22 GFF3 annotation](http://ftp.ensembl.org/pub/release-104/gff3/homo_sapiens/Homo_sapiens.GRCh38.104.chromosome.22.gff3.gz)
* [Chromossome 22 GTF annotation](http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz)
* OKT3  replica 1 = SRR6974025
* FvFcR replica 1 = SRR6974027

### Downloading the reads

`fastq-dump --accession SRR6974025 --outdir rawdata -v`

`fastq-dump --accession SRR6974027 --outdir rawdata -v`

### Filtering the reads quality using fastp

`fastqc rawdata/SRR6974025.fastq`

`fastqc rawdata/SRR6974027.fastq`

`mkdir filtered_data`

`fastp --thread 8 -p -q 30 -i rawdata/SRR6974025.fastq -o filtered_data/SRR6974025_FILTERED.fastq`

`fastp --thread 8 -p -q 30 -i rawdata/SRR6974027.fastq -o filtered_data/SRR6974027_FILTERED.fastq`

`mv fastp.* filtered_data/`

### Generating an index for the genome reference

`mkdir genome_reference && cd genome_reference`

`wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz`

`wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz`

`gunzip Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz && gunzip Homo_sapiens.GRCh38.104.chr.gtf.gz && cd ../`

`STAR --runThreadN 8 --runMode genomeGenerate --genomeDir genome_reference --genomeFastaFiles genome_reference/Homo_sapiens.GRCh38.dna.chromosome.22.fa --sjdbGTFfile genome_reference/Homo_sapiens.GRCh38.104.chr.gtf --sjdbOverhang 49`

### Mapping the filtered reads to the genome using STAR

`STAR --genomeDir genome_reference --runThreadN 8 --readFilesIn filtered_data/SRR6974025_FILTERED.fastq filtered_data/SRR6974027_FILTERED.fastq --outFileNamePrefix mapped_data --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --quantMode GeneCounts TranscriptomeSAM`

### Counting mapped genes with htseq-count

`htseq-count --nonunique all --format bam --stranded no --order pos --type exon --idattr gene_id <FILE_NAME>.out.bam genome_reference/Homo_sapiens.GRCh38.104.chr.gtf > table.counts
`

## References
Sousa IG, Simi KCR, do Almo MM, Bezerra MAG et al. Gene expression profile of human T cells following a single stimulation of peripheral blood mononuclear cells with anti-CD3 antibodies. BMC Genomics 2019 Jul 19;20(1):593. PMID: 31324145

## :notebook_with_decorative_cover: Practice 02.1 - EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival (Fu, 2019)

Análise de genes diferencialmente expressos em células de glândulas mamárias de camundongos fêmeas em duas situações: grávidas e lactantes.
Neste tutorial, a partir de dados já mapeados e contados, é realizada a análise de gene diferencialmente expressos (DEGs).

[Link par o Tutorial](https://github.com/waldeyr/DisciplinaBioinfo/blob/main/rna-seq.ipynb)


## :notebook_with_decorative_cover: Practice 03 - [COVID-19 variant calling from Illumina data](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)

This practice is an adaptation of 

* Nielsen R, Paul JS, Albrechtsen A, Song YS. Genotype and SNP calling from next-generation sequencing data. Nat Rev Genetics, 2011, 12:433-451.

* Olsen ND et al. Best practices for evaluating single nucleotide variant calling methods for microbial genomics. Front. Genet., 2015, 6:235.


### Enter in the conda environmet

`conda activate pipelines`


### Obtaining the raw material

* SARS-CoV-2 genome sequencing Rio Grande do Sul / Brazil, Dec 2020; Total RNA from SARS-CoV-2 positive samples was converted to cDNA. Viral whole-genome amplification was performed according to the Artic Network. Available in [SRR13510367](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13510367).


`fastq-dump --accession SRR13510367 --outdir rawdata -v`

* Download the Sars-Cov-2 reference genome

`mkdir genome && esearch -db nucleotide -query "NC_045512.2" | efetch -format fasta > genome/NC_045512.2.fasta`


### Alignment with Map with BWA-MEM

* create the genome index

`bwa index genome/NC_045512.2.fasta`

* Map reads to the reference genome

`bwa mem -t 4 genome/NC_045512.2.fasta rawdata/SRR13510367.fastq > SRR13510367.sam`

* Convert the SAM file to BAM format

`samtools view -S -b SRR13510367.sam > SRR13510367.bam`

* Sort BAM file by coordinates

`samtools sort -o SRR13510367_sorted.bam SRR13510367.bam`

* create an sam indexes for the genome and the mapped reads (bam)

`samtools faidx genome/NC_045512.2.fasta`

`samtools index SRR13510367_sorted.bam`

* Some stats about yout alignment

`samtools flagstat SRR13510367_sorted.bam`


### Call variants

* Create a folder to save the variant resuls

`mkdir variants`

* generating a file with the variants

`freebayes -p 1 -f genome/NC_045512.2.fasta SRR13510367_sorted.bam > variants/SRR13510367.vcf`

or (using samtools)

`bcftools mpileup -Ou -f genome/NC_045512.2.fasta SRR13510367_sorted.bam | bcftools call -mv -o variants/SRR13510367.vcf`

* Take a look in the first lines

`cat variants/SRR13510367.vcf | grep -v '##' | head -4`


## Some stats

* Compressing the file

`bgzip variants/SRR13510367.vcf`

* Creating an index for the variants file

`tabix -p vcf variants/SRR13510367.vcf.gz`

* Quick stats

`rtg vcfstats variants/SRR13510367.vcf.gz`

* More detailed stats

`bcftools stats -F genome/NC_045512.2.fasta -s - variants/SRR13510367.vcf.gz > variants/SRR13510367.vcf.gz.stats`

* Plots

`mkdir variants/plots`

`plot-vcfstats -p variants/plots/ variants/SRR13510367.vcf.gz.stats`
