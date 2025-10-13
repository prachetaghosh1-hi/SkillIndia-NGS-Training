# VARIANT CALLING  

## 🔹 Introduction  

**Variant calling** is the process of identifying genetic variants (such as SNPs and indels) by comparing sequencing reads from a sample to a reference genome.  

In this process, we require two main inputs:  
1. **Trimmed FASTQ files** → obtained after quality control and trimming.  
2. **Reference `.fna` file** → downloaded from the NCBI dataset for the target organism.  

> ⚠️ *Both these steps (trimming and reference download) were already discussed in previous files.*

---

## 🧫 Organism and Setup  

For this file, we used the organism **_Staphylococcus aureus_** with the **SRA ID: SRR35532843**.  

We created a main directory named **`staph`**, and inside it, another directory called **`variant_calling`**, where all data and outputs related to variant calling will be stored.  

Important tools required for variant calling were installed inside a **new environment** named **`snp`**.

---

## 🧰 Tools Installed in Environment `snp`

<pre>
mamba create -n snp
mamba activate snp
mamba install -c bioconda bwa samtools picard freebayes bcftools -y</pre>

---

*Purpose of tools*:

1. bwa → Align reads to reference genome.
2. samtools → File conversion, sorting, filtering, and indexing.
3. picard → Mark duplicate reads.
4. freebayes → Perform variant calling.
5. bcftools → Process and count variants.

---

## ⚙️ VARIANT CALLING WORKFLOW

***Step 1: Get trimmed fastq files***

***Step 2:Get reference fna file from NCBI***

***Step 3:Indexing of Reference Genome***

Before alignment, we must index the reference genome so that BWA can efficiently map reads.

<pre>bwa index GCF_000013425.1_ASM1342v1_genomic.fna</pre>

### 📝 Explanation:
This command creates several index files (.amb, .ann, .bwt, .pac, .sa) that allow the BWA aligner to quickly find matching positions for reads in the reference genome.

***Step 4: Alignment of Trimmed Reads to the Reference***

Use bwa mem to align the paired-end trimmed reads to the indexed reference genome.

<pre>bwa mem -t 4 -R '@RG\tSM:skill\tID:skill\tPL:illumina' GCF_000013425.1_ASM1342v1_genomic.fna ../trim_sa/SRR35532843_1_trimmed.fastq ../trim_sa/SRR35532843_2_trimmed.fastq > map.sam</pre>

### 📝 Explanation:
|flags | meaning |
|----|----|
| -t 4 | use 4 threads |
| -R | adds read group info (sample ID, platform) |

The output file map.sam contains alignment information for all reads.

***Step 5: Convert SAM to BAM***

SAM files are very large and text-based; converting to BAM reduces size and speeds up processing.

<pre>
samtools view -bS -@ 4 map.sam > map.bam</pre>

***Step 6: Filter Unmapped Reads***

Remove reads that did not align to the reference genome


<pre>samtools view -F 4 -@ 4 -o filtered.bam map.bam</pre>

### 📝 Explanation:

-F 4 excludes unmapped reads.

The output filtered.bam contains only mapped reads.

***Step 7: Sort BAM File***

Sorting arranges reads by their genomic coordinates, which is required for most downstream tools.

<pre>
samtools sort -@ 4 -o sorted.bam filtered.bam</pre>

***Step 8: Mark Duplicates***

During sequencing and PCR amplification, duplicate reads can appear.
picard marks such duplicates so that they are not counted multiple times in variant calling.

<pre>picard MarkDuplicates -I sorted.bam -O markdup.bam -M metrics.txt --CREATE_INDEX true</pre>

### 📝 Explanation:

markdup.bam → BAM file with duplicates flagged.

metrics.txt → report summarizing duplicate statistics.

The --CREATE_INDEX true option generates an index file (.bai).

***Step 9: Index the Reference FASTA File***

For variant calling, a FASTA index (.fai) of the reference is required.

<pre> samtools faidx GCF_000013425.1_ASM1342v1_genomic.fna</pre>

### 📝 Explanation:
This creates a .fai file allowing random access to specific genome regions during variant calling.

***Step 10: Perform Variant Calling (FreeBayes)***
Use freebayes to identify SNPs and Indels.

<pre>freebayes-parallel <(fasta_generate_regions.py GCF_000013425.1_ASM1342v1_genomic.fna.fai 100000) 4 -p 1 -f GCF_000013425.1_ASM1342v1_genomic.fna markdup.bam > variants.vcf</pre>

### 📝 Explanation:

The genome is divided into 100 kb regions for parallel processing.

The result is stored in variants.vcf, containing all identified variants.

***Step 11: Count the Total Number of Variants***

<pre>bcftools view variants.vcf | grep -v '^#' | wc -l</pre>

###📝 Explanation:
This counts the total number of variant entries in the VCF file, excluding header lines (which start with #).

**Summary**
After completing these steps:

The reference genome is indexed and aligned with trimmed reads.

A variant file (variants.vcf) is generated containing SNPs and indels.

The next step (in the following file) is variant annotation using snpEff.

📂 Final Output Files in variant_calling/
| File | Description |
|----|----|
| `map.sam` |	SAM file containing raw alignments |
| `map.bam`	| BAM format of SAM file |
| `filtered.bam` | Filtered BAM with mapped reads |
| `sorted.bam` | Sorted BAM by genomic coordinates |
| `markdup.bam` |	BAM file with duplicates marked |
| `metrics.txt` |	Duplicate read summary |
| `GCF_000013425.1_ASM1342v1_genomic.fna.fai` |	FASTA index file |
| `variants.vcf` |	Detected variants (SNPs + Indels) |

**Conclusion:**
This file documents the variant calling stage of the Staphylococcus aureus workflow.
Starting from the indexed reference genome and trimmed FASTQ files, the process aligns reads, processes BAM files, and produces *a final variant list ready for annotation*.

---

