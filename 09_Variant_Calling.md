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

## 🧬 VARIANT ANNOTATION USING snpEff

Once we have our variant file (`variants.vcf`), the next step is to **annotate** these variants to understand their biological significance — which genes they affect, whether they cause synonymous or nonsynonymous mutations, and where they occur in the genome.

We use **snpEff** (version: `snpEff_v4_3t_core`) for this purpose.  
snpEff uses a reference genome and its corresponding GFF annotation file to predict the effect of each variant.

---

***Step 12: Download snpEff***

Download snpEff (version: `snpEff_v4_3t_core`) from Google and move it to the working directory, `variant_calling`
<pre>mv /mnt/c/Users/Lenovo/Downloads/snpEff_v4_3t_core.zip </pre>

Activate the environment having sudo tool and unzip the file obtained from snpEff.
<pre> unzip snpEff_v4_3t_core.zip</pre>

Make several directories at once `variant_calling` --> `snpEff` --> `data` --> `staph`
<pre>cd snpEff
mkdir data
cd data
mkdir staph
cd staph</pre>

### 📝 Explanation:
We created a folder structure inside snpEff/data/ for the organism (Staphylococcus aureus) that will hold genome sequence (.fna) and annotation (.gff3) files.

***Step 13: Download Reference Files for Annotation***

We need both:

1. The genome sequence file (.fna)

2. The annotation file (.gff3)

Download them again using NCBI Datasets(don't forget to activate the environment hwere this toolis present):
<pre> datasets download genome accession GCF_000013425.1 --include genome,gff3
unzip ncbi_dataset.zip</pre>

since this reference fille is too far so , move this ref file to the `variant_calling` directory:
<pre>mv ncbi_dataset/data/GCF_000013425.1/* .</pre>


***Step 14: Rename Files for snpEff Compatibility***

snpEff requires standardized filenames to recognize genome and feature files properly.

<pre>mv GCF_000013425.1_ASM1342v1_genomic.fna sequences.fa
mv genomic.gff genes.gff3</pre>

  
### 📝 Explanation:

sequences.fa → the full genomic sequence (FASTA format)

genes.gff3 → gene annotation file describing genomic features (genes, exons, UTRs, etc.)

***Step 15: Configure snpEff***

snpEff uses a configuration file (snpEff.config) to locate available genomes and their data directories.
Since the configuration file is present 2 directories back in `snpEff` directory
<pre> cd ../../</pre>

Open the file and add a new entry:
*very important step*
<pre>nano snpEff.config</pre>

Scroll to the end (or anywhere safe in the file) and add:-  staph.genome : staph

### 📝 Explanation:

staph.genome → ID of the genome.

staph → directory name under snpEff/data/ that contains sequences.fa and genes.gff3.

⚠️ Important: Do not delete or modify any existing entries in snpEff.config.

***Step 16: Build snpEff Database***

Now we build the annotation database for Staphylococcus aureus using the provided GFF3 file.

<pre>java -jar snpEff.jar build staph -gff3</pre>

### 📝 Explanation:
This step parses genes.gff3 and sequences.fa, creating a local snpEff database for annotation.

***Step 17: Perform Variant Annotation***

Annotate the variants from the VCF file generated in the previous steps.
<pre>java -jar snpEff.jar ann staph ../variants.vcf > variants.ann.vcf</pre>

You can open the annotated VCF file to inspect the added information:
<pre> nano variants.ann.vcf </pre>

🧩 Compression and Indexing of Annotated File:
<pre>bgzip variants.ann.vcf

gunzip variants.ann.vcf.gz</pre>

### Explanation:

bgzip → compresses the VCF file.

tabix → creates an index for fast retrieval.

gunzip → decompresses when needed for manual inspection.

***Step 18: Count SNPs and Indels in Annotated File***

Count the total number of SNPs and indels detected in the annotated variant file.

<pre>bcftools view -v snps variants.ann.vcf | grep -v "^#" | wc -l
bcftools view -v indels variants.ann.vcf | grep -v "^#" | wc -l</pre>

### 📝 Explanation:

The first command counts all SNPs (single nucleotide polymorphisms).

The second command counts all indels (insertions/deletions).

***Step 19: Optional – Filtration of Annotated Variants***

Variant filtration is optional and is generally performed when the mapping percentage is below 70%.
In our case, the mapping percentage was 92.33%, so filtration is not mandatory — but shown here for reference.

<pre>bcftools filter -i 'TYPE="snp" && (FMT/DP)>10 && (FMT/GQ)>20' --threads 10 -o filtered.vcf.gz variants.ann.vcf.gz
bcftools filter -i 'TYPE="snp"' --threads 4 -o filtered.vcf.gz variants.ann.vcf.gz
gunzip filtered.vcf.gz
nano filtered.vcf</pre>

### 📝 Explanation:

Filters variants based on read depth (DP) and genotype quality (GQ).

Produces a new filtered VCF file (filtered.vcf).

***Step 20: Count SNPs and Indels in Filtered File***

If filtration was performed, count the number of SNPs and indels again.

<pre>bcftools view -v snps filtered.vcf.gz | grep -v '^#' | wc -l
bcftools view -v indels filtered.vcf.gz | grep -v '^#' | wc -l</pre>

### 📝 Explanation:
This provides counts for high-quality SNPs and indels after filtering.

**📊 Final Outputs in variant_calling/**
|File	| Description |
|-----|----|
|variants.vcf |	Raw variants identified by FreeBayes |
| variants.ann.vcf | Annotated variants after snpEff |
| filtered.vcf	(Optional) | Filtered variants based on quality metrics |
| sequences.fa |	Reference genome sequence (renamed) |
| genes.gff3 | Gene annotation file for snpEff |
| snpEff.config |	Configuration file with genome path |
| metrics.txt |	Picard duplicate read summary |

**✅ Summary**
In this section, we annotated the variants detected in the Staphylococcus aureus sample (SRR35532843) using snpEff.
This step enriches the variant list with biological context — identifying the gene affected, its location, and the predicted impact of each mutation.


---


## 🧾 RESULTS

### **1. Alignment and Sorting**

After converting the SAM file to a **sorted BAM**, the **percentage of properly paired reads** increased from **90.86% to 98.41%**.  
This improvement occurs because **coordinate sorting** groups mate pairs together, enabling tools such as `samtools flagstat` to correctly recognize properly paired reads.  

It was also observed that the **number of primary mapped reads** became **100%** after the **filtering step**.  
This is due to the use of the `-F 4` parameter in `samtools view`, which removes all unmapped reads.  

🧩 **Interpretation:**  
In BAM/SAM files, each read can have multiple alignments — primary, secondary, or supplementary.  
A **primary alignment** represents the best, main location for that read, which is expected after proper filtering and sorting.

---

### **2. Duplicate Marking (Picard Metrics)**

Using **Picard MarkDuplicates**, duplicate reads generated during PCR amplification were identified.  
The **metrics.txt** summary reported a **duplication rate of 0.095163**, indicating that approximately **9.5% of the reads** were duplicates.  

The remaining **90.5% unique reads** were retained for variant calling.  
Removing duplicate reads is crucial to prevent **false-positive SNPs**, which can arise from PCR amplification artifacts.

---

### **3. Variant Detection (FreeBayes)**

Variants were detected using **FreeBayes**, which compared the aligned reads (`markdup.bam`) to the reference genome in small chunks (100 kb windows).  
The resulting `variants.vcf` file listed all variant sites found in the sample.  

🧠 **Note:**  
When viewing the file (`nano variants.vcf`), we can see the variant positions and reference/alternate alleles.  
However, this raw variant file does **not** include biological information about affected genes or mutation effects.

---

### **4. Total Variants Identified**

A total of **45,503 variant sites** were detected in the *Staphylococcus aureus* sample (**SRR35532843**).  

📊 **Interpretation:**  
A higher number of variants indicates that the sample genome differs significantly from the reference genome at many positions, suggesting possible strain-specific mutations or evolutionary divergence.

---

### **5. Variant Annotation (snpEff)**

Annotation was performed using **snpEff**, which interprets the biological meaning of each variant.  
It identifies whether variants occur within **protein-coding genes**, and specifies the **gene name**, **functional impact**, and **mutation type**.  

🧬 **Example Interpretation:**  
A variant located on **NC_00795.1 at position 137** showed a change from **A (reference)** to **G (sample)**.  
The annotation revealed that this variant is part of a **protein-coding gene (dnaA region)**, located at the **380th base position**, and classified as a **modifier variant**.  

Thus, snpEff provides insights into the functional context of each variant beyond simple positional data.

---

### **6. Variant Type Summary**

Among the total **45,503 variants**, snpEff reported:

- **37,452 SNPs (Single Nucleotide Polymorphisms)**  
- **589 Indels (Insertions and Deletions)**  

This distribution reflects that most genomic differences between the sample and reference are single-base substitutions.

---

### **7. Filtered Variant Results (Optional Step)**

Variant filtration was optional in this workflow, as the mapping percentage was **92.33%** (above the 70% threshold).  
However, for demonstration, filtration was performed to retain only high-confidence SNPs based on depth (`DP > 10`) and genotype quality (`GQ > 20`).


