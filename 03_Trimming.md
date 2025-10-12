# Trimming 

## **Why We Do Trimming**
Trimming is an important preprocessing step for sequencing data. It removes:  
- Low-quality bases from the ends of reads  
- Adapter sequences  
- Poor-quality regions that may affect downstream analysis  

This ensures that subsequent analyses (like mapping, assembly, or variant calling) are more accurate and reliable.

---

**Note:**  
- For **short-read sequencing data** (Illumina), trimming can be done using **Trim Galore** or **fastp**. In this workflow, **fastp** is used.  
- For **long-read sequencing data** (Nanopore), trimming is typically done with **Porechop**.

---


## for SHORT READS
**Step 1: Create Environment and Install fastp**

- Create a new conda environment called `trimenv` and install **fastp**:
<pre>mamba create -n fastp </pre>
<pre>mamba activate fastp</pre>
<pre>mmaba install -c bioconda fastp</pre>

**Step 2:Create Directory for Trimmed Files**

Make a directory to store trimmed FASTQ files:
<pre>mkdir -p trim_pa</pre>

**Step 3: Run fastp on Paired-End FASTQ Files**

Run fastp on raw FASTQ files from rawdata_pa/ and save output to trim_pa/:

<pre>fastp -i rawdata_pa/SRR35114395_1.fastq -I rawdata_pa/SRR35114395_2.fastq -o trim_pa/SRR35114395_1_trimmed.fastq -O trim_pa/SRR35114395_2_trimmed.fastq --cut_front --cut_tail --cut_window_size 4 --thread 4 --detect_adapter_for_pe --qualified_quality_phred 20</pre>

###Explanation of key options:

`--cut_front / --cut_tail` → trims low-quality bases at the start/end of reads

`--cut_window_size 4` → sliding window of 4 bases for trimming

`--thread 4` → uses 4 CPU threads

`--detect_adapter_for_pe` → detects adapters automatically for paired-end reads

`--qualified_quality_phred 20` → trims bases below quality score Q20

***output:***
The trimmed files generated:

SRR35114395_1_trimmed.fastq (forward)

SRR35114395_2_trimmed.fastq (reverse)

**Step 4: Quality Check of Trimmed Files**

After trimming, it is important to check whether trimming improved read quality.
Run FastQC and MultiQC on trimmed files:
<pre>mkdir -p QC_trim_pa
fastqc -o QC_trim_pa trim_pa/SRR35114395_1_trimmed.fastq trim_pa/SRR35114395_2_trimmed.fastq
multiqc QC_trim_pa</pre>


***output:*** The MultiQC report summarizes the quality improvements after trimming.

---
