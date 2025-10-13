# Trimming 

## **Why We Do Trimming**
Trimming is an important preprocessing step for sequencing data. It removes:  
- Low-quality bases from the ends of reads  
- Adapter sequences  
- Poor-quality regions that may affect downstream analysis  

This ensures that subsequent analyses (like mapping, assembly, or variant calling) are more accurate and reliable.

---

*Prerequisites* Two fastq files from NCBI (discussed in last section)(in case the sequences were paired)

---

**Note:**  

- For **short-read sequencing data** (Illumina), trimming can be done using **Trim Galore** or **fastp**. In this workflow, **fastp** is used.  
- For **long-read sequencing data** (Nanopore), trimming is typically done with **Porechop** and QC for trimmed files is done with **NanoFilt**.

---


## for SHORT READS
**Step 1: Create Environment and Install fastp**

- Create a new conda environment called `trimenv` and install **fastp**:
<pre>mamba create -n fastp </pre>
<pre>mamba activate fastp</pre>
<pre>mmaba install -c bioconda fastp</pre>

**Step 2:Create Directory for Trimmed Files**

Make a directory to store trimmed FASTQ files:
<pre>mkdir -p trim_ab</pre>

**Step 3: Run fastp on Paired-End FASTQ Files**

Run fastp on raw FASTQ files from rawdata_pa/ and save output to trim_pa/:

<pre>fastp -i rawdata_ab/SRR35227129_1.fastq -I rawdata_ab/SRR35227129_2.fastq -o trim_ab/SRR35227129_1_trimmed.fastq -O trim_ab/SRR35227129_2_trimmed.fastq --cut_front --cut_tail --cut_window_size 4 --thread 4 --detect_adapter_for_pe --qualified_quality_phred 20</pre>

###Explanation of key options:

`--cut_front / --cut_tail` → trims low-quality bases at the start/end of reads

`--cut_window_size 4` → sliding window of 4 bases for trimming

`--thread 4` → uses 4 CPU threads

`--detect_adapter_for_pe` → detects adapters automatically for paired-end reads

`--qualified_quality_phred 20` → trims bases below quality score Q20

***output:***
The trimmed files generated:

SRR35227129_1_trimmed.fastq (forward)

SRR35227129_2_trimmed.fastq (reverse)

**Step 4: Quality Check of Trimmed Files**

After trimming, it is important to check whether trimming improved read quality.
Run FastQC and MultiQC on trimmed files:
<pre>mkdir QC_trim_ab
fastqc -o QC_trim_ab trim_pa/SRR35227129_1_trimmed.fastq trim_ab/SRR35227129_2_trimmed.fastq
multiqc QC_trim_ab</pre>


***output:*** The MultiQC report summarizes the quality improvements after trimming.

---

# LONG READS

organism used : *Mycobacterium tuberculosis* with SRA ID:SRR35504420 , Layout : Single

similar to short reads, the workflow is similar just the tools for trimming and its QC check is different

<pre> mamba create -n porechop
mamba activate porechop
mamba install -c bioconda porechop=0.2.4
porechop -i SRR.fastq -o  SRR35504420_trimmed.fastq --format fastq -t 2 –discard_middle
mamba deactivate
mamba create -n nanofilt
mamba install -c bioconda nanofilt -y
mamba activate nanofilt
cat SRR_trimmed.fastq | NanoFilt -l 500 -q 12  > SRR.filtered.fastq</pre>

---

*The only difference I observed in the trimming step of short and long reads is that in long reads only the middle adapters are discarded and long reads trimming takes alot of time.*

---

