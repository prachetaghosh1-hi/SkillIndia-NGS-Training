# **SRA Data Retrieval and Quality Assessment**

This file demonstrates how to:

1. Download SRA files from the NCBI Sequence Read Archive (SRA).  
2. Convert the SRA files into FASTQ format using SRA-Tools.  
3. Organize raw FASTQ files for analysis.  
4. Perform quality control:

   - For **short-read sequencing data** (Illumina), we use **FastQC**.  
   - For **long-read sequencing data** (Nanopore), we would use **NanoPlot / NanoComp** (or trimming with **NanoChop**).

5. Generate a MultiQC report to summarize FastQC results.

---

## short-read sequencing data of *Pseudomonas aeruginosa* with SRA ID: SRR35114395.
### Step-by-Step Workflow

**Step 1: Set Up Environment**
- Create a new mamba environment and install required tools: **SRA-Tools**, **FastQC**, and **MultiQC**.  

<pre>mamba create -n sra_tools_fastqc_multiqc</pre>
<pre> mamba activate sra_tools_fastqc_multiqc</pre>
<pre> mamba install -c bioconda sra-tools fastqc multiqc</pre>

---

**Step 2: Make a new directory for the organism**

<pre>mkdir pseudomonas</pre>
<pre> cd pseudomonas</pre>

**Step 3: Identify the SRA ID**

- Go to the **NCBI website** → **SRA database**.  
- Search for the organism **“Pseudomonas aeruginosa”**.  
- From the results, choose the appropriate dataset and copy the **SRA ID** (e.g., `SRR35114395`).

**Step 4: Download SRA File from NCBI using `prefetch`**

<pre>prefetch SRR35114395</pre>
This will create a folder named SRR35114395 containing the .sra file.

**Step 5: Convert SRA File to FASTQ**

*Note*: Check sequencing type on NCBI
-Before converting to FASTQ, always check on the NCBI SRA page whether the sequence is **paired-end** or **single-end**.
--Paired-end: use the `--split-files` option in fasterq-dump.
--Single-end: no need to add `--split-files`.

**For paired-end sequencing**
<pre>fasterq-dump --split-files SRR35114395/SRR35114395.sra</pre>

**For single-end sequencing, simply run**
<pre>fasterq-dump SRR35114395/SRR35114395.sra</pre>

This produces two FASTQ files for paired-end:

SRR35114395_1.fastq (forward reads)
SRR35114395_2.fastq (reverse reads)

For single-end, only one FASTQ file is generated.

**Step 6:  Organize FASTQ Files**

Create a directory `rawdata_pa` to store raw FASTQ files and move them there:
<pre>mkdir rawdata_pa</pre>
<pre>mv *.fastq rawdata_pa/</pre>

**Step 7:Run FastQC**

Create a directory to store FastQC results:
<pre>mkdir -p QC_pa</pre>

Run FastQC on the FASTQ files:
<pre>fastqc -o QC_pa rawdata_pa/SRR35114395_1.fastq rawdata_pa/SRR35114395_2.fastq </pre>


*Note: Using -o QC_pa ensures FastQC outputs go into the QC_pa directory*

**Step 8: Generate MultiQC Report**

Aggregate all FastQC results into a single report:
<pre>multiqc QC_pa</pre>

The report (multiqc_report.html) will be saved in the pseudomonas directory.

# Outputs

Raw FASTQ files: `rawdata_pa/`
FastQC reports: `QC_pa/`
MultiQC report: `pseudomonas/multiqc_report.html`




