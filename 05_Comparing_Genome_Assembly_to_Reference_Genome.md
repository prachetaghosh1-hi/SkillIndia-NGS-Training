# Comparing Genome Assembly to Reference Genome

*Objective:Compare the genome assembly generated from Pesudomonas aeruginosa SRA file to the reference genome from NCBI using QUAST.*

---
**Why Use a Reference Genome?**

A reference genome provides a template for evaluating assembly quality.
Before scaffolding, it helps identify misassemblies and gaps in contigs.
Before annotation, it allows mapping of genes and genomic features accurately by comparison.
Overall, it ensures that the assembled genome is structurally and functionally consistent with known genomes.

---

***Step1:  Install NCBI Datasets Cli***

<pre>mamba create -n ncbidatasets_sudo
mamba install -c bioconda -y ncbi-datasets-cli
sudo apt install unzip
mamba activate ncbidatasets_sudo</pre>

***Step 2: Identify Reference Genome***
Search Pseudomonas aeruginosa under “Genomes” in NCBI. Select the reference genome (green tick) and note the accession number (GCF_…).
Set Up Reference Directory

<pre>mkdir ref
cd ref</pre>


***Step 3: Download Reference Genome***

<pre>datasets download genome accession GCF_009035845.1 --include genome,protein,gff3</pre>

Files are downloaded as .zip.
Unzip Reference Genome

<pre>sudo apt install unzip   # if not already installed
unzip ncbi_dataset.zip</pre>


Locate FASTA File
Navigate through unzipped folders to locate .fna file. Copy it to the `ref` directory for easy access.
<pre>cp -r GCF_009035845.1.fna ../../..</pre>

***Step 4 : Quast run***
Activate QUAST Environment
<pre>mamba activate quast </pre>

Run QUAST Analysis
Compare reference genome with assembly outputs:

<pre>quast.py -r ref/GCF_009035845.fna -o quast_results contigs.fasta scaffolds.fasta</pre>


# *Outcome:*
QUAST generates a report including:

Number of contigs/scaffolds
N50 and L50 values
Genome fraction covered
Misassemblies and indels

---

Using a reference genome before scaffolding and annotation ensures accurate assembly evaluation and reliable gene prediction. Always verify the .fna file path before running QUAST.

