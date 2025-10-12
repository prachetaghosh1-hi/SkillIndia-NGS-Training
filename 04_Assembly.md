# Assembly

## **why we do Assembly?**
Trimming as seen in the last step improved the quality and accuracy of genome assembly, which is why assembly is always performed on trimmed reads rather than raw sequencing data

---

This file demonstrates the steps for **assembling a genome from short-read sequencing data**, evaluating the assembly, and exploring sequences.  
We use **SPAdes** for assembly, **QUAST** for evaluation, and **SeqKit** for sequence manipulation.

---

***Step1: Installation of SPAdes in new environment called `spades`***
<pre>mamba create -n spades
mamba activate spades
mamba install -c bioconda spades -y</pre>

***Step 2: Genome Assembly using SPAdes***
<pre>spades.py -1 trim_pa/SRR35114395_1_trimmed.fastq -2 trim_pa/SRR35114395_2_trimmed.fastq -o trim_pa/spades --threads 4 --memory 16 </pre>

- Check total scaffolds:
<pre>wc -l trim_pa/spades/scaffolds.fasta</pre>

***Step 3: Installation of QUAST and Run QUAST for Assembly Evaluation***

<pre>mamba create -n quast
mamba activate quast
mamba install -c bioconda quast -y</pre>

<pre>quast trim_pa/spades/scaffolds.fasta -o quast_output</pre>


*Output*: QUAST provides reports including metrics like total length, N50, GC content, and scaffold count, stored in `quast_output` directory.

***Step 4: Installation of Seqkit and Use of SeqKit for Sequence Analysis***
<pre>mamba create -n seqkit_sudo 
  mamba activate seqkit_sudo
  mamba install -c bioconda seqkit sudo -y</pre>

Purpose: SeqKit allows exploration, filtering, and summarizing of assembled sequences.

*Example Commands*

Get statistics of scaffolds:
<pre>seqkit stats trim_pa/spades/scaffolds.fasta </pre>


View the first 5 scaffolds:
<pre>seqkit head -n 5 trim_pa/spades/scaffolds.fasta </pre>




# *Outputs*

1. Assembled scaffolds: `trim_pa/spades/scaffolds.fasta`

2. Filtered scaffolds: `trim_pa/spades/scaffolds_filtered.fasta`

3. QUAST evaluation results: `quast_output/` containing reports and scaffold files

4. SeqKit summaries: printed in terminal or saved to file

---

