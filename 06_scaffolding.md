# Scaffolding 

*Objective:
Scaffold assembled contigs using a reference genome to generate longer, ordered sequences.*

---

**Prerequisites:
1. Reference genome in .fna format

2. Assembled contigs in contigs.fasta

---

***Step 1: Create  `ragtag` environment and install RagTag***

<pre>mamba create -n ragtag
mamba activate ragtag
mamba install -c bioconda ragtag -y </pre>

***Step 2: Prepare for Scaffolding***

- Verify the paths of your reference GCF_009035845.1.fna file and the contigs.fasta file.
- Run RagTag Scaffolding

<pre>ragtag.py scaffold reference.fna contigs.fasta -o ragtag_output</pre>


***Step 3: Inspect Scaffolds***

- Scaffolded sequences are saved in ragtag_output/ragtag.scaffold.fasta.

- Preview sequences with:

<pre>head ragtag_output/ragtag.scaffold.fasta</pre>


*Outcome:*

-Contigs are connected with NNNNN gaps.

---

The scaffolded FASTA is ready for downstream annotation.

