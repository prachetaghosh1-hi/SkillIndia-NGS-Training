# Antimicrobial Resistance (AMR) Analysis with NCBI AMRFinderPlus

*Objective:
Identify antimicrobial resistance (AMR) genes in scaffolded genome sequences using AMRFinderPlus.*

---

NOTE: We worked on multiple organisms for short reads, including *E. coli*, *Pseudomonas aeruginosa*, *Acinetobacter baumannii* and *Staphylococcus aureus*. For this part, we used *Pseudomonas aeruginosa* with SRA ID: SRR35514395 as the example organism.

---

*Prerequisites* : scaffolded file (ragtag.scaffold.fasta)

---

***Step1: Create  AMRFinder environment and install AMRFinderPlus***

<pre>mamba create -n amrfinder
conda activate amrfinder
mamba install -c bioconda ncbi-amrfinderplus -y</pre>


***step 2: Update AMRFinderPlus database***

<pre>amrfinder --update </pre>


***step 3: Check available commands and flags***

<pre>amrfinder --help</pre>


***Step 4: Nucleotide-based AMR search***

<pre>amrfinder --nucleotide pseudomonas/ragtag_output/ragtag.scaffold.fasta --organism Pseudomonas_aeruginosa --output amr.txt  --threads 4</pre>


***Step 5: View results:***

<pre>nano amr.txt</pre>


*Output*
-includes AMR genes and associated antibiotics.

***Step 6: Protein-based AMR search***

<pre> amrfinder --protein annotation/pseudomonas.faa --gff annotation/pseudomonas.gff --annotation_format prokka --organism Pseudomonas_aeruginosa --output amr_protein.txt --threads 10</pre>


***step 7: View results***

<pre>nano amr_protein.txt</pre>


Output lists protein sequences linked to AMR genes and resistance phenotypes.

*Outcome:*

- Identified nucleotide and protein AMR genes in the scaffolded genome.
- Results can be used for further analysis of antibiotic resistance profiles.
