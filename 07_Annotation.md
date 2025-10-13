# Annotation with Prokka

*Objective:
Annotate scaffolded genome sequences to identify genes, coding regions, RNAs, and proteins.*

---

**Prerequisites**:

Scaffolded genome FASTA: `ragtag.scaffold.fasta`


---

***Step 1: Create `prokka` environment and install Prokka***

<pre>mamba create -n prokka
conda activate prokka
mamba install -c bioconda prokka -y </pre>


***Step 2: Run Prokka Annotation***

Navigate to the directory containing ragtag.scaffold.fasta 
Run Prokka with appropriate flags:

<pre>prokka --outdir annotation --prefix acetinobacter --addgenes --locustag PROKKA --cpus 2 --norrna --notrna ragtag_output/ragtag.scaffold.fasta </pre>


***Step 3: Inspect Annotation Output***

Output directory prokka_annotation contains multiple files:

| file name | meaning |
|----|----|
| .gbk | GenBank format |
| .faa | Protein sequences |
| .ffn | Nucleotide sequences |
| .gff | Gene features |


Analyze sequences and annotation results.

# *Outcome:*

- Scaffolded genome is fully annotated.
- Protein-coding genes, RNAs, and other genomic features are identified and ready for downstream analysis.

---
