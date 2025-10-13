# 🧫 CLASSIFICATION

## 🔹 Introduction

**Classification** in metagenomic analysis helps determine **what types of organisms** are present in a biological sample based on sequencing reads.  
This is crucial for identifying contaminants, mixed infections, or microbial community composition.

In this workflow, we used **_Escherichia coli_ (E. coli)** with **SRA ID: `SRR35320464`** as the sample organism.  
The input files required for classification were the **assembled FASTA files** obtained from previous assembly steps, such as:

- `contigs.fasta`  
- `scaffolds.fasta`  

These files contain assembled genomic sequences that are analyzed to determine their taxonomic origin.


In this workflow, we use two primary tools:

- 🧬 **Kraken2** — for taxonomic classification of reads using k-mer matching.  
- 📊 **Krona** — for interactive visualization of classification results.

Because these tools are memory-intensive (requiring ≥ 4 GB RAM), we use a special configuration that prevents files from being loaded entirely into RAM, allowing them to run efficiently on limited-memory systems.

---

## ⚙️ Tools Installation

Activate the appropriate environment or create a new one for classification:

<pre>mamba create -n classify
mamba activate classify</pre>

Install Kraken2, Krona, and Aria2

<pre>mamba install -c bioconda kraken2 krona -y
mamba install -c bioconda aria2 -y</pre>

### 📝 Explanation:

`kraken2` → performs read classification based on k-mer mapping to a reference database.

`krona` → creates an interactive HTML chart for visualization.

`aria2` → helps in parallel downloading of large databases or files.

🧠 Handling Memory Limitations
Normally, running Kraken2 and Krona requires at least 8 GB of RAM, which can cause problems on systems with limited memory.

To overcome this, we use `memory mapping`, which allows Kraken2 to read database files directly from disk instead of loading them entirely into RAM.

### 💻 Step-by-Step Workflow

*Using the Advit Pipeline*

The Advit pipeline automates downloading, setup, and classification.
It installs the required tools, builds or accesses existing databases, runs classification using Kraken2, and creates visual output with Krona.

***Step 1: Create a Shell Script (Advit.sh)***

Use the nano editor to create a script that will run Kraken and Krona with memory mapping enabled.
<pre>nano Advit.sh</pre>

Inside the file, type (or add) the following line:
`--memory-mapping \`

*📝 Note:*
This command enables disk-based memory mapping so that large Kraken databases can be accessed efficiently even with limited RAM (4–8 GB).

Save and exit the file (Ctrl + O, then Ctrl + X).

***Step 2: Run the Script***
Execute the script with appropriate parameters.
Here:

`--threads 4` → uses 4 CPU threads.

`--have-db no` → skips database download if not needed.

`--skip-bracken` → skips abundance estimation.

`--fasta` → provides the path to the FASTA files to be classified.


<pre>./Advit.sh --threads 4 --have-db no --skip-bracken --fasta ../e.coli/SRR35320464_spades/FASTA</pre>

🧩 Explanation:
This command classifies the input FASTA sequences (assembled reads) and generates output files for taxonomic classification and visualization.

***Step 3: Viewing Results in Krona***
After the analysis is complete:

Go to the Advit directory in your files.

Open the results/ folder.

Locate and click the Krona output file (usually an .html file).

The Krona interactive chart will open in your web browser, displaying a circular visualization of all identified taxa.


*🧠 Tip:*
Hovering over each section of the Krona chart displays the name and relative abundance of each organism identified in the sample.


*📊 Example Interpretation*
The Krona chart provides an intuitive way to explore classification results:

The center represents the highest taxonomic level (e.g., Bacteria, Archaea, Eukaryota).

Moving outward, rings represent more specific classifications — Phylum → Class → Order → Family → Genus → Species.

The size of each wedge indicates the relative abundance of that taxon in the sample.
