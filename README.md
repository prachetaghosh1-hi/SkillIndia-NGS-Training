# SkillIndia-NGS-Training
Practical NGS pipeline and analysis workflow learned during Skill India Genomics and Molecular Technologies training


This repository documents my **1-month training under the Skill India Programme** on the topic:  
**"Genomics and Molecular Technologies for Industrial Applications"**, conducted by **CSIR-IGIB**.  

During this training, I learned and practiced the **entire workflow of Next-Generation Sequencing (NGS) analysis**, starting from basic Linux commands to genome classification and AMR detection.  

---

## Folder Structure

| Step | Folder | Description |
|------|--------|-------------|
| 1 | `01_Linux_basics.md` | Basic Linux commands and environment setup |
| 2 | `02_SRA_Data_Retrieval_and_Quality_Assessment.md` | SRA file download from NCBI (organism you wish to work on) Quality check for short and long reads |
| 3 | `03_Trimming.md` | trimming of sequencing reads for short and long reads  |
| 4 | `04_Assembly.md` | Genome assembly for short and long readsGenome annotation using Prokka |
| 5 | `05_Comparing_Genome_Assembly_to_Reference_Genome.md` | reference file needed for scaffolding and assembly |
| 6 | `06_scaffolding.md`  
| 7 | `07_Annotation.md` | see the biological significance of your organism |
| 8 | `08_AMR_Detection.md` 
| 9 | `09_Variant_calling.md` | Variant detection using reference genome |
| 10 | `10_Classification.md` | Genome classification and taxonomy analysis |


---

## Notes

When practicing genome analysis with ***short-read data***I chose a random Illumina sequencing dataset of an organism. I ensured the dataset was not too large, as processing very long or high-coverage data can be difficult on a personal computer.

The main goal of this analysis is to obtain biological information about a particular strain — for example, *Pseudomonas aeruginosa* (SRR ID: [example]). On the NCBI website, sequencing data for this organism is provided in FASTQ format, which is the raw output directly obtained from a sequencer.

A FASTQ file contains millions of short DNA sequences (reads), each separated by a header line starting with @. Each read entry has four lines:

@Read_ID

ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGC

+

IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

Line 1: Header starting with @ (read identifier)

Line 2: Nucleotide sequence (A, T, G, C)

Line 3: A plus sign (+), sometimes followed by the header again

Line 4: Quality scores for each nucleotide

However, a FASTQ file by itself cannot directly give us the information we need because these short reads are randomly arranged in the file. This randomness happens because:

1. The sequencer cannot read the entire genome at once.

2. During wet-lab sample preparation, the DNA is fragmented into small pieces, so the sequences recorded in FASTQ files are short fragments from random positions in the genome.

As a result, the FASTQ file contains short, non-sequential DNA fragments that must be assembled into longer stretches before any meaningful analysis, such as gene annotation, can be performed.

**Assembly: From FASTQ Reads to Contigs**

Before assembly, the reads are trimmed to remove adapters and low-quality bases. This ensures that only high-quality sequences are used for downstream analysis.

During the assembly step, the trimmed short reads are aligned based on overlapping regions and merged to form contigs — longer continuous sequences representing stretches of the original genome.

Read 1: ATGCGTACCTGATCGATCGGATCGTACGTA

                 ||||||||||||||||||||||
                 
Read 2:          CGATCGGATCGTACGTAGCTAGCTGATCGA

                          ||||||||||||||||||||||
                          
Read 3:                   GTACGTAGCTAGCTGATCGATCGTGGATCG



>contig_1

ATGCGTACCTGATCGATCGGATCGTACGTAGCTAGCTGATCGATCGTGGATCG


Not all reads overlap perfectly, so usually multiple contigs are generated rather than a single continuous sequence. Contigs are the foundation for further steps such as scaffolding and annotation, which allow us to identify genes and other genomic features, including AMR (antimicrobial resistance) genes.

---

- **Short-read workflows** were completed fully, including QC, trimming, assembly, annotation, variant calling, AMR detection, and classification.  
- **Long-read workflows** were completed up to assembly due to hardware limitations.  



---

## Tools used for short and long reads:

| Step                       | Short Reads Tools                  | Long Reads Tools                 | Tool Used |
|-----------------------------|----------------------------------|--------------------------------|-----------|
| Data Retrieval              | SRA Tools                         | SRA Tools                       | SRA Tools |
| Quality Check (FASTQ)       | FastQC, MultiQC                   | NanoPlot                        | FastQC, MultiQC, NanoPlot |
| Trimming                    | Trim Galore, Fastp                | Porechop                        | Fastp, Porechop |
| Assembly                    | SPAdes                            | Canu, Flye                      | SPAdes, Canu |
| Assembly Quality Check      | QUAST                             | QUAST                            | QUAST     |
| Scaffolding                 | RagTag                            | Not required if using Flye       | RagTag    |
| Annotation                  | Prokka                            | Prokka                           | Prokka    |

---

| Step                  | Tool        | Purpose                                           | Notes |
|-----------------------|------------|--------------------------------------------------|-------|
| Variant Calling       | BWA        | Align reads to reference genome                 | Generally used for short reads |
| Variant Calling       | SAMtools   | File conversion, sorting, filtering, and indexing | Works with BAM/SAM files |
| Variant Calling       | Picard     | Mark duplicate reads                             | Improves accuracy of variant calling |
| Variant Calling       | FreeBayes  | Perform variant calling                          | Detects SNPs and small indels |
| Variant Calling       | BCFtools   | Process and count variants                       | Works with VCF/BCF files |
| Variant Calling       | SnpEff     | Annotation of the variant file                   | Predicts effects of variants |
| Variant Calling       | Minimap    | Align long reads to reference genome without indexing | Used if long reads are available |
| Classification        | Kraken2    | Taxonomic classification of reads                | Same tool for short and long reads |
| Classification        | Krona      | Visualize taxonomic classification              | Generates interactive plots |
| AMR Detection         | AMRFinder  | Detect antimicrobial resistance genes            | Same tool for short and long reads |

---

This repository serves as a **complete record of my practical learning** during the Skill India training at CSIR-IGIB.  

