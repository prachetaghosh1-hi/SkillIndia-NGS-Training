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

This repository serves as a **complete record of my practical learning** during the Skill India training at CSIR-IGIB.  

