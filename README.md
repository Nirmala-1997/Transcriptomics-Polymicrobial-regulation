# Transcriptomics-Polymicrobial regulation
Transcriptomics-Polymicrobial regulation of Staphylococcus aureus
# RNA-seq Analysis Pipeline and Visualization
## Repository Structure
This repository contains two main components of the RNA-seq study:

---


## 1) RNA-seq Visualization Notebook
**File:** `GitHub_Copy_RNA_Seq_Visualization.ipynb`

- Jupyter/Colab notebook for visualizing RNA-seq data outputs.  
- Includes plotting functions (heatmaps, PCA, volcano plots, scatterplots, etc.) for interpreting differential expression results.  
- Designed for quick exploration and figure generation.

To open in Google Colab:  
[Open in Colab]
https://colab.research.google.com/github/Nirmala-1997/Transcriptomics-Polymicrobial-regulation/blob/main/GitHub_Copy_RNA_Seq_Visualization.ipynb

---

## 2) RNA-seq Pipeline Script
**File:** `pipeline.R`

This R script documents the complete RNA-seq workflow:

1. **Quality Control & Preprocessing**  
   - FastQC v0.12.1  
   - fastp v0.23.4  

2. **De novo Transcriptome Assembly**  
   - Trinity v2.15.1  
   - Assembly statistics via assembly-stats v1.0.1  

3. **Transcript Quantification & DEG Analysis**  
   - RSEM v1.3.3 for expression quantification  
   - edgeR v3.42.4 for differential expression  
   - Excludes uniquely expressed genes (expressed in only one sample)  
   - Groups compared:  
     - SA-C vs SA  
     - SA-P vs SA  
     - SA-E vs SA  

4. **Functional Annotation**  
   - TransDecoder v5.5.0 for coding region prediction  
   - BLAST+ v2.14.1 against UniProt Swiss-Prot  
   - Annotations filtered with ≥50% identity threshold  

---

#############
## License

This project is licensed under the [MIT License](LICENSE). Feel free to use, modify, and distribute—with attribution.

## Citation

If you use this pipeline in your research or publication, please cite it:

B. Nirmala (2025). RNA-seq QC, Preprocessing, Assembly, and DEG Analysis Pipeline. GitHub repository: https://github.com/Nirmala-1997/Transcriptomics-Polymicrobial-regulation
