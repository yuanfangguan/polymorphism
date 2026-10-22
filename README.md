# Sequence-only Deep Learning for Predicting Expected Neutral Polymorphism

This repository contains the official code for our paper: **"Detecting natural selection in both coding and noncoding sequences by learning features of genomic polymorphisms"**.

Our deep learning framework predicts the baseline sequence-driven mutability (expected neutral polymorphism) across the human genome using primary DNA sequence context alone. By comparing observed genomic variation against this neutral baseline, we compute the **OVEP (Observed vs. Expected Polymorphism)** metric to infer purifying and diversifying selection without relying on simulated evolutionary data.

## 🗂 Repository Structure

The codebase is organized into several key modules mirroring the analytical pipeline of the manuscript:

### 1. Data Processing (`data/`, `data_mask*/`)
Scripts to parse the reference genome (FASTA) and population variation data (VCFs like 1000 Genomes), apply strict genomic masking strategies, and encode them into binary arrays for model training.
* `data/`: Standard processing pipeline (excluding ambiguous 'N' bases).
* `data_maskfunctional/` & `data_maskrepeat/`: Processing scripts for functional and repeat-masked ablation studies.

### 2. Model Training & Prediction (`code_halftrain*/`)
Contains the 1D U-Net Convolutional Neural Network architecture and scripts for 5-fold cross-validation, training, and whole-genome prediction.
* `code_halftrainmaskmissing/`: The core model architecture (`cnn_batchnorm.py`) and training pipelines for the unannotated-region model. Includes ablation studies for different receptive fields (512bp, 1024bp, 2048bp, 4096bp).
* `code_halftrainmaskfunctional_maskrepeat/`: Models trained strictly on non-functional, non-repetitive sequences.
* `motif/` & `ratio/`: Code used to interpret the model, including extraction of highly mutable motifs and triplet analyses.

### 3. Evolutionary & Downstream Analysis (`post_analysis_maskmissing/`)
Scripts to calculate the **OVEP** metric and evaluate signatures of selection across functional genomic elements.
* `category_annotation_calibration/`: Computes observed vs. expected polymorphism ratios across diverse genomic annotations (CDS, enhancers, promoters, repeats).
* `exon_analysis/`: Computes the correlations between our sequence-based mutability predictions and inter-species conservation metrics ($d_N/d_S$, $K_a/K_s$) using Chimpanzee and Mouse/Rat outgroups.

### 4. PAML Integration (`PAML/`)
Pipeline to integrate predicted mutation rates with phylogenetic $d_N/d_S$ calculations using `hpcleap_dnds`.

### 5. Benchmarking & Revisions (`external_test/`, `revision/`)
Comprehensive benchmarking suite comparing our model against state-of-the-art tools and handling specific reviewer requests.
* `revision/JARVIS/`: Large-scale evaluation pipeline integrating and comparing against existing genome-wide scores (e.g., gnomAD, Roulette, JARVIS, CADD, DANN, LINSIGHT, ncRVIS, Orion).
* Includes stratification evaluation across protein-coding, annotated non-coding, and unannotated genomic regions.

---

## ⚙️ Installation & Requirements

The scripts are written in Python 3 and Bash. Core Python dependencies include:
```bash
pip install numpy pandas scipy scikit-learn tensorflow keras pyfaidx liftover


Use scripts in data/ to convert reference FASTA and VCF files into .npy or .pkl arrays.

python data/pickle_chr.py
Step 2: Train the U-Net Model
Navigate to the desired model configuration (e.g., base 2048 window) and run the training pipeline.

cd code_halftrainmaskmissing/base_excludeN_5fold
python train_batchnorm.py
Step 3: Genome-Wide Prediction
Generate calibrated expected polymorphism scores across all chromosomes.

python predict.py
Step 4: Calculate OVEP & Evaluate Selection
Evaluate predictions against functional GFF annotations.

cd post_analysis_maskmissing/category_annotation_calibration/
python calculate_snp_by_category.py
