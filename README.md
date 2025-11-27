# PBMC Multimodal Single-Cell Analysis

Comprehensive analysis of the GSE164378 dataset from the Hao et al. (2021) Cell paper, featuring both single-modality RNA-seq analysis and multimodal integration using Weighted-Nearest Neighbor (WNN) analysis.

## 📚 Dataset

**Source**: [GSE164378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378)  
**Paper**: Hao, Y., et al. (2021). "Integrated analysis of multimodal single-cell data." *Cell* 184, 3573-3587.e29.  
**DOI**: [10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048)

## 🧭 Background & Motivation

- **Previous experience**: I first built a single-cell RNA/ATAC multi-omics workflow for the SNARE-seq GSE126074 dataset ([bio2a-multiomics](https://github.com/stevehanstudio/bio2a-multiomics)). That project taught me how to align chromatin accessibility with transcription at the single-cell level.
- **Mentorship guidance**: During a Genentech mentorship meeting, Bryan Ocutt recommended the Hao et al. PBMC study (GSE164378) as the best next step to push my multimodal skills—especially WNN integration of RNA + protein data.
- **Career goal**: This repository documents that deep dive and serves as part of my portfolio for upcoming bioinformatics internships, highlighting deliberate skill growth and mentor-driven learning.

### Dataset Contents
- **Sample**: Human Peripheral Blood Mononuclear Cells (PBMCs)
- **Cells**: ~161,764 cells (3' data)
- **Modalities**:
  - **RNA-seq**: Gene expression (~33,000 genes)
  - **ADT (CITE-seq)**: Surface protein expression (~200 proteins)
  - **HTO**: Hashtag oligonucleotides for sample multiplexing
  - **TCR**: T cell receptor sequences

### Key Concepts

**What is ADT/CITE-seq?**
- **ADT** (Antibody-Derived Tags) measures surface protein expression in single cells
- **CITE-seq** (Cellular Indexing of Transcriptomes and Epitopes by Sequencing) is the technology that enables simultaneous measurement of RNA and proteins
- **How it works**: Antibodies conjugated to DNA oligonucleotides bind to surface proteins; the DNA tags are sequenced alongside RNA
- **Why it's valuable**: 
  - Direct protein measurement (not inferred from RNA)
  - Complements RNA data with different biological information
  - Particularly powerful for immune cells with well-characterized surface markers (CD4, CD8, CD19, etc.)
  - Validates RNA-based findings (e.g., CD4 RNA vs CD4 protein expression)

**What is WNN?**
- **Weighted-Nearest Neighbor (WNN)** analysis integrates multiple data modalities by learning modality-specific weights for each cell
- Determines which modality (RNA vs protein) is most informative for each cell's identity
- Creates an integrated similarity graph that leverages the strengths of both modalities
- Results in improved cell type resolution compared to single-modality analysis

## 🎯 Project Structure

```
pbmc-multimodal-analysis/
├── README.md                           # This file
├── environment.yml                     # Conda environment specification ✨ NEW
├── explore_GSE164378.ipynb            # Python: RNA-seq analysis (Scanpy)
├── explore_GSE164378_R.ipynb          # R: RNA-seq analysis (Seurat) ✨ NEW
├── wnn_multimodal_integration.ipynb   # Python: WNN multimodal integration
├── download_geo_dataset.py            # Utility module for data download
├── LANGUAGE_COMPARISON.md             # Python vs R detailed comparison ✨ NEW
├── R_SETUP_GUIDE.md                   # R installation guide ✨ NEW
├── .gitignore                         # Git ignore patterns
├── .gitattributes                     # Git file handling
├── data/
│   └── GSE164378/
│       ├── suppl/                      # Downloaded GEO supplementary files
│       ├── GSE164378_series_matrix.txt # Sample metadata
│       └── GSE164378_family.xml.tgz    # Dataset metadata
├── results/
│   └── GSE164378/
│       ├── figures/                    # RNA-seq visualizations
│       ├── wnn_figures/                # WNN visualizations
│       ├── GSE164378_rna_processed.h5ad         # Processed RNA data
│       ├── GSE164378_adt_processed.h5ad         # Processed ADT data
│       └── GSE164378_multimodal_wnn.h5mu        # Integrated multimodal data
```

### Key Files

**Python Notebooks:**
- **`explore_GSE164378.ipynb`**: Complete RNA-seq analysis (Scanpy)
- **`wnn_multimodal_integration.ipynb`**: Multimodal WNN integration

**R Notebooks:**
- **`explore_GSE164378_R.ipynb`**: Complete RNA-seq analysis (Seurat)

**Utilities:**
- **`download_geo_dataset.py`**: Data download module with functions:
  - `check_dataset_exists()`: Check if data is already downloaded
  - `download_geo_dataset()`: Download GEO dataset files
  - `ensure_dataset_available()`: Check and download if needed

**Documentation:**
- **`README.md`**: Main project documentation (this file)
- **`LANGUAGE_COMPARISON.md`**: Detailed Python vs R comparison ✨ NEW
- **`R_SETUP_GUIDE.md`**: R installation and setup guide ✨ NEW
- **`PROJECT_OVERVIEW.md`**: High-level project summary

## 📓 Analysis Notebooks

### 🐍 Python Analysis (Scanpy)

#### Notebook 1: `explore_GSE164378.ipynb`
**Single-modality RNA-seq analysis**

Performs standard single-cell RNA-seq analysis pipeline:
1. **Data Loading**: Import 10X Genomics format data and metadata
2. **Quality Control**: Filter low-quality cells and genes
   - Remove cells with < 200 genes
   - Remove genes in < 3 cells
   - Filter cells with > 20% mitochondrial content
3. **Normalization**: Total count normalization + log transformation
4. **Feature Selection**: Identify ~1,580 highly variable genes
5. **Dimensionality Reduction**: PCA → UMAP
6. **Clustering**: Leiden algorithm (29 clusters identified)
7. **Visualization**: Cell types and marker gene expression

**Output**: `GSE164378_rna_processed.h5ad` (used as input for Notebook 2)

#### Notebook 2: `wnn_multimodal_integration.ipynb`
**Multimodal WNN integration**

Integrates RNA-seq and protein (ADT) data using WNN:
1. **Load Pre-processed RNA**: Import from Notebook 1
2. **Load & Process ADT**: 
   - CLR (Centered Log-Ratio) normalization for proteins
   - Dimensionality reduction on protein data
3. **Create Multimodal Object**: Combine using MuData framework
4. **WNN Analysis**: 
   - Compute modality-specific weights per cell
   - Build integrated WNN graph
   - WNN-based UMAP and clustering
5. **Comparison**: Visualize RNA-only vs ADT-only vs WNN results
6. **Cell Type Refinement**: Identify rare populations

**Output**: `GSE164378_multimodal_wnn.h5mu` or combined AnnData files

### 📊 R Analysis (Seurat)

#### Notebook 3: `explore_GSE164378_R.ipynb`
**Single-modality RNA-seq analysis in R**

Equivalent analysis using R/Seurat (parallel to Python notebook):
1. **Load Data**: Read 10X format with Seurat
2. **QC & Filtering**: PercentageFeatureSet, subset
3. **Normalization**: NormalizeData (log-normalization)
4. **Feature Selection**: FindVariableFeatures (VST method)
5. **Scaling & PCA**: ScaleData, RunPCA
6. **Clustering**: FindNeighbors, FindClusters (Louvain)
7. **UMAP**: RunUMAP
8. **Visualization**: DimPlot, FeaturePlot (ggplot2-based)

**Output**: Seurat object (can be saved as .rds)

**Why R?** Demonstrates proficiency in both Python and R ecosystems. See `LANGUAGE_COMPARISON.md` for detailed comparison of the two approaches.

## 🚀 Getting Started

### Prerequisites

#### Python Environment

**Option 1: Using environment.yml (Recommended)**
```bash
# Create environment from file (ensures exact reproducibility)
conda env create -f environment.yml

# Activate environment
conda activate sc-multiomics
```

**Option 2: Manual installation**
```bash
# Create conda environment with Python 3.11
conda create -n sc-multiomics python=3.11

# Activate environment
conda activate sc-multiomics

# Install all packages via conda-forge (better dependency management)
conda install -c conda-forge scanpy pandas matplotlib seaborn jupyter muon
```

**Why conda instead of pip?** Conda provides better dependency management, pre-compiled binaries, and handles system libraries automatically.

**Why Python 3.11?** ~25% faster than 3.9, better type hints, longer support window (until 2027).

#### R Environment (Optional - for R notebooks)

```bash
# Install R (if not already installed)
# Ubuntu/Debian: sudo apt-get install r-base r-base-dev
# macOS: brew install r

# Install R packages (in R console)
R
```

```r
# In R console:
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork"))

# Install IRkernel for Jupyter
install.packages('IRkernel')
IRkernel::installspec(user = TRUE)

# Exit R
quit()
```

**Note**: R notebooks are optional. Python notebooks provide complete analysis.  
**See `R_SETUP_GUIDE.md` for detailed R installation instructions.**

**Package versions:** Using conda ensures compatible versions of all packages. The environment uses Python 3.11 for better performance (~25% faster than 3.9).

### Download Data

The notebooks will automatically check if the dataset is available. If not, you'll be prompted to download it.

**Option 1: Automatic check (Recommended)**
```bash
# Just open the notebook - it will check and guide you
jupyter notebook explore_GSE164378.ipynb
```

**Option 2: Pre-download the dataset**
```bash
# Download metadata and get instructions for supplementary files
python download_geo_dataset.py

# Then download supplementary files (count matrices, ~2-3 GB)
# Follow the wget command shown by the script
```

**Option 3: Manual download**
- Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378
- Download supplementary files to `data/GSE164378/suppl/`

### Run Analysis

**Step 1**: RNA-seq analysis
```bash
jupyter notebook explore_GSE164378.ipynb
```
- The notebook checks if data exists (downloads if needed)
- Execute all cells to perform RNA-seq analysis
- Saves `GSE164378_rna_processed.h5ad` for use in Step 2

**Step 2**: WNN multimodal integration
```bash
jupyter notebook wnn_multimodal_integration.ipynb
```
- Loads pre-processed RNA from Step 1
- Processes ADT (protein) data
- Performs WNN integration
- Compares single-modality vs multimodal results

## 📊 Key Results

### RNA-seq Analysis (Notebook 1)
- **Cells analyzed**: 161,764
- **Highly variable genes**: 1,580
- **Leiden clusters**: 29
- **Cell type categories**: 8 major types, 31 subtypes

### WNN Integration (Notebook 2)
- **Modalities integrated**: RNA + ADT
- **Improved resolution**: Better separation of closely related cell types
- **Rare populations**: Protein markers identify populations missed by RNA alone

## 🔬 Key Findings

1. **Complementary Information**: RNA and protein data capture different aspects of cell identity
2. **Improved Resolution**: Multimodal integration provides better cell type separation
3. **Validation**: Protein expression validates RNA-based predictions
4. **Rare Populations**: Surface protein markers reveal rare cell types

## 💡 Key Insights & Learnings

This section documents what I learned from this analysis and demonstrates my understanding of both the computational methods and biological interpretation.

### Why Different Normalization Methods?

**RNA-seq: Total count + log(x+1)**
- **Why**: RNA-seq data has high variance across genes and cells. Total count normalization accounts for sequencing depth differences, while log transformation stabilizes variance and makes the data more normally distributed for downstream PCA.
- **Learning**: The log transformation is crucial because gene expression follows a log-normal distribution, not a normal distribution. This is why we see better separation in PCA space after log transformation.

**ADT/Protein: CLR (Centered Log-Ratio)**
- **Why**: Protein counts are compositional—they represent proportions of total protein signal per cell. CLR normalizes relative to the geometric mean, preventing artifacts from cells with very high or low total protein counts.
- **Learning**: This was a key insight—protein data requires different normalization because it's fundamentally different from RNA data. Using log-normalization on protein data would incorrectly treat it as independent counts rather than compositional data.

### Why WNN Over Simple Concatenation?

**The Problem**: Simply concatenating RNA and protein features doesn't work well because:
- Different scales (RNA: ~33K features, ADT: ~200 features)
- Different noise profiles (RNA has dropout, ADT is more sparse but less noisy)
- Different biological information (RNA captures transcriptional state, ADT captures surface phenotype)

**The WNN Solution**: 
- Learns which modality is more informative for each cell's identity
- Some cells are better characterized by RNA (e.g., activated T cells with high cytokine expression)
- Others are better characterized by protein (e.g., resting T cells with stable surface markers)
- Creates a weighted graph that leverages the strengths of both modalities

**Learning**: This taught me that multimodal integration isn't just about combining data—it's about intelligently weighting different information sources based on their relevance to each cell's identity.

### Biological Insights from Multimodal Analysis

**1. Validation of RNA-based Annotations**
- Protein markers (e.g., CD4, CD8, CD19) directly validate cell type assignments made from RNA
- This is particularly important for immune cells where surface markers are well-characterized
- **Learning**: Multimodal data provides built-in validation—if RNA says "T cell" and CD3 protein is high, we have confidence in the annotation

**2. Identification of Rare Populations**
- Some cell types (e.g., plasmacytoid dendritic cells) have distinctive protein signatures but subtle RNA differences
- Protein markers can identify rare populations that might be missed in RNA-only analysis
- **Learning**: Different modalities excel at identifying different cell types—multimodal analysis captures the full diversity

**3. Resolution of Ambiguous Clusters**
- Some RNA clusters split into multiple distinct populations when protein data is considered
- Example: A single "T cell" cluster in RNA might separate into CD4+ and CD8+ T cells using protein markers
- **Learning**: Multimodal integration improves resolution, especially for closely related cell types

### Technical Challenges & Solutions

**Challenge 1: WNN Graph Computation**
- **Problem**: Initial attempts failed because `muon` requires modality-specific neighbor graphs to be computed first
- **Solution**: Learned that `mu.pp.neighbors()` expects pre-computed neighbors for each modality
- **Learning**: Understanding the underlying algorithm (WNN needs per-modality KNN graphs) helped debug the issue

**Challenge 2: Metadata Propagation**
- **Problem**: Cell type annotations from RNA data weren't available in the MuData object for plotting
- **Solution**: Explicitly copied annotations from RNA AnnData to MuData.obs
- **Learning**: MuData is a container—metadata must be explicitly propagated to the top level for visualization

**Challenge 3: Normalization Choice**
- **Problem**: Initially tried log-normalization for ADT data (like RNA)
- **Solution**: Researched and implemented CLR normalization, which is standard for compositional data
- **Learning**: Different data types require different preprocessing—understanding the data structure is crucial

### Methodological Choices & Reasoning

**Why Leiden over Louvain clustering?**
- Leiden algorithm is an improvement over Louvain that guarantees well-connected communities
- Better for large datasets (161K cells) where Louvain can produce disconnected clusters
- **Learning**: Algorithm choice matters for large-scale data—Leiden is now the standard for single-cell analysis

**Why PCA before UMAP?**
- UMAP is computationally expensive on high-dimensional data
- PCA reduces dimensions while preserving most variance (typically 50-100 PCs capture >90% variance)
- UMAP on PCA space is faster and often produces better results than UMAP on raw data
- **Learning**: Dimensionality reduction is a two-step process—PCA for linear reduction, UMAP for non-linear embedding

**Why Highly Variable Genes (HVG)?**
- Most genes are not informative for cell type identity (low variance across cells)
- Focusing on HVGs reduces noise and computational cost
- Typically ~1,500-3,000 HVGs capture most biological signal
- **Learning**: Feature selection is critical—more features isn't always better, especially with noisy single-cell data

### Comparison: Python vs R Implementation

**Why implement in both languages?**
- Demonstrates flexibility and understanding of both ecosystems
- Python/Scanpy: Better for large datasets and custom algorithms
- R/Seurat: Better for standard workflows and publication-quality plots
- **Learning**: Different tools have different strengths—being bilingual makes you more versatile

**Key Differences Learned:**
- Scanpy uses dispersion-based HVG selection; Seurat uses VST (Variance Stabilizing Transformation)
- Scanpy defaults to Leiden; Seurat defaults to Louvain (though Leiden is available)
- Scanpy uses matplotlib; Seurat uses ggplot2 (more publication-ready by default)
- **Learning**: The underlying algorithms are similar, but the implementations and defaults differ—understanding both helps choose the right tool

## 📖 Methods Summary

### RNA-seq Processing
- **Normalization**: Total count (10,000) + log(x+1)
- **Feature selection**: Highly variable genes (HVG)
- **Dimensionality reduction**: PCA (50 components) → UMAP
- **Clustering**: Leiden algorithm on KNN graph

### ADT Processing
- **Normalization**: CLR (Centered Log-Ratio) transformation
- **Dimensionality reduction**: PCA → UMAP
- **Clustering**: Leiden algorithm on KNN graph

### WNN Integration
- **Framework**: MuData for multimodal data management
- **Algorithm**: Weighted-Nearest Neighbor (WNN)
  - Learns per-cell modality weights
  - Constructs integrated cell similarity graph
  - Enables joint UMAP and clustering

## 🎓 Learning Resources

### Understanding the Analysis
- **Scanpy Tutorial**: https://scanpy-tutorials.readthedocs.io/
- **CITE-seq Protocol**: Stoeckius et al., Nat Methods (2017)
- **WNN Method**: Hao et al., Cell (2021)

### Tools & Frameworks
- **Scanpy**: https://scanpy.readthedocs.io/
- **MuData**: https://muon.readthedocs.io/
- **Seurat WNN**: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html

## 📝 Citation

If you use this analysis or the GSE164378 dataset, please cite:

```bibtex
@article{hao2021integrated,
  title={Integrated analysis of multimodal single-cell data},
  author={Hao, Yuhan and Hao, Stephanie and Andersen-Nissen, Erica and Mauck III, William M and Zheng, Shiwei and Butler, Andrew and Lee, Maddie J and Wilk, Aaron J and Darby, Charlotte and Zager, Michael and others},
  journal={Cell},
  volume={184},
  number={13},
  pages={3573--3587},
  year={2021},
  publisher={Elsevier},
  doi={10.1016/j.cell.2021.04.048}
}
```

## 🛠️ Troubleshooting

### Common Issues

**Issue**: `muon` not found
```bash
# Solution: Install muon via conda
conda install -c conda-forge muon
```

**Issue**: Processed RNA file not found in Notebook 2
```bash
# Solution: Run Notebook 1 first to generate GSE164378_rna_processed.h5ad
```

**Issue**: Memory errors with large dataset
```bash
# Solution: Use a subset of cells or increase available RAM
# Consider using sparse matrices and chunked processing
```

## 📄 License

This analysis code is provided for educational purposes. The GSE164378 dataset is publicly available through GEO and should be used according to the terms specified by the original authors.


