# PBMC Multimodal Single-Cell Analysis

Comprehensive analysis of the GSE164378 dataset from the Hao et al. (2021) Cell paper, featuring both single-modality RNA-seq analysis and multimodal integration using Weighted-Nearest Neighbor (WNN) analysis.

## 📚 Dataset

**Source**: [GSE164378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378)  
**Paper**: Hao, Y., et al. (2021). "Integrated analysis of multimodal single-cell data." *Cell* 184, 3573-3587.e29.  
**DOI**: [10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048)

### Dataset Contents
- **Sample**: Human Peripheral Blood Mononuclear Cells (PBMCs)
- **Cells**: ~161,764 cells (3' data)
- **Modalities**:
  - RNA-seq: Gene expression (~33,000 genes)
  - ADT (CITE-seq): Surface protein expression (~200 proteins)
  - HTO: Hashtag oligonucleotides for sample multiplexing
  - TCR: T cell receptor sequences

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
- **`explore_GSE164378_R.ipynb`**: Complete RNA-seq analysis (Seurat) ✨ NEW

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

### 📊 R Analysis (Seurat) ✨ NEW

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

## 🤝 Contributing

This is an educational project for learning single-cell multimodal analysis. Feel free to:
- Extend the analysis with additional modalities (HTO, TCR)
- Implement alternative integration methods
- Add differential expression analysis
- Explore trajectory analysis

## 📧 Contact

For questions about this analysis, please refer to:
- **Original paper**: Hao et al., Cell (2021)
- **Dataset**: GEO accession GSE164378
- **Scanpy documentation**: https://scanpy.readthedocs.io/

## 📄 License

This analysis code is provided for educational purposes. The GSE164378 dataset is publicly available through GEO and should be used according to the terms specified by the original authors.

---

**Last Updated**: November 2025  
**Analysis Pipeline**: Scanpy + MuData  
**Python Version**: 3.11+ (recommended)

