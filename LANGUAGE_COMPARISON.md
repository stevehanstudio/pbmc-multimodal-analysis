# Python vs R: Single-Cell Analysis Comparison

## 📊 Overview

This project demonstrates the **same single-cell RNA-seq analysis** in both Python (Scanpy) and R (Seurat), allowing direct comparison of the two most popular ecosystems for single-cell genomics.

## 📓 Parallel Notebooks

| Python/Scanpy | R/Seurat |
|---------------|----------|
| `explore_GSE164378.ipynb` | `explore_GSE164378_R.ipynb` |
| `wnn_multimodal_integration.ipynb` | *(future)* |

## 🔬 Key Differences

### Workflow Comparison

| Step | Python/Scanpy | R/Seurat |
|------|---------------|----------|
| **Data Loading** | `sc.read_10x_mtx()` | `Read10X()` |
| **Normalization** | Total count + log1p | `NormalizeData()` |
| **HVG Selection** | Dispersion-based | VST method |
| **Clustering** | Leiden (default) | Louvain (default) |
| **Visualization** | matplotlib/seaborn | ggplot2 |

## 🎯 When to Use Each

### Use Python/Scanpy When:
✅ Large datasets (>100k cells)  
✅ Custom algorithm development  
✅ Deep learning integration  
✅ Team prefers Python  

### Use R/Seurat When:
✅ Standard workflows  
✅ Publication-quality plots  
✅ Bioconductor integration  
✅ Team prefers R  

## 💡 This Project Shows Both!

Demonstrating proficiency in **both** Python and R is highly valuable for:
- Working with diverse teams
- Choosing the right tool for the job
- Bridging between ecosystems
- Maximum flexibility

See the notebooks for detailed implementation in each language.
