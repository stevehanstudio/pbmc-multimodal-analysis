# R Setup Guide for Jupyter Notebooks

## 📋 Quick Setup

### Step 1: Install R

**Ubuntu/Debian:**
```bash
sudo apt-get install r-base r-base-dev
```

**macOS:**
```bash
brew install r
```

**Windows:**  
Download from: https://cran.r-project.org/bin/windows/base/

### Step 2: Install R Packages

Open R console and run:
```r
# Install Seurat and dependencies
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork"))

# Install IRkernel for Jupyter
install.packages('IRkernel')
IRkernel::installspec(user = TRUE)

# Exit
quit(save = "no")
```

### Step 3: Verify Installation

```bash
# Check if R kernel is available
jupyter kernelspec list

# You should see 'ir' in the list
```

## 🚀 Running R Notebooks

```bash
# Start Jupyter
jupyter notebook

# Open explore_GSE164378_R.ipynb
# The kernel should automatically be set to "R"
```

## 🔍 Troubleshooting

### Issue: "IRkernel not found"
```r
# Reinstall IRkernel
install.packages("IRkernel")
IRkernel::installspec(user = TRUE)
```

### Issue: "Cannot install Seurat"
```bash
# Install system dependencies (Ubuntu)
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
```

## 📚 Resources

- **Seurat**: https://satijalab.org/seurat/
- **IRkernel**: https://irkernel.github.io/
- **Tutorials**: https://satijalab.org/seurat/articles/

## ✅ Quick Test

Test your setup:
```r
library(Seurat)
print(packageVersion("Seurat"))
```

If this works, you're ready to run the R notebook! 🎉
