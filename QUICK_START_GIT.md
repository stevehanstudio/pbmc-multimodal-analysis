# 🚀 Quick Start: Git Commands

## Initial Setup (First Time Only)

```bash
# Initialize repository
git init

# Add all files (respects .gitignore)
git add .

# First commit
git commit -m "Initial commit: PBMC multimodal analysis project"

# Connect to GitHub (create repo on GitHub first)
git remote add origin https://github.com/YOUR_USERNAME/pbmc-multimodal-analysis.git
git branch -M main
git push -u origin main
```

## Daily Workflow

```bash
# Check status
git status

# Add changes
git add explore_GSE164378.ipynb README.md
# or add everything: git add .

# Commit with message
git commit -m "Descriptive message about what changed"

# Push to GitHub
git push
```

## Useful Commands

```bash
# View commit history
git log --oneline

# See what changed
git diff

# Undo uncommitted changes
git checkout -- filename

# Create a branch
git checkout -b feature-name

# Switch branches
git checkout main
```

## For Collaborators

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/pbmc-multimodal-analysis.git
cd pbmc-multimodal-analysis

# Install dependencies (conda environment)
conda env create -f environment.yml
conda activate sc-multiomics

# Download data (not in git)
python download_geo_dataset.py
```

---

**See `GIT_GUIDE.md` for detailed documentation!**
