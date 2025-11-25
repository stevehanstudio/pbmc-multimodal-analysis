#!/usr/bin/env python3
"""
Utility functions for downloading and managing GEO datasets.

This module provides functions to download single-cell multimodal datasets
from the Gene Expression Omnibus (GEO) database.
"""

import os
import urllib.request
import gzip
import shutil
from pathlib import Path
from typing import Union, Dict


def check_dataset_exists(accession: str, data_dir: str = "data") -> bool:
    """
    Check if a GEO dataset has already been downloaded.
    
    Args:
        accession: GEO accession number (e.g., GSE164378)
        data_dir: Base directory where datasets are stored
        
    Returns:
        bool: True if dataset directory and key files exist, False otherwise
    """
    dataset_path = Path(data_dir) / accession / "suppl"
    
    # Check if the supplementary files directory exists and has files
    if dataset_path.exists() and any(dataset_path.iterdir()):
        return True
    return False


def download_geo_dataset(accession: str, output_dir: str = "data", 
                         download_suppl: bool = True) -> Path:
    """
    Download a GEO dataset by accession number.
    
    This function downloads metadata files and optionally the supplementary
    data files (which contain the actual count matrices for single-cell data).
    
    Args:
        accession: GEO accession number (e.g., GSE164378)
        output_dir: Directory to save the downloaded files
        download_suppl: Whether to download supplementary files (can be large)
        
    Returns:
        Path: Path to the downloaded dataset directory
        
    Example:
        >>> dataset_path = download_geo_dataset("GSE164378")
        >>> print(f"Dataset downloaded to: {dataset_path}")
    """
    # Create output directory
    output_path = Path(output_dir) / accession
    output_path.mkdir(parents=True, exist_ok=True)
    
    print(f"📥 Downloading dataset {accession}...")
    print(f"   Output directory: {output_path.absolute()}")
    
    # GEO FTP URLs
    base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series"
    
    # Extract series number (e.g., GSE164378 -> GSE164nnn)
    series_num = accession[::-1][3:][::-1] + "nnn"
    
    # Construct URLs
    urls = {
        "series_matrix": f"{base_url}/{series_num}/{accession}/matrix/{accession}_series_matrix.txt.gz",
        "family_xml": f"{base_url}/{series_num}/{accession}/miniml/{accession}_family.xml.tgz",
        "suppl": f"{base_url}/{series_num}/{accession}/suppl/",
    }
    
    # Download series matrix file
    print("\n1️⃣  Downloading series matrix file...")
    matrix_file = output_path / f"{accession}_series_matrix.txt.gz"
    try:
        urllib.request.urlretrieve(urls["series_matrix"], matrix_file)
        print(f"   ✅ Downloaded: {matrix_file.name}")
        
        # Extract the gzipped file
        print("   📦 Extracting...")
        with gzip.open(matrix_file, 'rb') as f_in:
            with open(matrix_file.with_suffix('').with_suffix('.txt'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"   ✅ Extracted: {matrix_file.with_suffix('').with_suffix('.txt').name}")
    except Exception as e:
        print(f"   ❌ Error downloading series matrix: {e}")
    
    # Download family XML file
    print("\n2️⃣  Downloading family XML file...")
    xml_file = output_path / f"{accession}_family.xml.tgz"
    try:
        urllib.request.urlretrieve(urls["family_xml"], xml_file)
        print(f"   ✅ Downloaded: {xml_file.name}")
    except Exception as e:
        print(f"   ❌ Error downloading family XML: {e}")
    
    # Supplementary files
    print("\n3️⃣  Supplementary files...")
    if download_suppl:
        print(f"   ⚠️  Note: Supplementary files can be very large (several GB)")
        print(f"   📍 URL: {urls['suppl']}")
        print(f"\n   To download supplementary files, use wget:")
        print(f"   wget -r -np -nH --cut-dirs=4 {urls['suppl']} -P {output_path}/")
        print(f"\n   Or visit the URL in a browser to download specific files.")
    else:
        print(f"   ⏭️  Skipping supplementary files (set download_suppl=True to download)")
        print(f"   📍 URL: {urls['suppl']}")
    
    print(f"\n✅ Dataset {accession} metadata downloaded!")
    print(f"   📁 Files saved to: {output_path.absolute()}")
    
    return output_path


def ensure_dataset_available(accession: str, data_dir: str = "data") -> Path:
    """
    Ensure a GEO dataset is available locally, downloading if necessary.
    
    This is a convenience function that checks if the dataset exists and
    downloads it if it doesn't. Perfect for use in notebooks.
    
    Args:
        accession: GEO accession number (e.g., GSE164378)
        data_dir: Base directory where datasets are stored
        
    Returns:
        Path: Path to the dataset directory
        
    Example:
        >>> # In a notebook, this will download only if needed
        >>> dataset_path = ensure_dataset_available("GSE164378")
    """
    if check_dataset_exists(accession, data_dir):
        dataset_path = Path(data_dir) / accession
        print(f"✅ Dataset {accession} already exists at: {dataset_path.absolute()}")
        return dataset_path
    else:
        print(f"📥 Dataset {accession} not found locally. Downloading...")
        return download_geo_dataset(accession, data_dir, download_suppl=False)

if __name__ == "__main__":
    """
    Command-line interface for downloading GEO datasets.
    
    Usage:
        python download_geo_dataset.py
    """
    accession = "GSE164378"
    
    print("="*70)
    print("GEO Dataset Downloader")
    print("="*70)
    
    # Check if dataset already exists
    if check_dataset_exists(accession):
        print(f"\n✅ Dataset {accession} already downloaded!")
        print(f"   Location: {Path('data') / accession}")
        print(f"\n   To re-download, delete the directory and run again.")
    else:
        # Download the dataset
        output_dir = download_geo_dataset(accession, download_suppl=False)
        
        print("\n" + "="*70)
        print("📥 METADATA DOWNLOAD COMPLETE")
        print("="*70)
        print(f"\n⚠️  To download the supplementary files (count matrices, ~2-3 GB):")
        print(f"\n   wget -r -np -nH --cut-dirs=4 \\")
        print(f"     https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164378/suppl/ \\")
        print(f"     -P {output_dir.absolute()}/")
        print(f"\n   Or download specific files from:")
        print(f"   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}")
        print("\n" + "="*70)


