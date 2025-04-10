## EnrichKit DB

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10257552.svg)](https://doi.org/10.5281/zenodo.10257552)

### Introduction

This repo contains codes for the data curation pipeline for [EnrichKit](https://github.com/liulihe954/EnrichKitWeb). We collected publicly available data set to form a portable sqlite object. Following is the database schema:

![Database Schema](_image/EnrichKitDB_diagram.png)

The sqlite object is available to download via [Zenodo](https://doi.org/10.5281/zenodo.10257552).

The DOI of this sqlite object is - doi.org/10.5281/zenodo.10257552

### Usage

This repository is designed for users with programming background. You may   
1. use the sqlite database object directly by downloading it from [Zenodo](https://doi.org/10.5281/zenodo.10257552). For example querys, please refer to the source code of our [python implementation](https://github.com/liulihe954/EnrichKitPy).


2. replicate our process of creating such sqlite database following:    
    1. Clone this repo and create an conda environment following `conda env create -f environment.yml`
    2. Adjust the content in `config.py` and `create_db.py` to fit your purpose.
    3. run `python3 create_db.py`

**Please Note**:
Please be precautious if you are tring to run it on your laptop. There is a step in the process where gene ID will be compare to online information in a web-scraping style, this is an I/O bound task limted by network factors. This step is time-consuming, for example, this current version used 12 hours with 64 cores. **There is a parameter `CHECK_WEB` in `config.py` file, if you set `CHECK_WEB = False`, the pipeline will skip this step. However, you may get a database with different records than the one posted.**

### Update Info & Citation

Last updated on **04/06/2025**.

Citation: `Liu, L., & Peñagaricano, F. (2023). EnrichKit: a multi-omics tool for livestock research (0.1.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10257552`
