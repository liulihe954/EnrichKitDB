## EnrichKit DB

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10257552.svg)](https://doi.org/10.5281/zenodo.10257552)

### Introduction

This repo contains codes for the data curation pipeline for EnrichKit ([Web App implementation](https://github.com/liulihe954/EnrichKitWeb) & [Python implementation](https://github.com/liulihe954/EnrichKitPy)). We collected publicly available data set to form a portable sqlite object. Following is the database schema:

![Database Schema](_image/EnrichKitDB_diagram.png)

The sqlite object is available to download via [Zenodo](https://doi.org/10.5281/zenodo.10257552).

The DOI of this sqlite object is - doi.org/10.5281/zenodo.10257552

### Usage

- Direct download from Zenodo and establish connections with `SQLite`

- In addition to download the sqlite object, you may replicate such a sqlite database by following:

1. Clone this repo and create an conda virtual environment using `conda env create -f environment.yml`
2. Adjust the content in `config.py` and `create_db.py` to fit your purpose.
3. run `python3 create_db.py`

Note:

1. The current version is created with python3.12
2. There is a step in the process where gene ID will be compare to online information in a web-scraping style, this is an I/O bound task which are potentially time-consuming, e.g., this current version used 12 hours with 64 cores.

### Update Info & Citation

Last updated on **01/15/2024**.

Citation: `Liu, L., & Peñagaricano, F. (2024). EnrichKit: a multi-omics tool for livestock research (0.1.2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10535657`
