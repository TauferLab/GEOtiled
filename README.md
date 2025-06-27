![GEOtiled Logo](docs/images/logo.png)

## About

Terrain parameters such as slope, aspect, and hillshading are essential in various applications, including agriculture, forestry, 
and hydrology. However, generating high-resolution terrain parameters is computationally intensive, making it challenging to 
provide these value-added products to communities in need. We present a scalable workflow called GEOtiled that leverages data 
partitioning to accelerate the computation of terrain parameters from digital elevation models, while preserving accuracy.

This repository contains the library for all functions used by GEOtiled, and includes a Jupyter Notebook walking through 
GEOtiled's workflow and core features.

## Dependencies

### Supported Operating Systems

1. All work has been developed in [Ubuntu 24.04](https://releases.ubuntu.com/noble/)

### Required Software
> Note: These have to be installed on your own (instructions for SAGA and Performance Copilot given)

1. [Git](https://git-scm.com/downloads)
2. [Python 3](https://www.python.org/downloads/)
3. [Anaconda 3](https://www.anaconda.com/download/)
4. [SAGA](https://sourceforge.net/p/saga-gis/wiki/Binary%20Packages/)
5. [Performance Copilot](https://jupyter.org/install)

### Required Libraries
> Note: These will be installed with GEOtiled

1. matplotlib
2. geopandas
3. notebook
4. tabulate
5. pandas
6. tqdm
7. GDAL

## Installation

### Install Conda
> If you already have Conda installed on your machine, skip to Install GEOtiled
1. Download Anaconda
```
wget https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh
```
2. Run the downloaded file and agree to all prompts
```
bash ./Anaconda3-2023.09-0-Linux-x86_64.sh
```
3. Restart the shell to complete the installation

### Install GEOtiled

1. Create a new conda environment
   > Note: This process might take some time
```
conda create -n geotiled -c conda-forge gdal=3.8.0
```
2. Change to the new environment
```
conda activate geotiled
```
3. Clone the repository in a desired working directory
```
git clone https://github.com/TauferLab/GEOtiled
```
4. Change to the geotiled directory
   > Note: `your_path` should be replaced with your working directory
```
cd your_path/GEOtiled/geotiled
```
5. Install editable library
```
pip install -e .
```

> Note: Installations can be verified with `conda list`

### Install System Packages

1. Update apt-get
```
sudo apt-get update
```
2. Install SAGA
```
sudo apt-get install saga
```
3. Install Performance Copilot
```
sudo apt-get install pcp-zeroconf
```

## How to Use the Library

1. Ensure you are in the correct conda environment
```
conda activate geotiled
```
2. Place the following code snippet towards the top of any Python code to use GEOtiled functions
```
import geotiled
```
> Note: Documentation on functions can be found under docs/build/html/index.html

## Publications

Camila Roa, Paula Olaya, Ricardo Llamas, Rodrigo Vargas, and Michela Taufer. 2023. **GEOtiled: A Scalable Workflow
for Generating Large Datasets of High-Resolution Terrain Parameters.** *In Proceedings of the 32nd International Symposium 
on High-Performance Parallel and Distributed Computing* (HPDC '23). Association for Computing Machinery, New York, NY, USA, 
311–312. [https://doi.org/10.1145/3588195.3595941](https://doi.org/10.1145/3588195.3595941)

## Copyright and License

Copyright (c) 2024, Global Computing Lab

GEOtiled is distributed under the 3-Clause BSD License.

See [LICENSE](./LICENSE) for more details.

## Acknowledgements

SENSORY is funded by the National Science Foundation (NSF) under grant numbers [#1724843](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1724843&HistoricalAwards=false), 
[#1854312](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1854312&HistoricalAwards=false), [#2103836](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2103836&HistoricalAwards=false), 
[#2103845](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2103845&HistoricalAwards=false), [#2138811](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2138811&HistoricalAwards=false), 
and [#2334945](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2334945&HistoricalAwards=false).
Any opinions, findings, and conclusions, or recommendations expressed in this material are those of the author(s) 
and do not necessarily reflect the views of the National Science Foundation. 

## Contact Info

Dr. Michela Taufer: mtaufer@utk.edu

Gabriel Laboy: glaboy@vols.utk.edu

Jay Ashworth: washwor1@vols.utk.edu