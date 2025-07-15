![GEOtiled Logo](docs/assets/logo.png)

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

> Note: These have to be installed on your own

1. [Git](https://git-scm.com/downloads)
2. [Python 3](https://www.python.org/downloads/)

### Required Libraries

> Note: These will be installed with GEOtiled

1. matplotlib
2. setuptools
3. geopandas
4. notebook
5. tabulate
6. pandas
7. scipy
8. wheel
9. tqdm
10. gdal

### Required System Packages

> Note: Instructions on how to download are given

1. Performance Copilot
2. libgdal-dev
3. SAGA

## Installation

### Preconfiguration

1. Create a new Python virtual environment in a desired directory
   > Note: `<your_path>` should be replaced with your desired working directory

```
python3 -m venv <your_path>/geotiled_env
```

2. Activate the new virtual environment

```
source <your_path>/geotiled_env/bin/activate
```

3. Update pip

```
pip install --upgrade pip
```

4. Update apt-get

```
sudo apt-get update
```

### Install System Packages

1. Install libgdal-dev

```
sudo apt install libgdal-dev=3.8.4+dfsg-3ubuntu3
```

2. Install SAGA

```
sudo apt-get install saga
```

3. Install Performance Copilot

```
sudo apt-get install pcp-zeroconf
```

### Install GEOtiled

1. Ensure virutal environment is activated

```
source <your_path>/geotiled_env/bin/activate
```

2. Clone the repository in a desired working directory

```
git clone https://github.com/TauferLab/GEOtiled
```

3. Change to the geotiled directory

```
cd <your_path>/GEOtiled/geotiled
```

4. Install editable library

```
pip install -e .
```

### Debugging

The `plot_raster()` function may throw the following error:

`ImportError: cannot import name '_gdal_array' from 'osgeo'`

If so, run the following in the virtual environment to correct the issue:

```
pip install --no-cache --force-reinstall gdal[numpy]==3.8.4
```

## How to Use the Library

1. Ensure you are in the virtual environment

```
source <your_path>/geotiled_env/bin/activate
```

2. Place the following code snippet towards the top of any Python script to use GEOtiled functions

```
import geotiled
```

## Publications

Camila Roa, Paula Olaya, Ricardo Llamas, Rodrigo Vargas, and Michela Taufer. 2023. **GEOtiled: A Scalable Workflow
for Generating Large Datasets of High-Resolution Terrain Parameters.** _In Proceedings of the 32nd International Symposium
on High-Performance Parallel and Distributed Computing_ (HPDC '23). Association for Computing Machinery, New York, NY, USA,
311–312. [https://doi.org/10.1145/3588195.3595941](https://doi.org/10.1145/3588195.3595941)

## Copyright and License

Copyright (c) 2025, Global Computing Lab

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
