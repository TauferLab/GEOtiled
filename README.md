# GEOtiled Package Installation Instructions

## Initial Setup

1. **IF USING JETSTREAM**
   - Otherwise, skip this step and start from step 2 on your local machine.
   - Use Ubuntu and desired settings.
   - Create a volume and attach it to `/media/volume/sdb`.

2. **Clone Repository**
   - Execute: `git clone https://github.com/TauferLab/Src_SOMOSPIE`.

3. **Switch Branch**
   - Run: `git checkout jay_dev`.

4. **GitHub Sign In**
   - Sign into GitHub if you are on a fresh installation.

## Conda Installation

1. **Download Conda**
   - Visit the Conda website: [Anaconda Download](https://www.anaconda.com/download/)
   - Locate the Linux download link: [Anaconda for Linux](https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh)
   - Copy the download link address.
   - Navigate to the desired download directory in the terminal.
   - Use `wget` to download: `wget [link]`, where `[link]` is the copied URL.
   - Wait for the download to complete.

2. **Install Conda**
   - Execute the script: `bash ./xyz.sh`, replacing `xyz` with the script name.
   - Follow the installation process, typically agreeing to all prompts.
   - Restart the shell to apply the installation.

## Conda Environment Setup

- Create environment: `conda create -n geotiled -c conda-forge gdal=3.8.0`. 
  > Note: This process might take a significant amount of time.

## Additional Installations

1. **Install Editable Library**
   - Navigate to `GEOtiled_Refactor/somospie_lib`.
   - Execute: `pip install -e .`.

2. **Install Grass Libraries**
   - Run: `sudo apt-get install grass grass-doc`.

3. **Obtain Shapefiles**
   - Navigate to `GEOtiled_Refactor/misc_data`.
   - Download zip file from: [Google Drive Link](https://drive.google.com/file/d/1ODrkOTMM1_szjuG6AczHpebNfo1IZs7l/view?usp=drive_link)
   - Move the zip file to the working directory.
   - Extract files: `unzip shp_files.zip`.
   - Remove unnecessary files: `rm -r __MACOSX/` and `rm shp_files.zip`.

> Final step: Verify installation with `conda list`.

## Running Notebooks in Visual Studio Code

1. **Install Extensions**
   - Install Jupyter and Python extensions in Visual Studio Code.

2. **Select Kernel**
   - Open a notebook and select "select kernel" in the top right.
   - Choose the kernel named after the previously created Conda environment.

3. **Run Code Block**
   - Attempt to run a code block.
   - If prompted, install `ipykernel`.
