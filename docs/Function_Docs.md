# GEOtiled Function Documentation


## Dictionary Codes


### USGS_DATASET_CODES

Codes specifying which dataset to download GEOTiff elevation files for off the [USGS webpage](https://apps.nationalmap.gov/downloader/#/).

Currently supported codes and their correlating dataset:
* 30m (National Elevation Dataset (NED) 1 arc-second Current)
* 10m (National Elevation Dataset (NED) 1/3 arc-second Current)


### TERRAIN_PARAMETER_CODES

Codes specifying terrain parameters to be computed by GEOtiled.

Currently supported codes and their correlating parameter:
* SLP (slope)
* ASP (aspect)
* HLSD (hillshade)


### SHAPE_FILE_CODES

Codes specifying a US state or territory for which to fetch the shape file for off the [USGS webpage](https://apps.nationalmap.gov/downloader/#/).

Codes specifying the region and their correlating parameter:
* AL (Alabama)
* AK (Alaska)
* AR (Arkansas)
* AZ (Arizona)
* CA (California)
* CO (Colorado)
* CT (Connecticut)
* DC (District of Columbia)
* DE (Delaware)
* FL (Florida)
* GA (Georgia)
* GU (Guam)
* HI (Hawaii)
* ID (Idahoo)
* IL (Illinois)
* IN (Indiana)
* IA (Iowa)
* KS (Kansas)
* KY (Kentucky)
* LA (Louisiana)
* MA (Massachusetts)
* MD (Maryland)
* ME (Maine)
* MI (Michigan)
* MN (Minnesota)
* MO (Missouri)
* MP (Northern Mariana Islands)
* MS (Mississippi)
* MT (Montana)
* NC (North Carolina)
* ND (North Dakota)
* NE (Nebraska)
* NH (New Hamshire)
* NJ (New Jersey)
* NM (New Mexico)
* NV (Nevada)
* NY (New York)
* OH (Ohio)
* OK (Oklahoma)
* OR (Oregon)
* PA (Pennsylvania)
* PR (Puerto Rico)
* RI (Rhode Island)
* SC (South Carolina)
* SD (South Dakota)
* TN (Tennessee)
* TX (Texas)
* UT (Utah)
* VA (Virginia)
* VI (Virgin Islands)
* VT (Vermont)
* WA (Washington)
* WI (Wisconsin)
* WV (West Virginia)
* WY (Wyoming)


## Functions


### `bash(argv)`

#### Executes a command in bash.
This function acts as a wrapper to execute bash commands using the subprocess.Popen() method. Commands are executed synchronously, and errors are caught and raised.

#### Required Parameters
argv : List
- List of arguments for a bash command.
- The list should be ordered the same way you would write the function in a command line (e.g., ["ls", "-lh", "~/"]).

#### Outputs
Command Outputs
- Any output(s) produced by the bash command.

#### Returns
None

#### Error States
RuntimeError
- Will raise a RuntimeError if Popen() returns an error and print out the error code, stdout, and stderr.


### `download_file(url, folder, pbar)`

#### Downloads a file found at a URL and stores it in the specified folder.
This function facilitates download of a single file and is used by the `download_files` function.

#### Required Parameters
url : str
* String specifying the URL to download the file from.
folder : str
* String specifying folder in 'data' directory where file will be stored.
pbar : tqdm object
* Reference to a tqdm progress bar to indiciate download progress.

#### Outputs
File
* Downloaded file in the specified folder.

#### Returns
int
* Integer specifying the number of bytes downloaded.

#### Notes
* This function is designed to be used specifically by `download_files`
* If the file being downloaded already exists, no download occurs, and the function returns 0.


### `get_file_size(url)`

#### Return the size of a file, in bytes, retrieved from a URL.
This function uses the requests.head() function to read and return the 'Content-Length' header specified by a file found at a URL. This function is intended to be used with `download_files` to calculate download size before file downloads begin.

#### Required Parameters
url : str
* String representing URL where file is found.

#### Outputs
None

#### Returns
int
* Size of the file specified at the URL in bytes. Returns 0 if file size can't be determined.


### `set_data_directory(path)`

#### Sets the path where data computed by GEOtiled will be stored.
This function sets the 'data' directory where data will be searched for and generated in by GEOtiled functions.

#### Required Parameters
path : str
* A string that specifies working directory.

#### Outputs
None

#### Returns
None

#### Notes
* It's good practice to run this function before executing anything else in the workflow.
* If not set, data will be searched for and stored in the working directory.


