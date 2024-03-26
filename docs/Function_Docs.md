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
  List of arguments for a bash command. The list should be ordered the same way you would write the function in a command line (e.g., ["ls", "-lh", "~/"]).

#### Outputs
None
  The function has no outputs, with the exception of any outputs produced by the bash command.

#### Returns
None
  The function does not return a value.

#### Error States
RuntimeError
  Will raise a RuntimeError if Popen() returns an error and print out the error code, stdout, and stderr.
