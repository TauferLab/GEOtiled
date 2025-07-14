from osgeo import gdal
import subprocess

gdal.UseExceptions()

#################
### FUNCTIONS ###
#################
 
def bash(argv):
    """
    Runs a valid CLI command as a subprocess. The subprocess is run syncrhonously, so
    Python script execution will not continue until the subprocess is done.
    
    Parameters
    ----------
    argv : List[str]
        A list with the elements of the CLI command to execute. Order is important.
    """
    
    try:
        arg_seq = [str(arg) for arg in argv] # Convert all arguments in list into a string
        proc = subprocess.Popen(arg_seq, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait() # Synchronize
        stdout, stderr = proc.communicate() # Get standard output and error of command

        # Print error message if execution returned error
        if proc.returncode != 0:
            raise RuntimeError("'%s' failed, error code: '%s', stdout: '%s', stderr: '%s'" % (
                ' '.join(arg_seq), proc.returncode, stdout.rstrip(), stderr.rstrip()))
    except Exception as e:
        print(f"Error executing subprocess: {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def convert_file_format(input_file, output_file, new_format):
    """
    Converts a raster file to a different format. Currently only supports
    conversion between either the SAGA format or GeoTIFF format.
    
    Parameters
    ----------
    input_file : str
        Name of file in working directory or full path to file to convert.
    output_file : str
        Name of file in working directory or full path to file to write converted data to.
    new_format: str
        GDAL supported code with file format to convert file to. Options are `SAGA` or `GTiff`.
    """
    
    try:
        # Assume format is SAGA by default, else update to GeoTIFF format options
        translate_options = gdal.TranslateOptions(format="SAGA")
        if new_format == "GTiff":
            translate_options = gdal.TranslateOptions(format="GTiff", creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])

        # Perform conversion
        gdal.Translate(output_file, input_file, options=translate_options)
    except Exception as e:
        print(f"Error converting file format for {input_file}: {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def sequential_compute(input_file, method, parameter):
    """
    Computes a specified terrain parameter using a specified GIS library.
    
    Parameters
    ----------
    input_file : str
        Name or full path to file to compute the terrain parameter from.
        The file should contain elevation data.
    method : str
        GIS library to compute terrain parameter with. Options are `GDAL` or `SAGA`.
    parameter : str
        The name of the terrain parameter to compute. Options are `slope`, `aspect`, or `hillshade`.
    """
    
    try:
        if method == 'GDAL':
            # Configure parameters based on desired terrain parameter
            dem_options = gdal.DEMProcessingOptions(format='GTiff', alg='ZevenbergenThorne', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
            if parameter == 'aspect':
                dem_options = gdal.DEMProcessingOptions(zeroForFlat=False, format='GTiff', alg='ZevenbergenThorne', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
                
            # Computer terrain parameter
            gdal.DEMProcessing(f"{parameter}.tif", input_file, processing=parameter, options=dem_options)
        elif method == 'SAGA':
            # Convert input file to SAGA-compatible SGRD
            convert_file_format(input_file, input_file.replace('.tif','.sdat'), "SAGA")
            
            cmd = []
            # Build command based on desired terrain parameter
            if parameter == 'slope':
                cmd = ['saga_cmd', '-c=1', 'ta_morphometry', '0', '-ELEVATION', input_file.replace('.tif','.sgrd'), '-SLOPE', 'slope.sgrd']
            elif parameter == 'aspect':
                cmd = ['saga_cmd', '-c=1', 'ta_morphometry', '0', '-ELEVATION', input_file.replace('.tif','.sgrd'), '-ASPECT', 'aspect.sgrd']
            elif parameter == 'hillshade':
                cmd = ['saga_cmd', '-c=1', 'ta_lighting', '0', '-ELEVATION', input_file.replace('.tif','.sgrd'), '-SHADE', 'hillshade.sgrd']
            
            # Run as a subprocess
            bash(cmd)
            
            # Convert output file to GeoTIFF
            convert_file_format(f"{parameter}.sdat", f"{parameter}.tif", "GTiff")
    except Exception as e:
        print(f"Error computing {parameter} using {method}: {e}")
        return