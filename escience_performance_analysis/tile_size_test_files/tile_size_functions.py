from osgeo import gdal
import subprocess

gdal.UseExceptions()

#################
### FUNCTIONS ###
#################

def bash(argv):
    arg_seq = [str(arg) for arg in argv] # Convert all arguments in list into a string
    proc = subprocess.Popen(arg_seq, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait() # Synchronize
    stdout, stderr = proc.communicate() # Get standard output and error of command

    # Print error message if execution returned error
    if proc.returncode != 0:
        raise RuntimeError("'%s' failed, error code: '%s', stdout: '%s', stderr: '%s'" % (
            ' '.join(arg_seq), proc.returncode, stdout.rstrip(), stderr.rstrip()))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def sequential_compute(input_file, method, parameter):
    if method == 'GDAL':
        if parameter == 'aspect':
            dem_options = gdal.DEMProcessingOptions(zeroForFlat=False, format='GTiff', alg='ZevenbergenThorne', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
            gdal.DEMProcessing(f"{parameter}.tif", input_file, processing=parameter, options=dem_options)
        else:
            dem_options = gdal.DEMProcessingOptions(format='GTiff', alg='ZevenbergenThorne', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
            gdal.DEMProcessing(f"{parameter}.tif", input_file, processing=parameter, options=dem_options)
    elif method == 'SAGA':
        if parameter == 'slope':
            cmd = ['saga_cmd', '-c=1', 'ta_morphometry', '0', '-ELEVATION', input_file, '-SLOPE', 'slope.tif']
            bash(cmd)
        elif parameter == 'aspect':
            cmd = ['saga_cmd', '-c=1', 'ta_morphometry', '0', '-ELEVATION', input_file, '-ASPECT', 'aspect.tif']
            bash(cmd)
        elif parameter == 'hillshade':
            cmd = ['saga_cmd', '-c=1', 'ta_lighting', '0', '-ELEVATION', input_file, '-SHADE', 'hillshade.tif']
            bash(cmd)