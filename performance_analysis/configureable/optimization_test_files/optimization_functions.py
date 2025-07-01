from pathlib import Path
from osgeo import gdal
import multiprocessing
import subprocess
import glob
import os

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

def crop_pixels(input_file, output_file, window):
    """
    Crops a GeoTIFF to desired bounds.
    
    Parameters
    ----------
    input_file : str
        Name of file in working directory or full path to file to crop.
    output_file : str
        Name of file in working directory or full path to file to write cropped GeoTIFF to.
    window : List[int]
        Bounds to crop GeoTIFF to. Should be formatted as [x0,y0,x1,y1] where x0 and y0 represent
        starting indices in a two-dimensional array, and x1 and y1 are the length to crop in the 
        respective dimensions.
    """
    
    try:
        # Configure options and crop GeoTIFF
        translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
        gdal.Translate(output_file, input_file, options=translate_options)
    except Exception as e:
        print(f"Error cropping {input_file}: {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def sequential_crop(input_file, output_folder, column_length, row_length, buffer=10):
    """
    Crops a GeoTIFF file into a series of smaller tiles with an applied "buffer" region
    containing copied pixels from neighboring tiles.
    
    Parameters
    ----------
    input_file : str
        Name of file in working directory or full path to file to crop into tiles.
    output_folder : str
        Name of folder in working directory or full path to store cropped tiles in.
    column_length : int
        The number of columns each cropped tile will have.
    row_length : int
        The number of rows each cropped tile will have.
    buffer : int, optional
        The number of pixels to pad each tile with (default is 10).
    """
    
    try:
        # Create the output folder if it does not yet exist
        Path(output_folder).mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"Unable to create new folder: {e}")
        return

    try: 
        # Get the total number of rows and columns of the input file
        ds = gdal.Open(input_file, 0)
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        ds = None # Close file
    except Exception as e:
        print(f"Unable to open {input_file}: {e}")
        return

    try:
        # Begin cropping process
        tile_count = 0
        for i in range(0, rows, row_length):
            # If the next iteration were to exceed the number of rows of the file, crop it off to the correct length
            nrows = row_length
            if i + row_length > rows:
                nrows = rows - i

            for j in range(0, cols, column_length):
                # If the next iteration were to exceed the number of columns of the file, crop it off to the correct length
                ncols = column_length
                if j + column_length > cols:
                    ncols = cols - j

                # Create path to new tile file and set initial crop window
                tile_file = os.path.join(output_folder, "tile_{0:04d}.tif".format(tile_count))
                window = [j, i, ncols, nrows]

                # Set coords of upper left corner with the buffer included
                window[0] = window[0] - buffer
                window[1] = window[1] - buffer

                # Set coords of bottom right corner with the buffer included
                window[2] = window[2] + buffer*2
                window[3] = window[3] + buffer*2
                
                # Crop the new tile
                crop_pixels(input_file, tile_file, window)
                tile_count += 1 
    except Exception as e:
        print(f"Error cropping tiles: {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_pixels_parallel(window, items):
    """
    Crops a GeoTIFF to desired bounds. Meant to be used in tandem with the
    `parallel_crop()` function.
    
    Parameters
    ----------
    window : List[int]
        Bounds to crop GeoTIFF to. Should be formatted as [x0,y0,x1,y1] where x0 and y0 represent
        starting indices in a two-dimensional array, and x1 and y1 are the length to crop in the 
        respective dimensions.
    items : List[str]
        Contains the name of the input file to crop at index 0 and the output file to store
        the cropped data to at index 1.
    """
    
    try:
        input_file = items[0]
        output_file =items[1]
        
        # Set options and perform crop
        translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
        gdal.Translate(output_file, input_file, options=translate_options)
    except Exception as e:
        print(f"Error cropping {input_file}: {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def parallel_crop(input_file, output_folder, column_length, row_length, num_processes=None, buffer=10):
    """
    Crops a GeoTIFF file into a series of smaller tiles with an applied "buffer" region
    containing copied pixels from neighboring tiles. The tiles are cropped in parallel.
    
    Parameters
    ----------
    input_file : str
        Name of file in working directory or full path to file to crop into tiles.
    output_folder : str
        Name of folder in working directory or full path to store cropped tiles in.
    column_length : int
        The number of columns each cropped tile will have.
    row_length : int
        The number of rows each cropped tile will have.
    num_processes : int, optional
        The number of concurrent processes to use for computation (default is None).
        If set to `None`, the maximum number of available cores on the machine are used.
    buffer : int, optional
        The number of pixels to pad each tile with (default is 10).
    """

    try:
        # Create the output folder if it does not yet exist
        Path(output_folder).mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"Unable to create new folder: {e}")
        return

    try:
        # Get the total number of rows and columns of the input file
        ds = gdal.Open(input_file, 0)
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        ds = None # close file
    except Exception as e:
        print(f"Unable to open {input_file}: {e}")
        return

    try:
        # Begin cropping process
        items = []
        tile_count = 0
        for i in range(0, rows, row_length):
            # If the next iteration were to exceed the number of rows of the file, crop it off to the correct length
            nrows = row_length
            if i + row_length > rows:
                nrows = rows - i

            for j in range(0, cols, column_length):
                # If the next iteration were to exceed the number of columns of the file, crop it off to the correct length
                ncols = column_length
                if j + column_length > cols:
                    ncols = cols - j

                # Create path to new tile file and set initial crop window
                window = [j, i, ncols, nrows]

                # Set coords of upper left corner with the buffer included
                window[0] = window[0] - buffer
                window[1] = window[1] - buffer

                # Set coords of bottom right corner with the buffer included
                window[2] = window[2] + buffer*2
                window[3] = window[3] + buffer*2
                
                # Add window and other relevant variables to items
                tile_file = os.path.join(output_folder, "tile_{0:04d}.tif".format(tile_count))
                items.append((window, [input_file,tile_file]))
                tile_count += 1 
    except Exception as e:
        print(f"Error configuring windows: {e}")
        return

    try:
        # Determine the number of concurrent processes to use if not set
        if num_processes is None:
            num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
            
        # Begin concurrent computation
        pool = multiprocessing.Pool(processes=num_processes)
        pool.starmap(crop_pixels_parallel, items)
    except Exception as e:
        print(f"Error cropping tiles: {e}")
        return
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def gdal_compute_parameters_parallel(input_file, items):
    """
    Computes a terrain parameter using algorithms implemented by the GDAL library.
    Meant to be used in tandem with `parallel_compute()` for parallel processing.

    Parameters
    ----------
    input_file : str
        Name or full path to file to compute the terrain parameter from.
        The file should contain elevation data.
    items : List[]
        A list that should only contain the name of the terrain parameter to compute.
    """

    try: 
        # Extract parameter to compute
        parameter = items[0]
        
        # Configure path to output file
        output_file = os.path.join(f"{parameter}_tiles", os.path.basename(input_file))
        
        # Configure compute options based on desired terrain parameter
        dem_options = gdal.DEMProcessingOptions(format='GTiff', alg='ZevenbergenThorne', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
        if parameter == 'aspect':
            dem_options = gdal.DEMProcessingOptions(zeroForFlat=False, format='GTiff', alg='ZevenbergenThorne', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])

        # Compute
        gdal.DEMProcessing(output_file, input_file, processing=parameter, options=dem_options)
    except Exception as e:
        print(f"Error computing terrain parameter: {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def saga_compute_parameters_parallel(input_file, items):  
    """
    Computes a terrain parameter using algorithms implemented by the SAGA library.
    Meant to be used in tandem with `parallel_compute()` for parallel processing.
    Handles reprojection to and from SAGA's file format, SGRD.

    Parameters
    ----------
    input_file : str
        Name or full path to file to compute the terrain parameter from.
        The file should contain elevation data.
    items : List[]
        A list that should only contain the name of the terrain parameter to compute.
    """

    try: 
        # Extract parameter to compute
        parameter = items[0]

        # Convert input file into the SAGA-compatible SGRD format
        convert_file_format(input_file, input_file.replace('.tif','.sdat'), 'SAGA')

        # Configure path to output file
        output_file = os.path.join(f"{parameter}_tiles", os.path.basename(input_file))
        
        # Configure command and compute in a subprocess
        if parameter == 'hillshade':
            cmd = ["saga_cmd", "-c=1", "ta_lighting", "0", "-ELEVATION", input_file.replace('.tif','.sgrd'), "-SHADE", output_file.replace('.tif','.sgrd')]
        elif parameter == 'slope':
            cmd = ["saga_cmd", "-c=1", "ta_morphometry", "0", "-ELEVATION", input_file.replace('.tif','.sgrd'), "-SLOPE", output_file.replace('.tif','.sgrd')]
        elif parameter == 'aspect':
            cmd = ["saga_cmd", "-c=1", "ta_morphometry", "0", "-ELEVATION", input_file.replace('.tif','.sgrd'), "-ASPECT", output_file.replace('.tif','.sgrd')]
        bash(cmd)
        
        # Convert output file back into a GeoTIFF
        convert_file_format(output_file.replace('.tif','.sdat'), output_file, 'GTiff')
    except Exception as e:
        print(f"Error computing terrain parameter: {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def parallel_compute(input_folder, method, parameter, num_processes=None):
    """
    Computes a desired terrain parameter from a collection of GeoTIFF files in parallel
    using a desired GIS library of choice for computation.
    
    Parameters
    ----------
    input_folder : str
        Name or full path to folder containing input files to compute the terrain parameter from.
        The files in the folder should contain elevation data.
    method : str
        GIS library to compute terrain parameter with. Options are `GDAL` or `SAGA`.
    parameter : str
        The name of the terrain parameter to compute. Options are `slope`, `aspect`, or `hillshade`.
    num_processes : int, optional
        The number of concurrent processes to use for computation (default is None).
        If set to `None`, the maximum number of available cores on the machine are used.
    """
    
    try:
        # Create output folder for the desired terrain parameter
        Path(f"{parameter}_tiles").mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"Error creating terrain parameter folder: {e}")
        return
    
    try:
        # Get all input files
        input_files = glob.glob(os.path.join(input_folder, "*.tif"))
        
        # Create items for multiprocessing pool
        items = []
        for input_file in input_files:
            items.append((input_file, [parameter]))
    except Exception as e:
        print(f"Error creating items for mutltiprocessing: {e}")
        return
    
    try:
        # Determine the number of concurrent processes to use for computation
        if num_processes is None:
            num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
            
        # Begin computation based off desired library to use for computation
        pool = multiprocessing.Pool(processes=num_processes)
        if method == 'GDAL':
            pool.starmap(gdal_compute_parameters_parallel, items)
        elif method == 'SAGA':
            pool.starmap(saga_compute_parameters_parallel, items)
    except Exception as e:
        print(f"Error computing terrain parameter for tiles {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def mosaic_average(input_folder, output_file):
    """
    Mosaics together GeoTIFF files into a single GeoTIFF. Overlapping pixels (based on spatial location)
    are averaged together.
    
    Parameters
    ----------
    input_folder : str
        Name or path to folder containing GeoTIFF files to mosaic.
    output_file : str
        Name of path to file to write contents of merged raster files tos.
    """
    
    try:
        # Get files from input folder to merge together
        input_files = glob.glob(os.path.join(input_folder, "*.tif"))
    except Exception as e:
        print(f"Error getting input files: {e}")
        return

    try:
        # Create VRT file to point to files to merge
        vrt_file = 'merged.vrt'
        vrt = gdal.BuildVRT(vrt_file, input_files)
        vrt = None  # closes file

        with open(vrt_file, "r") as f:
            contents = f.read()

        # Extract nodata value from VRT
        if "<NoDataValue>" in contents:
            nodata_value = contents[contents.index("<NoDataValue>") + len(
                "<NoDataValue>"): contents.index("</NoDataValue>")]  # To add averaging function
        else:
            nodata_value = 0

        # Create XML that will average overlapping pixels during mosaic
        code = '''band="1" subClass="VRTDerivedRasterBand">
<PixelFunctionType>average</PixelFunctionType>
<PixelFunctionLanguage>Python</PixelFunctionLanguage>
<PixelFunctionCode><![CDATA[
import numpy as np

def average(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize,raster_ysize, buf_radius, gt, **kwargs):
    data = np.ma.array(in_ar, mask=np.equal(in_ar, {}))
    np.mean(data, axis=0, out=out_ar, dtype="float32")
    mask = np.all(data.mask,axis = 0)
    out_ar[mask] = {}
]]>
</PixelFunctionCode>'''.format(nodata_value, nodata_value)

        sub1, sub2 = contents.split('band="1">', 1)
        contents = sub1 + code + sub2

        # Write averaging function to VRT file
        with open(vrt_file, "w") as f:
            f.write(contents)
    except Exception as e:
        print(f"Error building VRT: {e}")
        return

    try:
        # Merge together the files
        cmd = ["gdal_translate", "-co", "COMPRESS=LZW", "-co", "TILED=YES", "-co", 
            "BIGTIFF=YES", "--config", "GDAL_VRT_ENABLE_PYTHON", "YES", vrt_file, output_file]
        bash(cmd)
    except Exception as e:
        print(f"Error mosaicking files: {e}")
        return

    try:
        # Delete VRT file
        os.remove(vrt_file)
    except Exception as e:
        print(f"Error deleting VRT file: {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_buffer_regions(input_folder, output_folder, num_processes=None, buffer=10):
    """
    Crops buffer pixels off of all GeoTIFF files in a specified folder.
    The cropping is done in parallel.

    Parameters
    ----------
    input_folder : str
        Name or path of the folder where files to crop are stored.
    output_folder : str
        Name or path of the folder where cropped files are stored.
    num_processes : int
        Number of concurrent processes to use for cropping (default is None).
        It will be set to the maximum number of available cores if set to `None`.
    buffer : int
        Size of buffer region of files to crop off (default is 10).
    """

    try:
        # Create output folder if it does not exist
        Path(output_folder).mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"Error creating output folder: {e}")
        return

    try:
        # Get all input files
        input_files = glob.glob(os.path.join(input_folder, "*.tif"))
    except Exception as e:
        print(f"Error gathering input files: {e}")
        return
    
    try:
        items = []
        for input_file in input_files:
            file_name = os.path.basename(input_file)
            output_file = os.path.join(output_folder, file_name)

            # Create window and crop new file without buffer pixels
            ds = gdal.Open(input_file, 0)
            cols = ds.RasterXSize
            rows = ds.RasterYSize
            window = [buffer, buffer, cols-(buffer*2), rows-(buffer*2)]
            items.append((window, [input_file, output_file]))
            ds = None # Close file
    except Exception as e:
        print(f"Error creating windows: {e}")
        return

    try:
        # Determine number of concurrent processes to use
        if num_processes is None:
            num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
            
        # Concurrently compute tiles
        pool = multiprocessing.Pool(processes=num_processes)
        pool.starmap(crop_pixels_parallel, items)
    except Exception as e:
        print(f"Error cropping files: {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def mosaic_concat(input_folder, output_file, num_processes=None, buffer=10):
    """
    Mosaic together GeoTIFF files. It first removes buffer regions of the input files (places where GeoTIFF files
    overlap) in parallel, and then concatenates the unbuffered files to create the mosaic.
    
    Parameters
    ----------
    input_folder : str
        Name or path of the folder where files to crop are stored.
    output_file : str
        Name or path of the folder where cropped files are stored.
    num_processes : int
        Number of concurrent processes to use for cropping (default is None).
        It will be set to the maximum number of available cores if set to `None`.
    buffer : int
        Size of buffer region of files to crop off (default is 10).
    """

    try:
        # Remove buffer pixels from input files
        crop_buffer_regions(input_folder, "unbuffered_tiles", num_processes, buffer=buffer)
    except Exception as e:
        print(f"Error cropping buffered files: {e}")
        return
    
    try:
        # Check for valid folder and put input files into list
        input_files = glob.glob(os.path.join("unbuffered_tiles", "*.tif"))
    except Exception as e:
        print(f"Error gathering unbuffered files: {e}")
        return

    try:
        # Create the VRT file to point at the unbuffered tiles
        vrt_file = 'merged.vrt'
        vrt = gdal.BuildVRT(vrt_file, input_files)
        
        # Mosaic together tiles
        translate_options = gdal.TranslateOptions(creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"])
        gdal.Translate(output_file, vrt, options=translate_options)
        vrt = None  # close file
    except Exception as e:
        print(f"Error mosaicking files: {e}")
        return

    try:
        # Delete VRT file
        os.remove(vrt_file)
    except Exception as e:
        print(f"Error deleting VRT file: {e}")
        return