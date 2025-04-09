"""
SAGA Wrapper v0.0.1
GCLab 2024

Developed by Gabriel Laboy (@glaboy-vol)

Wrapper for SAGA 9.3.0 command line functions for use with GEOtiled.
Read more about SAGA command line functions at: https://saga-gis.sourceforge.io/saga_tool_doc/9.3.0/a2z.html
"""

import subprocess
import os

############################
### FUNCTION DEFINITIONS ###
############################

def bash(argv):
    """
    Executes a command in bash.

    This function acts as a wrapper to execute bash commands using the Python subprocess Popen method. 
    Commands are executed synchronously, stdout and stderr are captured, and errors can be raised.

    Parameters
    ----------
    argv : list
        List of arguments for a bash command. They should be in the order that you would arrange them in the command line (e.g., ["ls", "-lh", "~/"]).

    Raises
    ------
    RuntimeError
        Popen returns with an error if the passed bash function returns an error.
    """

    arg_seq = [str(arg) for arg in argv] # Convert all arguments in list into a string
    proc = subprocess.Popen(arg_seq, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait() # Synchronize
    stdout, stderr = proc.communicate() # Get standard output and error of command

    # Print error message if execution returned error
    if proc.returncode != 0:
        raise RuntimeError("'%s' failed, error code: '%s', stdout: '%s', stderr: '%s'" % (
            ' '.join(arg_seq), proc.returncode, stdout.rstrip(), stderr.rstrip()))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def build_param_path(input_file, param):
    """
    Create path to parameter file.

    Creates a path to a parameter file based off the path of the input file being computed from.
    It replaces the name of the folder the input file is in with the folder name 'param_tiles'.

    Parameters
    ----------
    input_file : str
        Full path to an input file.
    param : str
        The parameter name.
    """

    if (param == "channel_network") or (param == "drainage_basins"):
        return os.path.join(os.getcwd(), f"{param}_tiles", os.path.basename(input_file).replace(".tif",".shp"))
    else:
        return os.path.join(os.getcwd(), f"{param}_tiles", os.path.basename(input_file))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_curvature(cmd_prefix, input_file, parameter_list):
    """
    Compute all requested curvature parameters.

    SAGA command documentation: https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_morphometry_0.html

    Parameters
    ----------
    cmd_prefix : str
        First part of SAGA command to append rest of command to.
    input_file : str
        Full path to input elevation file to compute parameters from.
    parameter_list : List[str]
        List of paramters to compute.
    """
    
    # Add requested parameters to command
    cmd = cmd_prefix + ["ta_morphometry", "0", "-ELEVATION", input_file]
    if "slope" in parameter_list:
        param_path = build_param_path(input_file, "slope")
        cmd = cmd + ["-SLOPE", param_path]
    if "aspect" in parameter_list:
        param_path = build_param_path(input_file, "aspect")
        cmd = cmd + ["-ASPECT", param_path]
    if "profile_curvature" in parameter_list:
        param_path = build_param_path(input_file, "profile_curvature")
        cmd = cmd + ["-C_PROF", param_path]
    if "plan_curvature" in parameter_list:
        param_path = build_param_path(input_file, "plan_curvature")
        cmd = cmd + ["-C_PLAN", param_path]

    bash(cmd) # Run
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_hillshade(cmd_prefix, input_file):
    """
    Compute analytical hillshading.

    SAGA command documentation: https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_lighting_0.html

    Parameters
    ----------
    cmd_prefix : str
        First part of SAGA command to append rest of command to.
    input_file : str
        Full path to input elevation file to compute parameters from.
    """

    # Build rest of command
    cmd = cmd_prefix + ["ta_lighting", "0", "-ELEVATION", input_file]
    param_path = build_param_path(input_file, "hillshade")
    cmd = cmd + ["-SHADE", param_path]
    
    bash(cmd) # Run

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_convergence_index(cmd_prefix, input_file):
    """
    Compute convergence index.

    SAGA command documentation: https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_morphometry_1.html

    Parameters
    ----------
    cmd_prefix : str
        First part of SAGA command to append rest of command to.
    input_file : str
        Full path to input elevation file to compute parameters from.
    """

    # Build rest of command
    cmd = cmd_prefix + ["ta_morphometry", "1", "-ELEVATION", input_file]
    param_path = build_param_path(input_file, "convergence_index")
    cmd = cmd + ["-RESULT", param_path]

    bash(cmd) # Run

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_channel_network(cmd_prefix, input_file, parameter_list):
    """
    Compute channel network parameters.

    SAGA command documentation (CN & DB): https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_channels_5.html

    Parameters
    ----------
    cmd_prefix : str
        First part of SAGA command to append rest of command to.
    input_file : str
        Full path to input elevation file to compute parameters from.
    parameter_list : List[str]
        List of paramters to compute.
    """

    # Build rest of command
    cmd = cmd_prefix + ["ta_channels", "5", "-DEM", input_file]
    if "channel_network" in parameter_list:
        param_path = build_param_path(input_file, "channel_network")
        cmd = cmd + ["-SEGMENTS", param_path]
    if "drainage_basins" in parameter_list:
        param_path = build_param_path(input_file, "drainage_basins")
        cmd = cmd + ["-BASINS", param_path]

    bash(cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_channel_network_distance(cmd_prefix, input_file, parameter_list):
    """
    Compute channel network distance parameters.

    SAGA command documentation (CNBL & CND): https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_channels_3.html

    Parameters
    ----------
    cmd_prefix : str
        First part of SAGA command to append rest of command to.
    input_file : str
        Full path to input elevation file to compute parameters from.
    parameter_list : List[str]
        List of paramters to compute.
    """

    # Compute intermediary files if necessary
    cn_param_path = build_param_path(input_file, "channel_network")
    if not os.path.isfile(cn_param_path):
        compute_channel_network(cmd_prefix, input_file, ["channel_network"])
    
    # Build rest of command
    cmd = cmd_prefix + ["ta_channels", "3", "-ELEVATION", input_file, "-CHANNELS", cn_param_path]
    if "channel_network_distance" in parameter_list:
        param_path = build_param_path(input_file, "channel_network_distance")
        cmd = cmd + ["-DISTANCE", param_path]
    if "channel_network_base_level" in parameter_list:
        param_path = build_param_path(input_file, "channel_network_base_level")
        cmd = cmd + ["-BASELEVEL", param_path]

    bash(cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_valley_depth(cmd_prefix, input_file, parameter_list):
    """
    Compute valley depth and/or relative slope position parameters.

    SAGA command documentation: https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_morphometry_14.html

    Parameters
    ----------
    cmd_prefix : str
        First part of SAGA command to append rest of command to.
    input_file : str
        Full path to input elevation file to compute parameters from.
    parameter_list : List[str]
        List of paramters to compute.
    """
    
    # Build rest of command
    cmd = cmd_prefix + ["ta_morphometry", "14", "-DEM", input_file]
    if "valley_depth" in parameter_list:
        param_path = build_param_path(input_file, "valley_depth")
        cmd = cmd + ["-HU", param_path]
    if "relative_slope_position" in parameter_list:
        param_path = build_param_path(input_file, "relative_slope_position")
        cmd = cmd + ["-HO", param_path]

    bash(cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_total_catchment_area(cmd_prefix, input_file):
    """
    Compute total catchment area.

    SAGA command documentation: https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_hydrology_0.html

    Parameters
    ----------
    cmd_prefix : str
        First part of SAGA command to append rest of command to.
    input_file : str
        Full path to input elevation file to compute parameters from.
    """

    # Build rest of command
    param_path = build_param_path(input_file, "total_catchment_area")
    cmd = cmd_prefix + ["ta_hydrology", "0", "-ELEVATION", input_file, "-FLOW", param_path]
    
    bash(cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_specific_catchment_area(cmd_prefix, input_file):
    """
    Compute specific catchment area.

    SAGA command documentation: https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_hydrology_19.html
    This parameter is specifically computed for calculating topographic wetness index.

    Parameters
    ----------
    cmd_prefix : str
        First part of SAGA command to append rest of command to.
    input_file : str
        Full path to input elevation file to compute parameters from.
    """

    # Compute intermediary files if necessary
    tca_param_path = build_param_path(input_file, "total_catchment_area")
    if not os.path.isfile(tca_param_path):
        compute_total_catchment_area(cmd_prefix, input_file)
    
    # Build rest of command
    sca_param_path = build_param_path(input_file, "specific_catchment_area")
    cmd = cmd_prefix + ["ta_hydrology", "19", "-DEM", input_file, "-TCA", tca_param_path, "-SCA", sca_param_path]

    bash(cmd)
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_topographic_wetness_index(cmd_prefix, input_file):
    """
    Compute topographic wetness index.

    SAGA command documentation: https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_hydrology_20.html

    Parameters
    ----------
    cmd_prefix : str
        First part of SAGA command to append rest of command to.
    input_file : str
        Full path to input elevation file to compute parameters from.
    """

    # Compute intermediary files if necessary
    slope_param_path = build_param_path(input_file, "slope")
    if not os.path.isfile(slope_param_path):
        compute_curvature(cmd_prefix, input_file, ["slope"])
    
    sca_param_path = build_param_path(input_file, "specific_catchment_area")
    if not os.path.isfile(sca_param_path):
        compute_specific_catchment_area(cmd_prefix, input_file)

    # Build rest of command
    twi_param_path = build_param_path(input_file, "topographic_wetness_index")
    cmd = cmd_prefix + ["ta_hydrology", "20", "-SLOPE", slope_param_path, "-AREA", sca_param_path, "-TWI", twi_param_path]

    bash(cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_ls_factor(cmd_prefix, input_file):
    """
    Compute LS Factor.

    SAGA command documentation: https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_hydrology_22.html

    Parameters
    ----------
    cmd_prefix : str
        First part of SAGA command to append rest of command to.
    input_file : str
        Full path to input elevation file to compute parameters from.
    """

    # Compute intermediary files if necessary
    slope_param_path = build_param_path(input_file, "slope")
    if not os.path.isfile(slope_param_path):
        compute_curvature(cmd_prefix, input_file, ["slope"])
    
    sca_param_path = build_param_path(input_file, "specific_catchment_area")
    if not os.path.isfile(sca_param_path):
        compute_specific_catchment_area(cmd_prefix, input_file)

    # Build rest of command
    lsf_param_path = build_param_path(input_file, "ls_factor")
    cmd = cmd_prefix + ["ta_hydrology", "22", "-SLOPE", slope_param_path, "-AREA", sca_param_path, "-LS", lsf_param_path]

    bash(cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_parameters(input_file, parameter_list, saga_cores=1):
    """
    Computes user-specified parameters using SAGA command line API.

    Computes a user-specified list of terrain parameters from elevation data using SAGA.
    The SAGA calls are done using command line functions and called by the `bash()` function.

    Parameters
    ----------
    input_file : str
        Full path to input elevation file to compute parameters from.
    parameter_list : List[str]
        List of paramters to compute.
    saga_cores : int, optional
        Specify number of cores to use for multithreading with SAGA functions (default is 1).
        It is not recommended to set to more than 1 if using full multiprocessing of tiles.
    """
    
    # Build base of command line SAGA function
    cmd_base = ["saga_cmd", f"-c={saga_cores}"]
    
    # Slope, Aspect, and Plan & Profile Curvature
    if any(param in parameter_list for param in ["slope","aspect","plan_curvature","profile_curvature"]):
        compute_curvature(cmd_base, input_file, parameter_list)

    # Analytical Hillshading
    if "hillshade" in parameter_list:
        compute_hillshade(cmd_base, input_file)

    # Convergence Index
    if "convergence_index" in parameter_list:
        compute_convergence_index(cmd_base, input_file)

    # Channel Network & Drainage Basins
    if any(x in parameter_list for x in ["channel_network","drainage_basins"]):
        compute_channel_network(cmd_base, input_file, parameter_list)

    # Channel Network Distance & Base Level
    if any(x in parameter_list for x in ["channel_network_base_level","channel_network_distance"]):
        compute_channel_network_distance(cmd_base, input_file, parameter_list)

    # Valley Depth and Relative Slope Position
    if any(x in parameter_list for x in ["valley_depth","relative_slope_position"]):
        compute_valley_depth(cmd_base, input_file, parameter_list)

    # Total Catchment Area
    if "total_catchment_area" in parameter_list:
        compute_total_catchment_area(cmd_base, input_file)

    # Specific Catchment Area
    if "specific_catchment_area" in parameter_list:
        compute_specific_catchment_area(cmd_base, input_file)

    # Topographic Wetness Index
    if "topographic_wetness_index" in parameter_list:
        compute_topographic_wetness_index(cmd_base, input_file)

    # LS Factor
    if "ls_factor" in parameter_list:
        compute_ls_factor(cmd_base, input_file)
        
    # Closed Depressions
    # Unable to find computation - but can be computed by subtracting a "filled DEM" from the original DEM
    # Function to fill DEM: https://saga-gis.sourceforge.io/saga_tool_doc/7.8.2/ta_preprocessor_7.html

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_all_parameters(input_file):
    """
    Computes all parameters from the SAGA `ta_compound 0` command.

    SAGA command documentation: https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_compound_0.html
    It is important to note that this function runs substantially slower than individual functions.
    All computed terrain parameters are stored in the same directory as the input file.

    Parameters
    ----------
    input_file : str
        Full path to input elevation file to compute parameters from.
    """

    # Build and run command
    cmd = ["saga_cmd", "-c=1", "ta_compound", "0", "-ELEVATION", input_file, 
           "-SHADE", build_param_path(input_file, "hillshade"),
           "-SLOPE", build_param_path(input_file, "slope"),
           "-ASPECT", build_param_path(input_file, "aspect"),
           "-HCURV", build_param_path(input_file, "plan_curvature"),
           "-VCURV", build_param_path(input_file, "profile_curvature"),
           "-CONVERGENCE", build_param_path(input_file, "convergence_index"),
           "-SINKS", build_param_path(input_file, "closed_depressions"),
           "-FLOW", build_param_path(input_file, "total_catchment_area"),
           "-WETNESS", build_param_path(input_file, "topographic_wetness_index"),
           "-LSFACTOR", build_param_path(input_file, "ls_factor"),
           "-CHANNELS", build_param_path(input_file, "channel_network"),
           "-BASINS", build_param_path(input_file, "drainage_basins"),
           "-CHNL_BASE", build_param_path(input_file, "channel_network_base_level"),
           "-CHNL_DIST", build_param_path(input_file, "channel_network_distance"),
           "-VALL_DEPTH", build_param_path(input_file, "valley_depth"),
           "-RSP", build_param_path(input_file, "relative_slope_position")]
    bash(cmd)