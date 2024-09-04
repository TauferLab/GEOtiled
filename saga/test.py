import geotiledsaga as gts
import os

PARAMETERS = ['SLP','ASP','HLD','PLC','PFC','CI']

store_path = '/media/volume/geotiled-saga/tile_compute_test_ca'
elev_file = '/media/volume/geotiled-saga/tile_compute_test_ca/elevation.tif'

gts.set_working_directory(store_path)

individual_results_file = os.path.join(store_path, 'individual_results.csv')
f = open(individual_results_file, 'w')
f.write('parameter,file_name,execution_time,peak_mem_usage\n')
f.close()

param_codes = gts.__get_codes("parameter")
for parameter in PARAMETERS: 
    # Compute for each parameter
    param = param_codes[parameter]
    gts.__compute_params_saga(elev_file, [param, 1, 'ctt', False])

# gdal.UseExceptions()

# file1 = '/media/volume/geotiled-saga/tile_compute_test_de/mosaic.tif'
# file2 = '/media/volume/geotiled-saga/tile_compute_test_tn/elevation.tif'
# file3 = '/media/volume/geotiled-saga/tile_compute_test_ca/elevation.tif'

# def get_nodata_percentage(file):
#     matrix_data = gdal.Open(file)
#     array_data = matrix_data.ReadAsArray()
#     nodata_percentage = (len(array_data[array_data == -999999.0]) / (len(array_data) * len(array_data[0]))) * 100
#     print('%.2f' % nodata_percentage)

# get_nodata_percentage(file1)
# get_nodata_percentage(file2)
# get_nodata_percentage(file3)
