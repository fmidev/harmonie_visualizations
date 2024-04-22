# harmonie_visualizations
This repository comprises Python-based Harmonie-Arome visualizations. The primary structure of the repository consists three main directories: visualization tools for MEPS (Python/), corresponding tools for MNWC (Python_MNWC/), and essential eccodes dependencies (eccodes/). 

The scripts in both MEPS and MNWC directories are requiring data in GRIB-format as an input. With the eccodes libraries, parameter values are fetched from the input data and stored into arrays or binary format depending on the script. Following the decoding of parameter values, a variety of charts are generated in PNG format.

Example of required package installations using a conda environment.:
```
conda create -n env_name -c conda-forge numpy scipy pyproj matplotlib cartopy eccodes python-eccodes
conda activate env_name
export ECCODES_DEFINITION_PATH=path_to_repository_dir/eccodes/locals/gl/definitions/:path_to_repository_dir/eccodes/locals/const/definitions/:path_to_repository_dir/eccodes/definitions/
```
The guidelines for executing the scripts can be found within the main directories.
