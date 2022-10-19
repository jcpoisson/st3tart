# st3tart
This repository provide a tools for St3TART ESA project.

This tool allows to build WSH (Water Surface Height) over rivers from standard L2 Sentinel-3 STM products. The following actions are performed:
1. Download S3 data corresponding to a sqare area around a specified coordinate (lon / lat) from the Scihub
2. Unzip the product and read the NetCDF
3. Build th WSH for the concerned zone (the built WSH is dedicated for hydro, not for other surfaces)
4. The tool can delete the downloaded ptoduct

The tool is provided in 2 different formats:
- a jupyter notebook
- a simple python script

A account with username and password is required on the Copernicus Scihub platform to download the data.
You can create one here: https://scihub.copernicus.eu/dhus/#/home

