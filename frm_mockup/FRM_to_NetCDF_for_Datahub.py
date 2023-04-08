# -*- coding: utf-8 -*-


######################################################
#   save FRM into NetCDF compliant to the DataHub format
######################################################

# 08/04/2023.... JC.Poisson Creation


# Libraries
import datetime as dt
import numpy as np
import netCDF4 as nc4
import unidecode



def frm2netcdf(lat, lon, wse, wse_u, time, surface_type, geo_area, sensor_type, pf_type, pf_id, version) -> int:
    """
    Save FRM data to NetCDF format compliant with the DataHub format
    
    Input Parameters:
    - lat: array of latitude
    - lon: array of longitude
    - wse: array of water surface elevation w.r.t. the reference ellipsoid
    - wse_u: array of uncertainty related to water surface elevation
    - time: array of time in datetime format (UTC)
    - surface_type: 'IW' for InlandWaters, 'SI' for SeaIce, 'LI' for LandIce
    - geo_area: XXX_Yyy with XXX = 3 first letters of the country, yyy = 3 letters for the water body 
        (ex: FRA_Gar = Garonne River in France, FRA_Rhi = Rhine River in France, ITA_Por = Po River in Italy, etc...)
    - sensor_type: 'FIX' for fix sensor, 'MOV' for moving sensor
    - pf_type: 'ARB' for AirBorne, 'RIS' for river station, 'UAV' for drone, 'VES' for vessel, 'HUM' for human, 'MOO' for mooring, 'HLC' for Helicopter, 'DSB' for drifting surface buoy
    - pf_id: for vorteX-io micro-stations, the platform id is the name of the micro-station (as visible on the Maelstrom platform)
    - version: the dataset version (foramt: VX.Y)
    
    Output:
    - a integer representing the output status of the function 0 is OK and 1 is NOK
    
    More information: TD-9 FRM Data Hub Data Filename Convention and Format Specification Document (https://weboffice.noveltis.fr/Products/Files/DocEditor.aspx?fileid=31360)
    """
    
    # define the netcdf filename:
    date_first_meas = time[0].strftime('%Y%m%dT%H%M%S')
    date_last_meas = time[-1].strftime('%Y%m%dT%H%M%S')
    netcdf_filename = "{}_{}_{}_{}_{}_{}_{}_{}.nc".format(surface_type, geo_area, sensor_type, pf_type, pf_id, date_first_meas, date_last_meas, version)
    
    # create the netcdf file
    with nc4.Dataset(netcdf_filename, 'w') as nc:
        nc.title = "Fiducial Reference Measurement"
        nc.summary = "Fiducial Reference Measurement generated for the St3TART project"
        nc.institution = "To be filled"
        nc.contact = "To be filled"
        nc.project = "ESA St3TART Project"
        nc.date_update = dt.datetime.now().strftime('%Y%m%dT%H%M%S')
        nc.area = geo_area
        nc.platform_type = pf_type
        nc.platform_name = pf_id
        nc.sensor_type = sensor_type
        nc.sensor = "To be filled"
        nc.time_coverage_start = date_first_measurement
        nc.time_coverage_end = date_last_measurement
        nc.key_variable = "wse"
        nc.data_type = "FRM calculated"

        dim = nc.createDimension('index', len(time))
        
        # CRS
        crs_var = nc.createVariable('crs', 'i4')
        crs_var.grid_mapping_name = 'latitude_longitude'
        crs_var.longitude_of_prime_meridian = 0.0
        crs_var.semi_major_axis = 6378137.0
        crs_var.inverse_flattening = 298.257223563
        crs_var.crs_wkt = """GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]"""
        
        # time
        ref_2000 = dt.datetime(2000,1,1)
        time_2000 = [(t-ref_2000).total_seconds() for t in time]
        date = nc.createVariable('time', 'f8',('index'), fill_value=nc4.default_fillvals['f8'])
        date.standard_name = 'time'
        date.long_name = 'time (sec. since 2000-01-01)'
        date.units = 'seconds since 2000-01-01 00:00:00.0'
        date.calendar = 'gregorian'
        date.comments = 'The time is given in UTC'
        date[:] = time_2000
        
        # lat
        lat_var = nc.createVariable('lat', 'f8',('index') , fill_value=nc4.default_fillvals['f8'])
        lat_var.standard_name = 'latitude'
        lat_var.long_name = 'latitude'
        lat_var.units = 'degrees_east'
        lat_var.grid_mapping = 'crs'
        lat_var.comments = 'Positive latitude is North latitude, negative latitude is South latitude.'
        lat_var[:] = lat

        # lon
        lon_var = nc.createVariable('lon', 'f8',('index') , fill_value=nc4.default_fillvals['f8'])
        lon_var.standard_name = 'longitude'
        lon_var.long_name = 'longitude'
        lon_var.units = 'degrees_north'
        lon_var.grid_mapping = 'crs'
        lon_var.comments = 'East longitude relative to Greenwich meridian.'
        lon_var[:] = lon

        # wse
        wse_var = nc.createVariable('wse', 'f8',('index') , fill_value=nc4.default_fillvals['f8'])
        wse_var.standard_name = 'wse'
        wse_var.long_name = 'water surface elevation above reference ellipsoid'
        wse_var.units = 'm'
        wse_var.grid_mapping = 'crs'
        wse_var.comments = 'water surface elevation above reference ellipsoid provided as Fiducial Reference Measurement'
        wse_var[:] = wse

        # wse uncertainty
        wse_var_u = nc.createVariable('wse_uncertainty', 'f8',('index') , fill_value=nc4.default_fillvals['f8'])
        wse_var_u.standard_name = 'wse_uncertainty'
        wse_var_u.long_name = 'uncertainty related to water surface elevation above reference ellipsoid'
        wse_var_u.units = 'm'
        wse_var_u.grid_mapping = 'crs'
        wse_var_u.comments = 'uncertainty related to water surface elevation above reference ellipsoid provided as Fiducial Reference Measurement'
        wse_var_u[:] = wse_u



