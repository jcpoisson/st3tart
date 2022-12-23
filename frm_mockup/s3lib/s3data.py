# -*- coding: utf-8 -*-


##########################################
#   Lib for collecting Sentinel-3 data   #
##########################################

# 04/11/2022.... JC.Poisson Creation

# Libs
import os
import numpy as np
from sentinelsat import SentinelAPI
import datetime as dt
from shapely import geometry
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84
import zipfile
import netCDF4 as nc4
import pandas as pd
import geopy.distance



# S3 data object
class S3_HYDRO_DATA:
    """
        sentinel-3 data object
    """

    def __init__(self) -> None:
        """
            constructor
        """
        self.__user = ""
        self.__password = ""
        self.__api_url = 'https://apihub.copernicus.eu/apihub'
        self.relative_orbit_number = None
        self.product_file = None
        self.nc_file = None
        self.df = None
        self.area = None
        self.lat_target = None
        self.lon_target = None
        self.api = None

    # destructor
    def reset_data(self) -> None:
        """
            reset the object
        """
        self.relative_orbit_number = None
        self.product_file = None
        self.nc_file = None
        self.df = None
        self.lon_target = None
        self.lat_target = None

    # set user and password for the api
    def set_api_user(self, user, password) -> int:
        """
            set user and password for the compernicus apihub
        """
        self.__user = user
        self.__password = password

        return 0

    # connect to the api
    def api_connect(self) -> int:
        """
            connect to the API
        """

        # connect to the api
        self.api = SentinelAPI(user=self.__user, password=self.__password, api_url=self.__api_url)

        return 0


    # set the taget coordinates
    def set_target(self, lon, lat) -> int:
        """
            set the targetted coordinates (lon, lat)
        """
        self.lon_target = lon
        self.lat_target = lat

        return 0

    # set a box area of interest from lat/lon points and distance from the center
    def set_box_area_from_point(self, km_from_center=1.0) -> geometry.Polygon:
        """
            set area of interest from a point (lon/lat) and a radius in km
        """
        # draw the box
        lat_min = geod.Direct(self.lat_target, self.lon_target, 180, km_from_center*1000.0)['lat2']
        lat_max = geod.Direct(self.lat_target, self.lon_target, 0, km_from_center*1000.0)['lat2']
        lon_min = geod.Direct(self.lat_target, self.lon_target, 270, km_from_center*1000.0)['lon2']
        lon_max = geod.Direct(self.lat_target, self.lon_target, 90, km_from_center*1000.0)['lon2']
        bbox = geometry.box(lon_min, lat_min, lon_max, lat_max)

        return bbox


    # download sentinel-3 products
    def download_s3_l2_products(self, start_date, end_date, s3sat, relativeorbitnumber, timelineess='Non Time Critical', km_from_latlon=1.0) -> int:
        """
            download for sentinel-3 l2 prodcuts on the copernicus scihub
        """

        # define the area
        if km_from_latlon == 0.0:
            self.area = None
        else:
            self.area = self.set_box_area_from_point(km_from_center=km_from_latlon)

        # set satellite
        filename_pattern = s3sat+'_SR_2_LAN__*'

        # earch by polygon, time, and Hub query keywords
        if self.area is None:
            products = self.api.query(date = (start_date, end_date),
                                 platformname = 'Sentinel-3',
                                 filename = filename_pattern,
                                 timeliness = timelineess,
                                 instrumentshortname= 'SRAL',
                                 productlevel = 'L2',
                                 relativeorbitnumber = relativeorbitnumber)
        else:
            products = self.api.query(area=self.area.wkt,
                                 date = (start_date, end_date),
                                 platformname = 'Sentinel-3',
                                 filename = filename_pattern,
                                 timeliness = timelineess,
                                 instrumentshortname= 'SRAL',
                                 productlevel = 'L2',
                                 relativeorbitnumber = relativeorbitnumber)

        # download all results from the search
        try:
            l_product = self.api.download_all(products, fail_fast=True)
        except:
            l_product = []
        if len(l_product) != 0:
            self.product_file = l_product[0][list(l_product[0].keys())[0]]['path']

        return 0


    # unzip file
    def __get_netcdf_file_from_product(self) -> int:
        """
            unzip the product file and get the netcdf
        """
        with zipfile.ZipFile(self.product_file, "r") as zip_ref:
            l_ncfiles = zip_ref.namelist()
            for ncfile in l_ncfiles:
                if "enhanced_measurement.nc" in ncfile:
                    self.nc_file = ncfile
                    zip_ref.extract(ncfile)
        if '/' in self.nc_file:
            self.nc_file = os.path.join(self.nc_file.split('/')[0], self.nc_file.split('/')[1])
        return 0


    # remove ncfile
    def remove_ncfile(self) -> int:
        """
            remove netcdf file
        """
        os.remove(self.nc_file)
        os.rmdir(os.path.split(self.nc_file)[0])
        return 0


    # remove product file
    def remove_product(self) -> int:
        """
            remove product file
        """
        os.remove(self.product_file)
        return 0


    # load s3data into a pandas dataframe
    def load_product_into_dataframe(self) -> pd.DataFrame:
        """
            unzip product and load netCDF file into pandas dataframe
        """
        # unzip product
        self.__get_netcdf_file_from_product()

        # read netcdf
        d_values = dict()
        with nc4.Dataset(self.nc_file) as nc:
            d_values['time'] = nc.variables['time_20_ku'][:]
            time_1 = nc.variables['time_01'][:]
            d_values['lat'] = nc.variables['lat_20_ku'][:]
            d_values['lon'] = nc.variables['lon_20_ku'][:]
            d_values['orbit'] = nc.variables['alt_20_ku'][:]
            d_values['range_ocog'] = nc.variables['range_ocog_20_ku'][:]
            d_values['range_ice'] = nc.variables['range_ice_20_ku'][:]
            # interp 1Hz data to 20 Hz 
            d_values['iono_cor'] = np.interp(d_values['time'], time_1, nc.variables['iono_cor_gim_01_ku'][:]) # 1 Hz -> 20 Hz
            d_values['dry_tropo'] = np.interp(d_values['time'], time_1, nc.variables['mod_dry_tropo_cor_meas_altitude_01'][:]) # 1 Hz -> 20 Hz
            d_values['wet_tropo'] = np.interp(d_values['time'], time_1, nc.variables['mod_wet_tropo_cor_meas_altitude_01'][:]) # 1 Hz -> 20 Hz
            d_values['solid_earth_tide'] = np.interp(d_values['time'], time_1, nc.variables['solid_earth_tide_01'][:]) # 1 Hz -> 20 Hz
            d_values['pole_tide'] = np.interp(d_values['time'], time_1, nc.variables['pole_tide_01'][:]) # 1 Hz -> 20 Hz

        # store into pandas dataframe
        df = pd.DataFrame.from_dict(d_values)

        # select data in the area
        lat_min = np.unique(self.area.boundary.xy[1]).min()
        lat_max = np.unique(self.area.boundary.xy[1]).max()
        df = df[(df['lat'] >= lat_min) & (df['lat'] <= lat_max)]

        return df

    # compute WSH for all beams
    def compute_wsh(self, df) -> int:
        """
            compute wsh from s3data w.r.t to the reference ellipsoid (WGS84)
        """
        # comute wsh for the 2 retrackers: ice and ocog
        df['wsh_ocog'] = df['orbit'] - (df['range_ocog'] + df['iono_cor'] + df['dry_tropo'] + df['wet_tropo'] + df['solid_earth_tide'] + df['pole_tide'])
        df['wsh_ice'] = df['orbit'] - (df['range_ice'] + df['iono_cor'] + df['dry_tropo'] + df['wet_tropo'] + df['solid_earth_tide'] + df['pole_tide'])
        self.df = df

        return 0

    # select the closest WSH
    def select_wsh(self) -> dict:
        """
            select the  closest WSH from the targetted coordinates.
            return a dict with the following keys [date, lat, lon, wsh, dist_nearest]
        """

        # find the nearst point from the targetted corrdinates
        diff_lat = np.abs(self.df['lat'] - self.lat_target)
        i_nearest = np.where(diff_lat == diff_lat.min())[0][0]
        dist_nearest = geopy.distance.geodesic((self.lat_target, self.lon_target), (self.df['lat'].values[i_nearest], self.df['lon'].values[i_nearest])).km

        # build the output list
        d_meas = dict()
        d_meas['time_sec'] = self.df['time'].values[i_nearest]
        d_meas['lat'] = self.df['lat'].values[i_nearest]
        d_meas['lon'] = self.df['lon'].values[i_nearest]
        d_meas['wsh_wgs84'] = self.df['wsh_ocog'].values[i_nearest]
        d_meas['dist'] = dist_nearest
        d_meas['time'] = dt.datetime(2000, 1, 1) + dt.timedelta(seconds=d_meas['time_sec'])

        return d_meas


# class s3a_orf
class S3_ORF:
    """
        sentinel-3 ORF object
    """

    # init object
    def __init__(self) -> None:
        """
            constructor
        """
        self.orf_file = ""
        self.df = None

    # load orf into a pandas dataframe
    def load_orf(self, orf_file) -> int:
        """
            read orf file and store it in self.df
        """
        if '/' in orf_file:
            self.orf_file = os.path.join(orf_file.split('/')[0], orf_file.split('/')[1])
        else:
            self.orf_file = orf_file

        # read orf file
        df = pd.read_table(self.orf_file, sep="\t", index_col=None, skiprows=64, header=None, names=['date', 'ncycle', 'npass', 'nrev', 'lon', 'lat'])
        df['date'] = np.array([dt.datetime.strptime(d, "%Y/%m/%d %H:%M:%S.%f") for d in df['date'].values])
        self.df = df

        return 0

    # get_cycle_dates
    def get_cycle_dates(self, cycle_number) -> list:
        """
            get first date and last date of a cycle number
        """
        df_cycle = self.df[self.df['ncycle'] == cycle_number]
        start_date = df_cycle.iloc[0]['date'].to_pydatetime()
        df_cycle_next = self.df[self.df['ncycle'] == (cycle_number + 1)]
        end_date = df_cycle_next.iloc[0]['date'].to_pydatetime()

        l_dates = [start_date, end_date]

        return l_dates

    # get_last_cycle_number
    def get_last_cycle_number(self) -> int:
        """
            get the last cycle number
        """
        end_cycle = self.df.iloc[-1]['ncycle']

        return end_cycle