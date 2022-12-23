# -*- coding: utf-8 -*-


##############################################
#   Lib of tools for computing in-istu FRM   #
##############################################

# 04/11/2022.... JC.Poisson Creation

# libs
import datetime as dt
from shapely.geometry import Point, LineString
from fiona.drvsupport import supported_drivers
from scipy.spatial import cKDTree
import geopy.distance
import numpy as np
import pandas as pd
import geopandas as gpd
import netCDF4 as nc4
import scipy.optimize as opt
import scipy.signal as sgl

# ===============================================================
#    Objects
# ===============================================================


# centerline object
class CurvAbsc:
    """
        curvilinear abscissa object with different methods
    """

    def __init__(self):
        """
            constructor
        """
        self.__gdf = None
        self.start_loc = None
        self.end_loc = None
        self.coordinates = None
        self.distances = None
        self.indices = None


    def load_kml(self, kml_file):
        """
            read kml file to create 
        """
        supported_drivers['LIBKML'] = 'rw'
        self.__gdf =  gpd.read_file(kml_file, driver='LIBKML')


    def __set_start_loc(self, lon, lat):
        """
            set location of the starting point of the curvilinear abscissa
        """
        self.start_loc = Point(lon, lat)


    def __set_end_loc(self, lon, lat):
        """
            set location of the end point of the curvilinear abscissa
        """
        self.end_loc = Point(lon, lat)


    def compute_curvilinear_abscissa(self, lon1, lat1, lon2, lat2) -> int:
        """
            compute curvilinear abscissa from start and last point
        """
        
        # set starting and ending points
        self.__set_start_loc(lon1, lat1)
        self.__set_end_loc(lon2, lat2)

        # find the closest point in the kml from the starting point
        lon, lat = self.__gdf.geometry[0].coords.xy
        tmp_lon_lat = np.vstack((lon, lat))
        tree_path = cKDTree(tmp_lon_lat.T)
        self.__set_start_loc(lon1, lat1)
        self.__set_end_loc(lon2, lat2)
        first_point = tree_path.query([lon1, lat1], k=1)
        last_point = tree_path.query([lon2, lat2], k=1)
        if first_point[1] < last_point[1]:
            lon_lat = tmp_lon_lat[:, first_point[1]: last_point[1]+1]
        else:
            lon_lat = tmp_lon_lat[:, first_point[1]: last_point[1]-1: -1]
        self.coordinates = lon_lat

        # compute the distance along the coordinates
        len_coordinates = lon_lat.shape[1]
        dist = [geopy.distance.geodesic((lon_lat[1, i], lon_lat[0, i]), (lon_lat[1, i+1], lon_lat[0, i+1])).km for i in range(len_coordinates-1)]
        distances = np.array([0]+dist)
        self.distances = np.cumsum(distances)*1000.0
        self.indices = np.arange(len(self.distances))
        return 0


    def proj_to_curvilinear_abscissa(self, lons, lats) -> np.array:
        """
            proj lat/lon to the curvvilinear abscissa
        """

        # create Linestring from coordinates
        line = LineString(self.coordinates.T)

        # interpolate points to the Linestring
        lon_proj = []
        lat_proj = []
        distances_proj = []
        tree_path = cKDTree(self.coordinates.T)
        for longi, lati in zip(lons, lats):
            p_proj = line.interpolate(line.project(Point(longi, lati), normalized=False))
            # find the nearest point in the theoretical centerline from the projected point
            d, i_nearest = tree_path.query([p_proj.xy[0][0], p_proj.xy[1][0]], k=1)
            # find the second closest point to the projected point lying on the centerline
            ls_before = LineString([self.coordinates.T[i_nearest-1], self.coordinates.T[i_nearest]])
            ls_after = LineString([self.coordinates.T[i_nearest], self.coordinates.T[i_nearest+1]])
            if ls_before.distance(p_proj) < ls_after.distance(p_proj):
                dist_proj = self.distances[i_nearest] - (geopy.distance.geodesic((self.coordinates[0, i_nearest], self.coordinates[1, i_nearest]), (p_proj.xy[0][0], p_proj.xy[1][0])).km * 1000.0)
            else:
                dist_proj = self.distances[i_nearest] + (geopy.distance.geodesic((self.coordinates[0, i_nearest], self.coordinates[1, i_nearest]), (p_proj.xy[0][0], p_proj.xy[1][0])).km * 1000.0)

            # dist_proj = np.interp(p_proj.xy[0][0], np.array([self.coordinates[0, i_nearest[0]] , self.coordinates[0, i_nearest[1]]]),  np.array([self.distances[i_nearest[0]] , self.distances[i_nearest[1]]]))
            # store results
            lon_proj.append(p_proj.xy[0][0])
            lat_proj.append(p_proj.xy[1][0])
            distances_proj.append(dist_proj)
        
        coord_proj = np.vstack((np.array(lon_proj), np.array(lat_proj)))
        return coord_proj, np.array(distances_proj)


# class drone data
class DRONE_DATA:
    """
        drone data object with different methods
    """

    # init method
    def __init__(self):
        """
            constructor
        """
        self.df = None
        self.nc_file = None
        self.reference_time = None
        self.reference_time_sec = None
    
    # load drone data
    def load_drone_file(self, nc_file) -> int:
        """
            load netcdf file containing drone data
        """
        # read nc file
        self.nc_file = nc_file
        d_drone = dict()
        with nc4.Dataset(nc_file) as nc:
            d_drone['time_sec'] = nc.variables['time'][:]
            d_drone['lat'] = nc.variables['latitude'][:]
            d_drone['lon'] = nc.variables['longitude'][:]
            d_drone['wsh_wgs84'] = nc.variables['wsh_wgs84'][:]

        # add a time in datetime format
        time = np.array([dt.datetime(2000, 1, 1) + dt.timedelta(seconds=t) for t in d_drone['time_sec']])
        d_drone['time'] = time

        # store into pandas dataframe
        self.df = pd.DataFrame(data=d_drone)

        return 0


    def project_on_centerline(self, lon_proj, lat_proj, absc_curv) -> int:
        """
            modify the dataframe to account for the projection on the centerline
        """
        self.df['lon'] = lon_proj
        self.df['lat'] = lat_proj
        self.df['absc_curv'] = absc_curv
        return 0


    # correct from water elevation evolution
    def correct_from_water_evolution(self, avg_speed, station_time_sec, station_wsh, station_absc) -> int:
        """
            correct the drone profile from the water evolution during the flight
        """

        # find the nearest point from the station
        diff_drone_absc = station_absc - self.df['absc_curv'].values
        i_drone_nearest = np.where(np.abs(diff_drone_absc) == np.min(np.abs(diff_drone_absc)))[0][0]
        t_sec_drone_nearest = self.df.loc[i_drone_nearest, 'time_sec']

        # find the corresponding wsh of the station
        diff_tsec_station = station_time_sec - t_sec_drone_nearest
        i_station_nearest = np.where(np.abs(diff_tsec_station) == np.min(np.abs(diff_tsec_station)))[0][0]
        wsh_station_reference = station_wsh[i_station_nearest]

        # compute delta wsh in the station time series w.r.t. the wsh reference
        diff_wsh_station = station_wsh - wsh_station_reference

        # compute the propagation time for all points of the drone data
        time_drone_propagation = diff_drone_absc / avg_speed

        # transform this time in time at the station and interpolate the wsh evolution correction to apply
        time_station_corr = self.df['time_sec'].values + time_drone_propagation
        wsh_evolution_corr = np.interp(time_station_corr, station_time_sec, diff_wsh_station)
        self.df['wsh_evo_corrected'] = self.df['wsh_wgs84'] + wsh_evolution_corr
        self.reference_time = pd.to_datetime(self.df.loc[i_drone_nearest, 'time']).date()
        self.reference_time_sec = t_sec_drone_nearest
        self.reference_height_station = wsh_station_reference

        return 0
    
    # compute h difference from a reference point on the curvilinear abscissa
    def compute_h_diff_from_ref_absc(self, ref_absc) -> pd.DataFrame():
        """
            compute the height diffrerence on the drone profile w.r.t. a reference asbcssa point
        """

        # sort the dataframe using the curvilinear abscissa
        sorted_df = self.df.sort_values(by='absc_curv')

        # interpolate the reference wsh
        wsh_ref = np.interp(ref_absc, sorted_df['absc_curv'].values, sorted_df['wsh_wgs84'].values)

        # compute the height difference
        sorted_df['h_diff'] = sorted_df['wsh_wgs84'] - wsh_ref

        return sorted_df






# ===============================================================
#    functions
# ===============================================================


# estimate the time shift between two time series
def estimate_time_shift(time_1_sec, wsh_1, time_2_sec, wsh_2) -> float:
    '''
        compute time shift in seconds between 2 time series
    '''

    # remove the mean on both curves
    wsh_1 = wsh_1 - np.mean(wsh_1)
    wsh_2 = wsh_2 - np.mean(wsh_2)

    # 2 cases to manage
    if len(time_1_sec) >= len(time_2_sec):
        # interp wsh_2 on time_1_sec
        wsh_2_interp = np.interp(time_1_sec, time_2_sec, wsh_2)
        correlation  = sgl.correlate(wsh_1, wsh_2_interp, mode='full')
        lags = sgl.correlation_lags(len(wsh_1), len(wsh_2_interp), mode='full')
        lag = lags[np.argmax(correlation)]
        time_shift = lag * np.median(np.diff(time_1_sec))
    else:
        # interp wsh_1 on time_2_sec
        wsh_1_interp = np.interp(time_2_sec, time_1_sec, wsh_1)
        correlation  = sgl.correlate(wsh_1_interp, wsh_2, mode='full')
        lags = sgl.correlation_lags(len(wsh_1_interp), len(wsh_2), mode='full')
        lag = lags[np.argmax(correlation)]
        time_shift = lag * np.median(np.diff(time_2_sec))
    
    return abs(time_shift)



