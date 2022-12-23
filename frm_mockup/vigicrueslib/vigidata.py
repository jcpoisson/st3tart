# -*- coding: utf-8 -*-


##########################################
#   Lib for collecting Vigicrues data   #
##########################################

# 04/11/2022.... JC.Poisson Creation

# Libs
import requests
import json
import datetime as dt
import numpy as np
import pandas as pd


# Vigicrues data object
class VIGICRUES_DATA:
    """
        vigicrues data object
    """

    # init
    def __init__(self):
        """
            constructor
        """
        self.df = None
        self.station_id = None
        self.start_date = None
        self.end_date = None
        self.link = "https://www.hydro.eaufrance.fr/stationhydro/ajax/{}/series?hydro_series[startAt]={}&hydro_series[endAt]={}&hydro_series[variableType]=simple_and_interpolated_and_hourly_variable&hydro_series[simpleAndInterpolatedAndHourlyVariable]=H&hydro_series[statusData]=most_valid"


    # set params
    def set_params(self, station_id, start_date, end_date) -> int:
        """
            set param for requesting vigicrues data for a specific station id and set stat_date and end_date
        """
        # set params
        self.station_id = station_id
        self.start_date = start_date
        self.end_date =end_date
        
        return 0


    # request data for a station between start_date and end_date
    def request_data(self) -> int:
        """
            request vigicrues data for a station based on its id between a start_date and en_date
        """
        # prepare data
        start = self.start_date.strftime("%d/%m/%YT%H:%M:%S")
        end = self.end_date.strftime("%d/%m/%YT%H:%M:%S")

        # prepare url, payload and headers
        url = self.link.format(self.station_id, start, end)
        payload={}
        headers = {
        'Content-Type': 'application/json'
        }

        # request data
        response = requests.request("GET", url, headers=headers, data=payload)

        # convert to pandas dataframe
        #t = datetime
        #v = valeur => Hauteur d'eau en mm
        df_data = pd.DataFrame.from_dict(json.loads(response.text)["series"]["data"])
        df_data['v_m'] = df_data.v / 1000.0
        df_data.v_m_positive = df_data.v_m + np.abs(df_data.v_m.min())
        df_data['t'] = pd.to_datetime(df_data['t'], format='%Y-%m-%dT%H:%M:%SZ')
        self.df = df_data[['t', 'v_m']]
        self.df = self.df.rename(columns={'t': 'time', 'v_m': 'wsh'})
        self.df['time_sec'] = (pd.to_datetime(self.df['time']) - dt.datetime(2000,1,1)).dt.total_seconds().values

        return 0



