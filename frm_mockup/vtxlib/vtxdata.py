# -*- coding: utf-8 -*-


#########################################
#   Lib for collecting vorteX-io data   #
#########################################

# 04/10/2022.... JC.Poisson Creation


# Libraries
import requests
import datetime as dt
import numpy as np
import pandas as pd
import json



# Vigicrues data object
class VORTEXIO_DATA:
    """
        vorteX-io data object
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
        self.link = "https://maelstrom.vortex-io.fr/api/metrics?name={}&timezone=UTC&max_result=1000&date_before={}&date_after={}"


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
        start = self.start_date.strftime("%Y-%m-%dT%H:%M:%S")
        end = self.end_date.strftime("%Y-%m-%dT%H:%M:%S")

        # prepare url, payload and headers
        url = self.link.format(self.station_id, end, start)
        payload={}
        headers = {
            'Accept': 'application/json',
            'Authorization': 'Bearer your_token_here'
        }

        # request data
        response = requests.request("GET", url, headers=headers, data=payload)
        d_response = json.loads(response.text)
        # convert to pandas.DataFrame
        l_data = d_response['metrics']
        d_data = dict.fromkeys(['time', 'wsh', 'water_surface_speed'])
        for k in d_data.keys():
            d_data[k] = []
        for data in l_data:
            d_data['time'].append(dt.datetime.strptime(data["time_stamp"], "%Y-%m-%d %H:%M:%S"))
            if "height" in data:
                d_data['wsh'].append(float(data['height']))
            else:
                d_data['wsh'].append(np.nan)
            if "avgSpeed" in data:
                d_data['water_surface_speed'].append(float(data["avgSpeed"]))
            else:
                d_data['water_surface_speed'].append(np.nan)
        
        # build time_sec from 2000-1-1
        ref_time = dt.datetime(2000, 1, 1)
        d_data['time_sec'] = [(d - ref_time).total_seconds() for d in d_data['time']]
        
        self.df = pd.DataFrame.from_dict(d_data)

        return 0


