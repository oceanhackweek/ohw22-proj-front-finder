import pydap
from pydap.client import open_url
from pydap.cas.urs import setup_session
import xarray as xr
import gsw
import numpy as np


def import_cmems(variable, url, username, password, lon_min, lon_max, lat_min, lat_max, start_time, end_time):
    ''' Import CMEMS data
    Writers: Felipe Vilela da Silva and Alessio Arena
    ==============================================================================
    INPUT:
       variables = adt, sla, ugosa, ugos, etc (string)
       url       = copy the OPENDAP url from CMEMS (string)
       username and password are the credentials to access CMEMS data (string)
       lon_min and lon_max indicate the western and eastern boundaries, respectively (float)
       lat_min and lat_max indicate the southern and northern boundaries, respectively (float)
       start_time and end_time indicate the start and end time: (string)
    OUTPUT:
       CMEMS data in xarray
    ==============================================================================
    '''
    
    lon_slice = slice(lon_min, lon_max)
    lat_slice = slice(lat_min, lat_max)
    time_slice = slice(start_time, end_time)

    with setup_session(username, password, check_url=url) as session:
        pydap_ds = open_url(url, session=session)
        store = xr.backends.PydapDataStore(pydap_ds)
    
        ds = xr.open_dataset(store)[variable].sel(time=time_slice).sel(longitude=lon_slice, latitude=lat_slice)
    return ds



def distFromStart(latitude, longitude):
    '''Cumulative distance from the first data point (km)
     Writers: Maya Jakes
     ==============================================================================
     INPUT: 
     latitude = 1D array of latitude values
     longitude = 1D array of longitude values
     
     OUTPUT:
     Cumulative distance (km) from the first data point (1D array)
     ==============================================================================
     '''

    dist_between_points = np.concatenate((np.array([0]), gsw.distance(longitude, latitude)))
    dist_between_points_km = dist_between_points/1000
    dist_from_start = np.nancumsum(dist_between_points_km)
    
    return dist_from_start