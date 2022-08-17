try:
    # put import of packages that are not installed
    import pydap
    from pydap.client import open_url
    from pydap.cas.urs import setup_session
except ImportError:
    pass

import xarray as xr
import gsw
import shapely.geometry as shpgeom


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
    try:
        pydap
    except NameError:
        raise ImportError('Please install the pydap package')
    
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
    dist_between_points_km = dist_between_profiles/1000
    dist_from_start = np.nancumsum(dist_between_profiles_km)
    
    return dist_from_start


def isnullgeometry(geom):
    """Test for empty geometry, or generally null value
    
    Hint: on geodataframes do: gdf.geometry.apply(isnullgeometry)
    
    Parameters
    ----------
    geom : shapely.geometry or any other Python object
    
    Returns
    -------
    bool
        whether is an empty geometry, None or NaN
    
    """
    if isinstance(geom, (shpgeom.base.BaseGeometry, shpgeom.base.BaseMultipartGeometry)):
        return geom.is_empty
    else:
        return isnull(geom)
    
def isnull(val):
    """test for None or NaNs
    """
    return val != val or val is None