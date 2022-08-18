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
    dist_between_points_km = dist_between_points/1000
    dist_from_start = np.nancumsum(dist_between_points_km)
    
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



def rotate_vels(u, v, lons, lats):
    ''' Rotate velocities from eastward (u) and northward (v) to along-stream (u_rot) and cross-stream (v_rot), respectively.
     Writers: Maya Jakes
     ==============================================================================
     INPUT:
     u = eastward velocity component (could be 1D or 2D)
     v = northward velocity component (could be 1D or 2D)
     lons = 1D array of longitude values
     lats = 1D array of latitude values
     
     OUTPUT:
     Rotated velocities with one less observation in x than the input velocities (no forward azimuth from the last lat lon position)
     u_rot = along-stream velocity
     v_rot = cross-stream velocity 
     ==============================================================================
    '''
    # remove the last velocity observation (no forward azimuth from the last lat lon position)
    speed = np.sqrt(u**2 + v**2)[:-1]
    velocity_bearing = uv_bearing(u, v)[:-1]

    along_strm_brng = stream_bearing(lons, lats)
    # make 1D array into 2D in the same shape as velocity bearings
    along_strm_brng = np.tile(along_strm_brng,(len(u[0]), 1)).transpose()

    # find the angle between the velocity bearing and the along stream direction
    theta = along_strm_brng - velocity_bearing

    # calculate u and v using this new angle (converting degrees to radians)
    u_rot = np.cos(theta*np.pi/180)*speed
    v_rot = np.sin(theta*np.pi/180)*speed
    
    return u_rot, v_rot


def uv_bearing(u, v):
    '''Calculates the bearing (clockwise from True North) using eastward (u) and northward (v) components of velocity'''
    theta = np.rad2deg(np.arctan2(u, v))
    theta += 360
    theta = theta % 360
    return theta


def stream_bearing(lons, lats):
    ''' Calculates the bearing (clockwise from True North) between each lat and lon position.
    Writers: Maya Jakes
    ==============================================================================
    INPUTS:
    lons = 1D array of longitude values
    lats = 1D array of latitude values
    
    OUTPUT:
    1D array of bearings (in degrees)
    '''
    bearing = []
    for i in range(0, len(lats)-1):
        lat1, lat2 = lats[i], lats[i+1]
        lon1, lon2 = lons[i], lons[i+1]
        
        # WGS84 is the reference coordinate system (Earth's centre of mass) used by GPS.
        geodesic = pyproj.Geod(ellps='WGS84')  
        # Inverse computation to calculate the forward azimuths from two lat and lon coordinates. 
        fwd_azimuth = geodesic.inv(lon1, lat1, lon2, lat2)[0]

        # if the angle is negative (anticlockwise from N), add it to 360 to get the bearing clockwise from N.
        if fwd_azimuth < 0:
            fwd_azimuth = 360 + fwd_azimuth
            
        bearing.append(fwd_azimuth)
        
    return np.asarray(bearing)