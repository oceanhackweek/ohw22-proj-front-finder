import functions
import numpy as np
import xarray as xr
import geopandas as gpd

def get_closer_to_center(arr):
    """Iteratively check for valid elements of an array starting from the center
    The sequence of indices tried are (e.g. center is 6):
    6, then 7, then 5, then 8, then 4, etc...
    
    """
    
    data_len = arr.shape
    if len(data_len) != 1:
        raise ValueError("Only 1D array supported")
    data_len = data_len[0]
    center = round(data_len/2)
    
    if not functions.isnullgeometry(arr[center]):
        return arr[center]
    d_index = 1
    while True:
        try:
            el = arr[center + d_index]
            if not functions.isnullgeometry(el):
                return el
            
            el = arr[center - d_index]
            if not functions.isnullgeometry(el):
                return el
            d_index += 1
            
        except IndexError:
            return np.nan
        
def prepare_saildrone_data(ds, variables, resample_to='5Min', fill_up_to=15, aggregate_method='mean'):
    """Preprocess the saildrone data
    
    This function will do:
    - variable selection and conversion to dataframe
    - spatial interpolation to fill time gaps up to "fill_up_to" value
    - time resampling and aggregation using "resample_to" and "aggregate_method"
    - convert to geodataframe for spatial functionalities

    Parameters
    ----------
    ds : Dataframe or xarray object
        data to process (typically is an xarray object)
    variables : list of str
        variable list to preserve
    resample_to : str, optional
        time frequency string to use for the time resampling, by default '5Min'
    fill_up_to : int, optional
        maximum gap in time to fill using spatial interpolation (only applied to latitude and longitude columns), by default 15
    aggregate_method : str or dict, optional
        method to aggregate variables, by default 'mean'

    Returns
    -------
    gpd.GeoDataFrame
        output geodataframe of preprocessed variables
    """
    
    ds = ds[variables]
    if isinstance(ds, (xr.DataArray, xr.Dataset)):
        ds = ds.to_dataframe()
    
    #lat, lon and time are coordinates of the xarray, so should be included
    ds = ds.set_index('time', drop=False)
    
    ds['latitude'] = ds.latitude.interpolate(limit=fill_up_to)
    ds['longitude'] = ds.longitude.interpolate(limit=fill_up_to)
    
    if isinstance(aggregate_method, str):
        agg_dict = {}
        for v in variables:
            agg_dict[v] = aggregate_method
    elif isinstance(aggregate_method, dict):
        agg_dict = aggregate_method
        
    agg_dict['latitude'] = get_closer_to_center
    agg_dict['longitude'] = get_closer_to_center
    agg_dict['time'] = get_closer_to_center
    
    ds = ds.resample(resample_to).agg(
        agg_dict
    )
    
    gdf = gpd.GeoDataFrame(data=ds, geometry=gpd.points_from_xy(ds.longitude, ds.latitude), crs='4326')
    non_valid_geom = gdf.geometry.apply(functions.isnullgeometry)
    
    return gdf.loc[~non_valid_geom]
        
    