import pandas as pd
import numpy as np
import xarray as xr

def front_finder(altimeter, lon_sdrn, lat_sdrn, time_sdrn,
                 saild, var, criterion, 
                 threshold = 2*1e-6, x_bin = 500, min_obs = None, method = 'jet_front'):
    ''' Identifying fronts using altimeter and in situ observations. 
    First, the algorithm tracks the in situ measurements located within peaks in gradients of 
    altimeter observations (e.g., SSH). It gives the initial condition for the position of mesoscale
    fronts. Then, it computes level of strength in the gradients of the chosen value from the in situ 
    observation.    
    ================================================================================================
    INPUT:
       altimeter = adt, sla, etc [time, n, m]
       lon_sdrn  = longitude array of the in sity obs [n]
       lat_sdrn  = latitude array of the in sity obs [m]
       time_sdrn = time array of the in sity obs [t]
       saild     = the data of the in in situ observations. 
                   It must contain a column named 'distance_km' with 
                   the distance between each sampling
       var       = variable name (string), to compute the gradient (e.g. Temperature, Salinity, 
                   Chlorophyll, Density)
       criterion = gradient criterion to group fronts into different levels (list of 5 values)
       method    = 'jet_front' uses both the defition to identify jets and fronts
                   'front_grad' identifies frontal filaments.  
    OUTPUT:
       The in situ dataset with the variables at the position of the identified fronts
    ================================================================================================
    Writers: Alessio Arena, Felipe Vilela da Silva, Mackenzie Blanusa, Maya Jakes and Sophie Clayton
    '''    
    if method == 'jet_front':
        # When does the position of the in situ observation coincide with the peaks in SSH gradient 
        # from the altimeter?
        boolean_pos = np.ravel([np.array(match_obs_alt(altimeter.sel(time=str(time_sdrn[t]), method = 'nearest'), 
                                                     lon_sdrn[t], lat_sdrn[t]))
                                for t in range(len(saild))])

        # What is the strength the SST from the in situ observation at these locations?
        saild_gradSSH = saild.loc[boolean_pos].copy()
        front_level = detect_grad_1d(saild_gradSSH, var, criterion, x_bin = x_bin, min_obs = min_obs)
        return front_level
    elif method == 'front_grad':
        front_level = detect_grad_1d(saild, var, criterion, x_bin = x_bin, min_obs = min_obs)
        return front_level

def match_obs_alt(altimeter, lon_sdrn, lat_sdrn, threshold = 2*1e-6):
    ''' Match the position of the peaks in SSH grad with the position of the in situ observation
    Writers: Felipe Vilela da Silva and Alessio Arena
    ==============================================================================
    INPUTS:
    
    OUTPUTS:
    '''
    
    if ~np.isnan(lat_sdrn):
        ## Next, I extract the altimeter information nearby the in situ observation
        alt_at_obs = altimeter.sel(latitude = slice(lat_sdrn-1, lat_sdrn+1), 
                                   longitude = slice(lon_sdrn-1,lon_sdrn+1))

        # TODO: if interested, find a more accurate way to convert deg to m
        # we have first followed the latitudinal and longitudinal distances but it only worked
        # if we turned the xarray into a np.array. We preferred to keep it as a xarray atm.
        deg_to_m = 111195.
        var_x = alt_at_obs.differentiate('longitude')/deg_to_m
        var_y = alt_at_obs.differentiate('latitude')/deg_to_m
        ## Below, I compute the module of the gradient 
        gradient = np.sqrt(var_x**2 + var_y**2)

        grad_thld   = gradient.where(gradient > threshold)
        grad_at_obs = grad_thld.sel(latitude=lat_sdrn.data, longitude=lon_sdrn.data, 
                                    method = 'nearest').data.item()
    
        return ~np.isnan(grad_at_obs)
    else:
        return False


def calculate_2d_gradient(var, threshold = 2*1e-6, to_dataframe=True):
    """Calculates the 2D gradient of the variable, optionally thresholding it and converting to a dataframe
    
    Parameters
    ----------
    var : xr.DataArray
        variable to process
    threshold : float, optional
        minimum value to keep
    to_dataframe : bool, optional
        whether you want a masked xr.DataArray, or a pd.Dataframe of only values above threshold
    
    Returns
    -------
    xr.DataArray or pd.DataFrame
        thresholded gradient values
    
    """
    if not isinstance(var, xr.DataArray):
        raise TypeError("please provide a DataArray. This may be done by selecting the variable of interest")
    deg_to_m = 111195.
    var_x = var.differentiate('longitude')/deg_to_m
    var_y = var.differentiate('latitude')/deg_to_m
    ## Below, I compute the module of the gradient 
    gradient = np.sqrt(var_x**2 + var_y**2)
    grad_thld   = gradient.where(gradient > threshold)
    if not to_dataframe:
        return grad_thld
    
    else:
        grad_thld_df = grad_thld.to_dataframe().dropna()
        grad_thld_df = grad_thld_df.reset_index(level=[1,2], drop=False)
        return grad_thld_df

def match_saildrone_to_seasurface_data(df_saildrone, df_var, lat_res=0.25, lon_res=0.25, keep_only_matched=True):
    """Finds the nearest seasurface data to each drone observation in space and time
    
    Parameters
    ----------
    df_saildrone : pd.DataFrame
        saildrone data having a time index and 'latitude' and 'longitude' columns
    df_var : pd.DataFrame
        sea surface data having time index and 'latitude', 'longitude' and a measurement column
    lat_res : float, optional
        resolution in degrees
    lon_res : float, optional
        resolution in degrees
    keep_only_matched: bool, optional
        whether to remove saildrone data that did not align with any seasurface data
    
    Returns
    -------
    xr.DataArray or pd.DataFrame
        thresholded gradient values
    
    """
    # TODO this assumes that df_var has one data layer per day
    # it should be generalised at later stages
    df_saildrone['date'] = df_saildrone.index.date
    
    df_saildrone.apply(lambda row: _is_saildrone_within_front(row, df_var, lat_res, lon_res), axis=1)
    if keep_only_matched:
        return df_saildrone.dropna(subset='gradient_value')
    

def _is_saildrone_within_front(row, gradient, lat_res=0.25, lon_res=0.25):
    # select appropriate date in gradient
    # we do the datetime matching once as that is time consuming
    row['gradient_value'] = np.nan
    row['gradient_lat'] = np.nan
    row['gradient_lon'] = np.nan
    
    try:
        selected = gradient.loc[pd.Timestamp(row.date)]
    except KeyError:
        return row
    
    # simple nearest neighbour approach with minimum distance condition
    # if the nearest neighbour is further than the data resolution, we don't have a match
    selected = selected.set_index('latitude', drop=True)
    lats = selected.index.values
    dlat = np.abs(lats - row.latitude)
    if dlat.min() < lat_res:
        selected = selected.loc[lats[dlat.argmin()]]
        # lat_selected = lats[dlat.argmin()]
        
        # selected = selected.set_index('longitude', drop=True)
        lons = selected.index.values
        dlon = np.abs(lons - row.longitude)
        if dlon.min() < lon_res:
            # lon_selected = lons[dlon.argmin()]
            out = selected.iloc[dlon.argmin()].drop('longitude')
            row['gradient_value'] = out.item()
            row['gradient_lat'] = lats[dlat.argmin()]
            row['gradient_lon'] = lons[dlon.argmin()]
            return row
        else:
            return row
    else:
        return row    
    
def detect_grad_1d(df, var, criterion, x = 'distance_km', x_bin = 20, min_obs = 10):
    ''' Detect fronts from 1D data
    Writers: Maya Jakes and Alessio Arena
    ==============================================================================
    INPUT:
    df = pandas dataframe
    var = variable name (as string), to compute the gradient (e.g. Temperature, Salinity, Chlorophyll, Density)
    criterion = gradient criterion to group fronts into different levels
    
    Optional inputs:
    x = name of x variable to compute the gradient along. Default = 'distance_km'
    x_bin = size of the bins to average over. Default = 20 (km)
    min_obs = minimum number of observations in each bin. Dedault = 10.
    
    OUTPUT:
    New dataframe containing binned averages of var, latitude, longitude, x, dx, d_var and d_var_dx.
    Given the criterion input, each binned group is assigned a level using d_var_dx to characterise the fronts.
    
    '''
    df['dx'] = np.abs(np.gradient(df[x]))
    df['d_var'] = np.abs(np.gradient(df[var]))
    
    df_new = df[[var, 'latitude', 'longitude', x, 'dx', 'd_var']]
    
    # group into segments with size accoridng to x_bin
    df_grouped = df_new.groupby(df_new[x] // x_bin)
    
    # Deal with groups that cross the international date line using average_lat_lon function.
    df_grouped = df_grouped.apply(average_lat_lon)
    df_grouped = df_grouped.groupby(df_grouped[x] // x_bin)
    
    # apply minimum observations function
    if min_obs != None:
        df_group_mean = df_grouped.apply(minimum_obs, min_obs = min_obs)
    else:
        # don't apply to lat lon
        df_group_mean = df_grouped.mean()

    # calculate gradient of var with distance
    df_group_mean['d_var_dx'] = np.round(df_group_mean['d_var'], 3)/df_group_mean['dx']
    
    # front detection criteria
    criterion = [-np.inf] + criterion + [np.inf]
    
    df_group_mean['level'] = pd.cut(df_group_mean['d_var_dx'], criterion, labels = False)
    
    return df_group_mean

def minimum_obs(group, min_obs = 10):
    ''' Check for groups with less than a certain number of observations (min_obs).
    '''
    if len(group) < min_obs:
        return np.nan
    else:
        return group.mean()

def average_lat_lon(group):
    lat = group.latitude
    lon = group.longitude
    
    x = np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon));
    y = np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon));
    z = np.sin(np.deg2rad(lat))
    
    mean_x, mean_y, mean_z = np.nanmean(x), np.nanmean(y), np.nanmean(z)
    
    group['longitude'] = np.rad2deg(np.arctan2(mean_y, mean_x))
    hyp = np.sqrt(mean_x * mean_x + mean_y * mean_y)
    group['latitude'] = np.rad2deg(np.arctan2(mean_z, hyp))
    
    return group
    
    
    
    
    
    
    