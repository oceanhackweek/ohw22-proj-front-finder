import pandas as pd
import numpy as np

def match_obs_alt(altimeter, lon_sdrn, lat_sdrn, threshold = 2*1e-6):
    ''' Match the position of the peaks in SSH grad with the position of the in situ observation
    Writers: Felipe Vilela da Silva and Alessio Arena
    ==============================================================================
    INPUTS:
    
    OUTPUTS:
    '''
    
    if ~np.isnan(lat_sdrn.data.item()):
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

def detect_grad_1d(df, var, criterion, x = 'distance_km', x_bin = 20):
    ''' Detect fronts from 1D data
    ==============================================================================
    INPUT:
    df = pandas dataframe
    var = variable name (as string), to compute the gradient (e.g. Temperature, Salinity, Chlorophyll, Density)
    criterion = gradient criterion to group fronts into different levels
    
    Optional inputs:
    x = name of x variable to compute the gradient along. Default = 'distance_km'
    x_bin = size of the bins to average over. Default = 20 (km)
    
    OUTPUT:
    New dataframe containing binned averages of var, latitude, longitude, x, dx, d_var and d_var_dx.
    Given the criterion input, each binned group is assigned a level using d_var_dx to characterise the fronts.
    
    '''
    
    df['dx'] = np.abs(np.gradient(df[x]))
    df['d_var'] = np.abs(np.gradient(df[var]))
    
    df_new = df[[var, 'latitude', 'longitude', x, 'dx', 'd_var']]
    
    # group into 20 km segments
    df_grouped = df_new.groupby(df_new[x] // x_bin)
    
    # apply minimum observations function
    # returns nan is there are < 10 observations else computes the mean of the group
    df_group_mean = df_grouped.apply(minimum_obs)
    
    # calculate gradient of var with distance
    df_group_mean['d_var_dx'] = np.round(df_group_mean['d_var'], 3)/df_group_mean['dx']
    
    # front detection criteria
    criterion = [-np.inf] + criterion + [np.inf]
    
    df_group_mean['level'] = pd.cut(df_group_mean['d_var_dx'], criterion, labels = False)
    
    return df_group_mean


def minimum_obs(group):
    ''' Check for groups with less than 10 observations.
    '''
    if len(group) < 10:
        return np.nan
    else:
        return group.mean()
    
    