import pandas as pd
import numpy as np

def detect_fronts(df, var, criterion, x = 'distance_km', x_bin = 20):
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
    
    