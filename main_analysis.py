# -*- coding: utf-8 -*-
"""
@author: Huy Tran Duc - IUPWARE - KUL
Master thesis - Extreme Event Impact Attribution: Heat-related Mortality Attribution of the 2020 Heatwave in Belgium

Main calculation
"""
import os
import numpy as np
import pandas as pd
import datetime as dt
import xarray as xr
import rioxarray
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geopandas as gpd
import cartopy as crs
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
# from rasterstats import zonal_stats
from pyproj  import Transformer
from scipy.stats import linregress
from netCDF4 import Dataset
from shapely import Point
import calendar


plt.rcParams["font.family"] = "Times New Roman"
#%% file path
pth = "D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/"
result_pth = os.path.join(pth, 'results')
with xr.open_dataset(pth+"data/RMI/"+"RMI_processed_data.nc") as rf:
    time = rf['time'][:]
    lon = rf['longitude'][:]
    lat = rf['latitude'][:]
    
    tmax = rf['MaxTemperature'][:]
    tmax_median = np.nanmean(tmax, axis=0)
    # sdate = dt.datetime.strptime("1954-01-01 00:00:00",'%Y-%m-%d %H:%M:%S')         #1951-1953 blank data
    # time = np.array([sdate+dt.timedelta(seconds=int(t)*86400) for t in time])
#     print(rf.variables)

shapefile_Be_path = pth+"DEM/"+'provinces_L08_repro.shp'          
shapefile_crs = 'EPSG:31370'    #Belgium Lambert=31370   #Default = 4326 

# Assign coordinate to shapefile
gdf_Be = gpd.read_file(shapefile_Be_path).to_crs(shapefile_crs)
# gdf_Be = gpd.read_file(shapefile_Be_path)
# Transform the projection from EPSG:4326 to EPSG:31370
transformer = Transformer.from_crs("EPSG:4326", "EPSG:31370", always_xy=True)
X,Y = np.meshgrid(lon,lat)
X, Y = transformer.transform(X,Y)
#%% Calculate mean temperature for each region (1954-2023)
# Create GeoDataFrame for transformed coordinates
data = pd.DataFrame(columns=['lon','lat'])
data['lon'] = X.flatten()
data['lat'] = Y.flatten()
geometry = [Point(xy) for xy in zip(X.flatten(),Y.flatten())]
geometry = gpd.GeoDataFrame(data,geometry=geometry)
# Assign each point to a region
geometry['zone'] = None
for i, row in gdf_Be.iterrows():
    geometry.loc[geometry.within(row.geometry),'zone']=row['Name']
# Calculate mean temperature for each region !!Data run slow
for i, row in gdf_Be.iterrows():
    print(row['Name'])
    tmp = geometry.within(row.geometry)
    tmp = np.where(tmp)[0]
    
    i = tmp//X.shape[1]
    j = tmp%X.shape[1]
    t_region = tmax[:,i,j]
    t_region_mean = np.nanmean(t_region)
    print(f"Mean temperature for {row['Name']}: {t_region_mean:.2f}°C")
    
#%% Calculate mean temperature for each region heatwave 2020
#  Total mean temperature 
start_date = dt.datetime(2020, 8, 5)
end_date = dt.datetime(2020, 8, 16)

# Convert dates to numpy datetime64 format for comparison
start_date_idx = np.where(time == np.datetime64(start_date))[0][0]
end_date_idx = np.where(time == np.datetime64(end_date))[0][0]

# Slice tmax data for 5th August to 16th August 2020
tmax_period = tmax[start_date_idx:end_date_idx + 1, :, :]

# Assuming X and Y are transformed coordinates already defined
# Create GeoDataFrame for transformed coordinates
data_hw = pd.DataFrame({'lon': X.flatten(), 'lat': Y.flatten()})
geometry_hw = [Point(xy) for xy in zip(X.flatten(), Y.flatten())]
geometry_hw = gpd.GeoDataFrame(data_hw, geometry=geometry_hw)

# Spatial join to assign each point to a region
geometry_hw['zone'] = None
for i, row in gdf_Be.iterrows():
    geometry_hw.loc[geometry_hw.within(row.geometry), 'zone'] = row['Name']

# Calculate daily mean temperature for each day from 5th to 16th August 2020
for day_idx in range(start_date_idx, end_date_idx + 1):
    tmax_day = tmax[day_idx, :, :]
    
    for zone, group in geometry_hw.groupby('zone'):
        idx = group.index  # Get indices of points within the region
        idx_lon = idx % lon.shape[0]
        idx_lat = idx // lon.shape[0]
        
        # Ensure indices are within bounds of tmax_day
        valid_idx = (idx_lon < tmax_day.shape[1]) & (idx_lat < tmax_day.shape[0])
        
        if np.any(valid_idx):
            t_region_day_mean = np.nanmean(tmax_day[idx_lat[valid_idx], idx_lon[valid_idx]])
            print(f"Day {day_idx - start_date_idx + 1} (Index {day_idx}): Daily mean temperature for {zone}: {t_region_day_mean:.2f}°C")
        else:
            print(f"Day {day_idx - start_date_idx + 1} (Index {day_idx}): No valid points found for {zone}.")

#%% Find indices for the time dimension for 5th August to 16th August 2020
# Daily mean temperature for heatwaves event 2020
start_date = dt.datetime(1990, 8, 5)
end_date = dt.datetime(2015, 8, 16)

# Convert dates to numpy datetime64 format for comparison
start_date_idx = np.where(time == np.datetime64(start_date))[0][0]
end_date_idx = np.where(time == np.datetime64(end_date))[0][0]

# Slice tmax data for 5th August to 16th August 2020
tmax_period = tmax[start_date_idx:end_date_idx + 1, :, :]

# Create GeoDataFrame for transformed coordinates
data_hw = pd.DataFrame({'lon': X.flatten(), 'lat': Y.flatten()})
geometry_hw = [Point(xy) for xy in zip(X.flatten(), Y.flatten())]
geometry_hw = gpd.GeoDataFrame(data_hw, geometry=geometry_hw)

# Spatial join to assign each point to a region
geometry_hw['zone'] = None
for i, row in gdf_Be.iterrows():
    geometry_hw.loc[geometry_hw.within(row.geometry), 'zone'] = row['Name']

# Create an empty list to store DataFrames for each day
results = []

# Calculate daily mean temperature for each day from 5th to 16th August 2020
for day_idx in range(start_date_idx, end_date_idx + 1):
    tmax_day = tmax[day_idx, :, :]
    
    daily_results = []
    for zone, group in geometry_hw.groupby('zone'):
        idx = group.index  # Get indices of points within the region
        idx_lon = idx % lon.shape[0]
        idx_lat = idx // lon.shape[0]
        
        # Ensure indices are within bounds of tmax_day
        valid_idx = (idx_lon < tmax_day.shape[1]) & (idx_lat < tmax_day.shape[0])
        
        if np.any(valid_idx):
            t_region_day_mean = np.nanmean(tmax_day[idx_lat[valid_idx], idx_lon[valid_idx]])
            date = time[day_idx].astype('datetime64[D]').astype(str)
            daily_results.append({'Date': date, 'canton_name': zone, 'tmean': t_region_day_mean})
        else:
            print(f"Day {day_idx - start_date_idx + 1} (Index {day_idx}): No valid points found for {zone}.")
    
    # Convert daily results to DataFrame and append to results list
    daily_df = pd.DataFrame(daily_results)
    results.append(daily_df)

# Concatenate all daily DataFrames into one DataFrame
results_df = pd.concat(results, ignore_index=True)

# Define path to save CSV file
output_path = "D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/results/daily_mean_temperatures_regions.csv"

# Save results to CSV
results_df.to_csv(output_path, index=False)
print(f"Results saved to {output_path}")

#%% Daily mean temperature for heatwaves event from 1990-2015 for each regions
# Define start and end dates
start_date = dt.datetime(1990, 8, 5)
end_date = dt.datetime(2015, 8, 16)

# Convert dates to numpy datetime64 format for comparison
start_date_idx = np.where(time == np.datetime64(start_date))[0][0]
end_date_idx = np.where(time == np.datetime64(end_date))[0][0]

# Create an empty list to store DataFrames for each day
results = []

# Loop over each year from 1990 to 2015
for year in range(1990, 2016):
    # Calculate start and end date for the current year
    start_date_year = dt.datetime(year, 8, 5)
    end_date_year = dt.datetime(year, 8, 16)

    # Convert dates to numpy datetime64 format for comparison
    start_date_idx_year = np.where(time == np.datetime64(start_date_year))[0][0]
    end_date_idx_year = np.where(time == np.datetime64(end_date_year))[0][0]

    # Slice tmax data for 5th August to 16th August of the current year
    tmax_period = tmax[start_date_idx_year:end_date_idx_year + 1, :, :]

    # Create GeoDataFrame for transformed coordinates
    data_hw = pd.DataFrame({'lon': X.flatten(), 'lat': Y.flatten()})
    geometry_hw = [Point(xy) for xy in zip(X.flatten(), Y.flatten())]
    geometry_hw = gpd.GeoDataFrame(data_hw, geometry=geometry_hw)

    # Spatial join to assign each point to a region
    geometry_hw['zone'] = None
    for i, row in gdf_Be.iterrows():
        geometry_hw.loc[geometry_hw.within(row.geometry), 'zone'] = row['Name']

    # Calculate daily mean temperature for each day from 5th to 16th August of the current year
    for day_idx in range(start_date_idx_year, end_date_idx_year + 1):
        tmax_day = tmax[day_idx, :, :]

        daily_results = []
        for zone, group in geometry_hw.groupby('zone'):
            idx = group.index  # Get indices of points within the region
            idx_lon = idx % lon.shape[0]
            idx_lat = idx // lon.shape[0]

            # Ensure indices are within bounds of tmax_day
            valid_idx = (idx_lon < tmax_day.shape[1]) & (idx_lat < tmax_day.shape[0])

            if np.any(valid_idx):
                t_region_day_mean = np.nanmean(tmax_day[idx_lat[valid_idx], idx_lon[valid_idx]])
                date = time[day_idx].astype('datetime64[D]').astype(str)
                daily_results.append({'Date': date, 'canton_name': zone, 'tmean': t_region_day_mean})
            else:
                print(f"Day {day_idx - start_date_idx_year + 1} (Index {day_idx}): No valid points found for {zone}.")

        # Convert daily results to DataFrame and append to results list
        daily_df = pd.DataFrame(daily_results)
        results.append(daily_df)

# Concatenate all daily DataFrames into one DataFrame
results_df = pd.concat(results, ignore_index=True)

# Define path to save CSV file
output_path = "D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/results/daily_mean_temperatures_regions_1990_2015_test.csv"

# Save results to CSV
results_df.to_csv(output_path, index=False)
print(f"All results saved to {output_path}")

#%% Calculate tmean for Belgium during heatwave period from 1990-2015
start_date = dt.datetime(1990, 8, 5)
end_date = dt.datetime(1990, 8, 16)  # End date for the period of interest

# Convert dates to numpy datetime64 format for comparison
start_date_idx = np.where(tmax.time.values == np.datetime64(start_date))[0][0]
end_date_idx = np.where(tmax.time.values == np.datetime64(end_date))[0][0]

# Create an empty list to store dictionaries for each year
results = []

# Loop over each year from 1990 to 2015
for year in range(1990, 2016):
    # Calculate start and end date for the current year
    start_date_year = dt.datetime(year, 8, 5)
    end_date_year = dt.datetime(year, 8, 16)

    # Convert dates to numpy datetime64 format for comparison
    start_date_idx_year = np.where(tmax.time.values == np.datetime64(start_date_year))[0][0]
    end_date_idx_year = np.where(tmax.time.values == np.datetime64(end_date_year))[0][0]

    # Slice tmax data for 5th August to 16th August of the current year
    tmax_period = tmax.isel(time=slice(start_date_idx_year, end_date_idx_year + 1))

    # Calculate mean temperature for this period
    t_period_mean = np.nanmean(tmax_period)

    # Get the year and store the result
    year_str = str(year)
    results.append({'Year': year_str, 'tmean': t_period_mean})

# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Define path to save CSV file
output_path = "D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/results/daily_meantemp_all_be_1990-2015.csv"

# Save results to CSV
results_df.to_csv(output_path, index=False)
print(f"All results saved to {output_path}")

#%% Calculate tmean for Belgium period from 1990-2015
# Create an empty list to store dictionaries for each day
results = []

# Loop over each year from 1990 to 2022
for year in range(1990, 2023):
    # Loop over each month
    for month in range(1, 13):
        # Get the number of days in the current month
        num_days = calendar.monthrange(year, month)[1]
        # Loop over each day of the month
        for day in range(1, num_days + 1):
            # Calculate the current date
            current_date = dt.datetime(year, month, day)

            # Convert date to numpy datetime64 format for comparison
            current_date_idx = np.where(tmax.time.values == np.datetime64(current_date))[0][0]

            # Slice tmax data for the current date
            tmax_current_day = tmax.isel(time=current_date_idx)

            # Calculate mean temperature for this day
            t_current_day_mean = np.nanmean(tmax_current_day)

            # Get the year, month, day and store the result
            results.append({
                'Year': year,
                'Month': month,
                'Day': day,
                'tmean': t_current_day_mean
            })

# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Define path to save CSV file
output_path = "D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/results/daily_meantemp_be_1990-2023.csv"

# Save results to CSV
results_df.to_csv(output_path, index=False)
print(f"All results saved to {output_path}")

#%% Calculate yearly temperature in belgium (1954-2022)
# File path for the output
output_path = "D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/results/yearly_meantemp_be_1954-2022.csv"

# Initialize an empty list to store the results
yearly_results = []

# Loop over each year from 1990 to 2015
for year in range(1954, 2023):
    # Initialize a list to store the daily mean temperatures for the current year
    yearly_temps = []

    # Loop over each month
    for month in range(1, 13):
        # Get the number of days in the current month
        num_days = calendar.monthrange(year, month)[1]
        
        # Loop over each day of the month
        for day in range(1, num_days + 1):
            # Calculate the current date
            current_date = dt.datetime(year, month, day)

            # Convert date to numpy datetime64 format for comparison
            try:
                current_date_idx = np.where(tmax.time.values == np.datetime64(current_date))[0][0]
            except IndexError:
                # If the date is not found in tmax.time, continue to the next iteration
                continue

            # Slice tmax data for the current date
            tmax_current_day = tmax.isel(time=current_date_idx)

            # Calculate mean temperature for this day
            t_current_day_mean = np.nanmean(tmax_current_day)

            # Store the daily mean temperature
            yearly_temps.append(t_current_day_mean)

    # Calculate the yearly mean temperature if there are daily temperature data
    if yearly_temps:
        yearly_mean_temp = np.nanmean(yearly_temps)
    else:
        yearly_mean_temp = np.nan  # No data for this year

    # Store the result for the current year
    yearly_results.append({
        'Year': year,
        'YearlyMeanTemp': yearly_mean_temp
    })

# Convert the results to a DataFrame
yearly_results_df = pd.DataFrame(yearly_results)

# Save the results to a CSV file
yearly_results_df.to_csv(output_path, index=False)
print(f"All results saved to {output_path}")


#%% Calculate yearly temperature in belgium (1954-2022) for summer JJA
# Input for line regression (Warming level)
# File path for the output
output_path = "D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/results/jja_meantemp_be_1954-2022.csv"

# Initialize an empty list to store the results
jja_results = []

# Loop over each year from 1990 to 2015
for year in range(1954, 2023):
    # Initialize a list to store the daily mean temperatures for the current summer period
    summer_temps = []

    # Loop over the summer months: June, July, August
    for month in [6, 7, 8]:
        # Get the number of days in the current month
        num_days = calendar.monthrange(year, month)[1]
        
        # Loop over each day of the month
        for day in range(1, num_days + 1):
            # Calculate the current date
            current_date = dt.datetime(year, month, day)

            # Convert date to numpy datetime64 format for comparison
            try:
                current_date_idx = np.where(tmax.time.values == np.datetime64(current_date))[0][0]
            except IndexError:
                # If the date is not found in tmax.time, continue to the next iteration
                continue

            # Slice tmax data for the current date
            tmax_current_day = tmax.isel(time=current_date_idx)

            # Calculate mean temperature for this day
            t_current_day_mean = np.nanmean(tmax_current_day)

            # Store the daily mean temperature
            summer_temps.append(t_current_day_mean)

    # Calculate the JJA mean temperature if there are daily temperature data
    if summer_temps:
        jja_mean_temp = np.nanmean(summer_temps)
    else:
        jja_mean_temp = np.nan  # No data for this summer period

    # Store the result for the current year
    jja_results.append({
        'time': year,
        'jja': jja_mean_temp
    })

# Convert the results to a DataFrame
jja_results_df = pd.DataFrame(jja_results)

# Save the results to a CSV file
jja_results_df.to_csv(output_path, index=False)
print(f"All results saved to {output_path}")

#%% Daily mean temperature in Belgium summer 2020
# Path to save the output CSV file
output_path = "D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/results/daily_jja_meantemp_be_2020.csv"

# Initialize an empty list to store the results
daily_results = []

# Define the year for which we want to calculate the data
year = 2020

# Loop over the summer months: June, July, August
for month in [6, 7, 8]:
    # Get the number of days in the current month
    num_days = calendar.monthrange(year, month)[1]
    
    # Loop over each day of the month
    for day in range(1, num_days + 1):
        # Calculate the current date
        current_date = dt.datetime(year, month, day)

        # Convert date to numpy datetime64 format for comparison
        try:
            current_date_idx = np.where(tmax.time.values == np.datetime64(current_date))[0][0]
        except IndexError:
            # If the date is not found in tmax.time, continue to the next iteration
            continue

        # Slice tmax data for the current date
        tmax_current_day = tmax.isel(time=current_date_idx)

        # Calculate mean temperature for this day
        t_current_day_mean = np.nanmean(tmax_current_day)

        # Store the daily mean temperature
        daily_results.append({
            'date': current_date,
            'temperature': t_current_day_mean
        })

# Convert the results to a DataFrame
daily_results_df = pd.DataFrame(daily_results)

# Save the results to a CSV file
daily_results_df.to_csv(output_path, index=False)
print(f"Daily results for summer 2020 saved to {output_path}")


#%% Daily mean temperature in Belgium in 2020
# Path to save the output CSV file
output_path = "D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/results/daily_meantemp_be_2020.csv"

# Initialize an empty list to store the results
daily_results = []

# Define the year for which we want to calculate the data
year = 2020

# Loop over all months from January to December
for month in range(1, 13):
    # Get the number of days in the current month
    num_days = calendar.monthrange(year, month)[1]
    
    # Loop over each day of the month
    for day in range(1, num_days + 1):
        # Calculate the current date
        current_date = dt.datetime(year, month, day)

        # Convert date to numpy datetime64 format for comparison
        try:
            current_date_idx = np.where(tmax.time.values == np.datetime64(current_date))[0][0]
        except IndexError:
            # If the date is not found in tmax.time, continue to the next iteration
            continue

        # Slice tmax data for the current date
        tmax_current_day = tmax.isel(time=current_date_idx)

        # Calculate mean temperature for this day
        t_current_day_mean = np.nanmean(tmax_current_day)

        # Store the daily mean temperature
        daily_results.append({
            'date': current_date,
            'temperature': t_current_day_mean
        })

# Convert the results to a DataFrame
daily_results_df = pd.DataFrame(daily_results)

# Save the results to a CSV file
daily_results_df.to_csv(output_path, index=False)
print(f"Daily results for the year 2020 saved to {output_path}")
#%%Calculate anomaly summer 2020 vs 1990-2015 average (June, July, August) using axarray (different result)
# Assuming 'time' is a coordinate in your xarray Dataset or DataArray object
# Convert 'time' to a DataArray if it's not already
time_dataarray = xr.DataArray(time, dims='time')

# Create a boolean mask to select months 6, 7, and 8 (June, July, August)
m_summer = (time_dataarray.dt.month == 6) | (time_dataarray.dt.month == 7) | (time_dataarray.dt.month == 8)

#'m_summer' is a boolean DataArray with True values for the summer months
# Create a boolean mask to select year 2020
y_select = time_dataarray.dt.year == 2020
y_period = (time_dataarray.dt.year >= 1990) & (time_dataarray.dt.year <= 2015)
loc = np.logical_and(m_summer,y_select)
loc_summer = np.logical_and(m_summer,y_period) 
t_2020 = np.nanmean(tmax[loc,:,:],axis=0)
t_ave_summer = np.nanmean(tmax[loc_summer,:,:],axis=0)
anomaly_summer = t_2020-t_ave_summer
# Calculate the mean value of the anomaly_summer array
mean_anomaly_summer = np.nanmean(anomaly_summer)
min_anomaly_summer = np.nanmin(anomaly_summer)
max_anomaly_summer = np.nanmax(anomaly_summer)

# Print the results
print("Minimum value of anomaly_summer:", min_anomaly_summer)
print("Maximum value of anomaly_summer:", max_anomaly_summer)
print("Mean value of anomaly_summer:", mean_anomaly_summer)
#%% Set up heatwave event (5-16th August 2020)
# Create a boolean mask to select day 5-16th August 2020
hw_event = (time_dataarray.dt.month == 8) & (time_dataarray.dt.day >= 5) & (time_dataarray.dt.day <= 16)

#'hw_summer' is a boolean DataArray with True values for the heatwave period
# Create a boolean mask to select year 2020
y_select = time_dataarray.dt.year == 2020
y_period = (time_dataarray.dt.year >= 1990) & (time_dataarray.dt.year <= 2015)
loc = np.logical_and(hw_event,y_select)
loc_hw = np.logical_and(hw_event,y_period) 
t_2020_hw = np.nanmean(tmax[loc,:,:],axis=0)
t_ave_hw = np.nanmean(tmax[loc_hw,:,:],axis=0)
anomaly_hw = t_2020_hw-t_ave_hw  
# Calculate the mean value of the anomaly_summer array
mean_anomaly_hw = np.nanmean(anomaly_hw)
min_anomaly_hw = np.nanmin(anomaly_hw)
max_anomaly_hw = np.nanmax(anomaly_hw)

# Print the results
print("Minimum value of anomaly heatwaves:", min_anomaly_hw)
print("Maximum value of anomaly heatwaves:", max_anomaly_hw)
print("Mean value of anomaly heatwaves:", mean_anomaly_hw)  
#%% Ploting Daily Mean temperature 1954-2023
lonmin = 22000  #bounds[0]
lonmax = 300000 #bounds[2]
latmin = 21000  #bounds[1]
latmax = 245000 #bounds[3]

fig = plt.figure(figsize=(8,5),dpi=300)
ax = plt.subplot(projection=crs.crs.epsg(31370))
ax.add_geometries(gdf_Be['geometry'],crs=crs.crs.epsg(31370),edgecolor='k',facecolor='none')
bounds = gdf_Be.total_bounds  # [minx, miny, maxx, maxy]
ax.set_extent([bounds[0], bounds[2], bounds[1], bounds[3]], crs=crs.crs.epsg(31370))

cs = ax.contourf(X,Y,tmax_median,np.linspace(10,15,11),cmap='coolwarm')         #np.linspace(10,15,11)  11=(15-10)*2+1 --> equal the level
# cs = ax.imshow(tmax_median[::-1], extent=[lonmin, lonmax, latmin, latmax], cmap='coolwarm')
# cs = ax.pcolormesh(X,Y,tmax_median,cmap='coolwarm')
cb = plt.colorbar(cs,shrink = 0.8,ticks = range(0,16))
cb.ax.set_title('(℃)',rotation=0,fontsize=10,weight='bold')
cb.ax.tick_params(labelsize=10)

ax.set_xticks(np.linspace(lonmin,lonmax,6), crs=crs.crs.epsg(31370))
ax.set_yticks(np.linspace(latmin,latmax,6), crs=crs.crs.epsg(31370))
ax.tick_params(axis="x", labelsize=10, pad=10)
ax.tick_params(axis="y", labelsize=10, pad=10)
# Hide axis
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.axis('off')          #Turn of the axis
# Customizing the plot
plt.title('Daily Mean Temperature Over Time (1954-2023)')  # Title of the plot
# plt.grid(color='k',alpha=0.2,linestyle='--')  # Show grid lines
plt.show()

#%% Ploting anomaly for summer 2020
fig = plt.figure(figsize=(8,5),dpi=300)
ax = plt.subplot(projection=crs.crs.epsg(31370))
ax.add_geometries(gdf_Be['geometry'],crs=crs.crs.epsg(31370),edgecolor='k',facecolor='none')
bounds = gdf_Be.total_bounds  # [minx, miny, maxx, maxy]
ax.set_extent([bounds[0], bounds[2], bounds[1], bounds[3]], crs=crs.crs.epsg(31370))

cs = ax.contourf(X,Y,anomaly_summer,np.linspace(0,2.5,6), cmap='autumn_r')
cb = plt.colorbar(cs, shrink = 0.8,ticks = np.arange(0,3,0.5))
cb.ax.set_title('Anomaly (℃)',rotation=0,fontsize=10,weight='bold', pad=10)
cb.ax.tick_params(labelsize=10)

ax.set_xticks(np.linspace(lonmin,lonmax,6), crs=crs.crs.epsg(31370))
ax.set_yticks(np.linspace(latmin,latmax,6), crs=crs.crs.epsg(31370))
ax.tick_params(axis="x", labelsize=10, pad=10)
ax.tick_params(axis="y", labelsize=10, pad=10)
# Hide axis titles
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.axis('off')
# Customizing the plot
# plt.title('Daily Mean Temperature Anomaly Summer 2020 \nCompare to Summer 1990-2015 In Belgium')  # Title of the plot
plt.title('Daily Mean Temperature Anomaly in Summer 2020')
plt.show()

#%% Ploting anomaly heatwave event (5-16th August 2020)
fig = plt.figure(figsize=(8,5),dpi=300)
ax = plt.subplot(projection=crs.crs.epsg(31370))
ax.add_geometries(gdf_Be['geometry'],crs=crs.crs.epsg(31370),edgecolor='k',facecolor='none')
bounds = gdf_Be.total_bounds  # [minx, miny, maxx, maxy]
ax.set_extent([bounds[0], bounds[2], bounds[1], bounds[3]], crs=crs.crs.epsg(31370))

cs = ax.contourf(X,Y,anomaly_hw,np.linspace(6,10,9), cmap='autumn_r')
cb = plt.colorbar(cs, shrink = 0.8, ticks = np.arange(6,10.5))
cb.ax.set_title('Anomaly (℃)',rotation=0,fontsize=10,weight='bold',pad=10)
cb.ax.tick_params(labelsize=10)

ax.set_xticks(np.linspace(lonmin,lonmax,6), crs=crs.crs.epsg(31370))
ax.set_yticks(np.linspace(latmin,latmax,6), crs=crs.crs.epsg(31370))
ax.tick_params(axis="x", labelsize=10, pad=10)
ax.tick_params(axis="y", labelsize=10, pad=10)
# Hide axis titles
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.axis('off')
# Customizing the plot
# plt.title('Daily Mean Temperature Anomaly \nDuring Heatwaves 2020 Compare to 1990-2015 In Belgium')  # Title of the plot
plt.title('Daily Mean Temperature Anomaly of Heatwaves in 2020')
plt.show()



    
















