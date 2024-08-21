# -*- coding: utf-8 -*-
"""
@author: Huy Tran Duc - IUPWARE - KUL
Master thesis - Extreme Event Impact Attribution: Heat-related Mortality Attribution of the 2020 Heatwave in Belgium

Plotting script
"""

import os
import numpy as np
import pandas as pd
import datetime as dt
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geopandas as gpd
import cartopy as crs
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyproj  import Transformer
from scipy.stats import linregress
from netCDF4 import Dataset

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

shapefile_Be_path = pth+"DEM/"+'provinces_L08_repro.shp'          
shapefile_crs = 'EPSG:31370'    #Belgium Lambert=31370   #Default = 4326 
gdf_Be = gpd.read_file(shapefile_Be_path).to_crs(shapefile_crs)
#%% Plotting
# FIGURE 1
# Daily mean temperature over time (1954-2023)
# Transform the projection from EPSG:4326 to EPSG:31370
transformer = Transformer.from_crs("EPSG:4326", "EPSG:31370", always_xy=True)
X,Y = np.meshgrid(lon,lat)
X, Y = transformer.transform(X,Y)

lonmin = 22000  #bounds[0]
lonmax = 300000 #bounds[2]
latmin = 21000  #bounds[1]
latmax = 245000 #bounds[3]

fig = plt.figure(figsize=(8,5),dpi=300)
ax = plt.subplot(projection=crs.crs.epsg(31370))
ax.add_geometries(gdf_Be['geometry'],crs=crs.crs.epsg(31370),edgecolor='k',facecolor='none')
bounds = gdf_Be.total_bounds  # [minx, miny, maxx, maxy]
ax.set_extent([bounds[0], bounds[2], bounds[1], bounds[3]], crs=crs.crs.epsg(31370))

cs = ax.contourf(X,Y,tmax_median,cmap='coolwarm')
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

cs = ax.contourf(X,Y,anomaly_summer,cmap='autumn_r')
cb = plt.colorbar(cs, shrink = 0.8,ticks = range(0,3))
cb.ax.set_title('Anomaly (℃)',rotation=0,fontsize=10,weight='bold')
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
plt.title('Daily Mean Temperature Anomaly Summer 2020 \nCompare to Summer 1990-2015 In Belgium')  # Title of the plot
plt.show()
#%% Ploting anomaly heatwave event (5-16th August 2020)
fig = plt.figure(figsize=(8,5),dpi=300)
ax = plt.subplot(projection=crs.crs.epsg(31370))
ax.add_geometries(gdf_Be['geometry'],crs=crs.crs.epsg(31370),edgecolor='k',facecolor='none')
bounds = gdf_Be.total_bounds  # [minx, miny, maxx, maxy]
ax.set_extent([bounds[0], bounds[2], bounds[1], bounds[3]], crs=crs.crs.epsg(31370))

cs = ax.contourf(X,Y,anomaly_hw,cmap='autumn_r')
cb = plt.colorbar(cs, shrink = 0.8,ticks = range(5,10))
cb.ax.set_title('Anomaly (℃)',rotation=0,fontsize=10,weight='bold')
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
plt.title('Daily Mean Temperature Anomaly \nDuring Heatwaves 2020 Compare to 1990-2015 In Belgium')  # Title of the plot
plt.show()
#%% Plot daily mean data for summer 2020
# Slice the temperature data for the period June 1st to August 31st, 2020
tmax_period = tmax.sel(time=slice('2020-06-01', '2020-08-31'))

# Calculate the mean temperature for each day within the specified period
daily_average_temperature = tmax_period.mean(dim=['latitude', 'longitude'])

# Plot the daily temperature data with month names as x-axis labels
plt.figure(figsize=(10, 6))
daily_average_temperature.plot()

plt.title('Daily Average Temperature (June-August 2020)')
plt.xlabel('Month')
plt.ylabel('Temperature (°C)')
plt.grid(True)
plt.show()

# Insert mortality data
# heat - related all death
all_death_be_path = pth+ "data/" +'TF_DEATHS.xlsx'
all_death_be = pd.read_excel(all_death_be_path)

# Step 1: Load the data
# Ensure the date column is in datetime format
all_death_be['DATE_DEATH'] = pd.to_datetime(all_death_be['DATE_DEATH'], format= '%d/%m/%Y')

# Step 2: Filter data for dates between 2020-06-01 and 2020-08-31
start_date = '2020-06-01'
end_date = '2020-08-31'
mask = (all_death_be['DATE_DEATH'] >= start_date) & (all_death_be['DATE_DEATH'] <= end_date)
filtered_df = all_death_be.loc[mask]

# Step 3: Plotting
plt.figure(figsize=(10, 6))  # Set the figure size (optional)
plt.plot(filtered_df['DATE_DEATH'], filtered_df['CNT'], marker='o', linestyle='-', color='b')
plt.title('Total Death Number from 2020-06-01 to 2020-08-31')
plt.xlabel('Date')
plt.ylabel('Total Death Number')
plt.xticks(rotation=45)  # Rotate date labels for better readability
plt.tight_layout()  # Adjust layout to not cut off labels
plt.grid(True)  # Optional: Add grid for better readability
plt.show()

# Merge 2 graphs into 1
plt.figure(figsize=(15, 6), dpi=300)
plt.title('Daily Deaths Number and Daily Average Temperature (June-August 2020)',fontsize=16)

# Plot the total death number graph with markers and lines extending to the date axis
plt.plot(filtered_df['DATE_DEATH'], filtered_df['CNT'], markersize=8, marker='o', linestyle='', c='gray', mfc='r')
for DATE_DEATH, CNT in zip(filtered_df['DATE_DEATH'], filtered_df['CNT']):
    plt.vlines(DATE_DEATH, 0, CNT, colors='gray')

# plt.plot(filtered_df['DATE_DEATH'], filtered_df['CNT'], marker='o', linestyle='', mfc='r')
# plt.xlabel('Date')
plt.ylabel('Daily Deaths', fontsize=16)
plt.xticks(rotation=45)
plt.ylim(250,500)           #adjust limit of y-axis (normall from 200-500)
plt.tight_layout()
plt.grid(True)

# Create a second axis for the temperature data
ax2 = plt.twinx()

# Plot the daily temperature data
daily_average_temperature.plot(ax=ax2, color='k', linestyle='-', linewidth=2)
plt.ylabel('Daily Mean Temperature (°C)', fontsize=16, rotation= 270, labelpad=15)

# Add a horizontal dashed line at y=25 on the secondary y-axis
ax2.axhline(y=25, color='blue', linestyle='--')

plt.show()



#%% Zoom in heatwaves event (5-16th August 2020)
# Slice the temperature data for the period 5-16th August, 2020
tmax_period_hw = tmax.sel(time=slice('2020-08-05', '2020-08-16'))

# Calculate the mean temperature for each day within the specified period
daily_average_temperature_hw = tmax_period_hw.mean(dim=['latitude', 'longitude'])

# Plot the daily temperature data with month names as x-axis labels
plt.figure(figsize=(10, 6))
daily_average_temperature_hw.plot()

# Set the x-axis labels to the names of the months
# plt.xticks(range(len(month_names)), month_names)

plt.title('Daily Average Temperature (Heatwave event 2020)')
plt.xlabel('Month')
plt.ylabel('Temperature (°C)')
plt.grid(True)
plt.show()

# Insert mortality data
# heat - related all death
all_death_be_path = pth+ "data/" +'TF_DEATHS.xlsx'
all_death_be = pd.read_excel(all_death_be_path)

# Step 1: Load the data
# Ensure the date column is in datetime format
all_death_be['DATE_DEATH'] = pd.to_datetime(all_death_be['DATE_DEATH'], format= '%d/%m/%Y')

# Step 2: Filter data for dates between 2020-06-01 and 2020-08-31
start_date_hw = '2020-08-05'
end_date_hw = '2020-08-16'
mask_hw = (all_death_be['DATE_DEATH'] >= start_date_hw) & (all_death_be['DATE_DEATH'] <= end_date_hw)
filtered_df_hw = all_death_be.loc[mask_hw]

# Step 3: Plotting
plt.figure(figsize=(10, 6))  # Set the figure size (optional)
plt.plot(filtered_df_hw['DATE_DEATH'], filtered_df_hw['CNT'], marker='o', linestyle='-', color='b')
plt.title('Total Death Number from Heatwave Event')
plt.xlabel('Date')
plt.ylabel('Total Death Number')
plt.xticks(rotation=45)  # Rotate date labels for better readability
plt.tight_layout()  # Adjust layout to not cut off labels
plt.grid(True)  # Optional: Add grid for better readability
plt.show()

# Merge 2 graphs into 1
plt.figure(figsize=(15, 6), dpi=300)
plt.title('Daily Deaths Number and Daily Average Temperature in Heatwave event',fontsize=16)

# Plot the total death number graph with markers and lines extending to the date axis
plt.plot(filtered_df_hw['DATE_DEATH'], filtered_df_hw['CNT'], markersize = 8, marker='o', linestyle='', c='gray', mfc='r')
for DATE_DEATH, CNT in zip(filtered_df_hw['DATE_DEATH'], filtered_df_hw['CNT']):
    plt.vlines(DATE_DEATH, 0, CNT, colors='gray')

# plt.plot(filtered_df['DATE_DEATH'], filtered_df['CNT'], marker='o', linestyle='', mfc='r')
# plt.xlabel('Date')
plt.ylabel('Daily Deaths',fontsize=16)
plt.xticks(rotation=45)
plt.ylim(250,500)           #adjust limit of y-axis (normall from 200-500)
plt.tight_layout()
plt.grid(True)

# Create a second axis for the temperature data
ax2 = plt.twinx()

# Plot the daily temperature data
daily_average_temperature_hw.plot(ax=ax2, color='k', linestyle='-', linewidth=2)
plt.ylabel('Daily Mean Temperature (°C)',fontsize=16, rotation= 270, labelpad=15)

# Add a horizontal dashed line at y=25 on the secondary y-axis
ax2.axhline(y=25, color='blue', linestyle='--')

plt.show()

#%% Mapping mean daily temp values for each region during heatwave period
# mean_temp_reg = pd.read_csv(pth+'data/', 'Mean_temp_regions_be_values.csv')
mean_temp_reg = pd.read_csv("D:\Belgium\OneDrive - KU Leuven\KUL\Thesis\data\Mean_temp_regions_be_values.csv")
gdf_Be_region = gdf_Be.merge(mean_temp_reg, on = "provISO")

# fig, ax = plt.subplots(1,1)
# gdf_Be_region.plot(ax=ax, column="Mean Temperature", cmap='YlOrRd', legend = True)
fig, ax = plt.subplots(figsize=(8, 5), dpi=300)

# Plot the regions with their mean temperature values
gdf_Be_region.boundary.plot(ax=ax, color='k', linewidth=0.5)
cs = gdf_Be_region.plot(ax=ax, column="Mean Temperature 1954-2023", cmap='YlOrRd', legend=True)
# Customize the color bar (accessed through the figure's axes)
cbar = ax.get_figure().get_axes()[1]  # The second axis is the color bar
cbar.set_ylabel('Mean Temperature (°C)', fontsize=12)
cbar.tick_params(labelsize=10)
# Customize the axes
ax.set_xticks(np.linspace(lonmin, lonmax, 6))
ax.set_yticks(np.linspace(latmin, latmax, 6))
ax.tick_params(axis="x", labelsize=10, pad=10)
ax.tick_params(axis="y", labelsize=10, pad=10)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.grid(False)
ax.axis('off')

plt.title('Daily Mean Temperature During Heatwaves 2020 In Belgium (°C)', fontsize=12)
plt.show()









