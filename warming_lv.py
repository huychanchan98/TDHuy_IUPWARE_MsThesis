# -*- coding: utf-8 -*-
"""

@author: Huy Tran Duc - IUPWARE - KUL
Master thesis - Extreme Event Impact Attribution: Heat-related Mortality Attribution of the 2020 Heatwave in Belgium

Calculate warming level in Belgium during heatwaves period
# The resulting slopes from the linear regressions inform on the warming (in °C) per °C GMST increase
# convert these slopes to warming levels by multiplying with a GMST increase of 1.15 °C attributable 
to human influence (IPCC 2021 SPM or chapter 3)

"""
import os
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.colors as mcolors
import statsmodels.api as sm
import rioxarray
import seaborn as sns
from sklearn.linear_model import LinearRegression
from pyproj import Transformer
from scipy.stats import linregress
import geopandas as gpd
import cartopy.crs as ccrs 
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
plt.rcParams["font.family"] = "Times New Roman"
#%% file path
pth = "D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/"
result_pth = os.path.join(pth, 'results')
#%% For SWISS 1864-2022
#Load GMST data from HADCRUT5
gmst_df = pd.read_csv(pth + 'data/HadCRUT5/HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv')

# Filter GMST data to include only the years from 1864 to 2022
gmst_df_filtered = gmst_df[(gmst_df['time'] >= 1864) & (gmst_df['time'] <= 2022)]

# Load MeteoSwiss data
meteoswiss_df = pd.read_csv(pth + 'data/meteoswiss/meteoswiss.csv')

# Filter MeteoSwiss data to include only the years from 1864 to 2022
meteoswiss_df_filtered = meteoswiss_df[(meteoswiss_df['time'] >= 1864) & (meteoswiss_df['time'] <= 2022)]

# Drop rows with NaN values in 'time' or 'jja' columns
meteoswiss_df_filtered = meteoswiss_df_filtered.dropna(subset=['time', 'jja'])

# Merge the two dataframes on the 'time' column
merged_df = pd.merge(gmst_df_filtered, meteoswiss_df_filtered, on='time', suffixes=('_gmst', '_jja'))

# Extract independent variable (GMST Anomaly) and dependent variable (JJA Temperature Anomaly)
X = merged_df['Anomaly (deg C)'].values.reshape(-1, 1)  # Reshape X into a 2D array for scikit-learn
y = merged_df['jja'].values

# Create and fit linear regression model
model = LinearRegression()
model.fit(X, y)

# Make predictions
y_pred = model.predict(X)

# Plot the data and regression line
plt.figure(dpi=300)
plt.scatter(X, y, color='blue', label='Observed data')
plt.plot(X, y_pred, color='red', label='Fitted line')
plt.xlabel('GMST Anomaly (°C)')
plt.ylabel('JJA Temperature Anomaly (°C)')
plt.title('MeteoSwiss Data (1864-2022) vs. GMST Anomaly')
plt.legend()
plt.show()

# Print coefficients
intercept = model.intercept_
slope = model.coef_[0]
print('Intercept:', intercept)
print('Slope:', slope)

# Calculate warming level
gmst_increase = 1.25  # GMST increase attributable to human influence (IPCC 2021)
warming_level = slope * gmst_increase
print('Warming level (°C per °C GMST increase):', warming_level)


#%% For BELGIUM 
#Load GMST data from HADCRUT5
gmst_df = pd.read_csv(pth + 'data/HadCRUT5/HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv')

# Filter GMST data to include only the years from 1864 to 2022
gmst_df_filtered = gmst_df[(gmst_df['time'] >= 1954) & (gmst_df['time'] <= 2022)]

# Load RMI data
rmi_jja_be_df = pd.read_csv(pth + 'results/jja_meantemp_be_1954-2022.csv')

# Filter rmi_jja_be data to include only the years from 1864 to 2022
rmi_jja_be_df_filtered = rmi_jja_be_df[(rmi_jja_be_df['time'] >= 1954) & (rmi_jja_be_df['time'] <= 2022)]

# Drop rows with NaN values in 'time' or 'jja' columns
rmi_jja_be_df_filtered = rmi_jja_be_df_filtered.dropna(subset=['time', 'jja'])

# Merge the two dataframes on the 'time' column
merged_df = pd.merge(gmst_df_filtered, rmi_jja_be_df_filtered, on='time', suffixes=('_gmst', '_jja'))

# Extract independent variable (GMST Anomaly) and dependent variable (JJA Temperature Anomaly)
X = merged_df['Anomaly (deg C)'].values.reshape(-1, 1)  # Reshape X into a 2D array for scikit-learn
y = merged_df['jja'].values

# Create and fit linear regression model
model = LinearRegression()
model.fit(X, y)

# Make predictions
y_pred = model.predict(X)

# Plot the data and regression line
plt.figure(dpi=300)
plt.scatter(X, y, color='blue', label='Observed data')
plt.plot(X, y_pred, color='red', label='Fitted line')
plt.xlabel('GMST Anomaly (°C)')
plt.ylabel('JJA Temperature Anomaly (°C)')
plt.title('RMI Data (1954-2022) vs. GMST Anomaly')
plt.legend()
plt.show()

# Print coefficients
intercept = model.intercept_
slope = model.coef_[0]
print('Intercept:', intercept)
print('Slope:', slope)

# Calculate warming level
gmst_increase = 1.3  # GMST increase attributable to human influence (IPCC 2021)
warming_level = slope * gmst_increase
print('Warming level (°C per °C GMST increase):', warming_level)