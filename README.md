# climate_change_timeseries

In 2013, the topic of climate change was the focal point of much debate when reports claimed that global warming had stopped, while other research papers stated the opposite. A case study was conducted using data sources from NASA and the UK Met Office. 

The objective of this project was to repeat the analysis on the NASA and Met Office data but updated for 2021 and perform the following:

•	Forecast global average temperatures through to year 2100 and use the predictions to comment on the concerns about increasing global temperatures. 

•	Develop point predictions, including the 90% confidence intervals for global average temperatures for January and July 2030, 2050 and 2100. 

•	Repeat the analysis above using temperature data for Kingston, ON (postal code K7L3N6).

•	Repeat the analysis using the same data but for the time period covered in the infamous bet between University of Pennsylvania Professor J. Scott Armstrong, and former U.S. Vice President, Al Gore and comment. 

**Data Sources and Preparation**

_NASA Data_

Source: https://data.giss.nasa.gov/gistemp/

The NASA dataset used for this analysis is the combined land-surface air and sea-surface water monthly temperature anomalies from January 1880 to March 2021. The temperature anomalies are calculated as the deviations from the corresponding 1951-1980 means. This data was available in both text and CSV format. 

For the purpose of our analysis, the data was tidied such that the dates (by year, month, day) were in one column, followed by a column for temperature anomalies and a column for actual temperature. The actual temperature was converted by adding the average 30-year baseline temperature of 14 degrees Celsius. 

A check for missing values was also performed and none were identified. Lastly, the dataset was converted into a time series object starting from January 1950 to December 2020 with frequency equal to 12 to capture that there are 12 data points (i.e. months) in each year. 

_UK Met Data_

Source: https://www.metoffice.gov.uk/hadobs/hadcrut5/data/current/download.html

The UK MET dataset used for this analysis is the HadCRUT5 near surface temperature data set from January 1850 to January 2021. The HadCRUT5 is a blend of surface air temperature and sea-surface temperature anomalies in degrees’ Celsius relative to 1961 to 1990. The data was available in NetCDF format and included 10 variables. 

For the purpose of this analysis, the time and the tas_mean variable which is the average mean temperature anomalies, was used for this dataset. The data was tidied such that date (by year, month, day) was in one column, followed by a column for temperature anomalies and a column for actual temperature. The actual temperature was converted by adding the average 30-year baseline temperature of 14 degrees Celsius. 

A check for missing values was also performed and none were identified. Lastly, the dataset was converted into a time series object starting from January 1950 to December 2020 with frequency equal to 12 to capture that there are 12 data points (i.e. months) in each year. 

_Kingston ON Data_

Source: https://climate.weather.gc.ca/climate_data/monthly_data_e.html?hlyRange=%7C&dlyRange=1960-08-01%7C2007-11-30&mlyRange=1960-01-01%7C2006-12-01&StationID=4300&Prov=ON&urlExtension=_e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2021&selRowPerPage=25&Line=10&searchMethod=contains&Month=5&Day=11&txtStationName=kingston&timeframe=3&Year=2006

The Kingston dataset used for this analysis is the Kingston Pumping Station Ontario Monthly Data Report from 1960 to 2006. This data was available in CSV format and included 29 variables such as longitude, latitude, climate ID, mean temp, total rain, total snowfall, etc. 

For the purpose of this analysis, the following variables were extracted: year, month and mean temp variable which was in degrees Celsius. The data was tidied the such that date (by year, month, day) was in one column followed by the mean temperature (in degrees Celsius) in another column. 

As there were missing values until March 1978, the dataset was framed to start from January 1979 and end December 2006 to capture full years. 

**Exploratory Data Analysis**

NASA & UK MET DATA

_Time Plots_

As evidenced in the first time plot, the data appears to have an upward trend, but it is not clear whether there is seasonality. For time series analysis, the data must be stationary so transformations were investigated to see which could make the data flatter. Thus, the first difference of the data was taken to remove the trend. In this way, the change in data was viewed month to month, rather than the data itself. As evidenced in the second row of graphs, the data is relatively flat with large fluctuations. 

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Temp_Per_Month.PNG)
![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Change_Temp_Per_Mon.PNG)

_Seasonality Plots_

Now that the data appears trend stationary, seasonality was investidated (i.e. if the large fluctuations happen at the same month every year) or if there is no pattern. For this, a seasonality plot was created. As shown below, the data does not appear to have seasonality as temperatures can be higher or lower in the same month depending on the year. 

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Seasonality_NASA1.PNG)

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Seasonality_MET.PNG)

Another seasonal plot called the subseries plot was also explored. In this graph below, the temperature change in each month (i.e. all the Januarys, Februarys, etc.) is connected with a line which represents the averages. As shown below, there are slight temperature changes each month.

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Subseries_Plot.PNG)

_Decomposition_

Finally, the decomposition of the differenced NASA and MET data using STL (season and trend using Loess) was examined. The results are consistent with the analysis above in which the data became trend stationary once the data was transformed. In addition, it is noted that the confidence interval for seasonality is quite large. If there was strong evidence that seasonality was significant, the confidence interval would be much more narrow. 

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Diff_Data.PNG)

KINGSTON DATA

_Time Series Plot_

As per the time series plot below, there appears to be seasonality as there is a repeated pattern year over year. Intuitively, this makes sense given that the Kingston area experiences four different seasons every year. Thus, it is clear the predictive model will need to capture the seasonality feature. There does not appear to be a clear pattern of trend, however this can be confirmed with further analysis below. 

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Temp_Kingston.PNG)

_Seasonality Plot_

Similar to NASA and UK MET, the seasonality plot and subseries plot can be used to visualize seasonality within the data. As evidenced in this plots, temperature is lowest during the winter months, warmer in the spring, highest in the summer months and cools down in the fall.  

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Seasonality_Kingston.PNG)

_Decomposition_

Finally, the decomposition of the Kingston data using STL (season and trend using Loess) was examined. The results are consistent with the analysis above in which the small confidence interval for seasonal supports that there is strong evidence of seasonality. Trend also does not appear to be significant given the wide range of the confidence interval. 

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/decomposition_Kingston.PNG)




