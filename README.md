# climate_change_timeseries

In 2013, the topic of climate change was the focal point of much debate when reports claimed that global warming had stopped, while other research papers stated the opposite. A case study was conducted using data sources from NASA and the UK Met Office. 

The objective of this project was to repeat the analysis on the NASA and Met Office data but updated for 2021 and perform the following:

•	Forecast global average temperatures through to year 2100 and use the predictions to comment on the concerns about increasing global temperatures. 

•	Develop point predictions, including the 90% confidence intervals for global average temperatures for January and July 2030, 2050 and 2100. 

•	Repeat the analysis above using temperature data for Kingston, ON (postal code K7L3N6).

•	Repeat the analysis using the same data but for the time period covered in the infamous Climate Bet between University of Pennsylvania Professor J. Scott Armstrong, and former U.S. Vice President, Al Gore and comment. 

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

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Seasonality_MET1.PNG)

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

**Metric Selection**

Using the appropriate error metric is a critical decision in this analysis as it will impact which model is ultimately selected to create the predictions. To evaluate the accuracy of the model, the mean absolute error (MAE) and root mean squared error (RMSE) metrics were considered as they are reatively easy to calculate and interpret. To compute the MAE, the mean of the absolute differences between actual temperatures and predicted temperatures was taken. To compute the RMSE, the square root of the mean squared errors (the difference between actual temperatures and predicted temperatures) was taken. A difference between these metrics is that RMSE places greater weight to large errors compared to small errors.

**Forecast the global average temperatures through year 2100**

_Assumptions Made:_

•	The benchmark for average global temperature for the baseline 30-year period for NASA sourced data (1951 through 1980) is 14 deg Celsius.

•	For simplicity and consistency purpose, the benchmark for average global temperature for UK Met sourced data (1961 through 1990) was considered to be 14 deg Celsius as well.

•	The data for NASA and MET both start and end at different periods. To keep the dataset the same size for both and use full years, the data was framed from Jan-1950 to Dec-2020.

•	Removing the earlier periods in the data does not have significant impact in time series models, given the weights assigned to earlier years.

•	The anomalies are converted to temperatures with a reference of 14 deg Celsius for both the NASA and the UK Met datasets.

•	A core assumption in time series modelling is the stationarity of data. For the benchmarking models including mean, naïve and seasonal naïve, the non-stationary feature was addressed using the first difference data when evaluating model fit and forecast accuracy. For the ETS, TBATs and ARIMA models, the ets, tbats, auto.arima functions in the forecast package automatically corrected non-stationary mean and variance and returned the best fit models. 

_Modelling Process:_

In order to determine the best fit for our model, several different approaches were investigated including mean, naïve, seasonal naïve, exponential smoothing model (ETS), trigonometric box-cox autoregressive trend seasonality (TBATS) and autoregressive integrated moving average model (ARIMA). Residual diagnostics were performed on each method and the results were analyzed from the residuals time plot, ACF plot and the histogram of residuals. Inspecting residuals was useful for checking whether the model has adequately captured the information in the data. Specifically, when analyzing each of the plots, the following were inspected:

•	Residuals time plot: Residuals from each model were close to zero and the variation of residuals stays similar across the historical data. Thus, the residual mean and variance were constant which is a core assumption for time series models.

•	ACF plot: ideally, there would be no lines extending beyond the 95% confidence intervals (as indicated by the blue lines). If lines extended beyond these confidence intervals, this would indicate that there is information the model is not capturing. 

•	Histogram: the distribution of the residuals appeared to be normal. Thus, it was assumed that the distribution of possible future values follows a normal distribution as well. 

The lowest value of the metrics of MAE and RMSE were used to determine the model with the best fit. For both NASA and MET, the best model turned out to be ARIMA. Refer to the model fit results in the section below and the ARIMA residual diagnostic plots in the Appendix.  

_Forecast accuracy_

For the benchmarking models (mean, naïve and seasonal naïve), the forecast accuracy was determined by pulling a window (i.e. subset) of the data from January 1950 to December 2000 and creating a model based on the data in this period. Then the predictions for the next 20 years were calculated by setting the horizon equal to 240, compared those 20-year predictions to the actual temperatures from January 2001 to December 2020 and calculated the MAE and RMSE. 

For the ETS, TBATs and ARIMA models, time series cross validation was used to assess the forecast accuracy. The forecast horizon was set to 1 and the length of the rolling window to 600 in order to predict data points starting after the year 2000 (there are 50 years between 1950 and 2000 and 12 months in each year). The error of these predictions were computed and used to calculated the MAE and RMSE. 

Refer to the forecast error results for both NASA and MET in the section below. 

NASA Data: Point forecast of temperature from years 2021 through 2100

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/NASA_Predictions.PNG)

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/NASA%20Results.PNG)

Based on the lowest RMSE and MAE scores in the model analysis and the forecast errors, ARIMA was identified as the better suitable model.

UK MET Data: Point forecast of temperature from years 2021 through 2100

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/UKMET_Predictions.PNG)

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/UK%20MET%20Results.PNG)

ARIMA was identifed as the best model for UK MET because it has the lowest RMSE and MAE for model fit and lowest MAE for forecast errors. Although TBATs had a lower RMSE when evaluating for forecast errors, it was not significantly lower than the RMSE under ARIMA.  

NASA and UK MET Data: Point predications as well as the 90% confidence intervals for the global average temperatures for January and July 2030, 2050, and 2100.

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/NASA_UKMET_Pred_CI.PNG)

Kingston Data: Point predications as well as the 90% confidence intervals for the global average temperatures for January and July 2030, 2050, and 2100.

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Kingston_Results.PNG)

**The Climate Bet**

J.Scott Armstrong, a University of Pennsylvania Professor, wanted to make a bet with Al Gore, former U.S. Vice President, regarding future global average temperatures over a ten year period from 2008 to 2017. The purpose was to raise awareness about the application of scientific methods to forecast climate change which could eventually be used to influence public policy. Armstrong predicted that global mean temperature would not change and his forecasts were based on the naïve model. In this case, he used the most recent year’s average temperature as the forecast for each of the years in the future. 

In the end, Gore did not accept the bet so Armstrong used the U.N. Intergovernmental Panel on Climate Change’s Third Assessment Report in 2001 (IPCC) which predicted a 0.3-degree Celsius warming trend to represent Gore’s position and created theclimatebet.com to track the results. The data used for the bet was based on the satellite temperature data from the University of Alabama at Huntsville (UAH) which measures temperature anomalies of the atmospheric layer from the surface to approximately 10 km, centered in the Lower Troposphere. The cumulative absolute error was the metric used to assess forecasting accuracy. 

The 0.159 degrees Celsius is the average anomaly for the year 2007 and represents Armstrong’s baseline for his bet with Gore.

_Model Retraining and Forecasting_

From external research, the ten years of the bet started on January 2008. Thus, the period of forecast was assumed to be from January 2008 to December 2017. Based on the results, the ARIMA model fits the data better than the naïve model as the MAE and RMSE are lower for both the NASA and MET datasets.

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Climate_Bet_Results.PNG)

_Naive Model_

Using the naïve model applied by Armstrong, the forecasted temperatures for the 10-year period are predicted to be constant at 14.5 degrees Celsius for NASA and 14.419 degrees Celsius for MET. The graphs below use the actual time series data, however, when calculating the model fit, predictions and forecast errors, the first differenced data was used. As discussed in the exploratory data analysis, the data was transformed by taking the first difference to ensure the data was stationary as this is a core assumption in time series models. While non-stationary mean and variance is automatically corrected for ARIMA models using auto.arima, for naïve models this transformation must be applied seperately. 

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Naive_NASA_UKMET.PNG)

_ARIMA Model_

The model was re-trained to cover the period of January 1950 to Dec 2007 for both NASA and MET data. Then, the dataset was converted into a time series object with frequency equal to 12 to capture that there are 12 data points (i.e. months) in each year. Based on the previous analysis, the ARIMA model produced the best fit. Thus, the ARIMA model was fit on the re-shaped NASA and MET datasets. The result was ARIMA (2,1,3)(1,0,0)[12] with drift was the best fit for NASA and ARIMA (4,1,1)(2,0,1)[12] was the best fit for MET. 

In order to predict the following 10 years from 2008 to 2017, the horizon was set to 120 for both NASA and MET. The graphs below capture the predicted temperatures within a confidence interval of 80% to 90%.  

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Climate_BET_ARIMA.PNG)

A closer look at the forecasted temperatures for NASA during the 10-year period reveals that temperatures will drop to its lowest point of 14.55441 degrees Celsius in Jan 2008 and keep increasing until a high of 14.75385 in December 2017. On the other hand, MET forecasted temperatures starts with its lowest point of 14.4782 degrees Celsius in Jan 2008 and appears to have some seasonality in the beginning of the forecasted period but eventually flattens out to approximately 14.54 degrees Celsius closer to 2017. 

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Climate_Bet_NASA_MET.PNG)

_Testing and Forecast Accuracy_

Using the forecasted temperatures based on the naïve and ARIMA models for the period of January 2008 to December 2017, the predictions were compared to the actual temperatures and the errors for that same period were calculated. As cumulative absolute error was the metric to assess accuracy in the Climate Bet, these results were considered along with the chosen metrics of MAE and RMSE. The result is that the ARIMA models have better forecast accuracy than the naïve model for both NASA and MET data.

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Climate_Bet_Results2.PNG)

In the graphs below, the train period (Jan 1950 to Dec 2007) over which the moedel was trained is captured in black, the respective forecasted temperatures for naïve and ARIMA during the test period (Jan 2008 to Dec 2017) are in red and the actual temperatures for the test period are in blue. 

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/Naive_vs_Actual.PNG)

![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/ARIMA_vs_Actual.PNG)

Based on the analysis above, Armstrong would have lost the bet if the NASA and MET data sets were used. The forecast errors based on the naïve model are greater than the errors using an ARIMA model across multiple metrics. 

**Technical Appendix**

NASA Data: ARIMA Model Residual Plots
![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/ARIMA_Residual_Plot_NASA.PNG)

UK Met Data: ARIMA Model Residual Plots
![](https://github.com/joanneevangelista/climate_change_timeseries/blob/main/images/ARIMA_Residual_Plot_UKMET.PNG)







