install.packages("ncdf4")
library("ncdf4")

install.packages("raster")
library(raster)

install.packages("rgdal")
library(rgdal)

install.packages("tidyrverse")
library(tidyverse)

library(lubridate)

install.packages("forecast")
library(forecast)

library(ggplot2)

library(dplyr)

install.packages("smooth")
library(smooth)

#Load data

NASAdata <- read.csv(file.choose(), skip = 1, header = TRUE)

nc <- nc_open("HadCRUT.5.0.1.0.analysis.summary_series.global.monthly.nc")
print(nc)
#10 variables, temperature in Celsius and time is days since 1850-01-01
#2 dimensions - time (size 2053) and bnds (size 2)
#Global attributes - land and see water temperature expressed as monthly anomalies relative to 1961-1990 climatology 


########################################################################################################

#Tidy the data so that the columns represent variables and all rows contain observations

#Tidy NASA data
NASAdata %>% gather(Month, TempAnomaly,Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec ) %>% unite(text_date, Year, Month) %>% mutate(date=ym(text_date)) -> NASAdata

NASAdata <- NASAdata[1:1704,8:9] #pull only the date and temp column

NASAdata <- NASAdata[c("date", "TempAnomaly")] #switch columns Date and TempAnomaly

NASAdata <- arrange(NASAdata,(date)) #Arrange in chronological order

NASAdata <- NASAdata[841:1692, 1:2] #To pull data from Jan 1950 Dec 31 2020

NASAdata$TempAnomaly <- as.numeric(NASAdata$TempAnomaly)
NASAdata$date <- as.Date(NASAdata$date, format = "%Y-%m-%d")

NASAdata$Temp <- NASAdata$TempAnomaly + 14 #30-year base temp per case is 14 degrees C

dim(NASAdata)

#Tidy MET data and format column headers similar to NASA dataset 
attributes(nc$var) #this gives us the variables

t <- ncvar_get(nc,"time")
t <- as.Date(tm, origin = "1850-01-01")
t
dim(t)
head(t)
tail(t)

anom <- ncvar_get(nc,"tas_mean")
dim(anom)
anom

METdata <-data.frame(t,anom)

names(METdata)[1] <- "date"         
names(METdata)[2] <- "TempAnomaly"  

METdata$Temp <- METdata$TempAnomaly + 14 #30-year base temp per case is 14 degrees C
METdata <- METdata[1201:2052,]  #To pull data for time period that matches NASA dataset from Jan 1950 to Dec 2020 

dim(METdata)

########################################################################################################
#NASA DATASET
########################################################################################################

#Define the dataset as a timeseries starting January 1880 and having seasonality frequency 12 (monthly)
NASAdata_ts <- ts(NASAdata$Temp, start = c(1950,1), frequency = 12)
str(NASAdata)

#######################################################################################################

#Preliminary Analysis 

#Time Plot

autoplot(NASAdata_ts) + ggtitle("Time Plot: NASA Average Temperature Per Month") +  
  ylab("Temperature in degrees Celsius") + xlab("Year")
# Data appears to have an upward trend. Investigate transformations. Data should be stationary (i.e. flat). 
#Data also appears to be cyclic because it rises and falls not a fixed frequency
#If the trend is strong relative to seasonality, it make it hard to see if seasonal patterns in data. Transform the data so get rid of trend.
#Take the first difference of the data to remove the trend. This means look at the change in data month to month, rather than data itself

NASAdata_ts_diff <- diff(NASAdata_ts)

#Time Plot of differenced data.

autoplot(NASAdata_ts_diff) + ggtitle("Time Plot: Change in NASA Temperatures Per Month") +  
  ylab("Temperature in degrees Celsius")

#This is the change from month to month. The data is relatively flat with big fluctuations. 
#Series appears trend stationary, use to investigate seasonality. 

#Now we want to see if these big fluctuations happen at the same month every year (seasonality) or if they are irregular

ggseasonplot(NASAdata_ts_diff) + ggtitle("Seasonal Plot: Change in NASA Temps Per Month") + 
  ylab("Temperature in degrees Celsius")

#Let's look at another seasonal plot, the subseries plot
ggsubseriesplot(NASAdata_ts_diff) + ggtitle("Subseries Plot: Change in NASA Temps Per Month") +
  ylab("Temperature in degrees Celsius")
#This is plotting the change in January 1880,1881, etc and connecting with a line. The lines are the averages. 
#Looks like the average change each month does not vary significantly. 


decomp <- stl(NASAdata_ts_diff, t.window=12, s.window="periodic") #decompose using STL (Season and trend using Loess)
plot(decomp) + title("NASA Differenced Data")

####################################################################################################### 
#Our series, NASA Temp Anomalies, has trend and maybe seasonality 
#To remove the trend, we take the first difference.
#The first differenced series does not have seasonality. 

#Plot forecast with various methods
#######################################################################################################

#Use a simple forecasting method such as average method, naive method and seasonal naive
#Average method the forecasts of all future values are equal to the average of the historical data
#Naive method the forecasts of all future values are equal to the last observation 
#Seasonal naive method, we set each forecast to be equal to the last observed value from the same month of the previous year

NASAdata_ts_diff2 <- window(NASAdata_ts_diff, start = c(1950,1), end = c(2000,12)) #Create training data from Jan 1950 to Dec 2000


autoplot(NASAdata_ts_diff2) + autolayer(meanf(NASAdata_ts_diff2, h = 240), series = "Mean", PI = FALSE) +  
  autolayer(naive(NASAdata_ts_diff2, h=240), series = "Naive", PI = FALSE) +
  autolayer(snaive(NASAdata_ts_diff2, h=240), series = "Seasonal Naive", PI = FALSE) +
  ggtitle("Forecasts for changes in monthly NASA temperatures") +
  xlab("Year") + ylab("Degrees Celsius") +
  guides(colour = guide_legend(title = "Forecast"))

#######################################################################################################

#Residual diagnostics

#Residuals in a timeseries are what is left over after fitting the model. 
#Residuals are useful for checking whether a model has adequately captured the information in the data. 
#A good forecasting model will yield residuals with the following properties:

#1) Residuals are uncorrelated (Refer to top graph). 
#If there are correlations between residuals, then there is info left in the residuals which should be used in computing forecast

#2) Residuals have zero mean. If the residuals have a mean other than zero, then the forecasts are based
#If either of these properties is not satisfied, then the forecasting method can be modified to give better forecasts. 

#3)Residuals have constant variance

#4) The residuals are normally distributed (Refer to third graph)

fit <- meanf(NASAdata_ts_diff) #MAE is 0.09431116  RMSE is 0.1198332, std dev is 0.1199037
print(summary(fit))
checkresiduals(fit) 

#Mean of residuals is close to zero and there is no significant correlation in the residual series
#Time plot of residuals shows that the variation of residuals stays similar across the historical data, 
#therefore the residual variance can be treated as a constant.
#ACF of the residuals from the naive method shows lack of correlation suggesting the forecasts are good.
#ACF shows lines above and below the blue line which are the 95% confidence intervals. 
#This suggests there is information that the model is not capturing. 
#ACF plot also shows no seasonality because lines have no pattern
#Histogram of the residuals appears normal. Thus, we assume distribution of possible future values follows normal distribution
#Residual standard deviation is a measure of how well our model is fitting. Smaller #s closer to zero are better. 

fit <- naive(NASAdata_ts_diff) # MAE is 0.1575647, RMSE is 0.1992642, residual sd is 0.1993
print(summary(fit))
checkresiduals(fit)

fit <- snaive(NASAdata_ts_diff) #MAE is 0.1353516, RMSE is 0.1713002, residual sd is 0.1713
print(summary(fit))
checkresiduals(fit)

#Mean model provides the lowest MAE, RMSE and residual standard deviation

#######################################################################################################

#Evaluating forecast accuracy 

#Recall: forecast error is the difference between an observed value and its forecast (i.e. the unpredictable part of an observation)

fit1 <- meanf(NASAdata_ts_diff2,h=240)
fit2 <- naive(NASAdata_ts_diff2,h=240)
fit3 <- snaive(NASAdata_ts_diff2, h=240)

autoplot(window(NASAdata_ts_diff,start = 1950)) + autolayer(fit1, series = "Mean", PI = FALSE) +  
  autolayer(fit2, series = "Naive", PI = FALSE) +
  autolayer(fit3, series = "Seasonal Naive", PI = FALSE) +
  ggtitle("Forecasts for changes in monthly NASA temperature anomalies") +
  xlab("Year") + ylab("Degrees Celsius") +
  guides(colour = guide_legend(title = "Forecast"))

#Compute the forecast accuracy measures 
NASAdata_ts_diff3 <- window(NASAdata_ts_diff, start = 2001) #windom pulls a subset of the time series data
accuracy(fit1, NASAdata_ts_diff3)
accuracy(fit2, NASAdata_ts_diff3)
accuracy(fit3, NASAdata_ts_diff3)

# mean method has the lowest MAE and RMSE

####################################################################################################### 
#Exponential Smoothing Model (ETS)

#There is function tries every possible ETS model and return the one that is best. 
#We can use the data itself because it allows for trend. If computer thinks there is a trend, it will include this
#Thus, we can use the original data, rather than the differenced data

fit_ets <- ets(NASAdata_ts)
print(summary(fit_ets))
checkresiduals(fit_ets)

#The model picked is Multiplicative error, No trend, No Seasonality (MNN)
#Sigma (which is the residual standard deviation) is 0.0076. Smaller number represents a better fit.
#Residual standard deviation is the diff in std devs of observed values vs. predicted values) 
#MAE is 0.0868372
#RMSE is 0.1089844

#Create exponential smoothing model
NASAdata_MNN <- ets(NASAdata_ts, model = "MNN", damped = FALSE)

#Create prediction cone for 80 years to end of 2100 (960 months)
NASAdata_MNN_pred <- forecast(NASAdata_MNN, h=960, level = c(0.8, 0.90))
plot(NASAdata_MNN_pred, xlab = "Year", ylab = "Predicted Temperature")

#Look at what the models actually are -- ETS
NASAdata_MNN

####################################################################################################### 
#ARIMA Model

#Data in this model needs to be stationary! We can use differenced data to get rid of trend but there is also a built in function
#called auto.arima which will try out different arima models and return the one the fits the best
#We can use the regular data because we can account for the trend using d = 1
#Tells R to take first difference of the data before fit the ARIMA model. Does the same thing as differenced data
#D=1 tells R to take the seasonality difference 
#Approximation by default R tries to do things as fast as possible so approximates AIC instead of exact to save time. 
#However we do not need to save time since just trying out one model
#Stepwise - instead of trying out every possible ARIMA model, it will only try a few. 
#Set as false since only trying out one model so don't need to save too much time
#Trace = true will print out all the models it is trying out

fit_arima <- auto.arima(NASAdata_ts)
print(summary(fit_arima))
checkresiduals(fit_arima)
acf(residuals(fit_arima))
#Best model: ARIMA (2,1,3)(1,0,0)[12] with drift 
#Sigma^2 is 0.01124
#Sigma is 0.1060189
#MAE is 0.08326284
#RMSE is 0.105543



####################################################################################################### 
#TBATS

#Create a trigonometric box-cox autoregressive trend seasonality (TBATS) model

fit_tbats <- tbats(NASAdata_ts)
print(summary(fit_tbats))
checkresiduals(fit_tbats)

NASAdata_tbats <- tbats(NASAdata_ts)
NASAdata_tbats_pred <- forecast(NASAdata_tbats, h=960, level = c(0.8,0.90)) # h=960 to forecast to 2100 (inclusive)
par(mfrow=c(1,1))
plot(NASAdata_tbats_pred, xlab = "Year", ylab = "Predicted Temperatures")

# Lets look at what our models actually are -- TBATS
NASAdata_tbats #Sigma is 0.1054583

####################################################################################################### 
#Compare ETS, TBATS and ARIMA models -- Time series Cross Validation (Rolling Horizon Holdout)

f_MNN <- function(y,h) forecast(ets(y,model = "MNN"), h=h)
errors_MNN <- tsCV(NASAdata_ts, f_MNN, h=1, window = 600) #window = 600 because 50 years * 12 months per year

?tsCV

par(mfrow=c(1,1))
plot(errors_MNN, ylab = 'tsCV errors')
abline(0,0)

mean(abs(errors_MNN), na.rm = TRUE) #MAE is 0.09052108
sqrt(mean(errors_MNN^2,na.rm=TRUE)) #RMSE is 0.1118316

f_TBATS <- function(y,h) forecast(tbats(y), h=h)
errors_TBATS <- tsCV(NASAdata_ts, f_TBATS, h=1, window = 600) #window = 600 because 50 years * 12 months per year

plot(errors_MNN, ylab='tsCV errors')
abline(0,0)
lines(errors_TBATS, col="gray")
legend("left", legend = c("CV_error_MNN", "CV_error_TBATS"), col = c("black", "gray"), lty = 1:4)

mean(abs(errors_TBATS), na.rm = TRUE) #MAE is 0.08806389
sqrt(mean(errors_TBATS^2, na.rm = TRUE)) #RMSE is 0.1080439

f_ARIMA <- function(y,h) forecast(auto.arima(y), h=h)
error_ARIMA <- tsCV(NASAdata_ts, f_ARIMA, h=1, window = 600)
mean(abs(error_ARIMA), na.rm=TRUE) #MAE is 0.08759463
sqrt(mean(error_ARIMA^2, na.rm = TRUE)) #RMSE is 0.1098584

NASAdata_arima_pred <- forecast(fit_arima, h=960, level = c(0.8,0.9))
####################################################################################################### 
#For NASA dataset, ARIMA is the best model 

#Print the mean and confidence intervals for the model

NASAdata_arima_pred

#Graph the predictions
autoplot(NASAdata_arima_pred) + ggtitle("Time Plot: NASA Average Temperature Predictions Using ARIMA Model") +  
  ylab("Temperature in degrees Celsius") + xlab("Year")

#Export the selected mode's predictions into a CSV file

write.csv(NASAdata_arima_pred, "C:\\Users\\JO1\\OneDrive - Queen's University\\Desktop\\NASApredictions3.csv")

####################################################################################################### 

#MET DATASET

#######################################################################################################

#Define the dataset as a timeseries starting January 1880 and having seasonality frequency 12 (monthly)

METdata_ts <- ts(METdata$Temp, start = c(1950,1), frequency = 12)
str(METdata)

#######################################################################################################


#Preliminary Analysis

#Time Plot

autoplot(METdata_ts) + ggtitle("Time Plot: MET Average Temperature Per Month") +
  ylab("Temperature in degrees Celsius") + xlab("Year")

METdata_ts_diff <- diff(METdata_ts)

#Time Plot of differenced data

autoplot(METdata_ts_diff) + ggtitle("Time Plot: Change in MET Temperatures Per Month") +  
  ylab("Temperature in degrees Celsius")

ggseasonplot(METdata_ts_diff) + ggtitle("Seasonal Plot: Change in MET Temps Per Month") + 
  ylab("Temperature in degrees Celsius")

ggsubseriesplot(METdata_ts_diff) + ggtitle("Subseries Plot: Change in MET Temps Per Month") +
  ylab("Temperature in degrees Celsius")

decomp <- stl(METdata_ts_diff, t.window=12, s.window="periodic") #decompose using STL (Season and trend using Loess)
plot(decomp) + title("MET Differenced Data")

#Plot forecast with various methods
#######################################################################################################

METdata_ts_diff2 <- window(METdata_ts_diff, start = c(1950,1), end = c(2000,12))

autoplot(METdata_ts_diff2) + autolayer(meanf(METdata_ts_diff2, h = 480), series = "Mean", PI = FALSE) +  
  autolayer(naive(METdata_ts_diff2, h=480), series = "Naive", PI = FALSE) +
  autolayer(snaive(METdata_ts_diff2, h=480), series = "Seasonal Naive", PI = FALSE) +
  ggtitle("Forecasts for changes in monthly MET temperature") +
  xlab("Year") + ylab("Degrees Celsius") +
  guides(colour = guide_legend(title = "Forecast"))

#######################################################################################################

#Residual diagnostics

fit <- meanf(METdata_ts_diff) #MAE is 0.08690928, RMSE is 0.119458, sd is 0.1120116
print(summary(fit))
checkresiduals(fit)

fit <- naive(METdata_ts_diff) #MAE is 0.1447275, RMSE is 0.1853813, sd is 0.1854
print(summary(fit))
checkresiduals(fit)

fit <- snaive(METdata_ts_diff) #MAE is 0.1247288, RSME is 0.1597028, sd is 0.1597
print(summary(fit))
checkresiduals(fit)

#Mean model provides the lowest MAE, RMSE and residual standard deviation

####################################################################################################### 

#Evaluating forecast accuracy

fit1 <- meanf(METdata_ts_diff2,h=240)
fit2 <- naive(METdata_ts_diff2,h=240)
fit3 <- snaive(METdata_ts_diff2, h=240)

autoplot(window(METdata_ts_diff,start = 1950)) + autolayer(fit1, series = "Mean", PI = FALSE) +  
  autolayer(fit2, series = "Naive", PI = FALSE) +
  autolayer(fit3, series = "Seasonal Naive", PI = FALSE) +
  ggtitle("Forecasts for changes in monthly NASA temperature anomalies") +
  xlab("Year") + ylab("Degrees Celsius") +
  guides(colour = guide_legend(title = "Forecast"))

#Compute the forecast accuracy measures 
METdata_ts_diff3 <- window(METdata_ts_diff, start = 2001) #windom pulls a subset of the time series data
accuracy(fit1, METdata_ts_diff3)
accuracy(fit2, METdata_ts_diff3)
accuracy(fit3, METdata_ts_diff3)

# mean method has the lowest MAE and RMSE

####################################################################################################### 

#Exponential Smoothing Model (ETS)

fit_ets_MET <- ets(METdata_ts)
print(summary(fit_ets_MET))
checkresiduals(fit_ets_MET)

#The model picked is Additive error, No trend, No Seasonality (A,N,N)
#Sigma (which is the residual standard deviation) is 10.9405. Smaller number represents a better fit.
#MAE is 0.08026144
#RMSE is 0.102633
#sigma is 0.1028

#Create exponential smoothing model
METdata_ANN <- ets(METdata_ts, model = "ANN", damped = FALSE)

#Create prediction cone for 80 years to 2100 (960 months)
METdata_ANN_pred <- forecast(METdata_ANN, h=960, level = c(0.8, 0.90))
plot(METdata_ANN_pred, xlab = "Year", ylab = "Predicted Temperature")

#Look at what the models actually are -- ETS
METdata_ANN

####################################################################################################### 
#ARIMA Model

fit_arima_MET <- auto.arima(METdata_ts)
print(summary(fit_arima_MET))
checkresiduals(fit_arima_MET)

#Best model: ARIMA (4,1,1)(2,0,1)[12]
#Sigma^2 is 0.009905
#Simga is 0.099523
#MAE is 0.0780025
#RMSE is 0.09899628

####################################################################################################### 
#TBATS

#Create a trigonometric box-cox autoregressive trend seasonality (TBATS) model

fit_tbats_MET <- tbats(METdata_ts)
print(summary(fit_tbats_MET))
checkresiduals(fit_tbats_MET)

METdata_tbats <- tbats(METdata_ts)
METdata_tbats_pred <- forecast(METdata_tbats, h=960, level = c(0.8,0.90)) # h=960 to forecast to 2100
par(mfrow=c(1,1))
plot(METdata_tbats_pred, xlab = "Year", ylab = "Predicted Temperatures")

# Lets look at what our models actually are -- TBATS
METdata_tbats #Sigma is 0.101079

####################################################################################################### 
#Compare ETS, TBATS, and ARIMA models -- Time series Cross Validation (Rolling Horizon Holdout)

f_ANN_MET <- function(y,h) forecast(ets(y,model = "ANN"), h=h)
errors_ANN_MET <- tsCV(METdata_ts, f_ANN_MET, h=12, window = 600) #window = 600 because 50 years * 12 months per year

par(mfrow=c(1,1))
plot(errors_ANN, ylab = 'tsCV errors')
abline(0,0)

mean(abs(errors_ANN_MET), na.rm = TRUE)#MAE is 0.1133246
sqrt(mean(errors_ANN_MET^2, na.rm = TRUE)) #RMSE is 0.1452825

f_TBATS_MET <- function(y,h) forecast(tbats(y), h=h)
errors_TBATS_MET <- tsCV(METdata_ts, f_TBATS_MET, h=1, window = 600)

plot(errors_MNN_MET, ylab='tsCV errors')
abline(0,0)
lines(errors_TBATS_MET, col="gray")
legend("left", legend = c("CV_error_MNN", "CV_error_TBATS"), col = c("black", "gray"), lty = 1:4)

mean(abs(errors_TBATS_MET), na.rm = TRUE) #MAE is 0.07975525
sqrt(mean(errors_TBATS_MET^2, na.rm = TRUE)) #RMSE 0.09880734

f_ARIMA_MET <- function(y,h) forecast(auto.arima(y), h=h)
error_ARIMA_MET <- tsCV(METdata_ts, f_ARIMA, h=1, window = 600)

mean(abs(error_ARIMA_MET), na.rm=TRUE) #MAE is 0.07961633
sqrt(mean(error_ARIMA_MET^2,na.rm = TRUE)) #RMSE is 0.09935493

METdata_arima_pred <- forecast(fit_arima_MET, h=960, level = c(0.8,0.9))

####################################################################################################### 
#For MET dataset, ARIMA is the best model 

#Print the mean and confidence intervals for the model

METdata_arima_pred

#Graph the predictions
autoplot(METdata_arima_pred) + ggtitle("Time Plot: UK MET Average Temperature Predictions Using ARIMA Model") +  
  ylab("Temperature in degrees Celsius") + xlab("Year")

#Export the selected mode's predictions into a CSV file

write.csv(METdata_arima_pred, "C:\\Users\\JO1\\OneDrive - Queen's University\\Desktop\\METpredictions3.csv")

####################################################################################################### 
#KINGSTON DATASET

library(readxl) 
library(lubridate)

Kingstondata <- read.csv(file.choose(), header = TRUE)
str(Kingstondata)

Kingstondata %>% unite(text_date, Year, Month) %>% mutate(date=ym(text_date)) -> Kingstondata
Kingstondata <- Kingstondata[229:564,c(11,29)] #pull a full year starting Jan 1979 to Dec 2006
names(Kingstondata)[1] <-"Temp"
Kingstondata <- Kingstondata[c("date","Temp")]

#Define the dataset as a timeseries starting January 1979 and having seasonality frequency 12 (monthly)
Kingstondata_ts <- ts(Kingstondata$Temp, start = c(1979,1), frequency = 12)

#######################################################################################################

#Preliminary Analysis 

#Time Plot

autoplot(Kingstondata_ts) + ggtitle("Time Plot: Kingston Average Temperature per Month") + 
  ylab("Temperature in degrees Celsius")
#data appears stationary

ggseasonplot(Kingstondata_ts) + ggtitle("Seasonal Plot: Kingston Temps per Month") + 
  ylab("Temperature in degrees Celsius")

#Seasonality exists in the Kingston dataset

ggsubseriesplot(Kingstondata_ts) + ggtitle("Subseries Plot: Change in Kingston Temps per Month") +
  ylab("Temperature in degrees Celsius")

decomp <- stl(Kingstondata_ts, t.window=12, s.window="periodic") #decompose using STL (Season and trend using Loess)
plot(decomp) + title("Kingston Data")


Kingstondata_ts2 <- window(Kingstondata_ts, start = c(1979,1), end = c(1993,12)) #Create training data

autoplot(Kingstondata_ts2) + autolayer(meanf(Kingstondata_ts2, h = 156), series = "Mean", PI = FALSE) +  
  autolayer(naive(Kingstondata_ts2, h=156), series = "Naive", PI = FALSE) +
  autolayer(snaive(Kingstondata_ts2, h=156), series = "Seasonal Naive", PI = FALSE) +
  ggtitle("Forecasts for changes in monthly Kingston temperatures") +
  xlab("Year") + ylab("Degrees Celsius") +
  guides(colour = guide_legend(title = "Forecast"))

fit <- meanf(Kingstondata_ts) #MAE is 8.8035, RMSE is 10.0662, sd is 10.0812
print(summary(fit))
checkresiduals(fit) 

fit <- naive(Kingstondata_ts) #MAE is 5.073134, RMSE is 5.742195, sd is 5.7422
print(summary(fit))
checkresiduals(fit) 

fit <- snaive(Kingstondata_ts) #MAE is 2.045679, RMSE is 2.747951, sd is 2.748
print(summary(fit))
checkresiduals(fit) 

#Snaive provides the best benchmarking model'

########################################################################################################

#Exponential Smoothing Model (ETS)

fit_ets_kingston <- ets(Kingstondata_ts)
print(summary(fit_ets_kingston))
checkresiduals(fit_ets_kingston) 

#sigma is 1.9604
#MAE is 1.428686
#RMSE is 1.919168

#The model picked is Additive error, no trend and Additive seasonality 

#create exponential smoothing model
Kingstondata_ANA <- ets(Kingstondata_ts, model = "ANA", damped = FALSE)

Kingstondata_ANA_pred <- forecast(Kingstondata_ANA, h = 1128, c(0.8,0.9)) 
plot(Kingstondata_ANA_pred, xlab = "Year", ylab = "Predicted Temperature")
Kingstondata_ANA

####################################################################################################### 
#ARIMA Model

fit_arima_kingston <- auto.arima(Kingstondata_ts)
print(summary(fit_arima_kingston))
checkresiduals(fit_arima_kingston)

#Best model: ARIMA (2,0,2)(2,1,0)[12]
#Sigma^2 is 5.14
#Sigma is 2.267157
#MAE is 1.596699
#RMSE is 2.202145

####################################################################################################### 
#TBATS

fit_tbats_kingston <- tbats(Kingstondata_ts)
print(summary(fit_tbats_kingston))
checkresiduals(fit_tbats_kingston)

Kingstondata_tbats <- tbats(Kingstondata_ts)
Kingstondata_tbats_pred <- forecast(Kingstondata_tbats, h=1128, c(0.8,0.9))
par(mfrow=c(1,1))
plot(Kingstondata_tbats_pred, xlab = "Year", ylab = "Predicted Temperatures")

Kingstondata_tbats

#Sigma is 1.875174

####################################################################################################### 
#Compare ETS, TBATS and ARIMA models -- Time series Cross Validation (Rolling Horizon Holdout)

f_ANA_kingston <- function(y,h) forecast(ets(y,model = "ANA"),h=h)
errors_ANA_kingston <- tsCV(Kingstondata_ts, f_ANA, h=1, window = 162)

par(mfrow=c(1,1))
plot(errors_ANA_kingston, ylab = 'tsCV errors')
abline(0,0)

mean(abs(errors_ANA_kingston), na.rm = TRUE) #MAE is 1.5438
sqrt(mean(errors_ANA_kingston^2, na.rm = TRUE)) #RMSE is 1.985982

f_TBATS_kingston <- function(y,h) forecast(tbats(y),h=h)
errors_TBATS_kingston <- tsCV(Kingstondata_ts,f_TBATS_kingston,h=1,window = 162)


plot(errors_ANA_kingston, ylab='tsCV errors')
abline(0,0)
lines(errors_TBATS_kingston, col="gray")
legend("left", legend = c("CV_error_ANA", "CV_error_TBATS"), col = c("black", "gray"), lty = 1:4)

mean(abs(errors_TBATS_kingston), na.rm = TRUE) #MAE is 1.566274
sqrt(mean(errors_TBATS_kingston^2, na.rm = TRUE)) #RMSE is 1.991478

f_ARIMA_kingston <- function(y,h) forecast(auto.arima(y), h=h)
error_ARIMA_kingston <- tsCV(Kingstondata_ts, f_ARIMA_kingston, h=1, window = 162)

mean(abs(error_ARIMA_kingston), na.rm=TRUE) #MAE is 1.705851
sqrt(mean(error_ARIMA_kingston^2, na.rm = TRUE)) #RMSE is 2.2938

#Based on sigma of model fit, choose TBATS

Kingstondata_ANA_pred

write.csv(Kingstondata_ANA_pred, "C:\\Users\\JO1\\OneDrive - Queen's University\\Desktop\\Kingstonpredictions3.csv")

####################################################################################################### 
#The Climate Bet: Question 6

#Reshape the data and visuals
NASAdataQ6 <- NASAdata[1:816,1:3]

NASAdataQ6_ts <- ts(NASAdataQ6$Temp, start = c(1950,1), frequency = 12)

autoplot(NASAdataQ6_ts) + ggtitle("Time Plot: NASA Average Temperature Per Month") +
  ylab("Temperature in degrees Celsius") + xlab("Year")

NASAdataQ6_ts_diff <- diff(NASAdataQ6_ts)

autoplot(NASAdataQ6_ts_diff) + ggtitle("Time Plot: Chnage in NASA Temperatures per Month") +
  ylab("Temperature in degrees Celsius")


NASAdataQ6_ts_diff2 <- window(NASAdataQ6_ts_diff, start = c(1950,1), end = c(2007,12))

autoplot(NASAdataQ6_ts_diff2) + autolayer(naive(NASAdataQ6_ts_diff2, h=120), series = "Naive", PI = FALSE) +
  ggtitle("Naive forecasts for changes in monthly NASA temperatures") + xlab("Year") + 
  ylab("Degrees Celsius") + guides(colour=guide_legend(title = "Forecast"))


#Model Fit metrics
fit1Q6_NASA <- naive(NASAdataQ6_ts_diff2, h=120) #because 2007 to 2017 is 10 years * 12 months per year
print(summary(fit1Q6_NASA)) #MAE is 0.1599712 and RMSE is 0.202261


train_NASA <- window(NASAdataQ6_ts, start = c(1950,1), end = c(2007,12))
test_NASA <- window(NASAdataQ6_ts, start = c(2008,1), end = c(2017,12))
test_NASA

fit1Q6_NASA1 <- naive(train_NASA, h=120)
fit1Q6_NASA1

autoplot(window(train_NASA, start = 1950)) + autolayer(fit1Q6_NASA1, series = "Naive", PI = FALSE) + 
  ggtitle("Naive forecasts for NASA temperatures") + xlab("Year") + ylab("Degrees Celsius") 


autoplot(window(train_NASA, start = 1950)) + autolayer(fit1Q6_NASA1, series = "1. Naive", PI = FALSE) + 
  ggtitle("Naive Forecast vs. Actual for NASA temperatures") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(test_NASA, series = "2. Actual")

NASA_naive_accuracy <- data.frame(fit1Q6_NASA1, test_NASA)
NASA_naive_CAE <- sum(abs(NASA_naive_accuracy$Point.Forecast - NASA_accuracy$test_NASA))
NASA_naive_CAE
NASA_naive_MAE <- mean(abs(NASA_naive_accuracy$Point.Forecast - NASA_naive_accuracy$test_NASA))
NASA_naive_MAE
NASA_naive_RMSE <- sqrt(mean(NASA_naive_accuracy$Point.Forecast - NASA_naive_accuracy$test_NASA)^2)
NASA_naive_RMSE

#Reshape the data and visuals

METdataQ6 <- METdata[1:816, 1:3]
METdataQ6_ts <- ts(METdataQ6$Temp, start = c(1950,1), frequency = 12)

METdataQ6_ts_diff <- diff(METdataQ6_ts)
METdataQ6_ts_diff2 <- window(METdataQ6_ts_diff, start = c(1950,1), end = c(2007,12))


#Model Fit metrics
fit1Q6_MET <- naive(METdataQ6_ts_diff2, h=120)
print(summary(fit1Q6_MET)) #MAE is 0.1484013 and RMSE is 0.1902824


train_MET <- window(METdataQ6_ts, start = c(1950,1), end = c(2007,12))
test_MET <- window(METdataQ6_ts, start = c(2008,1), end = c(2017,12))

fit1Q6_MET1 <- naive(train_MET, h=120)
autoplot(window(train_MET, start = 1950)) + autolayer(fit1Q6_MET1, series = "Naive", PI = FALSE) +
  ggtitle("Naive forecasts for MET temperatures") + xlab("Year") + ylab("Degrees Celsius")

autoplot(window(train_MET, start = 1950)) + autolayer(fit1Q6_MET1, series = "1. Naive", PI = FALSE) +
  ggtitle("Naive Forecasts vs. Actual for MET temperatures") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(test_MET, series = "2. Actual")


fit1Q6_MET1

MET_naive_accuracy <- data.frame(fit1Q6_MET1, test_MET)
MET_naive_CAE <- sum(abs(MET_naive_accuracy$Point.Forecast - MET_accuracy$test_MET))
MET_naive_CAE
MET_naive_MAE <- mean(abs(MET_naive_accuracy$Point.Forecast - MET_naive_accuracy$test_MET))
MET_naive_MAE
MET_naive_RMSE <- sqrt(mean(MET_naive_accuracy$Point.Forecast - MET_naive_accuracy$test_MET)^2)
MET_naive_RMSE


#Compute the forecast accuracy measures
NASAdataQ6_ts_diff3 <- window(NASAdataQ6_ts_diff, start = 2007)
accuracy(fit1Q6_NASA, NASAdataQ6_ts_diff3) #Metrics on test set: MAE is 0.0982500, RMSE is 0.1308530


METdataQ6_ts_diff3 <- window(METdataQ6_ts_diff, start = 2007)
accuracy(fit1Q6_MET, METdataQ6_ts_diff3) #Metrics on test set: MAE is 0.09062094, RMSE is 0.1181082


#Climate Bet: ARIMA
NASAdataQ6_ARIMA <- auto.arima(NASAdataQ6_ts)
NASAdataQ6_ARIMA_pred <- forecast(NASAdataQ6_ARIMA,h=120,level = c(0.8,0.9))
autoplot(NASAdataQ6_ARIMA_pred, xlab = "Year", ylab = "Temperatures") + ggtitle("ARIMA forecasts for NASA temperatures")

fitQ6_arima_NASA <- auto.arima(train_NASA)
print(summary(fitQ6_arima_NASA)) 
#MAE is 0.08218317 and RMSE is 0.104383
#Best model: ARIMA(2,1,3)(1,0,0)

train_NASA_ARIMA_pred <- forecast(fitQ6_arima_NASA, h=120, level = c(0.8,0.9))
plot(train_NASA_ARIMA_pred, xlab = "Year", ylab = "NASA Predicted Temperatures")
train_NASA_ARIMA_pred
train_NASA_ARIMA_pred$mean

autoplot(train_NASA_ARIMA_pred, xlab = "Year", ylab = "Temperatures") + ggtitle("ARIMA forecasts for NASA temperatures")+
  autolayer(test_NASA)

autoplot(window(train_NASA, start = 1950)) + autolayer(train_NASA_ARIMA_pred, series = "1.ARIMA", PI = FALSE) +
  ggtitle("ARIMA Forecasts for NASA temperatures") + xlab("Year") + ylab("Degrees Celsius")

autoplot(window(train_NASA, start = 1950)) + autolayer(train_NASA_ARIMA_pred, series = "1.ARIMA", PI = FALSE) +
  ggtitle("ARIMA Forecasts vs. Actual for NASA temperatures") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(test_NASA, series = "2.Actual")

write.csv(train_NASA_ARIMA_pred, "C:\\Users\\JO1\\OneDrive - Queen's University\\Desktop\\NASAQ6v2.csv")


train_NASA <- window(NASAdataQ6_ts, start = c(1950,1), end = c(2007,12))
test_NASA <- window(NASAdataQ6_ts, start = c(2008,1), end = c(2017,12))


METdataQ6_ARIMA <- auto.arima(METdataQ6_ts)
METdataQ6_ARIMA_pred <- forecast(METdataQ6_ARIMA, h=120, level = c(0.8,0.9))
autoplot(METdataQ6_ARIMA_pred, xlab = "Year", ylab = "Temperatures") + ggtitle("ARIMA forecasts for MET temperatures")

fitQ6_arima_MET <- auto.arima(train_MET)
print(summary(fitQ6_arima_MET))

train_MET_ARIMA_pred <- forecast(fitQ6_arima_MET, h=120, level = c(0.8,0.9))

plot(train_MET_ARIMA_pred, xlab = "Year", ylab = "MET Predicted Temperatures")
train_MET_ARIMA_pred

autoplot(window(train_MET, start = 1950)) + autolayer(train_MET_ARIMA_pred, series = "1.ARIMA", PI = FALSE) +
  ggtitle("ARIMA Forecasts vs. Actual for MET temperatures") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(test_MET, series = "2.Actual")


write.csv(train_MET_ARIMA_pred, "C:\\Users\\JO1\\OneDrive - Queen's University\\Desktop\\METQ6v2.csv")


NASA_accuracy <- data.frame(train_NASA_ARIMA_pred, test_NASA)
NASA_CAE <- abs(NASA_accuracy$Point.Forecast - NASA_accuracy$test_NASA)
NASA_CAE
sum(NASA_CAE)
NASA_MAE <- mean(abs(NASA_accuracy$Point.Forecast - NASA_accuracy$test_NASA))
NASA_MAE
NASA_RMSE <- sqrt(mean(NASA_accuracy$Point.Forecast - NASA_accuracy$test_NASA)^2)
NASA_RMSE


MET_accuracy <- data.frame(train_MET_ARIMA_pred, test_MET)
MET_CAE <- sum(abs(MET_accuracy$Point.Forecast - MET_accuracy$test_MET))
MET_CAE
MET_MAE <- mean(abs(MET_accuracy$Point.Forecast - MET_accuracy$test_MET))
MET_MAE
MET_RMSE <- sqrt(mean(MET_accuracy$Point.Forecast - MET_accuracy$test_MET)^2)
MET_RMSE

write.csv(NASA_accuracy, "C:\\Users\\JO1\\OneDrive - Queen's University\\Desktop\\NASAQ6_Accuracy.csv")


####################################################################################################### 
#The Climate Bet: Question 7 
##Q7-----------------------------------------------------------------------------

# Repeat the same analyses starting in year 1999 with 10- and 20-year time-intervals. 
# Comment on what you observe. 

# re-train the models until 1999, have two models (NASA and MET) and four comparisons 
# (i.e. comparing 10 years from 1999, then 20 years from 1999 for both datasets). 

#NASA reshape
NASA_1 <- read.csv(choose.files(), sep = ",")

NASA_1 <- NASA_1%>%
  gather(month, Temp, -Year) %>%
  arrange(Year) %>%
  mutate(month = str_remove(month, "^X"))


NASA_Q7 <- NASA_1[-c(1441:1740,1:840),] # reshape from 1950-01 to 1999-12
NASA_Q7$Actual_Temp <- NASA_Q7$Temp + 14

NASA_ts.Q7 <- ts(NASA_Q7$Actual_Temp,start = 1950, frequency = 12)
plot(NASA_ts.Q7)

NASA_ts_diff.Q7 <- diff(NASA_ts.Q7)
plot(NASA_ts_diff.Q7)

#MET reshape
nc <- nc_open("HadCRUT.5.0.1.0.analysis.summary_series.global.monthly.nc")
attributes(nc$var)

t <- ncvar_get(nc,"time")
t <- as.Date(t,origin = "1850-01-01", tz = "UTC")
MET_anom <- ncvar_get(nc,"tas_mean")

MET_1 <- data.frame(t,MET_anom)
MET_Q7 <- MET_1[-c(1:1200,1801:2053),] # reshape from 1950-01 to 1999-12
MET_Q7$actual_Temp <- MET_Q7$MET_anom + 14

MET_ts.Q7 <- ts(MET_Q7$actual_Temp,start = 1950, frequency = 12)
plot(MET_ts.Q7)

MET_ts_diff.Q7 <- diff(MET_ts.Q7)
plot(MET_ts_diff.Q7)

#retrain on ARIMA model

fit_arima_NASA.Q7 <- auto.arima(NASA_ts.Q7)
print(summary(fit_arima_NASA.Q7))
checkresiduals(fit_arima_NASA.Q7)
# Best model: ARIMA(2,1,3)(1,0,0)[12]
# MAE is 0.08290307 and RMSE is 0.1054395

fit_arima_MET.Q7 <- auto.arima(MET_ts.Q7)
print(summary(fit_arima_MET.Q7))
checkresiduals(fit_arima_MET.Q7)
# Best model: ARIMA(4,1,1)(2,0,1)[12]
# MAE is 0.07843602 and RMSE is 0.1005348


#retrain naive model

fit_naive_NASA.Q7 <- naive(NASA_ts_diff.Q7)
print(summary(fit_naive_NASA.Q7))
checkresiduals(fit_naive_NASA.Q7)
# MAE is 0.1573579 and RMSE is 0.1993124

fit_naive_MET.Q7 <- naive(MET_ts_diff.Q7)
print(summary(fit_naive_MET.Q7))
checkresiduals(fit_naive_MET.Q7)
# MAE is 0.1466941 and RMSE is 0.1887985



#CV on ARIMA model

f_NASA_ARIMA.Q7  <- function(y, h) forecast(auto.arima(y), h = h)
errors_NASA_ARIMA.Q7 <- tsCV(NASA_ts.Q7, f_NASA_ARIMA.Q7, h=1)
mean(abs(errors_NASA_ARIMA.Q7),na.rm = TRUE) # MAE 0.08579599
sqrt(mean(errors_NASA_ARIMA.Q7^2,na.rm = TRUE)) # RMSE 0.1097137

f_MET_ARIMA.Q7  <- function(y, h) forecast(auto.arima(y), h = h)
errors_MET_ARIMA.Q7 <- tsCV(MET_ts.Q7, f_MET_ARIMA.Q7, h=1)
mean(abs(errors_MET_ARIMA.Q7),na.rm = TRUE) # MAE 0.08162268
sqrt(mean(errors_MET_ARIMA.Q7^2,na.rm = TRUE)) # RMSE 0.1054315

#reshape the test dataset, 10 years test data from 2000-2009, 20 years test data from 2000-2019
NASA_test_10yrs <- window(NASA_ts, start = 2000, end = c(2009,12))
NASA_test_20yrs <- window(NASA_ts, start = 2000, end = c(2019,12))
write.csv(NASA_test_20yrs, "E:\\Queens University\\MMA867-Predictive Modelling\\assignment 2\\NASA_test_20yrs.csv")


MET_test_10yrs <- window(MET_ts, start = 2000, end = c(2009,12))
MET_test_20yrs <- window(MET_ts, start = 2000, end = c(2019,12))
write.csv(MET_test_20yrs, "E:\\Queens University\\MMA867-Predictive Modelling\\assignment 2\\MET_test_20yrs.csv")


# prediction assessment, predict 10 yrs, 20 yrs future

#NASA---------------------------------------------------------------------------
#ARIMA
NASA_Q7_10yr_arima <- forecast(fit_arima_NASA.Q7, h=120, level = c(0.8,0.90)) # h=120 to forecast to 2009 (inclusive)
NASA_Q7_20yr_arima <- forecast(fit_arima_NASA.Q7, h=240, level = c(0.8,0.90)) # h=240 to forecast to 2019 (inclusive)
plot(NASA_Q7_10yr_arima, xlab = "Year", ylab = "Predicted Temperatures")
plot(NASA_Q7_20yr_arima, xlab = "Year", ylab = "Predicted Temperatures")

write.csv(NASA_Q7_10yr, "E:\\Queens University\\MMA867-Predictive Modelling\\assignment 2\\NASA_Q7_10yr_arima.csv")
write.csv(NASA_Q7_20yr, "E:\\Queens University\\MMA867-Predictive Modelling\\assignment 2\\NASA_Q7_20yr_arima.csv")

#NAIVE
NASA_Q7_10yr_naive <- naive(NASA_ts.Q7, h = 120)
autoplot(NASA_ts.Q7) + autolayer(NASA_Q7_10yr_naive, series = "NAIVE", PI = FALSE) +
  ggtitle("NAIVE forecasts for NASA temperature") + xlab("Year") + ylab("Degrees Celsius")

NASA_Q7_20yr_naive <- naive(NASA_ts.Q7, h = 240)
autoplot(NASA_ts.Q7) + autolayer(NASA_Q7_20yr_naive, series = "NAIVE", PI = FALSE) +
  ggtitle("NAIVE forecasts for NASA temperature") + xlab("Year") + ylab("Degrees Celsius")

write.csv(NASA_Q7_10yr_naive, "E:\\Queens University\\MMA867-Predictive Modelling\\assignment 2\\NASA_Q7_10yr_naive.csv")
write.csv(NASA_Q7_20yr_naive, "E:\\Queens University\\MMA867-Predictive Modelling\\assignment 2\\NASA_Q7_20yr_naive.csv")

#MET----------------------------------------------------------------------------
#ARIMA
MET_Q7_10yr_arima <- forecast(fit_arima_MET.Q7, h=120, level = c(0.8,0.90)) # h=120 to forecast to 2009 (inclusive)
MET_Q7_20yr_arima <- forecast(fit_arima_MET.Q7, h=240, level = c(0.8,0.90)) # h=240 to forecast to 2019 (inclusive)
plot(MET_Q7_10yr, xlab = "Year", ylab = "Predicted Temperatures")
plot(MET_Q7_20yr, xlab = "Year", ylab = "Predicted Temperatures")

write.csv(MET_Q7_10yr, "E:\\Queens University\\MMA867-Predictive Modelling\\assignment 2\\MET_Q7_10yr_arima.csv")
write.csv(MET_Q7_20yr, "E:\\Queens University\\MMA867-Predictive Modelling\\assignment 2\\MET_Q7_20yr_arima.csv")

#NAIVE
MET_Q7_10yr_naive <- naive(MET_ts.Q7, h = 120)
autoplot(MET_ts.Q7) + autolayer(MET_Q7_10yr_naive, series = "NAIVE", PI = FALSE) +
  ggtitle("NAIVE forecasts for MET temperature") + xlab("Year") + ylab("Degrees Celsius")

MET_Q7_20yr_naive <- naive(MET_ts.Q7, h = 240)
autoplot(MET_ts.Q7) + autolayer(MET_Q7_20yr_naive, series = "NAIVE", PI = FALSE) +
  ggtitle("NAIVE forecasts for MET temperature") + xlab("Year") + ylab("Degrees Celsius")

write.csv(MET_Q7_10yr_naive, "E:\\Queens University\\MMA867-Predictive Modelling\\assignment 2\\MET_Q7_10yr_naive.csv")
write.csv(MET_Q7_20yr_naive, "E:\\Queens University\\MMA867-Predictive Modelling\\assignment 2\\MET_Q7_20yr_naive.csv")







#Testing and foreacast accuracy
test_nasa <- read.csv(choose.files(),header = TRUE,sep = ",")
test_met <- read.csv(choose.files(),header = TRUE,sep = ",")
#NASA---------------------------------------------------------------------------

#NAIVE
autoplot(NASA_ts.Q7) + autolayer(NASA_Q7_10yr_naive, series = "1. NAIVE", PI = FALSE) +
  ggtitle("NAIVE forecasts for NASA temperature") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(NASA_test_10yrs, series = "2. Actual")
mean(abs(test_nasa$TEST_10-test_nasa$NAÏVE_10YR),na.rm = TRUE) #MAE 0.1969167
sum(abs(test_nasa$TEST_10-test_nasa$NAÏVE_10YR),na.rm = TRUE) #Sum abs ERROR 23.63
sqrt(mean((test_nasa$TEST_10-test_nasa$NAÏVE_10YR)^2,na.rm = TRUE)) #RMSE 0.2242785


autoplot(NASA_ts.Q7) + autolayer(NASA_Q7_20yr_naive, series = "1. NAIVE", PI = FALSE) +
  ggtitle("NAIVE forecasts for NASA temperature") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(NASA_test_20yrs, series = "2. Actual")
mean(abs(test_nasa$TEST_20-test_nasa$NAÏVE_20YR),na.rm = TRUE) #MAE 0.2980417
sum(abs(test_nasa$TEST_20-test_nasa$NAÏVE_20YR),na.rm = TRUE) #Sum abs ERROR 71.53
sqrt(mean((test_nasa$TEST_20-test_nasa$NAÏVE_20YR)^2,na.rm = TRUE)) #RMSE 0.3475492


#ARIMA
autoplot(NASA_ts.Q7) + autolayer(NASA_Q7_10yr_arima, series = "1. ARIMA", PI = FALSE) +
  ggtitle("ARIMA forecasts for NASA temperature") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(NASA_test_10yrs, series = "2. Actual")
mean(abs(test_nasa$TEST_10-test_nasa$ARIMA_10YR),na.rm = TRUE) #MAE 0.1900118
sum(abs(test_nasa$TEST_10-test_nasa$ARIMA_10YR),na.rm = TRUE) #Sum abs ERROR 22.80141
sqrt(mean((test_nasa$TEST_10-test_nasa$ARIMA_10YR)^2,na.rm = TRUE)) #RMSE 0.2173874


autoplot(NASA_ts.Q7) + autolayer(NASA_Q7_20yr_arima, series = "1. ARIMA", PI = FALSE) +
  ggtitle("ARIMA forecasts for NASA temperature") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(NASA_test_20yrs, series = "2. Actual")
mean(abs(test_nasa$TEST_20-test_nasa$ARIMA_20YR),na.rm = TRUE) #MAE 0.2902465
sum(abs(test_nasa$TEST_20-test_nasa$ARIMA_20YR),na.rm = TRUE) #Sum abs ERROR 69.65916
sqrt(mean((test_nasa$TEST_20-test_nasa$ARIMA_20YR)^2,na.rm = TRUE)) #RMSE 0.3403519


#MET----------------------------------------------------------------------------

#NAIVE
autoplot(MET_ts.Q7) + autolayer(MET_Q7_10yr_naive, series = "1. NAIVE", PI = FALSE) +
  ggtitle("NAIVE forecasts for MET temperature") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(MET_test_10yrs, series = "2. Actual")
mean(abs(test_met$TEST_10-test_met$NAÏVE_10YR),na.rm = TRUE) #MAE 0.218278
sum(abs(test_met$TEST_10-test_met$NAÏVE_10YR),na.rm = TRUE) #Sum abs ERROR 26.19336
sqrt(mean((test_met$TEST_10-test_met$NAÏVE_10YR)^2,na.rm = TRUE)) #RMSE 0.2426209


autoplot(MET_ts.Q7) + autolayer(MET_Q7_20yr_naive, series = "1. NAIVE", PI = FALSE) +
  ggtitle("NAIVE forecasts for MET temperature") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(MET_test_20yrs, series = "2. Actual")
mean(abs(test_met$TEST_20-test_met$NAÏVE_20YR),na.rm = TRUE) #MAE 0.3196976
sum(abs(test_met$TEST_20-test_met$NAÏVE_20YR),na.rm = TRUE) #Sum abs ERROR 76.72743
sqrt(mean((test_met$TEST_20-test_met$NAÏVE_20YR)^2,na.rm = TRUE)) #RMSE 0.3629082


#ARIMA
autoplot(MET_ts.Q7) + autolayer(MET_Q7_10yr_arima, series = "1. ARIMA", PI = FALSE) +
  ggtitle("ARIMA forecasts for MET temperature") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(MET_test_10yrs, series = "2. Actual")
mean(abs(test_met$TEST_10-test_met$ARIMA_10YR),na.rm = TRUE) #MAE 0.1778529
sum(abs(test_met$TEST_10-test_met$ARIMA_10YR),na.rm = TRUE) #Sum abs ERROR 21.34235
sqrt(mean((test_met$TEST_10-test_met$ARIMA_10YR)^2,na.rm = TRUE)) #RMSE 0.2019202


autoplot(MET_ts.Q7) + autolayer(MET_Q7_20yr_arima, series = "1. ARIMA", PI = FALSE) +
  ggtitle("ARIMA forecasts for MET temperature") + xlab("Year") + ylab("Degrees Celsius") +
  autolayer(MET_test_20yrs, series = "2. Actual")
mean(abs(test_met$TEST_20-test_met$ARIMA_20YR),na.rm = TRUE) #MAE 0.2747039
sum(abs(test_met$TEST_20-test_met$ARIMA_20YR),na.rm = TRUE) #Sum abs ERROR 65.92893
sqrt(mean((test_met$TEST_20-test_met$ARIMA_20YR)^2,na.rm = TRUE)) #RMSE 0.3209652







