# 2020-Energy-AAPOGFWT

# 2020 - Annual Averaged Power Output Generation for Wind Turbines

In this project, we develop a novel framework for annual averaged power output generation prediction for wind turbines based on different wind speed probability distributions and cubic spline interpolation or logistic regression of power curve modelling for wind turbines. To get the required wind speed data, we recommend to access the open data server of the German Weather Service (DWD - https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/wind/historical/). The data for power curve modelling can be obtained from https://www.wind-turbine-models.com/ for the respective manufacturers' wind turbines.

Since our approach bases on R and GNU Octave, we can easily replace probability density functions or power curve models and our framework is therefore relatively flexible.

Hint: The actual GNU Octave Version only works for a single wind speed data file. We suggest to use the R version of our algorithm with full loop of station data from station data directory.

GNU Octave actually only provides a proof of concept of our R version for one weather station only.
