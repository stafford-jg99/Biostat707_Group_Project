# Biostat707_Group_Project
The Repository for the 707 Project

# Project Goals:
## Project Overview & Predictive Objective
* Predict whether tomorrow’s PM2.5 levels will exceed a predefined health threshold by modeling hourly pollutant and weather readings from a fixed 24- or 48-hour observation window.
* Evaluate and quantify the optimistic performance bias caused by temporal leakage by actively comparing naive random cross-validation results against a strict temporal validation design.

## Dataset Characteristics
* Source: The Beijing Multi-Site Air Quality dataset, sourced from the UCI Machine Learning Repository.
* Structure: Multivariate, time-series data comprising 420,768 total instances. Feature types consist of real and integer values.
* Variables: The data includes 6 main air pollutants and 6 relevant meteorological variables.
* Collection Sites: Hourly data was collected from 12 nationally-controlled air-quality monitoring sites across Beijing. The pollutant readings were supplied by the Beijing Municipal Environmental Monitoring Center and matched with meteorological data from the nearest weather stations.

## Validation & Leakage Prevention Strategy
* Implement strict temporal validation by training models on earlier years (March 2013 to early 2017) and testing on later years to prevent temporal data leakage.
* Execute leave-one-site-out validation across the 12 monitoring stations to evaluate spatial generalizability across different regions of Beijing.

## Data Engineering Challenges
* Missingness: The dataset contains known missing values. The modeling pipeline must address informative missing data resulting from gaps in the monitoring site sensor records.
* Temporal Drift: Manage class imbalance and handle seasonal environmental drift across the longitudinal hourly observations.

# Group Members and Roles (Name, Email)
1. Andreas Tjavaras (andreas.tjavaras@duke.edu):
2. Tim Deng (tim.deng@duke.edu):
3. Joe Stafford (joseph.stafford@duke.edu): 
4. JT Rapp (jt.rapp@duke.edu):

# Tentative Assignment of Roles:
Data Wrangling: JS

Descriptive Analysis:

Inferential Stats: 

Project write-up/presentation: Everyone 

# Tasks to be completed: 
Data Cleaning, handling of missing data, descriptive analysis, predictor/Variable Selection, Cross Validation Process, Model fitting, etc.
