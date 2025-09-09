**What this code does**

This repo contains all the processing, analysis and plotting code for my INF6000 dissertation project, titled '**Evaluating the impact of the Sheffield Clean Air Zone on air quality and traffic volume'.** As per University of Sheffield assessment guidelines, the source data for this project has been uploaded separately in a Google Drive folder, which also contains a link from this repo.

In order of how they should be run (scripts may throw errors if run out of order as later scripts rely on outputs from the earlier ones):

-   0_master_script loads all reuqired packages and runs all the other scripts in order.

-   1_spatial_preprocessing loads in sensor spatial metadata, defines the boundary of the CAZ and spillover zones and plots some initial maps.

-   2_timeseries_preprocessing loads in pollutant/traffic flow data, and creates a master dataset for each, filtering for the correct data range, removing NA values and aggregating each sensor's time series to hourly

-   3_timeseries_normalisation loads in ERA5 meteorological data, and uses the rmweather package to normalise air quality/traffic timeseries for weather conditions and seasonality

-   4_EDA creates summary tables for the pre-/post- analysis and creates plots of the normalised time series for air quality/traffic as well as plotting rmweather partial dependency/variable importance/predicted vs actual plots

-   5_RDiT carries out the sharp and donut RDiT analyses using the rdrobust package and creates summary tables of pooled coefficients

-   6_sensitivity_testing carries out the bandwidth sensitivity and minimum detectable effects analysis on the RDiT models.

**How to run it**

Download all the scripts and the source data, open 0_master_script, press run. All the other scripts will then run in sequence automatically. Some scripts (particularly EDA and timeseries_normalisation) can take a while to run due to growing and fitting multiple rmweather random forest models.The scripts can also be run individually, as long as they are run in same order as in the master script.

Each script generates quite a few objects. To aid with navigation, key output tables/lists/objects are highlighted in the comments at the top of each script. In general, functions/objects/dataframes relating to the air quality analysis are prefixed or suffixed 'aq\_', and those relating to traffic 'tf\_/traffic\_', or have those keywords in the name somewhere. The exceptions to this rule are graphs included in the dissertation, which are assigned to ggplot objects named for their figure references (e.g. 'Figure_4' in the script corresponds to 'Figure 4' in the dissertation).
