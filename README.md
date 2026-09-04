# Variance Fixing for Small Area Estimation with DHS Data

This repository contains R code for evaluating and applying a variance-fixing approach for small area estimation with DHS survey data. The main empirical application uses the 2018 Zambia DHS child nutrition indicator `CN_NUTS_C_WH2`, compares direct estimates, fixed direct estimates, Fay-Herriot smoothing, GVF-based variance estimates, and Joint-style spatial variance smoothing, and includes simulation studies based on the Zambia administrative and sampling structure.

## Repository Structure

``` text
├── codes/
│   ├── directEST_1030_national.R
│   ├── fhModel_1030.R
│   ├── dhs/
│   │   ├── compare_survey.r
│   │   ├── interval_plot.r
│   │   ├── investigate_cluster_individual_sample_size.R
│   │   ├── map_script.r
│   │   ├── ranking_map.r
│   │   ├── summary_illegal.r
│   │   └── supplement.r
│   ├── k0_5_pop/
│   │   └── k0_5_pop.r
│   ├── sanity check/
│   │   └── large_sample.r
│   ├── simulation comparison/
│   │   └── scenario_comparison_fixed_illegal_unfixed_legal.R
│   ├── simulation_sensitivity/
│   │   ├── simulation_direct_GVF_Admin1MeanWeight.R
│   │   └── simulation_direct_GVF_benchmarks_phantomQ10.R
│   └── simultion/
│       └── simulation_direct_GVF_benchmarks_updated.R
├── data/
│   └── Zambia/
│   │   ├── country_shp_analysis.rds
│   │   ├── zmb_frame_ea.rds
│   │   ├── zmb_ppp_2010_1km_Aggregated_UNadj.tif
│   │   ├── zmb_ppp_2018_1km_Aggregated_UNadj.tif
│   │   └── zmb_sample_ea.rds
└── vairance fix.Rproj
```

## Code Overview

### Core Helper Scripts

-   `codes/directEST_1030_national.R`: Defines `directEST_1030()`, the central direct-estimation function. It computes admin-0, admin-1, and admin-2 direct estimates and implements the variance-fixing function when needed. It is installed in surveyPrev package.
-   `codes/fhModel_1030.R`: Defines `fhModel_1030()`, a Fay-Herriot smoothing function. It is installed in surveyPrev package.

### DHS Analysis Scripts

### `codes/dhs/`

-   `map_script.r`: Runs the main Zambia DHS direct and FH analyses, classifies unsampled and sparse/variance-fixed admin-2 areas, and creates map outputs.
-   `interval_plot.r`: Runs Zambia admin-2 direct, fixed direct, nested FH, GVF, GVF-FH, and Joint-BYM2 interval comparisons.
-   `ranking_map.r`: Builds ranking plots and maps for direct and FH outputs.
-   `summary_illegal.r`: Produces the table of admin-2 areas where the variance fix is triggered, including phantom means and weights used in the fixed estimator.
-   `supplement.r`: Produces supplementary analyses, including the Lavushimanda example, nested versus non-nested FH comparisons, admin-1 aggregation, and WorldPop under-5 population weighted aggregation.
-   `compare_survey.r`: Compares `directEST_1030()` results against survey-package calculations.
-   `investigate_cluster_individual_sample_size.R`: Summarizes observed Zambia DHS cluster-level individual sample sizes and weight distributions. This is exploratory/supporting analysis and can be run independently once DHS and boundary inputs are available.

### Population Preparation

-   `codes/k0_5_pop/k0_5_pop.r`: Downloads and generates the under-5 population raster for Zambia 2018 using WorldPop data. This is needed before running the WorldPop population-weighted aggregation part of `codes/dhs/supplement.r`.

### Simulation Scripts

-   `codes/simultion/simulation_direct_GVF_benchmarks_updated.R`: Main Zambia-based direct-level GVF default (national mean with average national sum of weight) simulation.
-   `codes/simulation_sensitivity/simulation_direct_GVF_benchmarks_phantomQ10.R`: Sensitivity simulation where phantom-cluster weights use the national 10% cluster-weight quantile by stratum.
-   `codes/simulation_sensitivity/simulation_direct_GVF_Admin1MeanWeight.R`: Sensitivity simulation using admin-1 mean/weight setting.
-   `codes/simulation comparison/scenario_comparison_fixed_illegal_unfixed_legal.R`: Compares the default, national-10%-weight, and admin-1-weight simulation scenarios. It should be run only after the corresponding simulation summary `.rds` outputs have been regenerated.

### Sanity Checks

-   `codes/sanity check/large_sample.r`: Independent validation/sanity-check simulation for a large-sample setting.

## Data

### `data/Zambia/`

-   **`country_shp_analysis.rds`**\
    Zambia administrative boundary file used for analysis. This object contains a list of administrative levels, where each level includes detailed information on the corresponding geographic areas.

    To generate this file, follow the workflow described in the [Stratification-Pipeline repository](https://github.com/jitongj/Stratification-Pipeline). Specifically:

    1.  Run `DataProcessing_helper.R`\
    2.  Run `step0_create_info.R`\
    3.  Run `step1_prepare_dat.R`

    Make sure to configure all required parameters according to the instructions provided in the repository README before executing the scripts.

-   **`zmb_frame_ea.rds`**\
    Enumeration Area (EA) frame data from the DHS Zambia 2018 survey report.

    To generate this file, follow the procedures described in Section 6 of the [Stratification-Pipeline repository](https://github.com/jitongj/Stratification-Pipeline).

-   **`zmb_sample_ea.rds`**\
    Sampled EA data used for survey estimation, based on the DHS Zambia 2018 survey report.

    To generate this file, follow the procedures described in Section 6 of the [Stratification-Pipeline repository](https://github.com/jitongj/Stratification-Pipeline).

## Recommended Workflow

The analyses are not one single linear pipeline. Use the following groups depending on what needs to be reproduced.

### 1. Prepare Data

1.  Open the project with `vairance fix.Rproj` or run scripts from the project root so `here::here()` resolves paths correctly.
2.  Configure access to DHS data for Zambia 2018.
3.  For supplementary WorldPop age-0/age-1 aggregation, run `codes/k0_5_pop/k0_5_pop.r` to generate `data/subpop/zmb_k0_5_2018_1km.tif`.

### 2. Core Zambia DHS Analysis

Run these scripts independently after the required DHS and boundary inputs are available:

1.  `codes/dhs/map_script.r`
2.  `codes/dhs/interval_plot.r`
3.  `codes/dhs/summary_illegal.r`
4.  `codes/dhs/ranking_map.r`
5.  `codes/dhs/supplement.r`

### 3. Simulation Study

The simulation scripts can be run independently of the DHS plotting scripts once DHS access, boundaries, and the WorldPop 2010 population raster are available.

1.  Run `codes/simultion/simulation_direct_GVF_benchmarks_updated.R`.
2.  Run sensitivity simulations as needed:
    -   `codes/simulation_sensitivity/simulation_direct_GVF_benchmarks_phantomQ10.R`
    -   `codes/simulation_sensitivity/simulation_direct_GVF_Admin1MeanWeight.R`
3.  Run `codes/simulation comparison/scenario_comparison_fixed_illegal_unfixed_legal.R` after the default and sensitivity simulation summary `.rds` files have been regenerated.

### 4. Supporting Validation

Run `codes/dhs/compare_survey.r` to check `directEST_1030()` against survey-package calculations.
