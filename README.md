Supplemental Code for AOAS1943
================

## Directory structure

The repo has the following directory structure.

    ## 
    ## project directory
    ##   |-- flode.Rproj
    ##   |  o-- object of type(s):file
    ##   |-- README.Rmd
    ##   |  o-- object of type(s):file
    ##   |-- README.md
    ##   |  o-- object of type(s):file
    ##   |-- analysis/
    ##   |  o-- object of type(s):dir
    ##   |-- data/
    ##   |  o-- object of type(s):dir
    ##   |-- output/
    ##   |  o-- object of type(s):dir
    ##   o-- preprocessing/
    ##      o-- object of type(s):dir

Each of the files and subdirectories in this repository are described
below.

- `README.Rmd/md`: Readme explaining the materials in this repository
- `analysis/`
  - `202311_dataAnalysis_steps.Rmd/.html`: Code/output for processing
    motivating data and generating Figures 1 and 5 in the manuscript.
  - `202311_dataAnalysis_prediction.Rmd/.html`: Code/output for
    obtaining confidence intervals for coefficient functions shown in
    Figure 5 and for performing 10-fold cross validation of flode,
    fconc, and fhist.
  - `202312_simulation_results.Rmd/.html`: Code/output for processing
    simulation results and generating Figures 3 and 4 in the manuscript.
  - `202310_simulation_results_supplement.Rmd/.html`: Code/output for
    processing simulation results in the supplement.
- `data/20231012_tidied_steps.RDA`: Tidied mouse gait data
- `output/`: Figures, data analysis, and simulation results
  - `*.png`: Figures from the manuscript and supplement
  - `pUp_results.RData`: Data analysis results from `flode` model.
  - `pUp_fhist.RData`: Data analysis results from `fhist` model.
  - `pUp_fconc.RData`: Data analysis results from `fonc` model.
- `preprocessing/`: Source files for executing core modeling code,
  generating simulated data, and running simulations
  - `estimate_flode.R`: Defines function for estimating `flode` model
    using EM algorithm
  - `make_pffr.R`: Defines function for running functional historical
    regression model and processing surfaces for visualization
  - `search_alpha.R`: Defines function to perform grid search over
    possible $\alpha_0$ values to choose best initial value.
  - `simulated_paw.R`: Defines function for simulating data.
  - `utils.R`: Defines helper functions called in the `estimate_flode.R`
    file.
  - `simulations/`: File for generating simulation results across all
    simulation scenarios.
    - `20240605_simulations_fromKalman.R`: Performs simulations from a
      Kalman filter model for Supplemental Figure 6.
    - `20240123_simulations_varyN.R`: Performs simulations for
      Supplemental Figure 4.
    - `20231212_simulations_fullComparison.R`: Generates model results
      for flode, fhist, and fconc when data is generated from a flode
      model
    - `20231212_simulations_fullComparison_genfhist.R`: Generates model
      results for flode, fhist, and fconc when data is generated from an
      fhist model
    - `20231115_simulations_confInt_boot.R`: Produces simulations for
      flode paper for getting coverage for bootstrap confidence
      intervals
    - `20231025_simulations_tuning.R`: Performs simulations for varying
      over values of flode tuning parameters for Supplemental Figure 3.
