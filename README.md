# Spring_Departure

Code for Tonelli et al., Nature Ecol. Evol.

This repository accompanies the manuscript: A temperate-tropical divide in climate-driven changes to avian spring migration departure timing.

Please  contact btonell [at] clemson [dot] edu with questions.

**Repository structure:**

  * `code/`
    * `A_data_processing/`
      * `0_phenometric_functions.R` Re-used functions for phenometric generation.
      * `1_ebird_filter_departure.txt` Awk code for filtering and selecting relevant data from eBird global dataset.
      * `2_migrant_species_sorting_palm.R` Search for migratory species departing non-breeding grounds.
      * `3_sp_cl_all.R` Fits Bayesian GAMs of probability of occurrence on eBird checklists across the year. (HPC required)
      * `4_phenometric_filtering.R` Filter out bad phenometrics that fail model checks.
      * `5_era5_processing.R` Convert raster climate data to cell-area level measures.
      * `6_checklist_counter.R` Figure out how many checklists are used in this study (90 million).
      * `7_migratory_dist_calc.R` Gets the migratory distance of each species (measured as distance between centroids).
      * `8_prepare_mdl_data.R` Processes all data needed for final Bayesian model.
    * `B_model_run/`
      * `9_run_dep_mdl_hm.R` Fit the stan model (cluster suggested).
      * `dep_clim.stan` Stan model.
    * `C_model_analysis/`
      * `10_main_results.R` Extracts estimate for model parameters of interest.
      * `11_PPO.R` Calculate prior-posterior overlap.
      * `12_variance_explained.R` Variance explained calculations for populations and species.
    * `D_depature_projections/`
      * `13_modeled_projections.R` Calculate departure trajectories given regional climate estimates.
      * `14_shift_by_lat.R` Get phenology shifts according to species latitudes.
    * `E_figures/`
      * `15_methods_fig.R` Create Fig. 1.
      * `16_pop_sensitivity_maps.R` Create Fig. 2, S1, S2
      * `17_sp_sens_climate.R` Create Fig. 3.
      * `18_projection_plots.R` Create Fig. 4.
      * `19_climate correlation.R` Get climatic connectivity metrics for Fig. S3
      * `20_var_exp_climate_connectivity.R` Fig. S3.
      * `21_temp_precip_change.R` Fig. S4.
      * `22_PPC.R` Fig. S5.
      * `23_cat_plots.R` Fig. S6.
  * `data/` (See Zenodo)
      * `input/`
        * `elevation/` - Location for GMTED 2010 Elevation model, .tif format. https://www.usgs.gov/coastal-changes-and-impacts/gmted2010
        * `environmental/` - Environmental data.
          * `climate_projections/` - Climate projections for 3 CMIP6 climate models, split by scenario and environmental covariate. Available via: https://cds.climate.copernicus.eu/
            * `CESM2/`
              * `precip/`
                * `SSP1_2p6/`
                * `SSP2_4p5/`
                * `SSP3_7p0/`
              * `temp/`
                * `SSP1_2p6/`
                * `SSP2_4p5/`
                * `SSP3_7p0/`
            * `CNRM-ESM2-1/`
              * `precip/`
                * `SSP1_2p6/`
                * `SSP2_4p5/`
                * `SSP3_7p0/`
              * `temp/`
                * `SSP1_2p6/`
                * `SSP2_4p5/`
                * `SSP3_7p0/`
            * `MIROC6/`
              * `precip/`
                * `SSP1_2p6/`
                * `SSP2_4p5/`
                * `SSP3_7p0/`
              * `temp/`
                * `SSP1_2p6/`
                * `SSP2_4p5/`
                * `SSP3_7p0/`
              * `temp_may_global.nc` - May temperature averages for climatic connectivity 2002-2024. From ERA5
          * `ERA5/` - Historical climate data, saved as monthly .nc files (e.g. temp_2024_12.nc)
          * `intermediate_precip/` - Storage for intermediate precipitation data created during processing
          * `intermediate_temp/` - Storage for intermediate temperature data created during processing
        * `mdl_data/` - Data going into the main Bayesian model
        * `mig_dist/` - Migratory distance information
        * `phenology/` - Location for zero-filled eBird files, migratory species lists, GAM fits, 
          * `phenometrics_cell_yr_files/` - GAM fits, by species, year
          * `report_cards/` - Migratory species lists, by year
          * `zf_files/` - Zero-filled eBird files
        * `spatial/` - Polygons of areas, polygons of breeding ranges
      * `output/` - Model outputs, downstream analyses outputs
        * `departure_model/` - Fit model object and diagnostics
        * `departure_projections/` - Departure and climate projections by population, model/scenario
        * `mig_dist/` - Migration distance analyses model fits
        * `spec_hm_pdfs/` - Used if user wants to visualize GAM outputs.
        * `variance_explained/`- Variance explained calculations
        
