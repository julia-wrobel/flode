filename = paste0("/Users/JWROBE8/onedrive/Data/flode/funreg_comparison")

# load results for flode and fconc
# For this run fhist did not have enough basis functions and was rerun later
files = list.files(filename, 
                   pattern = "20231212_",
                   full.names =  TRUE)

files_genFlode_fhist = list.files(filename, 
                                  pattern = "fhistOnly",
                                  full.names =  TRUE)

files_genFhist = list.files(paste0(filename, "/genFhist"), 
                            pattern = ".RDA",
                            full.names =  TRUE)



# get results when data is generated from flode model

filename = here::here("output",
                      "202312_simResults.RData")
if(file.exists(filename)){
  load(filename)
}else{
  df = map_dfr(files, get_ests)
  df2 = map_dfr(files_genFlode_fhist, get_ests_fhist_only)
  df_fhist = map_dfr(files_genFhist, get_ests)
  
  # remove fhist from the 1st set of flode generating data
  df = df %>%
    select(-surfaceErr_fhist, -sse_fhist, -time_fhist)
  
  df2 = df2 %>%
    select(iter, seed, alpha_true,time_fhist, surfaceErr_fhist, sse_fhist, genMod)
  
  df = left_join(df, df2)
  
  save(df, df_fhist, file = filename)
}


# get results for a couple of surfaces

filename = here::here("output",
                      "202312_surfaceResults.RData")
if(file.exists(filename)){
  load(filename)
}else{
  fl = files[grep("iter_19_", files)]
  fh = files_genFlode_fhist[grep("iter_19_", files_genFlode_fhist)]
  f_genfhist = files_genFhist[grep("iter_19_", files_genFhist)]
  
  surface_df = map2_dfr(fl, fh, get_surfaces)
  
  surfaceGenFhist_df = get_surfaces_genFhist(f_genfhist)
  save(surface_df, surfaceGenFhist_df, file = filename)
}


filename = here::here("output",
                      "202401_simResults_prediction.RData")
if(file.exists(filename)){
  load(filename)
}else{
  ## generate list of datasets
  pred_dats = lapply(c(0.1, 0.5, 1, 4, 8, 12), function(a){
    simulate_flode(I = 10000, D = 50, sigma_y0 = 5,
                   sigma = 0.1, alpha = a,
                   rand_int = TRUE, lambda0 = 20,
                   seed = 23423)$data %>% mutate(int = 1)
  })
  
  set.seed(23423)
  pred_df_genFhist = pffr_sim(scenario = c("int", "ff"),
                              n = 10000, nxgrid = 50,SNR = 10,
                              sigma_d = 1,
                              rand_int = TRUE)
  
  df_pred = map2_dfr(files, files_genFlode_fhist, 
                     get_predictions, new_data = pred_dats)
  
  df_pred_fhist = map_dfr(files_genFhist, get_predictions_genFhist, new_data = pred_df_genFhist)
  
  save(df_pred, df_pred_fhist, file = filename)
}


# Save all processed simulation data.
save(df, df_fhist, df_pred, df_pred_fhist, surface_df, surfaceGenFhist_df, 
     file = here::here("output", "sim_results.Rdata"))
