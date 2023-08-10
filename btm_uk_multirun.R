# remotes::install_github("ncsu-landscape-dynamics/rpops")
library(PoPS)
library(terra)
library(folderfun)
library(doParallel)

# Change this to the directory that you save the folder to and you will good to go
setff("In", "C:/Users/cmjone25/Desktop/workshop_data/boxtree_moth/UK_data/")

infected_file <- ffIn("btm_2008.tif")
host_file <- ffIn("buxus_host_corrected.tif")
total_populations_file <- ffIn("total_pop.tif")
temp <- FALSE
temperature_coefficient_file <- ""
precip <- FALSE
precipitation_coefficient_file <- ""
model_type <- "SI"
latency_period <- 0
time_step <- "month"
season_month_start <- 4
season_month_end <- 10
start_date <- "2009-01-01"
end_date <- "2021-12-31"
use_survival_rates <- FALSE
survival_rate_month <- 3
survival_rate_day <- 15
survival_rates_file <- ""
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -12.87
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
mortality_frequency <- "year"
mortality_frequency_n <- 1
management <- FALSE
treatment_dates <- c("")
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "network"
natural_dir <- "NONE"
natural_kappa <- 0
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
mask <- NULL
output_frequency <- "year"
output_frequency_n <- 1
movements_file <- ""
use_movements <- FALSE
start_exposed <- FALSE
generate_stochasticity <- TRUE
establishment_stochasticity <- TRUE
movement_stochasticity <- TRUE
dispersal_stochasticity <- FALSE
establishment_probability <- 0.5
dispersal_percentage <- 0.99
quarantine_areas_file <- ""
use_quarantine <- FALSE
use_spreadrates <- FALSE
use_overpopulation_movements <- FALSE
overpopulation_percentage <- 0
leaving_percentage <- 0
leaving_scale_coefficient <- 1
number_of_iterations <- 100000
exposed_file <- ""
verbose = TRUE
write_outputs <- "None"
output_folder_path <- ""
network_filename <- "H:/Shared drives/Data/Table/Regional/UK_road_segments.csv"
network_movement <- "walk"
use_distance <- FALSE
use_rmse <- FALSE
use_mcc <- FALSE

parameter_means <- read.csv(ffIn("posterior_means.csv"))
parameter_means <- as.numeric(parameter_means[2:9])
parameter_cov_matrix <- read.csv(ffIn("posterior_cov_matrix.csv"))


btm_run_uk <- pops_multirun(infected_file,
                            host_file,
                            total_populations_file,
                            parameter_means,
                            parameter_cov_matrix,
                            temp,
                            temperature_coefficient_file,
                            precip,
                            precipitation_coefficient_file,
                            model_type,
                            latency_period,
                            time_step,
                            season_month_start,
                            season_month_end,
                            start_date,
                            end_date,
                            use_survival_rates = FALSE,
                            survival_rate_month = 3,
                            survival_rate_day = 15,
                            survival_rates_file = "",
                            use_lethal_temperature = FALSE,
                            temperature_file = "",
                            lethal_temperature = -12.87,
                            lethal_temperature_month = 1,
                            mortality_on = FALSE,
                            mortality_rate = 0,
                            mortality_time_lag = 0,
                            mortality_frequency = "year",
                            mortality_frequency_n = 1,
                            management = FALSE,
                            treatment_dates = c(""),
                            treatments_file = "",
                            treatment_method = "ratio",
                            natural_kernel_type = "cauchy",
                            anthropogenic_kernel_type = "cauchy",
                            natural_dir = "NONE",
                            anthropogenic_dir = "NONE",
                            number_of_iterations = 10,
                            number_of_cores = 5,
                            pesticide_duration = 0,
                            pesticide_efficacy = 1.0,
                            random_seed = NULL,
                            output_frequency = "year",
                            output_frequency_n = 1,
                            movements_file = "",
                            use_movements = FALSE,
                            start_exposed = FALSE,
                            generate_stochasticity = TRUE,
                            establishment_stochasticity = TRUE,
                            movement_stochasticity = TRUE,
                            dispersal_stochasticity,
                            establishment_probability,
                            dispersal_percentage,
                            quarantine_areas_file,
                            use_quarantine,
                            use_spreadrates,
                            use_overpopulation_movements,
                            overpopulation_percentage,
                            leaving_percentage,
                            leaving_scale_coefficient,
                            exposed_file,
                            mask,
                            write_outputs,
                            output_folder_path,
                            network_filename,
                            network_movement)

hosts <- rast(host_file)
plot(hosts)
plot(btm_run_uk$probability[[10]])
plot(btm_run_uk$simulation_mean[[10]])
# plot(btm_run_uk$simulation_max[[10]])
# plot(btm_points2, add = TRUE)
plot(btm_21_points, add = TRUE)
plot(btm_run_uk$probability)
# writeRaster(btm$probability, "all_probs.tif")
# writeRaster(data_network$probability[[1]], "probs2021.tif")