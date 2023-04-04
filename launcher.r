if (!require("pacman")) install.packages("pacman")
pacman::p_load(RColorBrewer, mgcv, faraway, scam, lme4)

### Here I have the following structure for the code
#
# Within the 'project_folder' we have:
#
# a subfolder named 'code' that contains all the code
# a subfolder named 'data' that contains all the data
# a subfolder named 'output' where figures will go

# project_folder <- [Insert project folder name here]

setwd(folder)

reprep <- TRUE
Include_Negative_NA <- TRUE

source("code/prep_data.r")

source("code/create_models.r")

source("code/OR_figure.r")

mod_to_plot <- mod_gam_linear_bin

source("code/results_figure.r")
