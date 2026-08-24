# IMPORTANT: Thesis-specific distribution plots for quantitative renal MRI.
#
# Run this script from the ``Analysis/BRAINSFit_Registration`` folder after
# ``Descriptive_roi_values.py`` has created the two BRAINSFit descriptive
# master tables. It generates separate T1, MOLLI, and T2 distribution plots
# for healthy volunteers and CKD patients within this folder.
#
# This script requires the R package tidyverse. Do not commit the input tables
# or generated plots containing study data to a public repository.

library(tidyverse)

# ------------------------------------------------------------
# Distribution plots for T1/T2 kidney MRI biomarkers
# ------------------------------------------------------------
# Descriptive plots only.
# Uses every available whole-kidney ROI measurement from the descriptive master
# tables. Pairing and repeatability/reproducibility membership are irrelevant.
# Healthy volunteers and CKD patients are plotted separately.

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)

analysis_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()

healthy_master_file <- file.path(
  analysis_dir,
  "Descriptive_HealthyVolunteers_master.csv"
)
patient_master_file <- file.path(
  analysis_dir,
  "Descriptive_CKDPatients_master.csv"
)
output_dir <- file.path(analysis_dir, "DistributionPlots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_descriptive_master <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    transmute(
      ID = as.character(ID),
      population = as.character(population),
      site = as.character(site),
      vendor = as.character(vendor),
      raw_visit = as.character(raw_visit),
      parameter = as.character(parameter),
      value = as.numeric(mean_roi_value_ms),
      valid_voxel_count = as.integer(valid_voxel_count),
      status = as.character(status)
    ) %>%
    filter(status == "included", !is.na(value))
}

healthy_measurements <- read_descriptive_master(
  healthy_master_file
)

patient_measurements <- read_descriptive_master(
  patient_master_file
)

parameters <- c("T1", "MOLLI", "T2")

center_colors <- c(
  "01" = "#1b9e77",
  "02" = "#d95f02",
  "03" = "#7570b3",
  "04" = "#e7298a"
)

create_distribution_plot <- function(dataset, parameter_name, group_label) {
  
  plot_data <- dataset %>%
    filter(parameter == parameter_name) %>%
    mutate(
      site = factor(
        sprintf("%02d", as.integer(site)),
        levels = c("01", "02", "03", "04")
      )
    )
  
  x_label <- switch(
    parameter_name,
    "T1" = "T1 (ms)",
    "MOLLI" = "MOLLI T1 (ms)",
    "T2" = "T2 (ms)",
    paste(parameter_name, "value")
  )
  
  ggplot(
    plot_data,
    aes(
      x = value,
      colour = site
    )
  ) +
    
    geom_density(
      linewidth = 0.9
    ) +
    
    scale_colour_manual(
      values = c(
        "01" = "#0057e7",
        "02" = "#d62d20",
        "03" = "#2ca02c",
        "04" = "#800080"
      ),
      drop = TRUE
    ) +
    
    labs(
      title = paste(
        parameter_name,
        "Distribution -",
        group_label
      ),
      x = x_label,
      y = "Density",
      colour = "Center"
    ) +
    
    theme_bw() +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = rel(1.7)
      ),
      
      axis.title = element_text(
        size = rel(1.4)
      ),
      
      axis.text = element_text(
        size = rel(1.3)
      ),
      
      legend.title = element_text(
        size = rel(1.5)
      ),
      
      legend.text = element_text(
        size = rel(1.3)
      )
    )
  
}

groups <- list(
  HealthyVolunteers = list(
    data = healthy_measurements,
    label = "Healthy Volunteers"
  ),
  CKDPatients = list(
    data = patient_measurements,
    label = "CKD Patients"
  )
)

for (group_name in names(groups)) {
  for (parameter_name in parameters) {
    plot <- create_distribution_plot(
      groups[[group_name]]$data,
      parameter_name,
      groups[[group_name]]$label
    )

    output_file <- file.path(
      output_dir,
      paste0(group_name, "_", parameter_name, "_Distribution.pdf")
    )

    ggsave(
      filename = output_file,
      plot = plot,
      device = "pdf",
      width = 7.5,
      height = 5.5,
      units = "in"
    )
  }
}

cat("Distribution plots saved in:", output_dir, "\n")
