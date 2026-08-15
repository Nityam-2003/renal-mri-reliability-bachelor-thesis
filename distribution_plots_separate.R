library(tidyverse)

# ------------------------------------------------------------
# Distribution plots for T1/T2 kidney MRI biomarkers
# ------------------------------------------------------------
# Descriptive plots only.
# Uses all whole-kidney ROI measurements from the four finalized input tables.
# Healthy volunteers and CKD patients are plotted separately.

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)

analysis_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()

repeatability_healthy_file <- file.path(
  analysis_dir,
  "Repeatability_HealthyVolunteers_input.csv"
)
reproducibility_healthy_file <- file.path(
  analysis_dir,
  "Reproducibility_HealthyVolunteers_input.csv"
)
repeatability_patients_file <- file.path(
  analysis_dir,
  "Repeatability_Patients_input.csv"
)
reproducibility_patients_file <- file.path(
  analysis_dir,
  "Reproducibility_Patients_input.csv"
)
output_dir <- file.path(analysis_dir, "DistributionPlots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_input <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      ID = as.character(ID),
      site = as.character(site),
      vendor = as.character(vendor),
      visit = as.character(visit),
      anatomy = as.character(anatomy),
      parameter = as.character(parameter),
      value = as.numeric(value)
    ) %>%
    filter(anatomy == "whole")
}

prepare_group_data <- function(repeatability_file, reproducibility_file) {
  bind_rows(
    read_input(repeatability_file),
    read_input(reproducibility_file)
  ) %>%
    distinct(
      ID,
      site,
      vendor,
      visit,
      anatomy,
      parameter,
      value,
      .keep_all = TRUE
    ) %>%
    filter(!is.na(value))
}

healthy_measurements <- prepare_group_data(
  repeatability_healthy_file,
  reproducibility_healthy_file
)

patient_measurements <- prepare_group_data(
  repeatability_patients_file,
  reproducibility_patients_file
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
