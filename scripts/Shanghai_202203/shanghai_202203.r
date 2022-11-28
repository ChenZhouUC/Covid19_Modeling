setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../viscovid.R")

data_sh <- read.csv("viscovid_sh_test.csv")

# param setting
time_range <- c("2022-02-16", "2022-06-27")

color_shape_map <- dict()
color_shape_map[["label_map"]] <- c("1", "2", "3", "0")
color_shape_map[["color_map"]] <- setNames(c("coral", "gold3", "olivedrab4", "darkgrey"), color_shape_map[["label_map"]])
color_shape_map[["shape_map"]] <- setNames(c(1, 3, 2, 8), color_shape_map[["label_map"]])

axis_map <- dict()
axis_map[["line_yaxis_unit"]] <- "W"
axis_map[["line_ylimit"]] <- c(0, 700000)
axis_map[["line_ybreaks"]] <- 100000
axis_map[["line_size"]] <- 1.5
axis_map[["line_alpha"]] <- 0.5
axis_map[["point_size"]] <- 2.5

axis_map[["area_alpha"]] <- 0.1

axis_map[["bar_yaxis_unit"]] <- "K"
axis_map[["bar_ylimit"]] <- c(0, 30000)
axis_map[["bar_y_breaks"]] <- 5000
axis_map[["bar_width"]] <- 0.8
axis_map[["bar_alpha"]] <- 0.5

axis_map[["error_bar_width"]] <- 0.4
axis_map[["error_bar_alpha"]] <- 0.8
axis_map[["error_bar_size"]] <- 0.5

axis_map[["title"]] <- "Shanghai Accumulated Covid19 Affected"
axis_map[["ylab_line"]] <- "Accumulated Affected"
axis_map[["ylab_bar"]] <- "Accumulated Affected"

xlabel_size <- 10
legend_open <- TRUE

plot_sh <- DrawPredictPlot(
  df = data_sh,
  time_range,
  color_shape_map,
  axis_map,
  xlabel_size,
  legend_open
)
plot_sh