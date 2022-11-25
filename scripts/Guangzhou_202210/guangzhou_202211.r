setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../viscovid.R")

data_gz <- read.csv("../../data/covid19_guangzhou_221123.csv")
data_gz[is.na(data_gz)] <- 0
data_gz$confirmed
data_gz$asymptomatic

# time range def
time_range <- c("2022-10-15", "2022-11-22")
# color shape def
color_shape_map <- dict()
color_shape_map[["label_map"]] <- c("confirmed", "asymptomatic", "total")
color_shape_map[["color_map"]] <- setNames(c("red", "gold3", "olivedrab4"), color_shape_map[["label_map"]])
color_shape_map[["shape_map"]] <- setNames(c(3, 2, 1), color_shape_map[["label_map"]])
# axis line config def
axis_map <- dict()
axis_map[["line_yaxis_unit"]] <- "W"
axis_map[["line_ylimit"]] <- c(0, 120000)
axis_map[["line_y_breaks"]] <- 10000
axis_map[["line_size"]] <- 2.0
axis_map[["line_alpha"]] <- 0.5
# axis bar config def
axis_map[["bar_yaxis_unit"]] <- "K"
axis_map[["bar_ylimit"]] <- c(0, 10000)
axis_map[["bar_y_breaks"]] <- 2000
axis_map[["bar_width"]] <- 0.7
axis_map[["bar_alpha"]] <- 0.5
# ylab and title define
axis_map[["title"]] <- "Guangzhou Accumulated Covid19 Affected"
axis_map[["ylab_line"]] <- "Accumulated Affected Total"
axis_map[["ylab_bar"]] <- "Daily Affected Added"

xlabel_size <- 10
legend_open <- TRUE

plot_gz <- DrawLinePlot(
  df = data_gz,
  time_range,
  color_shape_map,
  axis_map,
  xlabel_size,
  legend_open
)
plot_gz
