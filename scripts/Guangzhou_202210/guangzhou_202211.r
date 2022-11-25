source("../viscovid.R")
data_gz <- read.csv("../../data/covid19_guangzhou_221123.csv")
data_gz[is.na(data_gz)] <- 0
data_gz$confirmed
data_gz$asymptomatic

time_range <- c("2022-10-01", "2022-11-22")

color_shape_map <- dict()
color_shape_map[["label_map"]] <- c("confirmed", "asymptomatic", "total")
color_shape_map[["color_map"]] <- setNames(c("coral", "gold3", "olivedrab4"), color_shape_map[["label_map"]])
color_shape_map[["shape_size_map"]] <- setNames(c(2, 3, 1), color_shape_map[["label_map"]])

axis_map <- dict()
axis_map[["line_yaxis_unit"]] <- "W"
axis_map[["line_yaxis_unit_scale"]] <- ifelse(axis_map[["line_yaxis_unit"]] == "K", 0.001,
  ifelse(axis_map[["line_yaxis_unit"]] == "W", 0.0001, 1)
)
axis_map[["line_ylimit"]] <- c(0, 110000 + 1000)
axis_map[["line_y_breaks"]] <- 20000
axis_map[["line_size"]] <- 2.5

axis_map[["bar_yaxis_unit"]] <- "W"
axis_map[["bar_yaxis_unit_scale"]] <- ifelse(axis_map[["bar_yaxis_unit"]] == "K", 0.001,
  ifelse(axis_map[["bar_yaxis_unit"]] == "W", 0.0001, 1)
)
axis_map[["bar_ylimit"]] <- c(0, 20000)
axis_map[["bar_y_breaks"]] <- 5000
axis_map[["bar_width"]] <- 0.8

axis_map[["title"]] <- "Guangzhou Accumulated Covid19 Affected"
axis_map[["ylab_top"]] <- "Accumulated Affected"
axis_map[["ylab_bottom"]] <- "Accumulated Affected"

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
