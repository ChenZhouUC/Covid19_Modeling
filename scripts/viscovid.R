library(tidyverse)
library(scales)
library(gridExtra)
library(dict)

ConvertDf <- function(df) {
  # single day row
  df <- rename(df, "confirmed_add" = "confirmed")
  df <- rename(df, "asymp_add" = "asymptomatic")
  df$total_add <- df$confirmed_add + df$asymp_add

  # accumulated row
  df$confirmed_acc <- cumsum(df$confirmed_add)
  df$asymp_acc <- cumsum(df$asymp_add)
  df$total_affected <- df$confirmed_acc + df$asymp_acc

  # confirmed data
  confirmed <- df[, c("date", "confirmed_add", "confirmed_acc")]
  confirmed$type <- "confirmed"
  colnames(confirmed) <- c("date", "add", "acc", "type")

  # asymp data
  asymp <- df[, c("date", "asymp_add", "asymp_acc")]
  asymp$type <- "asymptomatic"
  colnames(asymp) <- c("date", "add", "acc", "type")

  # total affected data
  ori <- df[, c("date", "total_affected")]
  ori$add <- 0
  ori$type <- "total"
  colnames(ori) <- c("date", "acc", "add", "type")

  # comb all
  df_acc <- rbind(confirmed, asymp, ori)

  return(df_acc)
}


DrawLinePlot <- function(df,
                         time_range,
                         color_shape_map,
                         axis_map,
                         xlabel_size = 10, legend_open = TRUE) {
  df$date <- as.Date(df$date, format = "%Y-%m-%d")
  # date chose
  begin_date_set <- as.Date(time_range[1], format = "%Y-%m-%d")
  term_date_set <- as.Date(time_range[2], format = "%Y-%m-%d")
  time_seq_set <- seq(begin_date_set, term_date_set, by = "days")

  # date default
  begin_date <- min(df$date)
  term_date <- max(df$date)
  time_seq <- seq(begin_date, term_date, by = "days")

  # compare the date range
  if (all(time_seq_set %in% time_seq)) {
    print("Dates chose were all in the range.")
    begin_date <- begin_date_set
    term_date <- term_date_set
  } else {
    return("Dates chose were NOT in the range!")
  }

  # convert original df
  df_acc <- ConvertDf(df[which((df$date >= begin_date) & (df$date <= term_date)), ])

  # check if legend needed
  if (legend_open == TRUE) {
    legend_position_top <- c(0.01, 0.99)
    legend_position_bottom <- c(0.01, 0.01)
  } else {
    legend_position_top <- "none"
    legend_position_bottom <- "none"
  }

  # assign params
  label_map <- color_shape_map[["label_map"]]
  color_map <- color_shape_map[["color_map"]]
  shape_size_map <- color_shape_map[["shape_size_map"]]

  line_yaxis_unit <- axis_map[["line_yaxis_unit"]]
  line_yaxis_unit_scale <- axis_map[["line_yaxis_unit_scale"]]
  line_ylimit <- axis_map[["line_ylimit"]]
  line_y_breaks <- axis_map[["line_y_breaks"]]
  line_size <- axis_map[["line_size"]]

  bar_yaxis_unit <- axis_map[["bar_yaxis_unitbar_yaxis_unit"]]
  bar_yaxis_unit_scale <- axis_map[["bar_yaxis_unit_scale"]]
  bar_ylimit <- axis_map[["bar_ylimit"]]
  bar_y_breaks <- axis_map[["bar_y_breaks"]]
  bar_width <- axis_map[["bar_width"]]

  title <- axis_map[["title"]]
  ylab_top <- axis_map[["ylab_top"]]
  ylab_bottom <- axis_map[["ylab_bottom"]]

  # plot accumulation
  total_trend <- ggplot() +
    geom_line(data = df_acc, aes(x = date, y = acc, color = type), size = line_size, alpha = .5) +
    geom_point(data = df_acc, aes(x = date, y = acc, color = type, shape = type), size = line_size) +
    scale_color_manual(name = "types", label = label_map, values = color_map) +
    scale_shape_manual(name = "types", label = label_map, values = shape_size_map) +
    scale_x_date(breaks = seq(begin_date, term_date, by = "1 days"), limits = c(begin_date - 1, term_date + 1)) +
    scale_y_continuous(
      breaks = seq(line_ylimit[1], line_ylimit[2], line_y_breaks),
      labels = label_number(suffix = paste("", line_yaxis_unit), scale = line_yaxis_unit_scale),
      limits = line_ylimit
    ) +
    theme(
      legend.position = legend_position_top,
      legend.justification = c("left", "top"),
      legend.direction = "vertical",
      legend.box = "horizontal",
      axis.text.x = element_text(size = xlabel_size, angle = 90, hjust = -2, vjust = 0),
      axis.title.x = element_blank()
    ) +
    ylab(ylab_top) +
    ggtitle(title)

  # plot new add
  sub_df <- filter(df_acc, type != "total")
  new_label_map <- label_map[!label_map == "total"]
  new_color_map <- color_map[names(color_map) != "total"]

  total_add <- ggplot(data = sub_df, aes(x = date, y = add)) +
    geom_bar(aes(fill = type), width = bar_width, position = "stack", stat = "identity") +
    scale_fill_manual(name = "types", label = new_label_map, values = new_color_map) +
    scale_x_date(breaks = seq(begin_date, term_date, by = "1 days"), position = "top", limits = c(begin_date - 1, term_date + 1)) +
    scale_y_reverse(
      breaks = seq(bar_ylimit[1], bar_ylimit[2], bar_y_breaks),
      labels = label_number(suffix = paste("", bar_yaxis_unit), scale = bar_yaxis_unit_scale),
      limits = rev(bar_ylimit)
    ) +
    ylab(ylab_bottom) +
    theme(
      legend.position = legend_position_bottom,
      legend.justification = c("left", "bottom"),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )

  return(grid.arrange(total_trend, total_add, ncol = 1, heights = c(3, 1)))
}
