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
    print(sprintf("Dates chosen are all in the range given: [%s - %s] . Start plotting...", begin_date, term_date))
    begin_date <- begin_date_set
    term_date <- term_date_set
  } else {
    print(sprintf("Dates chosen are NOT in the range given: [%s - %s] . Please re-specify!", begin_date, term_date))
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
  shape_map <- color_shape_map[["shape_map"]]

  line_yaxis_unit <- axis_map[["line_yaxis_unit"]]
  line_yaxis_unit_scale <- ifelse(line_yaxis_unit == "K", 0.001,
    ifelse(line_yaxis_unit == "W", 0.0001, 1)
  )
  line_ylimit <- axis_map[["line_ylimit"]]
  line_y_breaks <- axis_map[["line_y_breaks"]]
  line_size <- axis_map[["line_size"]]
  line_alpha <- axis_map[["line_alpha"]]

  bar_yaxis_unit <- axis_map[["bar_yaxis_unit"]]
  bar_yaxis_unit_scale <- ifelse(bar_yaxis_unit == "K", 0.001,
    ifelse(bar_yaxis_unit == "W", 0.0001, 1)
  )
  bar_ylimit <- axis_map[["bar_ylimit"]]
  bar_y_breaks <- axis_map[["bar_y_breaks"]]
  bar_width <- axis_map[["bar_width"]]
  bar_alpha <- axis_map[["bar_alpha"]]

  title <- axis_map[["title"]]
  ylab_line <- axis_map[["ylab_line"]]
  ylab_bar <- axis_map[["ylab_bar"]]

  # plot accumulation
  total_trend <- ggplot() +
    geom_line(data = df_acc, aes(x = date, y = acc, color = type), size = line_size, alpha = line_alpha) +
    geom_point(data = df_acc, aes(x = date, y = acc, color = type, shape = type), size = line_size) +
    scale_color_manual(name = "types", label = label_map, values = color_map) +
    scale_shape_manual(name = "types", label = label_map, values = shape_map) +
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
    ylab(ylab_line) +
    ggtitle(title)

  # plot new add
  sub_df <- filter(df_acc, type != "total")

  total_add <- ggplot(data = sub_df, aes(x = date, y = add)) +
    geom_bar(aes(fill = type), width = bar_width, position = "stack", stat = "identity", alpha = bar_alpha) +
    scale_fill_manual(name = "types", label = label_map, values = color_map) +
    scale_x_date(breaks = seq(begin_date, term_date, by = "1 days"), position = "top", limits = c(begin_date - 1, term_date + 1)) +
    scale_y_reverse(
      breaks = seq(bar_ylimit[1], bar_ylimit[2], bar_y_breaks),
      labels = label_number(suffix = paste("", bar_yaxis_unit), scale = bar_yaxis_unit_scale),
      limits = rev(bar_ylimit)
    ) +
    ylab(ylab_bar) +
    theme(
      legend.position = legend_position_bottom,
      legend.justification = c("left", "bottom"),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )

  return(grid.arrange(total_trend, total_add, ncol = 1, heights = c(3, 1)))
}

ConvertIncreDf <- function(df) {
  # col type check
  df$phase <- as.factor(df$phase)
  df$date <- as.Date(df$date, format = "%Y-%m-%d")
  df <- arrange(df, date)

  # daily real incre
  df$incr_add <- c(df$total_affected[1], diff(df$total_affected))

  # daily predicted incre (phase 0)
  df$incr_pred_min <- c(df$pred_min[1], diff(df$pred_min))
  df$incr_pred_max <- c(df$pred_max[1], diff(df$pred_max))
  df[df$phase != "0", "incr_pred_min"] <- 0
  df[df$phase != "0", "incr_pred_max"] <- 0
  return(df)
}

SegmentDf <- function(df) {
  phase_vec <- unique(df$phase)

  tail_miss_data <- data.frame()
  for (phase in phase_vec[phase_vec != "0"]) {
    tail_row <- tail(df[df$phase == phase, ], 1)
    tail_miss_data <- rbind(tail_miss_data, tail_row)
  }

  head_miss_data <- data.frame()
  for (phase in phase_vec[phase_vec != "1"]) {
    head_row <- head(df[df$phase == phase, ], 1)
    head_miss_data <- rbind(head_miss_data, head_row)
  }

  miss_df <- cbind(select(tail_miss_data, date, total_affected), select(head_miss_data, date, total_affected))
  names(miss_df) <- c("last_date", "last_total_affected", "next_date", "next_total_affected")
  return(miss_df)
}

DrawPredictPlot <- function(df,
                            time_range,
                            color_shape_map,
                            axis_map,
                            xlabel_size = 10,
                            legend_open = TRUE) {
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
    print(sprintf("Dates chosen are all in the range given: [%s - %s] . Start plotting...", begin_date, term_date))
    begin_date <- begin_date_set
    term_date <- term_date_set
  } else {
    print(sprintf("Dates chosen are NOT in the range given: [%s - %s] . Please re-specify!", begin_date, term_date))
  }

  # convert original df & generate segment df
  df_acc <- ConvertIncreDf(df[which((df$date >= begin_date) & (df$date <= term_date)), ])
  seg_df <- SegmentDf(df_acc)

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
  shape_map <- color_shape_map[["shape_map"]]

  line_yaxis_unit <- axis_map[["line_yaxis_unit"]]
  line_yaxis_unit_scale <- ifelse(
    line_yaxis_unit == "K", 0.001,
    ifelse(line_yaxis_unit == "W", 0.0001, 1)
  )
  line_ylimit <- axis_map[["line_ylimit"]]
  line_ybreaks <- axis_map[["line_ybreaks"]]
  line_size <- axis_map[["line_size"]]
  line_alpha <- axis_map[["line_alpha"]]
  point_size <- axis_map[["point_size"]]

  area_alpha <- axis_map[["area_alpha"]]

  bar_yaxis_unit <- axis_map[["bar_yaxis_unit"]]
  bar_yaxis_unit_scale <- ifelse(
    bar_yaxis_unit == "K", 0.001,
    ifelse(bar_yaxis_unit == "W", 0.0001, 1)
  )
  bar_ylimit <- axis_map[["bar_ylimit"]]
  bar_y_breaks <- axis_map[["bar_y_breaks"]]
  bar_width <- axis_map[["bar_width"]]
  bar_alpha <- axis_map[["bar_alpha"]]

  error_bar_width <- axis_map[["error_bar_width"]]
  error_bar_alpha <- axis_map[["error_bar_alpha"]]
  error_bar_size <- axis_map[["error_bar_size"]]

  title <- axis_map[["title"]]
  ylab_line <- axis_map[["ylab_line"]]
  ylab_bar <- axis_map[["ylab_bar"]]

  # plot accumulation
  total_trend <- ggplot(data = df_acc) +
    geom_line(aes(x = date, y = total_affected, group = phase, color = phase), size = line_size, alpha = line_alpha) +
    geom_point(aes(x = date, y = pred_affected, group = phase, color = phase, shape = phase), size = point_size) +

    # segment line
    geom_segment(
      data = seg_df, aes(x = last_date, y = last_total_affected, xend = next_date, yend = next_total_affected),
      colour = color_shape_map[["color_map"]][1:length(unique(df_acc$phase)) - 1],
      size = line_size,
      alpha = line_alpha
    ) +

    # confidence interval
    geom_ribbon(aes(x = date, ymin = pred_min, ymax = pred_max), alpha = area_alpha) +
    scale_color_manual(name = "phase", label = label_map, values = color_map) +
    scale_shape_manual(name = "phase", label = label_map, values = shape_map) +
    scale_x_date(
      breaks = seq(begin_date, term_date, by = "1 days"),
      limits = c(begin_date - 1, term_date + 1)
    ) +
    scale_y_continuous(
      breaks = seq(line_ylimit[1], line_ylimit[2], line_ybreaks),
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
    ylab(ylab_line) +
    ggtitle(title)

  # plot daily increment
  sub_df <- df_acc[df_acc$phase == "0", ]

  total_add <- ggplot(data = df_acc, aes(x = date, y = incr_add)) +
    geom_bar(aes(fill = phase), width = bar_width, position = "stack", stat = "identity", alpha = bar_alpha) +
    scale_fill_manual(name = "phase", label = label_map, values = color_map) +

    # error bar
    geom_errorbar(
      data = sub_df,
      aes(x = date, ymin = incr_pred_min, ymax = incr_pred_max),
      width = error_bar_width,
      alpha = error_bar_alpha,
      size = error_bar_size
    ) +
    scale_x_date(
      breaks = seq(begin_date, term_date, by = "1 days"),
      position = "top",
      limits = c(begin_date - 1, term_date + 1)
    ) +
    scale_y_reverse(
      breaks = seq(bar_ylimit[1], bar_ylimit[2], bar_y_breaks),
      labels = label_number(suffix = paste("", bar_yaxis_unit), scale = bar_yaxis_unit_scale),
      limits = rev(bar_ylimit)
    ) +
    ylab(ylab_bar) +
    theme(
      legend.position = legend_position_bottom,
      legend.justification = c("left", "bottom"),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )
  return(grid.arrange(total_trend, total_add, ncol = 1, heights = c(3, 1)))
}
