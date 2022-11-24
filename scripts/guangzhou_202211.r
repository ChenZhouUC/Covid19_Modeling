data_gz <- read.csv("covid19_guangzhou_221123.csv")
data_gz[is.na(data_gz)] <- 0
data_gz$confirmed_add <- data_gz$confirmed
data_gz$asymp_add <- data_gz$asymptomatic
data_gz$total_add <- data_gz$confirmed_add + data_gz$asymp_add

data_gz$confirmed_acc <- cumsum(data_gz$confirmed_add)
data_gz$asymp_acc <- cumsum(data_gz$asymp_add)
data_gz$total_affected <- data_gz$confirmed_acc + data_gz$asymp_acc

data_gz$date <- as.Date(data_gz$date, format = "%Y-%m-%d")
data_gz <- data_gz[which(data_gz$date >= "2022-10-15"), ]

confirmed <- data_gz[, c("date", "confirmed_add", "confirmed_acc")]
confirmed$type <- "confirmed"
colnames(confirmed) <- c("date", "add", "acc", "type")
asymp <- data_gz[, c("date", "asymp_add", "asymp_acc")]
asymp$type <- "asymptomatic"
colnames(asymp) <- c("date", "add", "acc", "type")
draw <- rbind(confirmed, asymp)

begin_date <- as.Date("2022-10-15", format = "%Y-%m-%d")
term_date <- as.Date("2022-11-22", format = "%Y-%m-%d")
total_trend <- ggplot() +
  geom_line(data = data_gz, aes(x = date, y = total_affected), size = 1.5, alpha = .2) +
  geom_point(data = data_gz, aes(x = date, y = total_affected), shape = 19) +
  geom_line(data = draw, aes(x = date, y = acc, color = type), size = 1, alpha = .5) +
  geom_point(data = draw, aes(x = date, y = acc, color = type), shape = 5, size = 2.0) +
  scale_x_date(breaks = seq(begin_date, term_date, by = "1 days"), limits = c(begin_date - 1, term_date + 1)) +
  scale_y_continuous(breaks = seq(0, max(data_gz$total_affected), 5e4)) +
  ylab("Accumulated Affected") +
  theme(
    legend.position = c(0.01, 0.99),
    legend.justification = c("left", "top"),
    legend.direction = "vertical",
    legend.box = "horizontal",
    axis.text.x = element_text(size = 10, angle = 90, hjust = -2, vjust = 0),
    axis.title.x = element_blank()
  ) +
  ggtitle("Guangzhou Accumulated Covid19 Affected")

total_add <- ggplot(data = draw, aes(x = date, y = add, color = type)) +
  geom_bar(aes(fill = type), size = 0.1, position = "stack", stat = "identity") +
  scale_x_date(
    breaks = seq(begin_date, term_date, by = "1 days"), position = "top",
    limits = c(begin_date - 1, term_date + 1)
  ) +
  scale_y_reverse(breaks = seq(0, 3e4, 5e3), labels = label_number(suffix = " K", scale = 1e-3), limits = c(3e4, 0)) +
  ylab("Daily Added Affected") +
  theme(
    axis.text.x = element_blank(),
    legend.position = c(0.01, 0.01),
    legend.justification = c("left", "bottom"),
    axis.title.x = element_blank()
  )

grid.arrange(total_trend, total_add, ncol = 1, heights = c(3, 1))
