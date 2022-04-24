library(ggplot2) 
library(fda)
library(gridExtra)
library(scales)
setwd("/home/chenzhou/Documents/DocX/WhaleDocs/Covid19/data")
options(digits = 4)
options(scipen = 0)

data_sh <- read.csv("data_0423_SH.csv")
# data_sh_whole <- data_sh[which(data_sh$map_name=='上海'), ]
# data_sh_whole$date <- as.Date(data_sh_whole$datetime, format="%Y-%m-%d")
# data_sh_whole$confirmed_accumulated <- data_sh_whole$confirmed_accumulated - data_sh_whole$confirmed_accumulated[1]
# data_sh_whole$asymptomatic_accumulated <- data_sh_whole$asymptomatic_accumulated - data_sh_whole$asymptomatic_accumulated[1]

linear_fix <- function(df, colname){
  for (i in 2:(nrow(df)-1)) {
    if (df[i, colname] == df[i-1, colname]) {
      df[i, colname] <- (df[i-1, colname] + df[i+1, colname])/2
    }
  }
  return(df)
}
# data_sh_whole <- linear_fix(data_sh_whole, "asymptomatic_accumulated")
# data_sh_whole$total_affected <- data_sh_whole$confirmed_accumulated + data_sh_whole$asymptomatic_accumulated
# data_sh_whole$total_diff <- c(0, diff(data_sh_whole$total_affected,1))

data_sh$confirmed_add <- data_sh$inbound_confirmed + data_sh$outbound_confirmed
data_sh$asymp_add <- data_sh$inbound_asymp + data_sh$outbound_asymp
data_sh$total_add <- data_sh$confirmed_add + data_sh$asymp_add

data_sh$confirmed_acc <- cumsum(data_sh$confirmed_add) - 380
data_sh$asymp_acc <- cumsum(data_sh$asymp_add) - 120
data_sh$total_affected <- data_sh$confirmed_acc + data_sh$asymp_acc

data_sh$date <- as.Date(data_sh$date, format="%Y-%m-%d")
data_sh <- data_sh[which(data_sh$date >= "2022-03-01"),]

confirmed <- data_sh[,c("date","confirmed_add","confirmed_acc")]
confirmed$type <- "confirmed"
colnames(confirmed)  <- c("date", "add", "acc", "type")
asymp <- data_sh[,c("date","asymp_add", "asymp_acc")]
asymp$type <- "asymptomatic"
colnames(asymp)  <- c("date", "add", "acc", "type")
draw <- rbind(confirmed,asymp)

begin_date <- as.Date("2022-03-01", format="%Y-%m-%d")
term_date <- as.Date("2022-04-23", format="%Y-%m-%d")
total_trend <- ggplot() +
  geom_line(data=data_sh, aes(x=date, y=total_affected), size=1.5, alpha=.2) + geom_point(data=data_sh, aes(x=date, y=total_affected), shape=19) +
  geom_line(data=draw, aes(x=date, y=acc, color=type), size=1, alpha=.5) + geom_point(data=draw, aes(x=date, y=acc, color=type), shape=5, size=2.0) +
  scale_x_date(breaks=seq(begin_date,term_date,by="1 days"),limits=c(begin_date-1,term_date+1)) +
  scale_y_continuous(breaks=seq(0,max(data_sh$total_affected), 5e4)) +
  ylab("Accumulated Affected") + theme(
    legend.position = c(0.01, 0.99),
    legend.justification = c("left", "top"),
    legend.direction = "vertical",
    legend.box = "horizontal",
    axis.text.x=element_text(size=10, angle=90, hjust=-2, vjust=0),
    axis.title.x=element_blank()
  ) + ggtitle("Shanghai Accumulated Covid19 Affected")

total_add <- ggplot(data=draw, aes(x=date, y=add, color=type)) +
  geom_bar(aes(fill=type),size=0.1,position="stack",stat="identity") +
  scale_x_date(breaks=seq(begin_date,term_date,by="1 days"),position="top",limits=c(begin_date-1,term_date+1)) +
  scale_y_reverse(breaks=seq(0, 3e4, 5e3),labels=label_number(suffix=" K",scale=1e-3),limits=c(3e4,0)) +
  ylab("Daily Added Affected") +
  theme(axis.text.x=element_blank(),
        legend.position = c(0.01, 0.01),
        legend.justification = c("left", "bottom"),
        axis.title.x=element_blank())

grid.arrange(total_trend, total_add, ncol=1, heights=c(3,1))

######################################################################################################

covid19_dde <- function(ser, ser_interp, interp, total_len, kappa, mu, tau){
  # initialization
  init_num <- ceiling(tau * interp) + 1
  new_ser <- numeric(0)
  for (n in 1:init_num) {
    new_val_index <- (n-1) / interp * ser_interp
    l_ratio <- new_val_index - floor(new_val_index)
    new_val <- ser[floor(new_val_index)+1] * (1-l_ratio) + ser[floor(new_val_index)+2] * l_ratio
    new_ser <- c(new_ser, new_val)
  }
  # initialization finished
  while(length(new_ser) <= (total_len-1) * interp){
    deriv_ser <- max(0, kappa * (new_ser[length(new_ser)] - mu * new_ser[round(length(new_ser)-tau*interp)]))
    new_ser <- c(new_ser, new_ser[length(new_ser)]+deriv_ser/interp)
  }
  return(new_ser)
}

# show <- covid19_dde(data_sh_whole$total_affected[28:length(data_sh_whole$total_affected)], 1, 100, 20, 1.3, 1.3, 1.1)
# qplot(seq(0,20,0.01), show)

######################################################################################################

opt_bs_coeff <- function(y_ser, bs_interp, kappa, mu, tau, lambda, interg_interp, lr, decay, iters, show_process, iter_opt){
  # log likelihood part
  bs_domain <- seq(0, length(y_ser)-1, 1)
  knots <- seq(0, length(y_ser)-1, length.out=round(bs_interp*(length(y_ser)-1)+1))
  bS <- bsplineS(x=bs_domain,breaks=knots,norder=4,nderiv=0,returnMatrix=FALSE)
  n_coeff <- ncol(bS)
  init_coeff <- rep(0,n_coeff)
  init_grad <- rep(0,n_coeff)
  
  # differential integration part
  integ_x <- seq(0, length(y_ser)-1, length.out=interg_interp*(length(y_ser)-1)+1)
  bs_integ <- bsplineS(x=integ_x,breaks=knots,norder=4,nderiv=0,returnMatrix=FALSE)
  bs_partial <- bsplineS(x=integ_x,breaks=knots,norder=4,nderiv=1,returnMatrix=FALSE)
  shift_ <- round(tau * interg_interp)
  now_bs_integ <- bs_integ[(1+shift_):nrow(bs_integ),]
  delay_bs_integ <- bs_integ[1:(nrow(bs_integ)-shift_),]
  now_bs_partial <- bs_partial[(1+shift_):nrow(bs_partial),]
  delay_bs_partial <- bs_partial[1:(nrow(bs_integ)-shift_),]
  # loss_ser <- numeric(0)
  for (i_ in 1:iters){
    # loss calculation part
    loglik_loss_vec <- bS %*% init_coeff - y_ser
    loglik_loss <- sqrt(t(loglik_loss_vec) %*% loglik_loss_vec)
    mat_ <- now_bs_partial - kappa * now_bs_integ + kappa * mu * delay_bs_integ
    integ_loss_vec <- mat_ %*% init_coeff
    integ_loss <- sqrt(t(integ_loss_vec) %*% integ_loss_vec / interg_interp * lambda)
    loss_ <- loglik_loss + integ_loss
    # loss_ser <- c(loss_ser, loss_)
    # gradient calculation part
    loglik_part <- t(bS %*% init_coeff - y_ser) %*% bS
    integ_part <- t(mat_ %*% init_coeff) %*% mat_ / interg_interp * lambda
    init_coeff <- init_coeff - lr * t(loglik_part + integ_part)
    lr <- lr * decay
    # print(c(loss_, loglik_loss, integ_loss))
  }
  coef_mat <- t(bS) %*% bS + t(mat_) %*% mat_ / interg_interp * lambda
  k_mat <- t(now_bs_integ - mu * delay_bs_integ) %*% (now_bs_integ - mu * delay_bs_integ - now_bs_partial) %*% init_coeff
  m_mat <- ((mu * kappa**2) * (t(delay_bs_integ) %*% delay_bs_integ) + kappa * t(delay_bs_integ) %*% (now_bs_partial - kappa * now_bs_integ)) %*% init_coeff
  t_mat <- -(2*kappa*mu * (t(now_bs_partial-kappa*now_bs_integ+kappa*mu*delay_bs_integ) %*% delay_bs_partial) + 
               (now_bs_partial-kappa*now_bs_integ+kappa*mu*delay_bs_integ)[1,]**2) %*% init_coeff
  partial_resid_theta <- bS %*% solve(coef_mat) %*% cbind(k_mat,m_mat,t_mat)*2*lambda/interg_interp
  theta_grad <- solve(t(partial_resid_theta) %*% partial_resid_theta) %*% t(partial_resid_theta) %*% (y_ser - bS %*% init_coeff)
  if (show_process) {
    print(c(loss_, loglik_loss, integ_loss))
    print(bS[1,] %*% init_coeff)
    plot(0:(length(y_ser)-1), y_ser, cex=0.5,)
    lines(integ_x, bs_integ %*% init_coeff)
    Sys.sleep(1)
  }
  if (iter_opt) {
    return(c(theta_grad, loss_))
  } else {
    return(c(bS %*% init_coeff, (loglik_loss + integ_loss / sqrt(lambda)) / sqrt(length(y_ser))))
  }
}

opt_dde_theta <- function(y_ser, bs_interp, kappa_init, mu_init, tau_init, lambda, interg_interp, lr, decay, iters, lr_theta, iters_theta, show_process){
  losses <- c(Inf)
  theta <- c(kappa_init, mu_init, tau_init)
  for (i in 1:iters_theta) {
    theta_grad <- opt_bs_coeff(y_ser, bs_interp, kappa_init, mu_init, tau_init, lambda, interg_interp, lr, decay, iters, show_process, TRUE)
    if (theta_grad[4] == Inf | theta_grad[4] / min(losses) > 1.1) {
      break
    } else {
      if (min(losses) >= theta_grad[4]) {
        theta <- c(kappa_init, mu_init, tau_init)
        if (show_process) {
          print(theta)
        }
      }
      kappa_init <- kappa_init - lr_theta * theta_grad[1]
      mu_init <- mu_init - lr_theta * theta_grad[2]
      tau_init <- tau_init - lr_theta * theta_grad[3]
      lr_theta <- lr_theta * decay
      losses <- c(losses, theta_grad[4])
    }
  }
  return(theta)
}

opt_bs_coeff(data_sh$total_affected[1:28], 0.5, 1.60, 1.17, 1.42, 25., 500, 0.1, 0.99, 2e2, TRUE, FALSE) # 1.60, 1.17, 1.42
theta <- opt_dde_theta(data_sh$total_affected[1:28], 0.5, 1.60, 1.17, 1.42, 25., 500, 0.1, 0.99, 2e2, 0.05, 3e2, FALSE)
# show <- covid19_dde(data_sh$total_affected[1:28], 1, 10, 28, 1.6, 1.2, 1.4)
# qplot(seq(1,28,0.1), show)

opt_bs_coeff(data_sh$total_affected[28:nrow(data_sh)], 0.5, 0.72, 1.05, 2.18, 25., 500, 0.1, 0.99, 2e2, TRUE, FALSE) # 0.72, 1.05, 2.18
theta <- opt_dde_theta(data_sh$total_affected[28:nrow(data_sh)], 0.5, 0.72, 1.05, 2.18, 25., 500, 0.1, 0.99, 2e2, 0.05, 1e2, TRUE)
# show <- covid19_dde(data_sh$total_affected[29:nrow(data_sh)], 1, 10, nrow(data_sh), 0.93, 1.04, 1.11)
# qplot(seq(29,nrow(data_sh),0.1), show)

######################################################################################################


comp_bootstrap_std <- function(y_ser, bs_interp, kappa, mu, tau, lambda, interg_interp, lr, decay, iters, lr_theta, iters_theta, bootstrap_round){
  rst <- opt_bs_coeff(y_ser, bs_interp, kappa, mu, tau, lambda, interg_interp, lr, decay, iters, FALSE, FALSE)
  std_ <- rst[length(rst)]
  mean_ <- rst[1:(length(rst)-1)]
  theta_mat <- t(matrix(c(kappa, mu, tau)))
  for (br in 1:bootstrap_round) {
    y_ser_turb <- mean_ + rnorm(length(mean_), 0, std_)
    theta <- opt_dde_theta(y_ser_turb, bs_interp, kappa, mu, tau, lambda, interg_interp, lr, decay, iters, lr_theta, iters_theta, FALSE)
    theta_mat <- rbind(theta_mat, t(matrix(theta)))
  }
  return(theta_mat)
}


theta_mat_prev <- comp_bootstrap_std(data_sh$total_affected[1:28], 0.5, 1.60, 1.17, 1.42, 25., 500, 0.1, 0.99, 2e2, 0.05, 2e1, 100)
theta_mat_lockdown <- comp_bootstrap_std(data_sh$total_affected[28:nrow(data_sh)], 0.5, 0.72, 1.05, 2.18, 25., 500, 0.1, 0.99, 2e2, 0.05, 2e1, 100)

######################################################################################################

covid19_dde_predict <- function(ser, interp, total_len, kappa, mu, tau, sd_maxifier, err_maxifier){
  # initialization
  init_num <- floor((length(ser) - 1) * interp)
  init_ser <- c(ser[length(ser)])
  for (n in 1:(init_num-1)) {
    tail_left_index <- length(ser) - n / interp
    l_ratio <- tail_left_index - floor(tail_left_index)
    inter_val <- ser[floor(tail_left_index)] * (1-l_ratio) + ser[ceiling(tail_left_index)] * l_ratio
    init_ser <- c(inter_val, init_ser)
  }
  minus_steps <- round(tau*interp)
  lift_bias <- numeric(0)
  for (j in 1:(init_num - minus_steps)) {
    lift_bias <- c(lift_bias, (init_ser[j+minus_steps]-init_ser[j+minus_steps-1])*interp - kappa * (init_ser[j+minus_steps] - mu * init_ser[j]))
  }
  fit <- smooth.spline(1:length(lift_bias),lift_bias,df=3)
  # plot(fit)
  # initialization finished
  while(length(init_ser)-init_num < total_len * interp){
    deriv_ser <- max(0, sd_maxifier * sd(fit$y - fit$yin) + err_maxifier * predict(fit, length(init_ser)-minus_steps+1)$y + kappa * (init_ser[length(init_ser)] - mu * init_ser[length(init_ser)-minus_steps]))
    init_ser <- c(init_ser, init_ser[length(init_ser)]+deriv_ser/interp)
  }
  return(init_ser[init_num:length(init_ser)])
}

gen_covid19_dde_cis <- function(ser, interp, total_len, theta_mat, end_date, sd_maxifier, err_maxifier){
  pred_cis <- matrix(nrow=0, ncol=total_len*interp+1)
  for (r in 1:nrow(theta_mat)) {
    # pred_ser <- covid19_dde_predict(ser, interp, total_len, theta_mat[r,1], theta_mat[r,2], theta_mat[r,3], 0, err_maxifier)
    # pred_cis <- rbind(pred_cis, t(matrix(pred_ser)))
    pred_ser <- covid19_dde_predict(ser, interp, total_len, theta_mat[r,1], theta_mat[r,2], theta_mat[r,3], sd_maxifier, err_maxifier)
    pred_cis <- rbind(pred_cis, t(matrix(pred_ser)))
  }
  pred_ci_min <- apply(pred_cis,2,min)[seq(1,interp*total_len+1,interp)]
  pred_ci_max <- apply(pred_cis,2,max)[seq(1,interp*total_len+1,interp)]
  pred_ci_avg <- apply(pred_cis,2,mean)[seq(1,interp*total_len+1,interp)]
  return(data.frame(date=seq(end_date, end_date+total_len, 1), ci_min=pred_ci_min, ci_max=pred_ci_max, pred=(pred_ci_avg+pred_ci_min+pred_ci_max)/3.))
}

prev_df <- data_sh[1:28,c("date", "total_add", "total_affected")]
prev_df$period <- "before lockdown"
prev_pred <- opt_bs_coeff(data_sh$total_affected[1:28], 0.5, 1.60, 1.17, 1.42, 25., 500, 0.1, 0.99, 2e2, FALSE, FALSE)
prev_df$pred <- prev_pred[1:(length(prev_pred)-1)]
prev_df$total_add[nrow(prev_df)] <- 0

lock_df <- data_sh[28:nrow(data_sh),c("date", "total_add", "total_affected")]
lock_df$period <- "after lockdown"
lock_pred <- opt_bs_coeff(data_sh$total_affected[28:nrow(data_sh)], 0.5, 0.72, 1.05, 2.18, 25., 500, 0.1, 0.99, 2e2, FALSE, FALSE)
lock_df$pred <- lock_pred[1:(length(lock_pred)-1)]
fit_len <- 7
lock_tail_fit <- smooth.spline(1:fit_len,seq(0,1,length.out=fit_len)*(lock_df$total_affected-lock_df$pred)[(nrow(lock_df)-fit_len+1):nrow(lock_df)],df=3)
lock_df$pred[(nrow(lock_df)-fit_len+1):nrow(lock_df)] <- lock_df$pred[(nrow(lock_df)-fit_len+1):nrow(lock_df)] + lock_tail_fit$y


clip <- function(x, vmin, vmax){
  return(min(vmax,max(vmin,x)))
}

pred_cis <- gen_covid19_dde_cis(lock_df$pred, 10, 14, theta_mat_lockdown, data_sh$date[nrow(data_sh)],4.85, 1.0)
pred_std <- lock_pred[length(lock_pred)] * c(0, exp(seq(-1,-1, length.out = nrow(pred_cis)-1)))
pred_cis$total_add <- c(NaN,diff(pred_cis$pred))
pred_cis$period <- "forecasted"
pred_cis$total_affected <- NaN
pred_cis$ci_min <- sapply(pred_cis$ci_min-pred_std/2,clip,vmin=0,vmax=Inf)
pred_cis$ci_max <- pred_cis$ci_max + pred_std/2
pred_cis$add_min <- sapply(c(NaN,diff(pred_cis$ci_min))-pred_std/2,clip,vmin=0,vmax=Inf)
pred_cis$add_max <- c(NaN,diff(pred_cis$ci_max)) + pred_std/2


final_df <- rbind(prev_df, lock_df)
final_df$ci_min <- final_df$total_affected
final_df$ci_max <- final_df$total_affected
final_df$add_min <- NaN
final_df$add_max <- NaN
final_df <- rbind(final_df, pred_cis)

base_plot <- ggplot(data=final_df, aes(x=date, y=total_affected, color=period, linetype=period, shape=period, alpha=period, size=period)) + 
  scale_size_manual(values=c("before lockdown"=2,"after lockdown"=2,"forecasted"=0)) +
  scale_alpha_manual(values=c("before lockdown"=1,"after lockdown"=1,"forecasted"=1)) +
  scale_shape_manual(values=c("before lockdown"=16,"after lockdown"=16,"forecasted"=4)) +
  scale_color_manual(values=c("before lockdown"="gold3","after lockdown"="red3","forecasted"="blue3")) +
  scale_fill_manual(values=c("before lockdown"="gold3","after lockdown"="red3","forecasted"="blue3")) +
  scale_linetype_manual(values=c("before lockdown"="solid","after lockdown" = "solid","forecasted"="dashed")) +
  theme(
    legend.direction = "vertical",
    legend.box = "horizontal",
    axis.title.x=element_blank(),
  )

start_date <- as.Date("2022-03-01", format="%Y-%m-%d")
end_date <- as.Date("2022-05-07", format="%Y-%m-%d")

trend <- base_plot + geom_line(aes(y=pred), size=1, alpha=0.3) + geom_point() + 
  geom_ribbon(aes(ymin=ci_min, ymax=ci_max), linetype=3, alpha=0.2, size=0) +
  scale_x_date(breaks=seq(start_date,end_date,by="1 days"),limits=c(start_date-1,end_date+1)) +
  scale_y_continuous(breaks=seq(0, 6e5, 5e4),labels=label_number(suffix=" W",scale=1e-4)) + 
  ylab("Accumulated Confirmed Affected") +
  theme(axis.text.x=element_text(size=10, angle=90, hjust=-2, vjust=0),
        legend.position = c(0.01, 0.99),
        legend.justification = c("left", "top"))

growth <- base_plot +  
  geom_bar(aes(y=total_add, fill=period),size=0.2,position="stack",stat="identity") +
  geom_errorbar(aes(ymin=add_min, ymax=add_max),width=.4,alpha=0.5,size=.4,linetype="solid") +
  scale_x_date(breaks=seq(start_date,end_date,by="1 days"),position="top",limits=c(start_date-1,end_date+1)) +
  scale_y_reverse(breaks=seq(0, 3e4, 5e3),labels=label_number(suffix=" K",scale=1e-3),limits=c(3e4,0)) +
  scale_alpha_manual(values=c("before lockdown"=.5,"after lockdown"=.5,"forecasted"=.2)) +
  ylab("Daily Added Confirmed") +
  theme(axis.text.x=element_blank(),
        legend.position = c(0.01, 0.01),
        legend.justification = c("left", "bottom"),)

grid.arrange(trend, growth, ncol=1, heights=c(3,1))

######################################################################################################

base_plot <- ggplot(data=final_df[which(final_df$period!="before lockdown"),], aes(x=date, y=total_affected, color=period, linetype=period, shape=period, alpha=period, size=period)) + 
  scale_size_manual(values=c("after lockdown"=2,"forecasted"=0)) +
  scale_alpha_manual(values=c("after lockdown"=1,"forecasted"=1)) +
  scale_shape_manual(values=c("after lockdown"=16,"forecasted"=4)) +
  scale_color_manual(values=c("after lockdown"="red3","forecasted"="blue3")) +
  scale_fill_manual(values=c("after lockdown"="red3","forecasted"="blue3")) +
  scale_linetype_manual(values=c("after lockdown" = "solid","forecasted"="dashed")) +
  theme(
    legend.direction = "vertical",
    legend.box = "horizontal",
    axis.title.x=element_blank(),
  )

start_date <- as.Date("2022-03-28", format="%Y-%m-%d")
end_date <- as.Date("2022-05-07", format="%Y-%m-%d")

trend <- base_plot + geom_line(aes(y=pred), size=1, alpha=0.3) + geom_point() + 
  geom_ribbon(aes(ymin=ci_min, ymax=ci_max), linetype=3, alpha=0.2, size=0) +
  scale_x_date(breaks=seq(start_date,end_date,by="1 days"),limits=c(start_date-1,end_date+1)) +
  scale_y_continuous(breaks=seq(0, 6e5, 5e4),labels=label_number(suffix=" W",scale=1e-4)) + 
  ylab("Accumulated Confirmed Affected") +
  theme(axis.text.x=element_text(size=10, angle=90, hjust=-2, vjust=0),
        legend.position = c(0.01, 0.99),
        legend.justification = c("left", "top"))

growth <- base_plot +  
  geom_bar(aes(y=total_add, fill=period),size=0.2,position="stack",stat="identity") +
  geom_errorbar(aes(ymin=add_min, ymax=add_max),width=.4,alpha=0.5,size=.4,linetype="solid") +
  scale_x_date(breaks=seq(start_date,end_date,by="1 days"),position="top",limits=c(start_date-1,end_date+1)) +
  scale_y_reverse(breaks=seq(0, 3e4, 5e3),labels=label_number(suffix=" K",scale=1e-3),limits=c(3e4,0)) +
  scale_alpha_manual(values=c("after lockdown"=.5,"forecasted"=.2)) +
  ylab("Daily Added Confirmed") +
  theme(axis.text.x=element_blank(),
        legend.position = c(0.01, 0.01),
        legend.justification = c("left", "bottom"),)

grid.arrange(trend, growth, ncol=1, heights=c(3,1))


library(scatterplot3d)
source('addgrids3d.r')
theta_df <- data.frame(theta_mat_lockdown)
colnames(theta_df) <- c("kappa", "mu", "tau")
scatterplot3d(theta_df, type="h", pch = 1, color="darkgreen", grid=FALSE, box=FALSE)
addgrids3d(theta_df, grid = c("xy", "xz", "yz"))
