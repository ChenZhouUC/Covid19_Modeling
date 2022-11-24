import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skfda.representation.basis import BSpline
from copy import deepcopy

pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)
pd.set_option("display.float_format", lambda x: "%.3f" % x)
pd.set_option("chained_assignment", None)

plt.rcParams["font.family"] = ["Arial Unicode MS"]  # Chinese Labels
plt.rcParams["axes.unicode_minus"] = False  # Minus Sign
plt.rcParams["figure.max_open_warning"] = 0

sns.set(
    style="darkgrid",
    rc={
        "figure.figsize": (12, 6),
        "font.sans-serif": ["Arial Unicode MS", "Arial"],
    },
)

def revise_ser(ser, intercept=None):
    ser = np.array(ser)
    revised_zero_diff = np.clip(
        np.concatenate([[1], (ser[1:] - ser[:-1])]), 0, np.inf
    )
    check_zeros = (revised_zero_diff.cumprod() != 0)
    revised_zero = (revised_zero_diff * check_zeros).cumsum()
    if intercept is None:
        intercept = ser.mean() - revised_zero.mean()
    else:
        intercept = intercept - 1
    return revised_zero + intercept

def bs_coeff_local_opt(
    ser,
    kappa,
    mu,
    tau,
    bs_interp_prop,
    interg_interp_prop,
    labd,
    lr,
    decay,
    iters,
    early_stop_ratio,
    early_stop_steps,
    show_process,
    ser_loss_expansion=0,
):
    # define interp
    x = range(len(ser))
    knots = np.linspace(0, len(ser) - 1, round(bs_interp_prop * (len(ser) - 1) + 1))
    integ_x = np.linspace(
        0, len(ser) - 1, round(interg_interp_prop * (len(ser) - 1) + 1)
    )
    loss_weights = np.linspace(1 - ser_loss_expansion, 1 + ser_loss_expansion, len(ser))
    # define bs funcs
    bs_func = BSpline(domain_range=(min(x), max(x)), order=4, knots=knots)
    bs_func_d = bs_func.derivative()
    # obtain interpped bs rst
    bS = bs_func(x).squeeze(2)
    bs_integ = bs_func(integ_x).squeeze(2)
    bs_partial = bs_func_d(integ_x).squeeze(2)
    # construct diff features
    shift_ = round(tau * interg_interp_prop)
    now_bs_integ = bs_integ[:, shift_:]
    delay_bs_integ = bs_integ[:, :-shift_]
    now_bs_partial = bs_partial[:, shift_:]
    delay_bs_partial = bs_partial[:, :-shift_]
    loss_weights_shift = np.linspace(
        1 - ser_loss_expansion,
        1 + ser_loss_expansion,
        round(interg_interp_prop * (len(ser) - 1) + 1),
    )
    loss_weights_shift = loss_weights_shift[
        int(shift_ / 2) : (len(loss_weights_shift) - shift_ + int(shift_ / 2))
    ]
    # initialize parameters
    n_coeff = bS.shape[0]
    init_coeff = np.repeat(0, n_coeff)
    init_grad = np.repeat(0, n_coeff)

    loss_update = []
    losses = [np.inf]
    continuous_remain = 0
    for _iter in range(iters):
        # loss calculation
        loglik_loss_vec = (np.dot(init_coeff, bS) - ser) * loss_weights
        loglik_loss = np.linalg.norm(loglik_loss_vec, 2)
        mat_ = now_bs_partial - kappa * now_bs_integ + kappa * mu * delay_bs_integ
        integ_loss_vec = np.dot(init_coeff, mat_ * loss_weights_shift).astype(np.float)
        integ_loss = (
            np.linalg.norm(integ_loss_vec, 2) / np.sqrt(interg_interp_prop) * labd
        )
        loss_ = loglik_loss + integ_loss

        if (
            loss_ == np.inf
            or loss_ / min(losses) > early_stop_ratio[0]
            or continuous_remain > early_stop_steps
        ):
            break
        if loss_ / min(losses) >= early_stop_ratio[1]:
            continuous_remain = np.clip(continuous_remain + 1, 0, np.inf)
        else:
            continuous_remain = 0

        # gradient descent
        loglik_part = np.dot(bS * loss_weights, loglik_loss_vec).astype(np.float)
        integ_part = (
            np.dot(
                mat_ * loss_weights_shift, np.dot(init_coeff, mat_ * loss_weights_shift).astype(np.float)
            ).astype(np.float)
            / interg_interp_prop
            * labd**2
        )
        init_coeff = init_coeff - lr * (loglik_part + integ_part)
        lr *= decay
        losses.append(loss_)
        if show_process:
            # print(
            #     "loglik: {} integ: {}".format(
            #         loglik_loss / np.sqrt(len(ser)),
            #         integ_loss / labd / np.sqrt(len(ser)),
            #     )
            # )
            loss_update.append(
                [
                    _iter,
                    loglik_loss / np.sqrt(len(ser)),
                    integ_loss / labd / np.sqrt(len(ser)),
                ]
            )

    coef_mat = np.dot(bS, bS.T).astype(np.float) + np.dot(mat_, mat_.T) / interg_interp_prop * labd**2
    k_mat = np.dot(
        np.dot(
            now_bs_integ - mu * delay_bs_integ,
            (now_bs_integ - mu * delay_bs_integ - now_bs_partial).T,
        ).astype(np.float),
        init_coeff,
    ).astype(np.float)
    m_mat = np.dot(
        (mu * kappa**2) * np.dot(delay_bs_integ, delay_bs_integ.T)
        + kappa * np.dot(delay_bs_integ, (now_bs_partial - kappa * now_bs_integ).T),
        init_coeff,
    ).astype(np.float)
    t_mat = -np.dot(
        (
            2
            * kappa
            * mu
            * np.dot(
                now_bs_partial - kappa * now_bs_integ + kappa * mu * delay_bs_integ,
                delay_bs_partial.T,
            )
            + (now_bs_partial - kappa * now_bs_integ + kappa * mu * delay_bs_integ)[
                :, [0]
            ]
            ** 2
        ),
        init_coeff,
    ).astype(np.float)
    partial_resid_theta = (
        np.dot(np.dot(bS.T, np.linalg.inv(coef_mat)).astype(np.float), np.array([k_mat, m_mat, t_mat]).T)
        * 2
        * labd**2
        / interg_interp_prop
    )
    theta_grad = np.dot(
        np.dot(
            np.linalg.inv(
                np.dot(partial_resid_theta.T, partial_resid_theta).astype(np.float)
            ),
            partial_resid_theta.T,
        ).astype(np.float),
        ser - np.dot(bS.T, init_coeff),
    ).astype(np.float)
    fitted_ser = np.dot(bS.T, init_coeff).astype(np.float)
    if show_process:
        f, axes = plt.subplots(1, 2, figsize=(12, 4))
        sns.scatterplot(x=x, y=ser, alpha=1.0, color="blue", s=15, ax=axes[0])
        sns.lineplot(
            x=integ_x,
            y=np.dot(bs_integ.T, init_coeff),
            alpha=0.8,
            color="yellow",
            ax=axes[0],
        )

        loss_df = pd.DataFrame(
            loss_update, columns=["iter", "loglik_loss", "integ_loss"]
        )
        ax = sns.lineplot(
            data=loss_df, x="iter", y="loglik_loss", color="green", ax=axes[1]
        )
        ax2 = ax.twinx()
        sns.lineplot(data=loss_df, x="iter", y="integ_loss", color="red", ax=ax2)
    return theta_grad, loss_, fitted_ser, (init_coeff, bs_func)

def model_theta_global_opt(
    ser,
    kappa_init,
    mu_init,
    tau_init,
    bs_interp_prop,
    interg_interp_prop,
    labd,
    bs_lr_schedule,
    theta_lr_schedule,
    early_stop_ratio,
    early_stop_steps,
    show_process="outer",
    ser_loss_expansion=0.5,
):
    losses = [np.inf]
    theta = np.array([kappa_init, mu_init, tau_init])
    theta_selected = theta.copy()
    iters_bs, lr_bs, decay_bs = bs_lr_schedule
    iters_theta, lr_theta, decay_theta = theta_lr_schedule

    if show_process == "outer":
        show_process_bs, show_process_theta = False, True
    elif show_process == "inner":
        show_process_bs, show_process_theta = True, False
    elif show_process == "both":
        show_process_bs, show_process_theta = True, True
    else:
        show_process_bs, show_process_theta = False, False

    continuous_remain = 0
    for iter_theta_ in range(iters_theta):
        theta_grad, loss_, fitted_ser, bs_funcs = bs_coeff_local_opt(
            ser,
            theta[0],
            theta[1],
            theta[2],
            bs_interp_prop,
            interg_interp_prop,
            labd,
            lr_bs,
            decay_bs,
            iters_bs,
            early_stop_ratio,
            early_stop_steps,
            show_process_bs,
            ser_loss_expansion,
        )
        if (
            loss_ == np.inf
            or loss_ / min(losses) > early_stop_ratio[0]
            or continuous_remain > early_stop_steps
        ):
            break
        if loss_ / min(losses) >= early_stop_ratio[1]:
            continuous_remain = np.clip(continuous_remain + 1, 0, np.inf)
        else:
            continuous_remain = 0

        if min(losses) > loss_:
            theta_selected = theta.copy()
            fitted_ser_selected = fitted_ser.copy()
            bs_funcs_selected = deepcopy(bs_funcs)

        theta -= lr_theta * theta_grad
        lr_theta = lr_theta * decay_theta
        losses.append(loss_)
        if show_process_theta:
            print("inter {}: Theta: {} Loss: {}".format(iter_theta_, theta, loss_))
    if show_process_theta:
        plt.figure(figsize=(12, 6))
        ax = sns.lineplot(x=range(len(losses)), y=losses)
        ax = sns.scatterplot(x=range(len(losses)), y=losses)
    return theta_selected, losses, fitted_ser_selected, bs_funcs_selected


def next_pred(ser_for_now, pred_interp, theta):
    kappa, mu, tau = theta
    delta_t = 1 / pred_interp
    tau_shift_floor = int(tau * pred_interp)
    tau_shift_float = tau * pred_interp - tau_shift_floor
    ser_tau = (
        tau_shift_float * ser_for_now[-(tau_shift_floor + 2)]
        + (1 - tau_shift_float) * ser_for_now[-(tau_shift_floor + 1)]
    )
    return ser_for_now[-2] + np.clip(
        2 * delta_t * kappa * (ser_for_now[-1] - mu * ser_tau), 0, np.inf
    )


def pred_ser(ser, pred_length, pred_interp, pred_prev_len, theta):
    pred_prepare = np.array([ser[-1]])
    for _i in range(pred_prev_len):
        pred_prepare_append = np.linspace(
            ser[-2 - _i], ser[-1 - _i], pred_interp, endpoint=False
        )
        pred_prepare = np.concatenate([pred_prepare_append, pred_prepare])
    for _i in range(pred_length * pred_interp):
        next_tic = next_pred(pred_prepare, pred_interp, theta)
        pred_prepare = np.append(pred_prepare, next_tic)
    return pred_prepare

def detect_explosion_point(ser, top=1, daily_add_thresh=20):
    ser_diff = ser[1:] - ser[:-1]
    thresh_explosion = np.argmax(ser_diff >= daily_add_thresh) + 1
    return (1 + np.argsort((ser[2:] + ser[:-2]) - 2 * ser[1:-1])[::-1])[:top], thresh_explosion

def growth_decision(fitted_ser, theta_opt):
    return (theta_opt.prod() - 1) * (fitted_ser[-1] - fitted_ser[-2]) >= theta_opt[
        0
    ] * (theta_opt[1] - 1) * fitted_ser[-1]


def covid19_seer(cumulative_affected, seer_type, theta_city, five_day_filter=25, lock_down_hard_thresh=200):
    assert seer_type in [
        "WARNING",
        "LOCKDOWN",
    ], "seer_type should be in [WARNING,LOCKDOWN]!!!"
    
    if seer_type == "WARNING":
        last_five_day = cumulative_affected[-5:]
        if (last_five_day[-1] - last_five_day[0]) / len(last_five_day) <= five_day_filter / 5:
            return False, None, None, None

    kappa_init, mu_init, tau_init = theta_city
    theta_opt, losses, fitted_opt, bs_funcs = model_theta_global_opt(
        cumulative_affected,
        kappa_init=kappa_init,
        mu_init=mu_init,
        tau_init=tau_init,
        bs_interp_prop=0.5,
        interg_interp_prop=200.0,
        labd=5.0,
        bs_lr_schedule=(int(5e2), 0.02, 0.999),
        theta_lr_schedule=(int(5e2), 0.05, 0.99),
        early_stop_ratio=[1.5, 0.999],
        early_stop_steps=10,
        show_process="neither",
        ser_loss_expansion=0.5,
    )
    print("fitted theta: {}".format(theta_opt))
    kappa_fit, mu_fit, tau_fit = theta_opt
    fitted_ser = revise_ser(fitted_opt)

    if growth_decision(fitted_ser, theta_opt):
        growth_flag = True
        theta_selection = theta_city
    else:
        growth_flag = False
        theta_selection = (theta_opt + theta_city) / 2.0

    if seer_type == "WARNING":
        pred_interp = 10
        pred_len = 30
        pred_prev_len = int(np.ceil(tau_fit))

        pred_ser_combo = pred_ser(
            cumulative_affected, pred_len, pred_interp, pred_prev_len, theta_opt
        )
        pred_series = pred_ser_combo[pred_interp * pred_prev_len :][::pred_interp]
        pred_series = revise_ser(pred_series, intercept=pred_series[0])
        whole_ser = np.concatenate([cumulative_affected, pred_series[1:]])
        
        cc_selection, thresh_selection = detect_explosion_point(whole_ser, 1)
        turning_point = max(min(cc_selection[0], thresh_selection), len(cumulative_affected))
        before_lockdown = whole_ser[: (turning_point + 1)]
        before_lockdown = before_lockdown[before_lockdown <= lock_down_hard_thresh]
        turning_point = min(turning_point, len(before_lockdown) - 1)
    else:
        turning_point = None
        before_lockdown = cumulative_affected
    print(turning_point)
    pred_interp = 10
    pred_len = 60
    pred_prev_len = int(np.ceil(tau_fit))

    pred_ser_combo = pred_ser(
        before_lockdown, pred_len, pred_interp, pred_prev_len, theta_selection
    )
    pred_series = pred_ser_combo[pred_interp * pred_prev_len :][::pred_interp]
    pred_series = revise_ser(pred_series, intercept=pred_series[0])
    whole_ser = np.concatenate([before_lockdown, pred_series[1:]])
    stop_point = (pred_series[1:] - pred_series[:-1]).argmin() + len(before_lockdown)

    return (
        growth_flag,
        turning_point,
        stop_point,
        whole_ser,
    )