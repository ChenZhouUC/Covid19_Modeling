import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skfda.representation.basis import BSpline

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
    revised_zero = np.clip(
        np.concatenate([[0], (ser[1:] - ser[:-1])]), 0, np.inf
    ).cumsum()
    if intercept is None:
        intercept = ser.mean() - revised_zero.mean()
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
    return theta_grad, loss_, fitted_ser

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
        theta_grad, loss_, fitted_ser = bs_coeff_local_opt(
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

        theta -= lr_theta * theta_grad
        lr_theta = lr_theta * decay_theta
        losses.append(loss_)
        if show_process_theta:
            print("inter {}: Theta: {} Loss: {}".format(iter_theta_, theta, loss_))
    if show_process_theta:
        plt.figure(figsize=(12, 6))
        ax = sns.lineplot(x=range(len(losses)), y=losses)
        ax = sns.scatterplot(x=range(len(losses)), y=losses)
    return theta_selected, losses, fitted_ser_selected