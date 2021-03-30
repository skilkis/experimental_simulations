import os
from contextlib import contextmanager
from functools import cached_property
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Sequence, Tuple, Union

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from numpy.polynomial import Polynomial
from scipy.io import loadmat

ROOT_DIR: Path = Path(__file__).parent.parent
FIGURE_PATH: Path = ROOT_DIR / "figures"

DATA = loadmat(
    ROOT_DIR / "data" / "processed.mat",
    squeeze_me=True,
    struct_as_record=False,
)

PLOT: bool = True


@contextmanager
def figure(
    name: str,
    directory: Path = FIGURE_PATH,
    style: str = "bmh",
    figsize: Tuple[float, float] = (12, 6),
    **kwargs: Dict[str, Any],
):
    """Creates figure handles and automatically shows/saves on exit.

    Args:
        name: Name of the figure which becomes the filename
        directory: Location to save the current figure
        style: Style used for the matplotlib plot
        kwargs: Keyword arguments to pass to plt.subplots
    """
    try:
        with plt.style.context(style):
            # Modify rcParams locally
            plt.style.use({"lines.linewidth": 1})
            fig, ax = plt.subplots(num=name, figsize=figsize, **kwargs)
            yield fig, ax
    finally:
        if os.getenv("CI") is None and PLOT:
            plt.show()
        fig.savefig(fname=directory / f"{name}.pdf", bbox_inches="tight")


class Regression:
    """Runs a regression on (x, y) with a polynomial of ``deg``."""

    def __init__(self, x: np.ndarray, y: np.ndarray, deg: int = 1):
        self.x = x
        self.y = y
        self.deg = deg

    @cached_property
    def polynomial(self) -> Polynomial:
        """Fitted polynomial."""
        return Polynomial.fit(x=self.x, y=self.y, deg=self.deg)

    @cached_property
    def y_hat(self) -> np.ndarray:
        """Predicted data values."""
        return self.polynomial(self.x)

    @cached_property
    def y_bar(self) -> np.ndarray:
        """Mean of the observed data."""
        return np.sum(self.y) / len(self.y)

    @cached_property
    def rss(self) -> float:
        """Residual sum of squares."""
        return np.sum((self.y - self.y_hat) ** 2)

    @cached_property
    def tss(self) -> float:
        """Total sum of squares."""
        return np.sum((self.y - self.y_bar) ** 2)

    @cached_property
    def r_squared(self) -> float:
        """Coefficient of Determination (R^2)."""
        return 1 - (self.rss / self.tss)

    def get_plot_label(
        self,
        x_label: str,
        y_label: str,
        power_fmt: str,
        coef_fmt: str,
        show_stats: bool = True,
    ) -> str:
        rhs = ""
        for power, coef in reversed(tuple(enumerate(self.polynomial.coef))):
            c_str = f"{coef_fmt.format(coef)}{x_label}"
            x_str = f"{x_label}" if power > 1 else ""
            p_str = f"^{{{power_fmt.format(power)}}}" if power > 1 else ""
            if power == self.polynomial.coef.size - 1:  # 1st coefficient?
                rhs += f"{c_str}{x_str}{p_str}"
            else:
                oper = " - " if c_str.startswith("-") else " + "
                rhs += f"{oper}{c_str}{x_str}{p_str}"
        equation = f"${y_label} = {rhs}$"
        stats = f"$R^2 = {self.r_squared:.4f}$"
        return f"{equation}\n{stats}" if show_stats else equation


# TODO fully document
def regplot(
    lines: Union[Line2D, Sequence[Line2D]],
    deg: int = 1,
    ci: Optional[float] = None,
    x_label: str = "x",
    y_label: str = "y",
    coef_fmt: str = "{:.4f}",
    power_fmt: str = "{:.0f}",
) -> Union[Regression, Sequence[Regression]]:
    """Adds a regression of ``degree`` onto ``lines``.

    Args:
        depvar: Dependent variable label"""
    regs = []
    for line in lines if isinstance(lines, Iterable) else [lines]:
        reg = Regression(x=line.get_xdata(), y=line.get_ydata())
        label = reg.get_plot_label(
            x_label=x_label,
            y_label=y_label,
            coef_fmt=coef_fmt,
            power_fmt=power_fmt,
        )
        sns.regplot(
            x=line.get_xdata(),
            y=line.get_ydata(),
            order=deg,
            ci=ci,
            color=line.get_color(),
            line_kws=dict(linestyle="--"),
            label=label,
            scatter=False,
        )
        regs.append(reg)
    return regs if len(regs) > 1 else regs[0]


# Quantify effect of the slipstream on directional stability at Beta = 0
with figure("SlipstreamJSweep") as (fig, ax):
    d = DATA["BalData"].windOn.j_sweep
    unc = ax.plot(d.TC1, -d.CMy, label="Total", marker="+", linewidth=1)
    reg = regplot(
        lines=unc,
        deg=1,
        ci=95,
        x_label=r"\beta",
        y_label=r"C_n",
    )
    unc_aero = ax.plot(
        d.TC1, -d.CMya, label="Aerodynamic", marker="x", linewidth=1
    )
    ax.set_xlabel("Thrust Coefficient $T_{c_{L}}$")
    ax.set_ylabel("Yaw Moment Coefficient $C_n$")
    ax.legend(loc="best")

# Plot of Cn beta for velocity 20 and 40, zero thrust
# TODO find out why negative Cn_beta does not work?
with figure("CnBetaRudder0Thrust0") as (fig, ax):
    d = DATA["BalData"].windOn.rudder0  # Select relevant data
    d_corr = DATA["BalDataCorr"].rudder0.Total
    for velocity, zero_thrust_rpm in (points := {20: 41.6, 40: 78.5}).items():
        # Filtering data at the current velocity, AoA, and zero thrust
        idx_v = np.isclose(d.V, velocity, atol=1)
        idx_aoa = np.isclose(d.AoA, 0, atol=1)
        idx_m1 = np.isclose(d.rpsM1, zero_thrust_rpm, atol=1)
        idx_m2 = np.isclose(d.rpsM2, zero_thrust_rpm, atol=1)
        idx = idx_v & idx_aoa & idx_m1 & idx_m2

        # Creating main line plot
        sorted_idx = np.argsort(d.AoS[idx])
        corrected = ax.plot(
            d_corr.AoS[idx][sorted_idx],
            -d_corr.CMy[idx][sorted_idx],
            label=f"V={velocity} m/s",
            marker="+",
            zorder=2 if velocity == 20 else 1,
        )

        # Adding uncorrected data
        uncorrected = ax.plot(
            d.AoS[idx][sorted_idx],
            -d.CMy[idx][sorted_idx],
            label=f"V={velocity} m/s Uncorrected",
            marker="x",
            color=corrected[-1].get_color(),
            zorder=0,
            alpha=0.5,
        )

        # Adding Regression
        reg = regplot(
            lines=corrected,
            deg=1,
            ci=95,
            x_label=r"\beta",
            y_label=r"C_n",
        )
    ax.set_xlabel("Angle of Sideslip $\\beta$")
    ax.set_ylabel("Yaw Moment Coefficient $C_n$")
    ax.legend(loc="best", ncol=len(points))

# Cn vs beta, Cndelta vs. beta, and w/ AoA

# Only V40 for AoA effect, this is limited to potentially not capturing
# non-linearties at low reynolds due to sepration.

# TODO show all the different correction contributions on Cn B
# TODO Noise measurements (idea use noise to estimate separation)
# Non-dimensionalize w/ pressure jump

# TODO Find scaling for noise