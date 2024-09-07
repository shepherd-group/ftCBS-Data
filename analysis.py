#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from numpy.typing import NDArray as Array
from scipy.optimize import curve_fit as cfit
from matplotlib.backends.backend_pdf import PdfPages

mpl.rc('text', usetex=True)
mpl.rc('savefig', dpi=100)
mpl.rc('hatch', lw=0.25)
mpl.rc('lines', lw=2, markersize=6)
mpl.rc('legend', fontsize=8, numpoints=1)
mpl.rc(('axes', 'xtick', 'ytick'), labelsize=8)
mpl.rc('figure', dpi=200, figsize=(3.37, 3.37*(np.sqrt(5)-1)/2))
mpl.rc('font', **{'family': 'serif', 'sans-serif': 'Computer Modern Roman'})

# GLOBAL DATA
Eh = r'$E_\mathrm{h}$'
# GLOBAL DATA


class Getter:
    def __init__(self, csvs: str | list[str], minM: int = None) -> None:
        self.read_csv(csvs, minM)
        self.M = np.unique(self.raw.M.values)
        self.theta = np.unique(self.raw.theta.values)
        self._Mloc = self.raw.M.values
        self._Tloc = self.raw.theta.values

    def read_csv(self, csvs: str | list[str], minM: int) -> None:
        if isinstance(csvs, str):
            csvs = [csvs]

        raw = pd.concat([pd.read_csv(csv) for csv in csvs])

        if minM is not None:
            raw = raw[raw.M.values >= minM].reset_index(drop=True)

        if 'theta' not in raw.columns:
            raw['theta'] = raw.temp.values

        self.raw = raw
        self.raw.reset_index(drop=True, inplace=True)
        return

    def __getitem__(self, vlocs: tuple) -> pd.DataFrame:
        M, T = vlocs

        if isinstance(M, slice):
            M = self._Mloc[M].round(12)
        elif isinstance(M, float):
            M = round(M, 12)

        if isinstance(T, slice):
            T = self._Tloc[T].round(12)
        elif isinstance(T, float):
            T = round(T, 12)

        mask = np.isin(self._Mloc.round(12), M)
        mask = mask & np.isin(self._Tloc.round(12), T)

        sf = self.raw.copy()[mask].reset_index(drop=True)

        sf.sort_values(by=['M', 'beta'], inplace=True, ignore_index=True)

        return sf


class Data(Getter):
    def __init__(
            self,
            csvs: str | list[str],
            minM: int = None,
        ) -> None:

        Getter.__init__(self, csvs, minM)

        self.raw['-TS'] = -(1.0/self.raw.beta.values)*self.raw.S.values


class nqData(Getter):
    def __init__(
            self,
            csvs: str | list[str],
            rs: float = 1.0,
            N: int = 2,
        ) -> None:

        Getter.__init__(self, csvs)

        self.rs = rs
        self.N = N
        self.V = (4*np.pi*N*rs*rs*rs/3)
        self.L = self.V**(1/3)


class SqData(Getter):
    def __init__(
            self,
            csvs: str | list[str],
            rs: float = 1.0,
            N: int = 2,
        ) -> None:

        Getter.__init__(self, csvs)

        self.rs = rs
        self.N = N
        self.V = (4*np.pi*N*rs*rs*rs/3)
        self.L = self.V**(1/3)

        self.raw['V_q'] = 4*np.pi/(self.V*(self.raw.q.values**2))
        self.raw['Ec_q'] = self.raw['Sc_q']*self.raw['V_q']*self.raw['D_q']
        self.raw['Eex_q'] = self.raw['Sx_q']*self.raw['V_q']*self.raw['D_q']
        self.raw['Eud_q'] = self.raw['Sud_q']*self.raw['V_q']*self.raw['D_q']

    @staticmethod
    def _vlocs_check(
            M: None | int | list[int],
            T: None | float | list[float],
        ) -> tuple:
        if M is None:
            M = slice(None, None, None)
        if T is None:
            T = slice(None, None, None)
        return M, T

    def getitem_wrapper(
            self,
            M: None | int | list[int],
            T: None | float | list[float],
        ) -> tuple[pd.DataFrame, list[int], list[float]]:
        M, T = self._vlocs_check(M, T)
        sf = self[M,T]
        uM = np.unique(sf.M.values)
        uT = np.unique(sf.theta.values)
        return sf, uM, uT

    def parallel_exchange(
            self,
            M: int | list[int] = None,
            T: float | list[float] = None,
        ) -> pd.DataFrame:

        _, uM, uT = self.getitem_wrapper(M, T)

        data = {
            'M': [],
            'theta': [],
            'Eex': [],
        }

        for m in uM:
            for t in uT:
                sf = self[m,t]
                data['M'].append(m)
                data['theta'].append(t)
                data['Eex'].append(sf.Eex_q.values.sum())

        return pd.DataFrame(data)

    def parallel_correlation(
            self,
            M: int | list[int] = None,
            T: float | list[float] = None,
        ) -> pd.DataFrame:

        _, uM, uT = self.getitem_wrapper(M, T)

        data = {
            'M': [],
            'theta': [],
            'Ec': [],
        }

        for m in uM:
            for t in uT:
                sf = self[m,t]
                data['M'].append(m)
                data['theta'].append(t)
                data['Ec'].append(sf.Ec_q.values.sum())

        return pd.DataFrame(data)

    def antiparallel_correlation(
            self,
            M: int | list[int] = None,
            T: float | list[float] = None,
        ) -> pd.DataFrame:

        _, uM, uT = self.getitem_wrapper(M, T)

        data = {
            'M': [],
            'theta': [],
            'Ec': [],
        }

        for m in uM:
            for t in uT:
                sf = self[m,t]
                data['M'].append(m)
                data['theta'].append(t)
                data['Ec'].append(sf.Eud_q.values.sum())

        return pd.DataFrame(data)


def scinot(vals: list[str]) -> str | list[str]:
    # Assume v is form aEb, where a and b are an integer.
    svals = []
    for v in vals:
        a, b = v.split('E')
        svals.append(fr'${a}\times 10^{{ {b} }}$')
    return svals


def theta_label(theta: float) -> str:
    sflookup = {str(i): None for i in range(1,10)}

    for i, c in enumerate(f'{theta:.12f}'):
        if c in sflookup:
            ilast = i
        elif c == '.':
            idecimal = i

    dp = max(ilast - idecimal, 1)

    dpfmt = f'.{dp}f'
    label = fr'$\theta_{{\xi=0}}={theta:{dpfmt}}$'

    return label


def temp_label(temp: float) -> str:
    sflookup = {str(i): None for i in range(1,10)}

    for i, c in enumerate(f'{temp:.12f}'):
        if c in sflookup:
            ilast = i
        elif c == '.':
            idecimal = i

    dp = max(ilast - idecimal, 1)

    dpfmt = f'.{dp}f'
    label = fr'$T={temp:{dpfmt}}$'

    return label


def linear(x: Array, a: float, b: float) -> Array:
    return a*x + b


def npoint(
        n: int,
        x: Array,
        y: Array,
        return_error: bool = False,
    ) -> tuple[float | tuple[float]]:

    if len(x) != n or len(y) != n:
        raise RuntimeError(f'lengths of x and/or y are not {n}!')

    info = cfit(
        linear,
        x,
        y,
        method='lm',
        maxfev=100000,
        ftol=1E-14,
        xtol=1E-14,
        full_output=True,
    )

    ab, ab_cov = info[:2]

    if return_error:
        return ab, np.sqrt(np.diag(ab_cov))

    return ab


def plot_internal_energy(
        d: Data,
        thetas: list[float],
        plotname: str,
    ) -> None:

    report = {'theta': [], 'CBS': [], 'ERR': []}

    cmap = mpl.colormaps['plasma'].resampled(len(thetas) + 3)

    with PdfPages(plotname) as pdf:

        plt.close()

        fig, axes = plt.subplots(
            nrows=1,
            ncols=2,
            sharex=False,
            sharey=False,
        )

        ax0, ax1 = axes

        for itheta, theta in enumerate(thetas):
            sf = d[:,theta]

            x = sf.M.values**(-1.0)
            y = sf.U.values

            n = 3
            xn = x[-n:]
            yn = y[-n:]

            (_, cbs), (_, cbs_err) = npoint(3, xn, yn, return_error=True)

            report['theta'].append(theta)
            report['CBS'].append(cbs)
            report['ERR'].append(cbs_err)

            dy = y - cbs

            ax0.plot(
                x,
                y,
                color=cmap(itheta+2),
                label=theta_label(theta),
                zorder=len(thetas) - itheta,
            )

            ax1.plot(
                x,
                dy,
                color=cmap(itheta+2),
                zorder=len(thetas) - itheta,
            )

        ax0.axhline(
            y=0,
            ls='--',
            color='k'
        )

        ax1.axhline(
            y=0,
            ls='--',
            color='k'
        )

        ax0.set_xlabel('$M^{-1}$')
        ax0.set_ylabel(f'$U$ / {Eh}')

        ax1.set_xlabel('$M^{-1}$')
        ax1.set_ylim(-0.0005, 0.00025)
        ax1.set_xlim(right=0.010)
        ax1.yaxis.tick_right()
        ax1.yaxis.set_label_position('right')
        ax1.set_ylabel(f'$U - U(\mathrm{{CBS}})$ / {Eh}')

        fig.legend(
            ncol=2,
            loc='upper center',
            bbox_to_anchor=(0.5, -0.04),
        )

        pdf.savefig(bbox_inches='tight')

    report = pd.DataFrame(report)
    print(report)

    return


def plot_all_ueg_quantities_seperated_delta_and_potential_single(
        d: Data,
        Sq: SqData,
        thetas: list[float],
        plotname: str,
        rs: float = 1.0,
    ) -> None:

    EudM = r'$E_{\uparrow\downarrow}^\mathrm{c}(M)$'
    EuucM = r'$E_{\uparrow\uparrow}^\mathrm{c}(M)$'
    EuuxM = r'$E_{\uparrow\uparrow}^\mathrm{ex}(M)$'
    EudCBS = r'$E_{\uparrow\downarrow}^\mathrm{c}(\mathrm{CBS})$'
    EuucCBS = r'$E_{\uparrow\uparrow}^\mathrm{c}(\mathrm{CBS})$'
    EuuxCBS = r'$E_{\uparrow\uparrow}^\mathrm{ex}(\mathrm{CBS})$'

    cmap = {
        0: 'C0',
        1: 'C1',
        2: 'C2',
        3: 'C3',
        4: 'C4',
        5: 'C6',
        6: 'k',
        7: 'C9',
    }

    lsmap = {
        0: '-',
        1: ':',
        2: '-',
        3: '-',
        4: '--',
        5: '--',
        6: ':',
        7: ':',
    }

    if rs == 1.0:
        ttloc = [
            [0.0025, 0.002],
            [0.0025, 0.002],
            [0.0025, 0.014],
            [0.0025, 0.004],
            [0.01, -20],
            [0.0025, -0.02],
        ]
    elif rs == 10.0:
        ttloc = [
            [0.0025, 0.00012],
            [0.0025, -0.00008],
            [0.0025, 0.0018],
            [0.0025, 0.0009],
            [0.01, -0.23],
            [0.01, -0.0],
        ]
    else:
        raise ValueError('The rs value could not be identified!')

    def _callable_getitem(T: float | list[float] = None) -> pd.DataFrame:
        if T is None:
            raise ValueError('Default value of T is not allowed!')
        return d[:,T]

    labels = {
        fr'$U(M) - U(\mathrm{{CBS}})$': ('U', _callable_getitem, 0),
        fr'$K(M) - K(\mathrm{{CBS}})$': ('TrT', _callable_getitem, 0),
        fr'$V(M) - V(\mathrm{{CBS}})$': ('TrV', _callable_getitem, 1),
        fr'{EuuxM}$ - ${EuuxCBS}': ('Eex', Sq.parallel_exchange, 1),
        fr'{EuucM}$ - ${EuucCBS}': ('Ec', Sq.parallel_correlation, 1),
        fr'{EudM}$ - ${EudCBS}': ('Ec', Sq.antiparallel_correlation, 1),
        fr'$-TS(M) + TS(\mathrm{{CBS}})$': ('-TS', _callable_getitem, 0),
        fr'$F(M) - F(\mathrm{{CBS}})$': ('F', _callable_getitem, 0),
    }

    with PdfPages(plotname) as pdf:
        plt.close()

        nrows = 3
        ncols = 2

        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharex=False,
            sharey=False,
        )

        iplt = 0

        for itemp, T in enumerate(thetas):
            ax0, ax1 = axes[itemp]

            ax1.yaxis.tick_right()
            ax1.yaxis.set_label_position("right")

            report = {'label': [], 'CBS': [], 'ERR': []}

            for axi in axes[itemp]:
                axi.text(
                    *ttloc[iplt],
                    theta_label(T),
                    ha='left',
                    va='center',
                )
                iplt += 1

            for i, (label, (key, method, iax)) in enumerate(labels.items()):
                sf = method(T=T)
                x = sf.M.values
                y = sf[key].values

                n = 3
                pwr = -1.0
                xn = x[-n:]**pwr
                yn = y[-n:]

                (_, cbs), (_, cbs_err) = npoint(3, xn, yn, return_error=True)

                report['label'].append(key)
                report['CBS'].append(cbs)
                report['ERR'].append(cbs_err)

                y -= cbs

                if iax == 0:
                    axi = ax0
                else:
                    axi = ax1

                axi.plot(
                    x**pwr,
                    y,
                    ls=lsmap[i],
                    color=cmap[i],
                    label=label if itemp == 0 else '',
                )

                if key == 'TrV':
                    ax0.plot(
                        x**pwr,
                        y,
                        ls=lsmap[i],
                        color=cmap[i],
                    )

            for axi in [ax0, ax1]:
                axi.tick_params(
                    top=True,
                    bottom=True,
                )

                axi.set_xlabel('$M^{-1}$')
                axi.tick_params(axis='x', which='minor', top=True)

            if rs == 1.0 and T == 1.0:
                ax0.set_ylim(-0.02, 0.02)
            ax0.set_ylabel(f'$\Delta$ / {Eh}')
            ax0.ticklabel_format(axis='y', useOffset=False, style='plain')
            ax1.set_ylabel(f'$\Delta$ / {Eh} (potential)')
            ax1.ticklabel_format(axis='y', useOffset=False, style='plain')

            print(f'--- theta = {T} ---')
            report = pd.DataFrame(report)
            print(report.round(10).to_string(
                index=False,
                float_format='% .2E',
            ))


        fig.subplots_adjust(
            wspace=0.2,
            hspace=0.4,
        )

        fig.set_size_inches(
            ncols*3.37/2,
            nrows*3.37*(5**0.5-1)/2,
        )

        fig.legend(
            ncol=2,
            loc='upper center',
            bbox_to_anchor=(0.5, 0.04),
        )

        pdf.savefig(bbox_inches='tight')

    return


def plot_H0_free_energy_loglog(
        d: Data,
        thetas: list[float],
        plotname: str,
    ) -> None:

    F0M = r'$F^{(0)}(M)$'
    F0CBS = r'$F^{(0)}(\mathrm{CBS})$'

    pwr = -1.0

    report = {
        'theta': [],
        'CBS': [],
        'ERR': [],
    }

    cmap = mpl.colormaps['plasma'].resampled(len(thetas) + 3)

    with PdfPages(plotname) as pdf:

        plt.close()

        fig, axes = plt.subplots(
            nrows=1,
            ncols=1,
            sharex=False,
            sharey=False,
        )

        ax0 = axes

        for itheta, theta in enumerate(thetas):
            sf = d[:,theta]

            x = sf.M.values
            y = sf.F.values

            n = 3
            xn = x[-n:]**pwr
            yn = y[-n:]

            (_, cbs), (_, cbs_err) = npoint(
                3,
                xn,
                yn,
                return_error=True,
            )

            report['theta'].append(theta)
            report['CBS'].append(cbs)
            report['ERR'].append(cbs_err)

            y -= cbs

            ax0.plot(
                x,
                y,
                color=cmap(itheta+2),
                label=theta_label(theta),
                zorder=len(thetas) - itheta,
            )

        xline = np.linspace(30, 1100, 100)

        yline = 0.1*xline*np.exp(-0.3168*(xline**(2/3)))
        lline = '$0.1\, M\, e^{-0.3168 M^{2/3}}$'

        ax0.plot(
            xline,
            yline,
            ls='--',
            color='C9',
            label=lline,
            zorder=100,
        )

        ax0.set_xscale('log')
        ax0.set_xlabel(f'$M$')

        ax0.set_yscale('log')
        ax0.set_ylim(bottom=1E-6)
        ax0.set_ylabel(f'{F0M}$ - ${F0CBS} / {Eh}')

        fig.legend(
            ncol=1,
            loc='upper center',
            bbox_to_anchor=(1.13, 0.9),
        )

        pdf.savefig(bbox_inches='tight')

    report = pd.DataFrame(report)
    print(report.to_string(float_format='% .12f'))

    return


def plot_free_energy_corr_loglog(
        d: Data,
        d0: Data,
        thetas: list[float],
        plotname: str,
    ) -> None:

    FcorrM = r'$F_\mathrm{corr.}(M)$'
    FcorrCBS = r'$F_\mathrm{corr.}(\mathrm{CBS})$'

    pwr = -1.0

    report = {'theta': [], 'CBS': [], 'ERR': []}

    cmap = mpl.colormaps['plasma'].resampled(len(thetas) + 3)

    with PdfPages(plotname) as pdf:

        plt.close()

        fig, axes = plt.subplots(
            nrows=1,
            ncols=1,
            sharex=False,
            sharey=False,
        )

        ax0 = axes

        for itheta, theta in enumerate(thetas):
            sf = d[:,theta]
            sf0 = d0[:,theta]

            x = sf.M.values
            y = sf.F.values - sf0.F.values

            n = 3
            xn = x[-n:]**pwr
            yn = y[-n:]

            (_, cbs), (_, cbs_err) = npoint(3, xn, yn, return_error=True)

            report['theta'].append(theta)
            report['CBS'].append(cbs)
            report['ERR'].append(cbs_err)


            y -= cbs

            ax0.plot(
                x,
                y,
                color=cmap(itheta+2),
                label=theta_label(theta),
                zorder=len(thetas) - itheta,
            )

        for pol in [False, True]:
            xline = np.linspace(34, 1100, 100)

            if pol:
                yline = 3*(xline**(-5/3))
                lline = '$3\, M^{-5/3}$'
            else:
                yline = 0.05*(xline**(-1))
                lline = '$0.05\, M^{-1}$'

            ax0.plot(
                xline,
                yline,
                ls=':' if pol else '--',
                color='lime' if pol else 'C9',
                label=lline,
                zorder=100,
            )

        ax0.set_xscale('log')
        ax0.set_xlabel(f'$M$')

        ax0.set_yscale('log')
        ax0.set_ylabel(f'{FcorrM}$ - ${FcorrCBS} / {Eh}')

        fig.legend(
            ncol=1,
            loc='upper center',
            bbox_to_anchor=(1.07, 0.89),
        )

        pdf.savefig(bbox_inches='tight')

    report = pd.DataFrame(report)
    print(report)

    return


def plot_h_analytical_nmax(
        h: Data,
        plotname: str,
        temp: float = 0.1,
    ) -> None:

    M = h.raw.M.values
    nmax = [np.roots([2/3, 1, 1/3, -m]).real.max() for m in M]
    h.raw['nmax'] = nmax

    if temp is None:
        temps = h.theta.copy()
        temps = temps[temps <= 0.5]
    else:
        temps = [temp]

    cmap = mpl.colormaps['plasma_r'].resampled(len(temps) + 2)

    with PdfPages(plotname) as pdf:
        plt.close()

        nrows = 2
        ncols = 1

        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharex=True,
            sharey=False,
        )

        ax0, ax1 = axes

        ax0.axhline(
            y=0,
            ls='--',
            lw=1,
            color='k',
        )

        for itemp, temp in enumerate(temps):
            sf = h[:,temp]
            x = sf.nmax.values
            y0 = sf.U.values
            y1 = sf.F.values

            ax0.plot(
                x,
                y0,
                color=cmap(itemp+1),
                label=temp_label(temp),
            )

            ax1.plot(
                x,
                y1,
                color=cmap(itemp+1),
            )

        for axi in axes:
            axi.tick_params(
                top=True,
                bottom=True,
                left=True,
                right=True,
            )

            axi.set_xscale('log')
            axi.set_xticks([2E0, 5E0, 1E1, 2E1, 5E1])
            axi.set_xticklabels(scinot(['2E0', '5E0', '1E1', '2E1', '5E1']))
            axi.set_xlabel('$n_\mathrm{max}$')
            axi.tick_params(axis='x', which='minor', top=True)

        ax0.set_ylabel(f'$U$ / {Eh}')
        ax0.ticklabel_format(axis='y', useOffset=False)
        ax1.set_ylabel(f'$F$ / {Eh}')
        ax1.ticklabel_format(axis='y', useOffset=False)


        fig.subplots_adjust(
            wspace=0.0,
            hspace=0.0,
        )

        fig.set_size_inches(
            ncols*3.37,
            nrows*3.37*(5**0.5-1)/2,
        )

        fig.legend(
            ncol=2,
            loc='upper center',
            bbox_to_anchor=(0.5, 0.02),
        )

        pdf.savefig(bbox_inches='tight')

    return


def plot_helium_free_energy_loglog(
        he: Data,
        d: Data,
        temps: list[float],
        plotname: str,
        theta: float = 0.9,
    ) -> None:

    pwr = -1.0

    cmap = mpl.colormaps['plasma'].resampled(len(temps) + 3)

    with PdfPages(plotname) as pdf:

        plt.close()

        fig, axes = plt.subplots(
            nrows=1,
            ncols=1,
            sharex=False,
            sharey=False,
        )

        ax0 = axes

        # UEG
        ax1 = ax0.twinx()

        sf = d[:,theta]

        x = sf.M.values
        y = sf.F.values

        n = 3
        xn = x[-n:]**pwr
        yn = y[-n:]

        _, cbs = npoint(3, xn, yn)

        y -= cbs

        ax1.plot(
            x,
            y,
            color='gray',
            label=theta_label(theta) + ' (UEG)',
            zorder=1,
            #alpha=0.01,
        )

        # Helium PBC
        for itemp, temp in enumerate(temps):
            sf = he[:,temp]
            x = sf.M.values
            y = sf.F.values

            n = 3
            xn = x[-n:]**pwr
            yn = y[-n:]

            _, cbs = npoint(3, xn, yn)

            y -= cbs

            ax0.plot(
                x,
                y,
                color=cmap(itemp+2),
                label=temp_label(temp),
                zorder=len(temps) - itemp + 1,
                marker='d',
                markerfacecolor='None',
            )

        if he.M.max() > 500:
            xline = np.linspace(16, 700, 100)
            yline = 2.3*(xline**pwr)
            lline = '$2.3\, M^{-1}$'
        else:
            xline = np.linspace(16, 700, 100)
            yline = 1.0*(xline**pwr)
            lline = '$1.0\, M^{-1}$'

        ax0.plot(
            xline,
            yline,
            ls='--',
            color='C9',
            label=lline,
            zorder=100,
        )

        ax0.set_xscale('log')
        ax0.set_xlabel(f'$M$')

        ax0.set_yscale('log')
        ax0.set_ylabel(f'$F(M) - F(\mathrm{{CBS}})$ / {Eh}')

        ax1.set_yscale('log')
        ax1.set_ylabel(f'$F(M) - F(\mathrm{{CBS}})$ / {Eh} (UEG)')

        fig.legend(
            ncol=2,
            loc='upper center',
            bbox_to_anchor=(0.5, -0.04),
        )

        pdf.savefig(bbox_inches='tight')

    return


def plot_nq(
        d: nqData,
        theta: float,
        plotname: str,
    ) -> None:

    nbasis = d.M.copy()

    cmap = mpl.colormaps['plasma'].resampled(len(nbasis) + 6)

    with PdfPages(plotname) as pdf:
        plt.close()

        nrows = 1
        ncols = 1

        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharex=True,
            sharey=False,
        )

        axi = axes

        for iM, M in enumerate(nbasis):
            sf = d[M,theta]
            x = sf.q.values
            y = sf.n_q.values

            axi.plot(
                x,
                y,
                color=cmap(iM),
                zorder=len(nbasis) - iM,
            )

        axi.tick_params(
            top=True,
            bottom=True,
            left=True,
            right=True,
        )

        # colorbar after subplots_adjust because the other way
        # around locks the position of the colorbar for some reason >:(
        cb = fig.colorbar(
            mpl.cm.ScalarMappable(
                cmap=mpl.colors.ListedColormap([
                    cmap(iM) for iM, _ in enumerate(nbasis)
                ]),
                norm=mpl.colors.BoundaryNorm(
                    np.arange(nbasis.shape[0] + 1),
                    nbasis.shape[0],
                ),
            ),
            ax=axes,
            orientation='vertical',
            label='$M$',
            ticks=nbasis,
            pad=0.03,
        )

        cb.set_ticks(np.arange(nbasis.shape[0]) + 0.5)
        cb.set_ticklabels(nbasis, fontdict={'fontsize': 6})
        cb.minorticks_off()

        axi.set_xlabel('$q$')

        axi.set_yscale('log')
        axi.set_ylabel(r'$n(q)$')

        fig.set_size_inches(
            1.2*ncols*3.37,
            nrows*3.37*(5**0.5-1)/2,
        )

        pdf.savefig(bbox_inches='tight')

    return


def plot_Sq(
        d: SqData,
        theta: float,
        plotname: str,
    ) -> None:

    nbasis = d.M.copy()

    cmap = mpl.colormaps['plasma'].resampled(len(nbasis) + 6)

    with PdfPages(plotname) as pdf:
        plt.close()

        nrows = 3
        ncols = 1

        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharex=True,
            sharey=False,
        )

        ax0, ax1, ax2 = axes

        for iM, M in enumerate(nbasis):
            sf = d[M,theta]
            x = sf.q.values

            ax0.plot(
                x,
                sf.Sx_q.values,
                color=cmap(iM),
                zorder=iM,
            )

            ax0.plot(
                x[-1],
                sf.Sx_q.values[-1],
                color=cmap(iM),
                zorder=iM,
                ls='None',
                marker='d',
                markerfacecolor='None',
            )

            ax1.plot(
                x,
                sf.Sc_q.values,
                color=cmap(iM),
                zorder=iM,
            )

            ax1.plot(
                x[-1],
                sf.Sc_q.values[-1],
                color=cmap(iM),
                zorder=iM,
                ls='None',
                marker='d',
                markerfacecolor='None',
            )

            ax2.plot(
                x,
                sf.Sud_q.values,
                color=cmap(iM),
                zorder=iM,
            )

            ax2.plot(
                x[-1],
                sf.Sud_q.values[-1],
                color=cmap(iM),
                zorder=iM,
                ls='None',
                marker='d',
                markerfacecolor='None',
            )

        for axi in axes:
            axi.axhline(
                y=0,
                lw=1,
                ls=':',
                color='k',
                zorder=100,
            )

            axi.tick_params(
                top=True,
                bottom=True,
                left=True,
                right=True,
            )

        fig.subplots_adjust(
            wspace=0.0,
            hspace=0.0,
        )

        # colorbar after subplots_adjust because the other way
        # around locks the position of the colorbar for some reason >:(
        cb = fig.colorbar(
            mpl.cm.ScalarMappable(
                cmap=mpl.colors.ListedColormap([
                    cmap(iM) for iM, _ in enumerate(nbasis)
                ]),
                norm=mpl.colors.BoundaryNorm(
                    np.arange(nbasis.shape[0] + 1),
                    nbasis.shape[0],
                ),
            ),
            ax=axes,
            orientation='vertical',
            label='$M$',
            ticks=nbasis,
            pad=0.03,
        )

        cb.set_ticks(np.arange(nbasis.shape[0]) + 0.5)
        cb.set_ticklabels(nbasis, fontdict={'fontsize': 6})
        cb.minorticks_off()

        ax0.set_xlabel('')
        ax1.set_xlabel('')
        ax2.set_xlabel('$q$')

        ax0.set_ylabel(r'$S_{\uparrow\uparrow}^\mathrm{ex}(q)$')
        ax1.set_ylabel(r'$S_{\uparrow\uparrow}^\mathrm{c}(q)$')
        ax2.set_ylabel(r'$S_{\uparrow\downarrow}(q)$')

        fig.set_size_inches(
            1.2*ncols*3.37,
            nrows*3.37*(5**0.5-1)/2,
        )

        pdf.savefig(bbox_inches='tight')

    return
