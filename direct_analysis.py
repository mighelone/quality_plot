"""
Plot results of the direct chemistry solution for the 3 reference
cases
"""

import matplotlib.pyplot as plt
import numpy as np
import pyFLUT
# import pyFLUT.ulf as ulf
import coal_calculation
import glob
import os
import yaml
import cantera

# plt.style.use(['mystyle', 'mystyle-paper'])
cycler = plt.rcParams['axes.prop_cycle']
colors = [i['color'] for i in cycler]

arr_th = 0.0025


def plot_arrow(ax, pos=(0, 0), color=None, direction='right'):
    # axt.arrow(X[index], hMean[index], arr_th, 0,
    #          color=colors[2], length_includes_head=True)
    if direction == 'left':
        shift = -arr_th
    else:
        shift = arr_th
    ax.annotate("", xytext=pos, xy=(pos[0] + shift, pos[1]),
                arrowprops=dict(arrowstyle='-|>', color=color, linewidth=1))


def read_coal1D(dir, pv, flut):
    """
    Read direct solutions
    """
    files = glob.glob(os.path.join(dir, 'coal1D_alpha*_U*.ulf'))
    results = [coal_calculation.Coal1D(f) for f in files]
    results.sort()
    print('Read results in {}'.format(dir))
    for r in results:
        print(r)
        r.calc_progress_variable(pv)
        r['alpha0'] = r['M_p'] / r['rho']
        r['rhoWPV'] = r['WPV'] * r['rho']
        # r['CH2O*100'] = r['CH2O'] * 100
    calc_Hnorm(results=results, flut=flut)
    return results


def calc_zstoich(yml_file):
    yml_dict = read_yml(yml_file)
    fuel = yml_dict['volatiles']
    fuel['T'] = 300
    ox = yml_dict['oxidizer']
    ox['T'] = 300
    path = os.path.dirname(yml_file)
    mech = os.path.join(path, yml_dict['mechanism'])
    gas = cantera.Solution(mech)
    ox = pyFLUT.utilities.fix_composition(ox, gas)
    fuel = pyFLUT.utilities.fix_composition(fuel, gas)
    return pyFLUT.utilities.calc_zstoich(fuel=fuel, oxidizer=ox, gas=gas)


def read_yml(yml_input):
    with open(yml_file, 'r') as f:
        yml_dict = yaml.load(f)
    return yml_dict


def calc_Hnorm(results, flut):
    H_o = np.array([flut.extract_values('hMean', Z=0, Hnorm=h)
                    for h in (0, 1)])
    H_f = np.array([flut.extract_values('hMean', Z=1, Hnorm=h)
                    for h in (0, 1)])
    [r.calc_Hnorm(H_o=H_o, H_f=H_f) for r in results]


def clean_axes(ax):
    ax.set_xlabel('')
    ax.set_ylabel('')


def plot_temperature(res, n=1, ax=None):
    if not ax:
        _, ax = plt.subplots()
    res.plotdata('X', 'T', ax=ax, legend=False)
    X = res['X']
    X_stg = res.stagnation_X()
    cond = X < X_stg
    ax.plot(X[cond], res['T_p'][cond])
    ax.locator_params(nbins=6)
    clean_axes(ax)
    X_stg = plot_stagnation(res, ax)
    if n == 0:
        ax.legend(['$T_g$', '$T_p$'], loc='best')
        ax.set_ylabel('Temperature, K')

    # axt = ax.twinx()
    # hMean = res['hMean'] * 1e-6
    # axt.plot(X, hMean, label='$h_g$', color=colors[2])
    # index = hMean.argmin()

    # axt.arrow(X[index], hMean[index], arr_th, 0,
    #          color=colors[2], length_includes_head=True)
    # axt.annotate("", xytext=(X[index], hMean[index]), xy=(X[index] + arr_th, hMean[index]),
    # arrowprops=dict(arrowstyle='-|>', color=colors[2], linewidth=2))
    # plot_arrow(axt, pos=(X[index], hMean[index]), color=colors[2])

    # Z = res['Z']
    # hMean_max = 1e-6 * (res.H_o[1] * (1 - Z) + res.H_f[1] * Z)
    # axt.plot(X, hMean_max, color=colors[2],
    #         linestyle='dashed', label='$h_g^{max}$')
    # axt.set_ylim([-0.5, 0.5])
    # axt.locator_params(axis='y', nbins=3)
    # if n == 2:
    #    axt.legend(loc='best')
    #    axt.set_ylabel('Total enthalpy, MJ/kg')
    # elif n != 2:
    #    axt.set_yticklabels([])
    return ax


def plot_stagnation(res, ax):
    """
    Plot stagnation line and return its value
    """
    X_stg = res.stagnation_X()
    ax.axvline(X_stg, linestyle='dashed', color='red')
    return X_stg


def plot_mixture_fraction(res, n=1, ax=None):
    """
    Plot mixture fraction and SDR
    """
    ax = set_ax(ax)
    res.plotdata('X', 'Z', legend=False, ax=ax,
                 label='Z', color=colors[0])
    ax.axhline(Z_st, linestyle='dashed', color='red')
    if n == 0:
        ax.annotate("$Z_{st}$", xytext=(0.005, Z_st + 0.03),
                    xy=(0.001, Z_st), arrowprops=dict(arrowstyle='-|>',
                                                      linewidth=1,
                                                      color='black'))
    plot_stagnation(res, ax)

    # shift the arrow 40 points left to the maximum
    index = res['Z'].argmax() - 40
    plot_arrow(ax, pos=(res['X'][index], res['Z'][index]), color=colors[0],
               direction='left')
    clean_axes(ax)
    axt = ax.twinx()
    res.plotdata('X', 'chi', legend=False, ax=axt,
                 label='$\chi$', color=colors[1])
    axt.set_ylim([0, 25])

    if n != 0:
        index = res['chi'].argmax()
        plot_arrow(axt, pos=(res['X'][index], res[
                   'chi'][index]), color=colors[1], direction='right')
    clean_axes(axt)
    if n != 2:
        axt.set_yticklabels([])
    if n == 0:
        ax.set_ylabel('Z')
        ax.legend(loc='upper left')
    elif n == 2:
        axt.legend(loc='upper right')
        axt.set_ylabel('$\chi$, 1/s')
    return ax


def plot_coal_status(res, n=1, ax=None):
    ax = set_ax(ax)
    X_st = res.stagnation_X()
    X = res['X']
    cond = X < X_st
    ax.plot(X[cond], res['alpha0'][cond],
            label='$\\alpha$', color=colors[0])
    clean_axes(ax)

    ax.plot(X[cond], res['m_volatiles'][cond] / res['m_volatiles'][0],
            label='$Y_v/Y_{v,0}$', color=colors[1])
    ax.set_ylim([0, 1.5])
    ax.locator_params(axis='y', nbins=3)

    X_st = plot_stagnation(res, ax)
    if n == 0:
        ax.set_ylabel('')
        ax.legend(loc='upper left')
        ax.annotate("Stagnation\nplane", xytext=(X_st + 0.001, 0.4),
                    xy=(X_st, 0.2), arrowprops=dict(arrowstyle='-|>',
                                                    linewidth=1,
                                                    color='black'))
    return ax


def set_ax(ax):
    if not ax:
        ax = plt.subplots()
    return ax


def plot_species(res, n=1, ax=None):
    ax = set_ax(ax)
    for sp in species:
        res.plotdata('X', sp, legend=False, ax=ax)
    clean_axes(ax)
    if n == 0:
        ax.set_ylabel('Mass fractions')
    if n == 2:
        ax.legend(species_label, loc='upper left',
                  bbox_to_anchor=(1.02, 1))
    plot_stagnation(res, ax)
    return ax


def plot_cema(res, n=1, ax=None):
    cem = '$\gamma_{e}$'
    ax = set_ax(ax)
    res.plotdata('X', 'CEMA', label=cem, ax=ax,
                 color=colors[0], legend=False)
    ax.axhline(0, color=colors[0], linestyle='dashed')
    cema = res['CEMA']
    X = res['X']
    X_cemapeak = X[cema.argmax()]
    X_st = plot_stagnation(res, ax)
    index = (X > X_cemapeak) & (X < X_st)
    X = X[index]
    cema = cema[index]
    X_cross = X[cema >= 0][-1]
    ax.plot([X_cross], [0], marker='o', color=colors[0],
            markeredgecolor=colors[0], markersize=8)

    index = res['CEMA'].argmax()
    pos = (res['X'][index], res['CEMA'][index])
    plot_arrow(ax, pos=pos, color=colors[0], direction='left')

    if n == 0:
        ax.annotate("zero-crossing", xytext=(X_cross + 0.005, 5),
                    xy=(X_cross, 0), arrowprops=dict(arrowstyle='-|>',
                                                     linewidth=1,
                                                     color='black'))

    axt = ax.twinx()
    hrr = -res['HRR'] * 1e-6
    axt.plot(res['X'], hrr, label='$\Omega_T$', color=colors[1])
    clean_axes(ax)
    clean_axes(axt)
    ax.locator_params(axis='y', nbins=4)
    axt.locator_params(axis='y', nbins=4)

    index = hrr.argmax()
    pos = (res['X'][index], hrr[index])
    plot_arrow(axt, pos=pos, color=colors[1], direction='right')

    if n != 2:
        axt.set_yticklabels([])
    if n == 0:
        ax.set_ylabel(cem + ' , 1/s')
        # ax.legend(loc='upper right')
    elif n == 2:
        axt.set_ylabel('$\Omega_T$, W')
        # axt.legend(loc='upper left')
    return ax

results_dir = 'results'
dc_dir = os.path.join(results_dir, 'FTC')
fpv_dir = os.path.join(results_dir, 'DFLUT')
flut_dir = fpv_dir
species = ['C6H6', 'CO', 'CO2', 'H2O', 'O2']
species_label = ['$C_6 H_6$', '$CO$', '$CO_2$',
                 '$H_2O$', '$O_2$']

yml_file = os.path.join(flut_dir, 'input.yml')
Z_st = calc_zstoich(yml_file=yml_file)

print('Z_st={:4.3f}'.format(Z_st))
pv = {'CO': 1, 'CO2': 1}

L = 11  # inches
H = 3   # inches

titles = ['$\\alpha_{in}$=0.2 - $U_{in}$=0.2 m/s',
          '$\\alpha_{in}$=0.5 - $U_{in}$=0.3 m/s',
          '$\\alpha_{in}$=0.5 - $U_{in}$=0.7 m/s']


if __name__ == '__main__':
    main()


def main():
    flut = ulf.read_ulf_bin(os.path.join(flut_dir, 'coalFLUT-cc.h5'))
    coal1D = read_coal1D(dc_dir, pv=pv, flut=flut)
    plot_fcts = [plot_coal_status,
                 plot_temperature,
                 plot_mixture_fraction,
                 plot_species,
                 plot_cema
                 ]

    n_rows = len(plot_fcts)
    fig, axes = plt.subplots(ncols=len(coal1D), nrows=n_rows, sharey='row',
                             sharex='col', figsize=(L, H * n_rows))

    for n_row, (fct, axes_row) in enumerate(zip(plot_fcts, axes)):
        for n_col, (c, ax) in enumerate(zip(coal1D, axes_row)):
            fct(c, n_col, ax=ax)
            ax.set_xlim(xmax=0.02)

    for i, t in enumerate(titles):
        axes[0][i].set_title(t)

    for ax in axes[-1]:
        ax.set_xlabel('x, m')

    # plt.show()
    fig.savefig('direct_analysis.pdf', dpi=80, bbox_inches='tight')
