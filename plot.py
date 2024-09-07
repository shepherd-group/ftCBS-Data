#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from analysis import (
    Data,
    nqData,
    SqData,
    plot_internal_energy,
    plot_H0_free_energy_loglog,
    plot_helium_free_energy_loglog,
    plot_all_ueg_quantities_seperated_delta_and_potential_single,
    plot_h_analytical_nmax,
    plot_Sq,
    plot_nq,
    plot_free_energy_corr_loglog,
)


def main() -> None:

    # r_s = 1.0 UEG
    d = Data([
        './rs001_data/ftfci.csv',
        './rs001_data/ftfci_finer_theta.csv',
        './rs001_data/ftfci_finer_theta_v2.csv',
    ])

    d0 = Data([
        './rs001_data/ftH0.csv',
    ])

    Sq = SqData([
        './rs001_data/Sq.csv',
        './rs001_data/Sq_finer_theta.csv',
        './rs001_data/Sq_finer_theta_v2.csv',
    ])

    nq = nqData([
        './rs001_data/nq.csv',
        './rs001_data/nq_finer_theta.csv',
        './rs001_data/nq_finer_theta_v2.csv',
    ])

    # r_s = 10.0 UEG
    d10 = Data([
        './rs010_data/ftfci.csv'
    ])

    Sq10 = SqData(
        ['./rs010_data/Sq.csv'],
        rs=10.0,
    )

    # Atoms
    h = Data('./atoms_data/ftfci_h_analytical_v2.csv')

    he333 = Data('./atoms_data/ftfci_he_pbc_3x3x3_v2.csv')

    he888 = Data('./atoms_data/ftfci_he_pbc_8x8x8.csv')

    # Plots
    thetas = np.array([
        5.0,
        2.0,
        1.0,
        0.5,
        0.2,
        0.1,
        0.05,
        0.02,
    ])

    plot_internal_energy(
        d,
        thetas,
        'Fig001.pdf',
    )

    plot_all_ueg_quantities_seperated_delta_and_potential_single(
        d,
        Sq,
        np.array([0.1, 1.0, 10.0]),
        'Fig002.pdf',
    )

    plot_all_ueg_quantities_seperated_delta_and_potential_single(
        d10,
        Sq10,
        np.array([0.1, 1.0, 10.0]),
        'Fig003.pdf',
        rs=10.0,
    )

    plot_H0_free_energy_loglog(
        d0,
        thetas,
        'Fig004.pdf',
    )

    plot_free_energy_corr_loglog(
        d,
        d0,
        thetas,
        'Fig005.pdf',
    )

    plot_h_analytical_nmax(
        h,
        'Fig006.pdf',
        temp=None, 
    )

    plot_helium_free_energy_loglog(
        he333,
        d,
        he333.theta[he333.theta <= 1.0][::-1],
        'Fig007.pdf',
    )

    plot_helium_free_energy_loglog(
        he888,
        d,
        he888.theta[he888.theta <= 1.0][::-1],
        'Fig008.pdf',
    )

    plot_nq(
        nq,
        2.5,
        'Fig009.pdf',
    )

    plot_Sq(
        Sq,
        2.5,
        'Fig010.pdf',
    )

    return


if __name__ == '__main__':
    main()
