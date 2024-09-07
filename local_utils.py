#!/usr/bin/env python

import numpy as np
import pandas as pd

from numpy.typing import NDArray as Array


def parse_determinant_file(output: str, ndets: int, nel: int) -> Array:
    dets = np.zeros((ndets, nel), dtype=int)

    with open(output, 'rt') as stream:
        for line in stream:
            line = line.replace('|', '').replace('>', '').split()
            dets[int(line[0])-1,:] = [int(k) for k in line[1:]]

    return dets


def parse_hamiltonian_file(output: str, ndets: int, eigs: Array) -> Array:
    ham = np.zeros((ndets, ndets), dtype=float)

    # Special case, no file exists since the matrix element is the eigenvalue.
    if ndets == 1:
        ham[0,0] = eigs[0]
        return ham

    with open(output, 'rt') as stream:
        for line in stream:
            i, j, hij = line.split()
            i, j, hij = int(i) - 1, int(j) - 1, float(hij)
            ham[i,j] = hij
            ham[j,i] = hij

    return ham


def parse_wavefunction_file(output: str, ndets: int) -> Array:
    wfns = np.zeros(ndets*ndets, dtype=float)

    # Special case, no file exists since the wave function is unity.
    if ndets == 1:
        wfns[0] = 1.0
        return wfns

    iloc = 0

    with open(output, 'rt') as stream:
        for line in stream:
            Ci = float(line.split()[-1])
            wfns[iloc] = Ci
            iloc += 1

    if iloc != wfns.shape[0]:
        raise ValueError(f'The number of determinants was wrong in {output}!')

    wfns = wfns.reshape((ndets, ndets))

    return wfns


def parse_basis_set(data: dict, line: str) -> dict:
    ''' Parse lines for basis set data.
    '''
    if len(data['basis_set']) == 7:
        # Vacuum basis
        header = 'index  spatial symmetry sym_index lz     ms       <i|f|i>'
    elif len(data['basis_set']) == 6:
        # PBC basis
        header = 'index   k-point            ms       <i|f|i>'
    else:
        raise ValueError('The basis_set dict has too many headers!')

    ld = line.replace(',', '').replace('(', '').replace(')', '').split()

    if data['stored']:
        pass
    elif not data['found'] and header in line:
        data['found'] = True
    elif data['found'] and len(ld) == 7:
        # Vacuum basis
        data['basis_set']['index'].append(int(ld[0]))
        data['basis_set']['spatial'].append(int(ld[1]))
        data['basis_set']['symmetry'].append(int(ld[2]))
        data['basis_set']['sym_index'].append(int(ld[3]))
        data['basis_set']['lz'].append(int(ld[4]))
        data['basis_set']['ms'].append(int(ld[5]))
        data['basis_set']['<i|f|i>'].append(float(ld[6]))
    elif data['found'] and len(ld) == 6:
        # PBC basis
        data['basis_set']['index'].append(int(ld[0]))
        data['basis_set']['kx'].append(int(ld[1]))
        data['basis_set']['ky'].append(int(ld[2]))
        data['basis_set']['kz'].append(int(ld[3]))
        data['basis_set']['ms'].append(int(ld[4]))
        data['basis_set']['<i|f|i>'].append(float(ld[5]))
    elif data['found']:
        data['stored'] = True

    return data


def parse_spectrum(data: dict, line: str) -> dict:
    ''' Parse the spectrum and related information from a HANDE
    FCI output.
    '''
    header = 'State     Energy'

    ld = line.replace(',', '').split()

    if '"Ms":' in line:
        data['mS'].append(int(ld[-1]))
        data['energy'].append([])
    elif '"symmetry":' in line:
        data['symmetry'].append(int(ld[-1]))
    elif not data['store'] and header in line:
        data['store'] = True
    elif data['store'] and len(ld) == 2:
        data['energy'][-1].append(float(ld[-1]))
    elif data['store']:
        data['store'] = False

    return data


def average_momentum_distribution(
        nq: pd.DataFrame,
        theta: float,
        beta: float,
        M: int,
    ) -> pd.DataFrame:
    
    u = nq.loc[nq.index % 2 == 0].reset_index(drop=True)
    d = nq.loc[nq.index % 2 == 1].reset_index(drop=True)

    ave = pd.DataFrame({
        'q': u.q.values.round(6),
        'n_q_alpha': u.n_q.values,
        'n_q_beta': d.n_q.values,
        'n_q': u.n_q.values + d.n_q.values,
    })

    cnt = ave.groupby('q').count().reset_index(drop=False)
    ave = ave.groupby('q').mean().reset_index(drop=False)

    ave['D_q'] = cnt['n_q'].values.flatten()
    ave['beta'] = [beta for _ in range(ave.shape[0])]
    ave['theta'] = [theta for _ in range(ave.shape[0])]
    ave['M'] = [M for _ in range(ave.shape[0])]
    ave['T_q'] = ave['q']*ave['q']*ave['D_q']*ave['n_q']/2.0

    return ave


def total_static_structure_factor(
        Sq: pd.DataFrame,
        theta: float,
        beta: float,
        M: int,
        N: int,
        rs: float,
    ) -> pd.DataFrame:

    cnt = Sq.groupby('qdot').count().reset_index(drop=False)
    tot = Sq.groupby('qdot').mean().reset_index(drop=False)

    tot['Suu_q'] += 1.0
    tot['S_q'] = tot['Suu_q'] + tot['Sud_q']

    L = rs*(((4.0*np.pi*N)/3)**(1/3))
    tot['q'] = 2.0*np.pi*(tot['qdot'].values**0.5)/L
    tot['D_q'] = cnt['Suu_q']
    tot['beta'] = [beta for _ in range(tot.shape[0])]
    tot['theta'] = [theta for _ in range(tot.shape[0])]
    tot['M'] = [M for _ in range(tot.shape[0])]

    return tot
