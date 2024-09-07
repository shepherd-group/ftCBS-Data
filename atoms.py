#!/usr/bin/env python

import numpy as np
import pandas as pd

from time import time
from typing import Tuple
from numpy.typing import NDArray as Array
from statesum import SumOverStatesAtomHighPrecision
from local_utils import (
    parse_spectrum,
    parse_basis_set,
)


class Atom(SumOverStatesAtomHighPrecision):
    def __init__(
            self,
            output: str,
            path: str,
            pbc: bool,
        ) -> None:
        ''' A wrapper to store somewhat templated read_in FCI data.
        See <class 'ueg.UEG'> for more information.
        '''
        ti = time()

        self.path = path

        basis, fci, Sfci, Emin = self.parse_hande_output(path, output, pbc)
        print('Read in eigenenergies!')

        self.M = basis.shape[0]
        self.basis = basis
        self.fci = fci
        self.Sfci = Sfci
        self.Emin = Emin

        self.ndets = {k: fci.shape[0] for k, fci in self.fci.items()}

        self.totdets = np.sum([ndets for k, ndets in self.ndets.items()])

        SumOverStatesAtomHighPrecision.__init__(self)

        tf = time()
        print(f'Total time to initialize system: {tf-ti:>.8f}')

    @staticmethod
    def appender(A: Array, B: Array) -> Array:
        return np.append(A, B)

    def flattened_fci(self) -> Array:
        flatfci = reduce(self.appender, [fci for _, fci in self.fci.items()])
        flatfci = np.sort(flatfci)
        return flatfci

    @staticmethod
    def parse_hande_output(
                path: str,
                output: str,
                pbc: bool,
            ) -> Tuple[pd.DataFrame, dict, dict, np.longdouble]:
        ''' Collect the relevant information from a HANDE FCI output.
        '''
        basis_keys = ['index']

        if pbc:
            basis_keys += ['kx', 'ky', 'kz']
        else:
            basis_keys += ['spatial', 'symmetry', 'sym_index', 'lz']

        basis_keys += ['ms', '<i|f|i>']

        basis_set_data = {
            'basis_set': {bk: [] for bk in basis_keys},
            'found': False,
            'stored': False,
        }

        spectrum_data = {
            'mS': [],
            'symmetry' : [],
            'energy': [],
            'store': False,
        }

        with open(f'{path}/{output}', 'rt') as stream:
            for line in stream:
                parse_basis_set(basis_set_data, line)
                parse_spectrum(spectrum_data, line)

        basis_set = pd.DataFrame(basis_set_data['basis_set'])

        zipper = zip(
            spectrum_data['mS'],
            spectrum_data['symmetry'],
            spectrum_data['energy'],
        )
        spectrum = {
            (mS, isym): np.longdouble(fci)
             for mS, isym, fci in zipper
        }

        empty = {k: v.shape[0] == 0 for k, v in spectrum.items()}

        for (mS, isym), remove in empty.items():
            if remove:
                print(
                    'WARNING: Removing empty symmetry block for: '
                    f'{mS = }, {isym = } !'
                )
                del spectrum[mS,isym]

        Emin = np.min([fci.min() for (_, _), fci in spectrum.items()])

        del spectrum_data

        shifted = {
            (mS, isym): fci.copy() - Emin
            for (mS, isym), fci in spectrum.items()
        }

        return basis_set, spectrum, shifted, Emin


class AnalyticalH(SumOverStatesAtomHighPrecision):
    def __init__(self, nmax: int, mS: int = None) -> None:
        ''' A wrapper to store the analytical eigenvalues for the H atom.
        See <class '__ueg__.UEG'> for more information.
        '''
        ti = time()

        basis, fci, Sfci, Emin = self.generate_spectrum(nmax, mS)

        self.M = basis.Dn.values.sum()
        self.basis = basis
        self.fci = fci
        self.Sfci = Sfci
        self.Emin = Emin
        print('Generated basis and eigenenergies!')

        self.ndets = {k: fci.shape[0] for k, fci in self.fci.items()}

        self.totdets = np.sum([ndets for k, ndets in self.ndets.items()])

        SumOverStatesAtomHighPrecision.__init__(self)

        tf = time()
        print(f'Total time to initialize system: {tf-ti:>.8f}')

    def generate_spectrum(
            self,
            nmax: int,
            mS: int = None,
        ) -> Tuple[pd.DataFrame, dict, dict, np.longdouble]:
        ''' Generate an analytical spectrum for the H atom.
        '''
        Emin = 0.0
        fci = {}
        basis = {
            'index': [],
            'n': [],
            'Dn': [],
            'ms': [],
            '<i|f|i>': [],
        }

        for n in range(1, nmax + 1):
            en = self.En(n)
            dn = self.Dn(n)

            if en < Emin:
                Emin = en

            for spin in [-1, 1]:
                basis['index'].append(len(basis['index']) + 1)
                basis['n'].append(n)
                basis['Dn'].append(dn)
                basis['ms'].append(spin)
                basis['<i|f|i>'].append(en)

                if mS is not None and mS != spin:
                    continue

                if (spin, n) not in fci:
                    fci[spin,n] = []

                fci[spin,n] = np.array([en for _ in range(dn)])

        basis = pd.DataFrame(basis)

        Sfci = {k: v.copy() - Emin for k, v in fci.items()}

        return basis, fci, Sfci, Emin

    @staticmethod
    def En(n: int) -> np.longdouble:
        return np.longdouble(-0.5/(n**2))

    @staticmethod
    def Dn(n: int) -> int:
        return n**2


def main() -> None:

    h = AnalyticalH(nmax=3)

    beta = '0.2'

    data = h[beta]

    for k, v in data.items():
        print(k)
        print(v)
        print()

    ref = pd.read_csv('./atoms_data/ftfci_h_analytical_v2.csv')
    ref = ref.iloc[(ref.beta.values == 0.2) & (ref.M.values == 28)]

    for c in ref:
        if c in ['theta', 'M']:
            continue
        assert round(ref[c].values[0], 8) == round(data[c], 8)
        print(f'Found {c} agreement to 1e-8 precision.')

    he01 = Atom(
        '001-HeATOM-M0018-FCI.out',
        './he_pbc_3x3x3',
        True,
    )

    data = he01[beta]

    for k, v in data.items():
        print(k)
        print(v)
        print()

    ref = pd.read_csv('./atoms_data/ftfci_he_pbc_3x3x3_v2.csv')
    ref = ref.iloc[(ref.beta.values == 0.2) & (ref.M.values == 18)]

    for c in ref:
        if c in ['theta', 'M']:
            continue
        assert round(ref[c].values[0], 8) == round(data[c], 8)
        print(f'Found {c} agreement to 1e-8 precision.')

    he02 = Atom(
        '001-HeATOM-M0018-FCI.out',
        './he_pbc_8x8x8',
        True,
    )

    data = he02[beta]

    for k, v in data.items():
        print(k)
        print(v)
        print()

    ref = pd.read_csv('./atoms_data/ftfci_he_pbc_8x8x8.csv')
    ref = ref.iloc[(ref.beta.values == 0.2) & (ref.M.values == 18)]

    for c in ref:
        if c in ['theta', 'M']:
            continue
        assert round(ref[c].values[0], 8) == round(data[c], 8)
        print(f'Found {c} agreement to 1e-8 precision.')

    return


if __name__ == '__main__':
    main()
