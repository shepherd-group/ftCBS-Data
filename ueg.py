#!/usr/bin/env python

import numpy as np
import pandas as pd

from time import time
from typing import Tuple
from statesum import SumOverStates
from numpy.typing import NDArray as Array
from local_utils import (
    parse_spectrum,
    parse_basis_set,
    parse_determinant_file,
    parse_hamiltonian_file,
    parse_wavefunction_file,
)


class UEG(SumOverStates):
    def __init__(
            self,
            output: str,
            rs: float = 1.0,
            path: str = './',
            get_wfnfile: callable = None,
            get_detfile: callable = None,
            get_hamfile: callable = None,
            zeroth_order_hamiltonian: bool = False,
        ) -> None:
        ''' A wrapper to store somewhat templated UEG FCI data. There are
        some exceptions that can be controlled such as the file names for
        the determinants, wavefunctions, and hamiltonians. For more control
        over other things it is recommended to write a new wrapper.
        '''
        ti = time()
        SumOverStates.__init__(self)

        self.path = path
        self.zeroth_order_hamiltonian = zeroth_order_hamiltonian

        basis, fci, Sfci, Emin = self.parse_hande_output(path, output)
        print('Read in eigenenergies!')

        self.rs = rs
        self.N = 2
        self.M = basis.shape[0]
        self.basis = basis
        self.kvectors = {
            'orbital': np.arange(self.M) + 1,
            'k': self.basis.loc[:, ['kx', 'ky', 'kz']].values,
        }
        self.kvectors['q'] = [np.dot(k, k) for k in self.kvectors['k']]
        self.fci = fci
        self.Sfci = Sfci
        self.Emin = Emin

        self.ndets = {k: fci.shape[0] for k, fci in self.fci.items()}

        self.totdets = np.sum([ndets for k, ndets in self.ndets.items()])

        self.set_file_getters(get_wfnfile, get_detfile, get_hamfile)

        self.wfn = self.parse_wavefunction_outputs()
        print('Read in wavefunctions!')
        self.det = self.parse_determinant_outputs()
        print('Read in determinants!')
        self.ham = self.parse_hamiltonian_outputs()
        print('Read in Hamiltonian!')
        self.kin = self.build_kinetic_matrix()
        print('Built kinetic matrices!')
        self.pot = self.build_potential_matrix()
        print('Built potential matrices!')

        if self.zeroth_order_hamiltonian:
            self.zeroth_order_downsample()
            print('Downsampled for H^(0) approximation!')

        tf = time()
        print(f'Total time to initialize system: {tf-ti:>.8f}')

    def set_file_getters(
            self,
            get_wfnfile: callable,
            get_detfile: callable,
            get_hamfile: callable,
        ) -> None:
        if callable(get_wfnfile):
            self.get_wfnfile = get_wfnfile
        if callable(get_wfnfile):
            self.get_wfnfile = get_wfnfile
        if callable(get_wfnfile):
            self.get_wfnfile = get_wfnfile
        return

    @staticmethod
    def get_detfile(path: str, mS: int, isym: int) -> str:
        return f'{path}/SYM{isym:06}-mS{mS}.det'

    @staticmethod
    def get_wfnfile(path: str, mS: int, isym: int) -> str:
        return f'{path}/SYM{isym:06}-mS{mS}.wfn'

    @staticmethod
    def get_hamfile(path: str, mS: int, isym: int) -> str:
        return f'{path}/SYM{isym:06}-mS{mS}.ham'

    def parse_determinant_outputs(self) -> dict:
        det = {}

        for (mS, isym), fci in self.fci.items():
            ndets = len(fci)
            detfile = self.get_detfile(self.path, mS, isym)
            det[mS,isym] = parse_determinant_file(detfile, ndets, 2)

        return det

    def parse_hamiltonian_outputs(self) -> dict:
        ham = {}

        for (mS, isym), fci in self.fci.items():
            ndets = len(fci)
            hamfile = self.get_hamfile(self.path, mS, isym)
            ham[mS,isym] = parse_hamiltonian_file(hamfile, ndets, fci)

        return ham

    def parse_wavefunction_outputs(self) -> dict:
        wfn = {}

        for (mS, isym), fci in self.fci.items():
            ndets = len(fci)
            wfnfile = self.get_wfnfile(self.path, mS, isym)
            wfn[mS,isym] = parse_wavefunction_file(wfnfile, ndets)

        return wfn

    def build_potential_matrix(self) -> dict:
        pot = {}

        for (mS, isym), ham in self.ham.items():
            pot[mS,isym] = ham - np.diag(self.kin[mS,isym])

        return pot

    def build_kinetic_matrix(self) -> dict:
        kin = {}

        L = (4*np.pi*self.N*self.rs*self.rs*self.rs/3)**(1/3)
        kinetic = 0.5*((2*np.pi/L)**2)*np.array(self.kvectors['q'])

        for (mS, isym), det in self.det.items():
            ndets = len(self.fci[mS,isym])
            v = np.zeros(ndets, dtype=float)

            for idet in range(ndets):
                v[idet] += kinetic[det[idet][0]-1]
                v[idet] += kinetic[det[idet][1]-1]

            kin[mS,isym] = v

        return kin

    def zeroth_order_downsample(self) -> None:
        for (mS, isym) in self.fci:
            self.pot[mS,isym] = np.einsum('ii->i', self.pot[mS,isym])
            self.ham[mS,isym] = np.einsum('ii->i', self.ham[mS,isym])
            self.wfn[mS,isym] = np.eye(self.ham[mS,isym].shape[0])
            self.fci[mS,isym] = self.ham[mS,isym].copy()

        self.Emin = np.min([Ei.min() for _, Ei in self.fci.items()])

        for (mS, isym) in self.Sfci:
            self.Sfci[mS,isym] = self.fci[mS,isym].copy() - self.Emin

        return

    def qvector(self, i: int, a: int) -> Array:
        return self.kvectors['k'][i] - self.kvectors['k'][a]

    def theta_to_beta(self, theta: float) -> float:
        if theta == np.inf:
            return 0.0
        E_F = 0.5*((((9*np.pi/4)**(1/3))/self.rs)**2)
        return 1.0/(theta*E_F)

    @staticmethod
    def parse_hande_output(
                path: str,
                output: str,
            ) -> Tuple[pd.DataFrame, dict, dict, np.longdouble]:
        ''' Collect the relevant information from a HANDE FCI output.
        '''

        basis_set_data = {
            'basis_set': {
                'index': [],
                'kx': [],
                'ky': [],
                'kz': [],
                'ms': [],
                '<i|f|i>': [],
            },
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

        Emin = np.min([fci.min() for (_, _), fci in spectrum.items()])

        del spectrum_data

        shifted = {
            (mS, isym): fci.copy() - Emin
            for (mS, isym), fci in spectrum.items()
        }

        return basis_set, spectrum, shifted, Emin


def main() -> None:

    sys = UEG(
        'FCI.out',
        rs=1.0,
        path='./rs001_M00038',
    )

    theta = 0.5
    beta = sys.theta_to_beta(theta)

    data = sys[beta]

    from local_utils import average_momentum_distribution
    nq = average_momentum_distribution(
        data['n_q'],
        theta,
        beta,
        sys.M,
    )

    from local_utils import total_static_structure_factor
    Sq = total_static_structure_factor(
        data['S_q'],
        theta,
        beta,
        sys.M,
        sys.N,
        sys.rs,
    )

    for k, v in data.items():
        print(k)
        print(v)
        print()

    print('nq')
    print(nq)
    print()

    print('Sq')
    print(Sq)
    print()

    ref = pd.read_csv('./rs001_data/ftfci.csv')
    ref = ref.iloc[(ref.theta.values == theta) & (ref.M.values == 38)]

    for c in ref:
        if c in ['theta', 'M']:
            continue
        assert round(ref[c].values[0], 8) == round(data[c], 8)
        print(f'Found {c} agreement to 1e-8 precision.')

    return


if __name__ == '__main__':
    main()
