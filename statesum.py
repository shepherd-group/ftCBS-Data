#!/usr/bin/env python

import numpy as np
import pandas as pd
import mpmath as mp

from time import time
from typing import Tuple
from functools import reduce
from numpy.typing import NDArray as Array


class SumOverStates:
    def __init__(self) -> None:
        pass

    @staticmethod
    def get_rho(ndets: int, z: Array, wfn: Array, zeroth_order: bool) -> Array:
        if zeroth_order:
            rho = np.diag(z)
        else:
            rho = np.zeros((ndets, ndets), dtype=np.longdouble)

            for i in range(ndets):
                rho += z[i]*np.outer(np.transpose(wfn[i]),  wfn[i])

        return rho

    @staticmethod
    def update_n_q(iocc: int, rho_ii: np.longdouble, traces: dict) -> dict:
        traces['n_q']['n_q'][iocc] += rho_ii
        return traces

    @staticmethod
    def update_Sq_keys(q: Array, traces: dict) -> dict:
        dx, dy, dz = q
        if (dx, dy, dz) not in traces['S_q_keys']:
            traces['S_q_keys'].append((dx, dy, dz))
        return traces

    @staticmethod
    def update_Suu_q(q: Array, rho_ij: np.longdouble, traces: dict) -> dict:
        dx, dy, dz = q
        if (dx, dy, dz) not in traces['Suu_q']:
            traces['Suu_q'][dx,dy,dz] = 0.0
        traces['Suu_q'][dx,dy,dz] += rho_ij
        return traces

    @staticmethod
    def update_Sx_q(q: Array, rho_ij: np.longdouble, traces: dict) -> dict:
        dx, dy, dz = q
        if (dx, dy, dz) not in traces['Sx_q']:
            traces['Sx_q'][dx,dy,dz] = 0.0
        traces['Sx_q'][dx,dy,dz] += rho_ij
        return traces

    @staticmethod
    def update_Sc_q(q: Array, rho_ij: np.longdouble, traces: dict) -> dict:
        dx, dy, dz = q
        if (dx, dy, dz) not in traces['Sc_q']:
            traces['Sc_q'][dx,dy,dz] = 0.0
        traces['Sc_q'][dx,dy,dz] += rho_ij
        return traces

    @staticmethod
    def update_Sud_q(q: Array, rho_ij: np.longdouble, traces: dict) -> dict:
        dx, dy, dz = q
        if (dx, dy, dz) not in traces['Sud_q']:
            traces['Sud_q'][dx,dy,dz] = 0.0
        traces['Sud_q'][dx,dy,dz] += rho_ij
        return traces

    def update_nq_Sq(
            self,
            ndets: int,
            traces: dict,
            det: Array,
            rho: Array,
        ) -> dict:
        ''' Loop through the density matrix elements and update the
        corresponding momentum distribution n(q) and static structure
        factor S(q) elements. S(q) is broken into the up-up (eq. down-down)
        and up-down (eq. down-up) contributions.
        '''
        for idet in range(ndets):
            rho_ii = rho[idet,idet]

            # Update n(q) and S(q) from diagonal elements.
            if abs(rho_ii) >= 1E-12:
                for i in range(self.N):
                    iocc = det[idet,i] - 1

                    self.update_n_q(iocc, rho_ii, traces)

                    for j in range(i+1, self.N):
                        jocc = det[idet,j] - 1

                        if (iocc % 2) == (jocc % 2):
                            q = self.qvector(iocc, jocc)
                            self.update_Sq_keys(q, traces)
                            self.update_Suu_q(q, -2.0*rho_ii, traces)
                            self.update_Sx_q(q, -2.0*rho_ii, traces)

            if self.zeroth_order_hamiltonian:
                continue

            # Update S(q) from off-diagonal elements.
            for jdet in range(ndets):
                rho_ij = rho[idet,jdet]

                if idet == jdet or abs(rho_ij) < 1E-12:
                    continue

                if any(orb in det[jdet] for orb in det[idet]):
                    raise RuntimeError('Single excitation encountered!')

                # Warning, there is a permutation factor implicitly here, but
                # we only have double excitations so it's always positive unity.
                i, j = det[idet] - 1
                a, b = det[jdet] - 1

                if (i % 2) == (j % 2):
                    q = self.qvector(i, a)
                    self.update_Sq_keys(q, traces)
                    self.update_Suu_q(q, 2.0*rho_ij, traces)
                    self.update_Sc_q(q, 2.0*rho_ij, traces)

                    q = self.qvector(i, b)
                    self.update_Sq_keys(q, traces)
                    self.update_Suu_q(q, -2.0*rho_ij, traces)
                    self.update_Sc_q(q, -2.0*rho_ij, traces)
                elif (i % 2) == (a % 2):
                    q = self.qvector(i, a)
                    self.update_Sq_keys(q, traces)
                    self.update_Sud_q(q, 2.0*rho_ij, traces)
                elif (i % 2) == (b % 2):
                    q = self.qvector(i, b)
                    self.update_Sq_keys(q, traces)
                    self.update_Sud_q(q, -2.0*rho_ij, traces)

        return traces

    def __getitem__(self, beta: float) -> dict:
        ti = time()

        b = np.longdouble(beta)
        t = 1.0/b if b > 0.0 else np.inf

        if b == 0.0:
            raise ValueError('Infinite temperature logic not in place!')

        zi = {
            (mS, isym): np.exp(-b*Sfci)
            for (mS, isym), Sfci in self.Sfci.items()
        }

        zsum = np.sum([z.sum() for _, z in zi.items()])

        zi = {k: z/zsum for k, z in zi.items()}

        S = np.sum([-np.log(np.power(z, z)).sum() for _, z in zi.items()])

        U = np.sum([np.dot(z, self.fci[k]) for k, z in zi.items()])

        F = U - t*S if b > 0.0 else -np.inf

        Z = np.exp(-b*F)

        traces = {
            'rho': 0.0,
            'rho_H': 0.0,
            'rho_T': 0.0,
            'rho_V': 0.0,
            'n_q': {
                'orbital': np.arange(self.M) + 1,
                'q': (2*self.basis['<i|f|i>'].values)**0.5,
                'n_q': [0.0 for _ in range(self.M)],
            },
            'S_q_keys': [],
            'Suu_q': {},
            'Sud_q': {},
            'Sx_q': {},
            'Sc_q': {},
        }

        for (mS, isym), z in zi.items():
            det = self.det[mS,isym]
            wfn = self.wfn[mS,isym]
            ndets = self.ndets[mS,isym]

            rho = self.get_rho(ndets, z, wfn, self.zeroth_order_hamiltonian)

            traces = self.update_nq_Sq(ndets, traces, det, rho)

            if self.zeroth_order_hamiltonian:
                subscripts = 'ii,i->'
            else:
                subscripts = 'ij,ji->'

            traces['rho'] += rho.trace()
            traces['rho_H'] += np.einsum(subscripts, rho, self.ham[mS,isym])
            traces['rho_T'] += np.einsum('ii,i->', rho, self.kin[mS,isym])
            traces['rho_V'] += np.einsum(subscripts, rho, self.pot[mS,isym])

        Sq = {
            'qdot': [],
            'S_q': [],
            'Suu_q': [],
            'Sud_q': [],
            'Sx_q': [],
            'Sc_q': [],
        }

        for key in traces['S_q_keys']:
            dx, dy, dz = key

            Sq['qdot'].append(np.dot([dx, dy, dz], [dx, dy, dz]))

            if key in traces['Suu_q']:
                Sq['Suu_q'].append(traces['Suu_q'][key]/self.N)
            else:
                Sq['Suu_q'].append(0.0)

            if key in traces['Sud_q']:
                Sq['Sud_q'].append(traces['Sud_q'][key]/self.N)
            else:
                Sq['Sud_q'].append(0.0)

            if key in traces['Sx_q']:
                Sq['Sx_q'].append(traces['Sx_q'][key]/self.N)
            else:
                Sq['Sx_q'].append(0.0)

            if key in traces['Sc_q']:
                Sq['Sc_q'].append(traces['Sc_q'][key]/self.N)
            else:
                Sq['Sc_q'].append(0.0)

        Sq['S_q'] = np.array(Sq['Suu_q']) + np.array(Sq['Sud_q'])

        quantities = {
            'beta': beta,
            'Z': Z,
            'U': U,
            'S': S,
            'F': F,
            'TrU': traces['rho_H']/traces['rho'],
            'TrT': traces['rho_T']/traces['rho'],
            'TrV': traces['rho_V']/traces['rho'],
            'n_q': pd.DataFrame(traces['n_q']),
            'S_q': pd.DataFrame(Sq),
        }

        tf = time()

        quantities['time'] = tf - ti

        return quantities


class SumOverStatesAtomHighPrecision:
    def __init__(self, dps: int = 1000) -> None:
        mp.mp.dps = dps

        self.E = []

        print(f'Storing the spectrum data with precision:')
        print(mp.mp)
        for (mS, isym), fci in self.fci.items():
            for Ei in fci:
                self.E.append(mp.mpf(str(Ei)))

    @staticmethod
    def exp(C: mp.mpf, E: list) -> list:
        ''' Evaluate exp(C * Ei) for a list of Ei.
        '''
        return [mp.exp(C*Ei) for Ei in E]

    @staticmethod
    def sum(A: list) -> mp.mpf:
        return reduce(lambda x, y: x + y, A)

    @staticmethod
    def mult(A: list, B: list) -> list:
        return [Ai*Bi for Ai, Bi in zip(A, B)]

    def __getitem__(self, beta: str) -> dict:
        ti = time()

        if not isinstance(beta, str):
            raise ValueError('beta must be type string!')

        b = mp.mpf(beta)
        t = mp.mpf('1.0')/b

        # Z = Tr[exp(-\beta H)]
        e = self.exp(-b, self.E)
        Z = self.sum(e)

        # F = -1/\beta ln(Z)
        F = -t*mp.ln(Z)

        # U = -d/d\beta ln(Z)
        #   = -Z'/Z
        # where
        # Z' = -Tr[H exp(-\beta H)]
        Zprime = -self.sum(self.mult(e, self.E))
        U = -Zprime/Z

        # S = -Tr[w ln(w)]
        # where
        # w = exp(-\beta H)/Z
        w = [v/Z for v in e]
        lnw = [mp.ln(v) for v in w]
        S = -self.sum(self.mult(w, lnw))

        data = {
            'beta': float(b),
            'temp': float(t),
            'Z': float(Z),
            'F': float(F),
            'U': float(U),
            'S': float(S),
            'time': time() - ti,
        }

        return data
