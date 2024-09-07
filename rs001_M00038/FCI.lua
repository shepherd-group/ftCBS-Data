function ipad(i)
    local p
    local s = '' .. i .. ''

    for p = string.len(s), 5, 1 do
        s = '0' .. s
    end

    return s
end

for isym = 1, 93, 1 do
    sym_label = 'SYM' .. ipad(isym)

    if isym <= 81 then
        fci {
            sys = ueg {
                nel = 2,
                ms = -2,
                sym = isym,
                dim = 3,
                cutoff = 1.0,
                rs = 1.0,
                use_mom_vec = true,
            },
            fci = {
                write_hamiltonian = true,
                hamiltonian_file = sym_label .. '-mS-2.ham',
                write_nwfns = -1,
                wfn_file = sym_label .. '-mS-2.wfn',
                write_determinants = true,
                determinant_file = sym_label .. '-mS-2.det',
            },
        }
    end

    fci {
        sys = ueg {
            nel = 2,
            ms = 0,
            sym = isym,
            dim = 3,
            cutoff = 1.0,
            rs = 1.0,
            use_mom_vec = true,
        },
        fci = {
            write_hamiltonian = true,
            hamiltonian_file = sym_label .. '-mS0.ham',
            write_nwfns = -1,
            wfn_file = sym_label .. '-mS0.wfn',
            write_determinants = true,
            determinant_file = sym_label .. '-mS0.det',
        },
    }

    if isym <= 81 then
        fci {
            sys = ueg {
                nel = 2,
                ms = 2,
                sym = isym,
                dim = 3,
                cutoff = 1.0,
                rs = 1.0,
                use_mom_vec = true,
            },
            fci = {
                write_hamiltonian = true,
                hamiltonian_file = sym_label .. '-mS2.ham',
                write_nwfns = -1,
                wfn_file = sym_label .. '-mS2.wfn',
                write_determinants = true,
                determinant_file = sym_label .. '-mS2.det',
            },
        }
    end

end
