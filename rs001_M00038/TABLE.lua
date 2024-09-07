for ms = 0, 2, 2 do
    fci {
        sys = ueg {
            nel = 2,
            ms = ms,
            sym = 1,
            dim = 3,
            cutoff = 1.0,
            rs = 1.0,
            use_mom_vec = true,
        },
        fci = {
            determinant_table_only = true,
        },
    }
end
