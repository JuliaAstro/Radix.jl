"""
calczone()

    Calculates all relevant quantities for one spatial zone. First calls dsec for
    thermal equilibrium calculation (if needed). Then calls calc_hmc_all and
    calc_emis_all.
"""
#xstarcalc
function calczone()

    bremsmap()

    lprisv, lpri = lpri, 0
    if nlimdt != 0
        dsec()
    end
    calc_hmc_all()
    calc_emisab_all()
    calc_emis_all()
    nry = nbinc(13.6, epi, ncn2) + 2
    lpri = lprisv
end