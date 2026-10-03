# Tables of the superlevel rates (ucalc types 71 and 77): a log₁₀ quantity on a (log₁₀ T, log₁₀ n) grid,
# with the temperature running fastest in the stored table.

const super_T_margin = 1.0            # log₁₀ T may exceed the table by this

"""
    interpolate_super_table(logne, logT, table, T, den; guard=0.0)

Linear interpolation of `table[it, in]` in log₁₀ T and log₁₀ n at the temperature `T` (K) and the
density `den`, as XSTAR's `calt71` and `calt77` do: the density is limited to the table, the
temperature to one dex outside it (the end segments are extrapolated), and `guard` is added to the
grid spacings (`calt77` uses 1e-36).
"""
function interpolate_super_table(logne, logT, table, T, den; guard=0.0)
    nden, ntem = length(logne), length(logT)
    rne = min(log10(den), logne[nden])
    rte = clamp(log10(T), logT[1] - super_T_margin, logT[ntem] + super_T_margin)
    in = 1
    if rne > logne[1]
        in = 0
        while true
            in += 1
            in < nden && rne >= logne[in + 1] && continue
            break
        end
    end
    it = 1
    if rte >= logT[1]
        it = 0
        while true
            it += 1
            if it >= ntem
                it = ntem - 1
            elseif rte >= logT[it + 1]
                continue
            end
            break
        end
    end
    at(j) = table[it, j] + (table[it + 1, j] - table[it, j])/(logT[it + 1] - logT[it] + guard)*(rte - logT[it])
    rec = at(in)
    in < nden && (rec += (at(in + 1) - rec)/(logne[in + 1] - logne[in] + guard)*(rne - logne[in]))   # (at the last density rne equals logne[in])
    rec
end
