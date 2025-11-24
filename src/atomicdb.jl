####    Atomic Database    ####

using FITSFiles

function load(io::IO)

    atomdb = fits(io)

    ptr  = reshape(atomdb[2].data[:pointers][1], 10, :)
    rdat = atomdb[3].data[:reals][1]
    idat = atomdb[4].data[:integers][1]
    sdat = atomdb[5].data[:char][1]

    N = length(ptr[1,:])
    rates = Array{Union{Atom, Ion, AbstractRate, Missing, Nothing}}(undef, N)
    for j=1:N
        if ptr[2,j] == 0
            break
        end
        rate = ratemap[ptr[2,j]]
        if rate in (Missing, Nothing)
            println("Rate $(ptr[2,j]) not implemented.")
            continue
        end
        type  = ptr[3,j]
        label = String(sdat[ptr[10,j]:ptr[10,j]+ptr[7,j]-1])
        ivec  = idat[ptr[ 9,j]:ptr[ 9,j]+ptr[6,j]-1]
        rvec  = rdat[ptr[ 8,j]:ptr[ 8,j]+ptr[5,j]-1]
        # if ptr[2,j] in (2, 9, 10, 13, 14)
        #     println("$(ptr[2,j]) => $rate: $ivec, $rvec, $label")
        # end
        rates[j] = rate(type, label, ivec, rvec)
    end
    rates
end
