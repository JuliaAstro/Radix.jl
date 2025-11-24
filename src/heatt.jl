"""
heatt

    Calculates heating, cooling, and performs radiative transfer.
"""
function heatt()
    
    xnx = xpx*xee
    r19 = 1e-19*r
    if lpri > 1 lpri=2 end
    fpr2 = 12.56*r19^2
    # if lpri > 1 write end
    # if lpri > 1 write end
    numcon = ncn2m

    # comment the following out to implment scattering
    clbrems, lskp, tmp2 = 0.0, 1, 0.0

    for k1 = 1:numcon

        for jkk = 1:nlsvn
        end

        # calculate rrc luminosities
        # first look for element data (jk is element index)
        # if lpri >= 2 write end
        klel = 11
        mlel = derivedpointer.npfirst(klel)
        jk, jkk = 0, 0

        while mlel != 0
            jk += 1
            drd()
            if nidt > 0
                mllel = masterdata.idat1[np1i+nidt+1]
                xeltp = masterdata.rdat1[np1r]
                xeltp = abel[mllel]
                nnz = masterdata.idat1[np1i]
                # ignore if abundance is small
                if xeltp < 1e-10
                    jkk += nnz
                else
                    # step through ions (jkk is ion index)
                    klion = 12
                    mlion = derivedpointer.nfpirst[lkion]
                    jkk, kl = 0, 0
                    while mlion != 0 && lk < nnz
                        jkk += 1
                        # retrieve ion name from kdati
                        mlm = mlion
                        drd()
                        
                    end
                end
            end
        end
    end



end