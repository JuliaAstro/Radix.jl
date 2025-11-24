module Radix

"""
    Radix(; trad=1.0, xlum=1.0, lpri=0, lwri=0, r=0.0, r0=r, t=1.0,
    xpx=1000.0, xpx0=xpx, p=0.03, lcdd=1, numrec=1, npass=0, nnmax=1,
    nlimd=0, rmax=1.0, xpxcol=1e16, zeta=0.0, xi=10.0^zeta, lfix=0, ncn2=999,
    cfrac=0.0, emult=0.75, taumax=5.0, xeemin=0.1, nloopctl=0, critf=1.0e-8,
    vturbi=1.0, xcol=0.0, radexp=0.0)
"""
function radix(; trad=1.0, xlum=1.0, lpri=0, lwri=0, r=0.0, r0=r, t=1.0,
    xpx=1000.0, xpx0=xpx, p=0.03, lcdd=1, numrec=1, npass=0, nnmax=1,
    nlimd=0, rmax=1.0, xpxcol=1e16, zeta=0.0, xi=10.0^zeta, lfix=0, ncn2=999,
    cfrac=0.0, emult=0.75, taumax=5.0, xeemin=0.1, nloopctl=0, critf=1.0e-8,
    vturbi=1.0, xcol=0.0, radexp=0.0)

    #  Allocate log files

    #  Display opening message

    #  Read parameter values

    ncn2 = 999
    zremsz = zeros(ncn2)
    ener(epi, ncn2)
    ncn2m = 999
    ener(epim, ncn2m)
    nry = nbinc(13.6, epi, ncn2) + 2
    #  Write nry
    if radexp < -99
        #   open density file
        r, xpx = rnew, dennew
    end

    #  Calculate absolute value of line wavelengths
    llinabs=1
    xstarsetup(lpri, llinabs, abndtbl, ababs, abel, atcredate)

    stotg,  atotg  = zeros(nni), zeros(nni)
    fstotg, fatotg = zeros(3, nni), zeros(3, nni)

    xee = 1.0
    #  nb special for pure fe nlte problem
    xi = 10^zeta

    #  Write the parameter list to FITS file

    #  Initialize
    rdel = 0.0

    #  Main loop
    for kk = 1:npass
        if lwri > 0 || npass > 1
            #  Initialize FITS file for detailed ionic data
        end

        if spectype == "file"
            xlum2 = 1.0
            ispecg()
        elseif spectype == "bbody"
            starf()
        elseif spectype == "pow"
            ispec4()
        elseif spectype == "brems"
            ispec()
        end
        ispecgg()
        ispcg2()

        init()

        lpri = 0
        ldir = (-1)^kk
        #   Xwrite

        #  Print heading for current iteration

        zrems[1,:]  .= zremsz[:]
        zremso[1,:] .= zremsz[:]

        #  Iterate through radius zones
        jkp, rdel0, rdel, ierr = 0, 0.0, 0.0, 0

        while (((kk == 1 && xcol < xpxcol && xee > xeemin && t > tinf*0.99) &&
            numrec > 0) || (kk > 1 && jkp < numrec)) && ierr == 0
            
            #  Save for variable density
            if kk == 1 && lcdd == 1 && jkp == 0
                xpx0, r0 = xpx, r
            end

            jkp += 1
            if jkp > 3999
                # stop
            end

            #  Increment step
            if 1 < kk <= npass
                jk = numrec + 1 - jkp
                lpriu = 0
                rdelo = rdel
                nlimdt = ldir > 0 ? 0 : nlimd
                unsavd()
            end

            if kk == 1
                nlimtd, delr, ectt, jk, lpris = nlimd, 0.0, 1.0, jkp, 0
                if jk > 1
                    step()
                end
            end

            #  Calculate flux
            r19, lpri3 = 1e-19*r, 0
            xi = xlum/r19^2/xpx
            zeta = log10(xi)

            trnfrc()

            if t < 1.02*tinf t = 1.01*tinf end

            #  Calculate temperature, ionization, etc.
            lpri2, lprid = 0, 0
            
            xstarcalc()

            if vturbi > 1e-34
                gsmooth()
            end
            #  Do transfer. assumes comp2 and bremem have called already
            heatt()

            #  Print 

            #  If this is the final pass over all zones
            if kk == npass
                #  Add the abundance information to the output array zrtmp
                #  Print
            end
            #  Write radial profile
            if lwri > 0 && npass > 1
                #  Save
            end
            #  New position etc.
            if radexp < -99
                delr = rnew - r
                if delr < 0 # stop
                end
                r, xpx = rnew, dennew
                #  write
            else
                r += delr
                if lcdd == 1 xpx = xpx0*(r/r0)^radexp end
            end
            rdel += delr
            xcol += xpx*delr

            #  Transfer

            elim1, elimh, lpril2 = 1.0, 1e6, 0
            stpcut()
            trnfrn()

        end

        if kk == 1 numrec = jkp + 1 end

        # Another printout to get the last step
        # Print

        # If final pass over all zones
        if kk == npass
            # Print
        end
        if lwri > 0 || npass > 1
            # Save
        end
        if mod(kk,2) == 0
            # Close FITS files
        end
        lpriu = 0
        # End main iteration loop

    end

    xstarcalc()
    heatt()
    lpri2 = 0
    sppcut()

    # Print

end

end
