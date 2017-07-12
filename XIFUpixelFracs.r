    XIFUpixelFracs <- function(ctrate, deltat, t0, t1, t2, tolsep){
        #
        # Calculate (Quality) Fractions of photons for XIFU pixels, using photon simulation AND
        # Poisson statistics
        #
        # Poissonian Distribution of Photons
        #   tp = Tprevious
        #   tn = Tnext
        #   L : event with "Low" Quality
        #   M : event with "Medium" Quality
        #   H : event with "High" Quality
        #   I : event "Invalid"
        #   P2 : secondary events piled-up (same arrival time) (preceeded by close pulses)
        #   P1: first pulse in pile-up (followed by close pulses but no pulse close before)
        #   P : Pile-up = P1+P2
        #===================================================
        #            | tolsep>tp  tolsep<tp<t0   tp>t0     =
        # -------------------------------------------------=
        #   tn<t1    |      P2           I        L        =
        # t1<tn<t2   |      P2           I        M        =
        #   tn>t2    |      P2           I        H        =
        #  tn<tolsep |      P2           P1       P1       =
        #===================================================
        #
        # ctrate (in) : count rate (phs/s)
        # deltat (in) : time interval for simulations
        # t0, t1, t2, tolsep as in table above
        # fractions (out): list with H(igh)Q, M(edium)Q, L(ow)Q, I(nvalid)Q, P(i)L(eup) quality
        
        N.mean.photons <- deltat*ctrate
        stopifnot(N.mean.photons>0)
        N.real.photons <- rpois(1,N.mean.photons)
        #cat("psfrac=",psf, "N.mean.photons=", N.mean.photons,"\n")
        points <- runif(N.real.photons,0.,deltat)
        points <- sort(points)
        
        # Calculate fractions by photon counting (points in plots)
        # ----------------------------------------------------------
        Frac.HQ <- 0.
        Frac.MQ <- 0.
        Frac.LQ <- 0.
        Frac.IQ <- 0.
        Frac.PL <- 0.
        Frac.P1L <- 0.
        Frac.P2L <- 0.
        
        N.HQ.photons <- 0.
        N.MQ.photons <- 0.
        N.LQ.photons <- 0.
        N.IQ.photons <- 0.
        N.PL.photons <- 0.
        N.P1L.photons <- 0.
        N.P2L.photons <- 0.
        
        if(N.real.photons < 2){  # if less than 2 photons in total
            Frac.HQ <- 1.
            Frac.MQ <- 0.
            Frac.LQ <- 0.
            Frac.IQ <- 0.
            Frac.PL <- 0.
            Frac.P1L <- 0.
            Frac.P2L <- 0.
        }else{
            for (i in 1:N.real.photons){
                if(i==1){            # first photon case
                    if((points[i+1]-points[i])>t2){  # HQ
                        N.HQ.photons <- N.HQ.photons +1
                    }else if((points[i+1]-points[i])>t1 && (points[i+1]-points[i])<t2){  # MQ
                        N.MQ.photons <- N.MQ.photons +1  
                    }else if((points[i+1]-points[i])<tolsep){
                        #N.PL.photons <- N.PL.photons +1  
                        N.P1L.photons <- N.P1L.photons +1  
                    }else{                           # LQ
                        N.LQ.photons <- N.LQ.photons +1
                    }
                }else if(i==N.real.photons){ # last photon case
                    if((points[i]-points[i-1])>t0){   # HQ
                        N.HQ.photons <- N.HQ.photons +1
                    }else if((points[i]-points[i-1])<tolsep){
                        #N.PL.photons <- N.PL.photons +1  
                        N.P2L.photons <- N.P2L.photons +1  
                    }else{                           # LQ
                        N.I.photons <- N.IQ.photons +1
                    }
                }else{
                    if((points[i]-points[i-1])>t0 && (points[i+1]-points[i])>t2){  # HQ
                        N.HQ.photons <- N.HQ.photons +1 
                    }else if((points[i]-points[i-1])>t0 && (points[i+1]-points[i])>t1 && 
                             (points[i+1]-points[i])<t2){  # MQ
                        N.MQ.photons <- N.MQ.photons +1 
                    }else if((points[i]-points[i-1])>t0 && ((points[i+1]-points[i])<t1) &&
                             (points[i+1]-points[i])>tolsep){ # LQ
                        N.LQ.photons <- N.LQ.photons +1
                    }else if((points[i+1]-points[i])<tolsep && (points[i]-points[i-1])>tolsep){
                        N.P1L.photons <- N.P1L.photons +1     # P1L
                    }else if((points[i]-points[i-1])<tolsep){
                        N.P2L.photons <- N.P2L.photons +1     # P2L
                    #}else if((points[i+1]-points[i])<tolsep || (points[i]-points[i-1])<tolsep){
                    #    N.PL.photons <- N.PL.photons +1  
                    }else{                          # IQ
                        N.IQ.photons <- N.IQ.photons +1
                    }
                } # if for first/last/interm. photons
            }# foreach photon
            N.PL.photons <- N.P1L.photons + N.P2L.photons
            Frac.HQ <- N.HQ.photons/N.real.photons
            Frac.MQ <- N.MQ.photons/N.real.photons
            Frac.LQ <- N.LQ.photons/N.real.photons
            Frac.IQ <- N.IQ.photons/N.real.photons
            Frac.PL <- N.PL.photons/N.real.photons
            Frac.P1L <- N.P1L.photons/N.real.photons
            Frac.P2L <- N.P2L.photons/N.real.photons
        }# if >2 photons
        
        # Calculate fractions from Expression by Poissonian theory: given a photon...
        #----------------------------------------------------------------------------
        # HQ => pH: exp(-ctrate*t0) * exp(-ctrate*t2)
        # Prob of having 0 prev. photons in t0  and     0 next photons in t2
        pH <- dpois(0,lambda*t0)*dpois(0,lambda*t2)
        
        # MQ => pM: exp(-ctrate*t1) *  exp(-ctrate*t1)- pH
        # Prob of having 0 prev. photons in t1  and     0 next photons in t1    - pH
        pM <- dpois(0,lambda*t0)*dpois(0,lambda*t1)-pH
        
        # LQ => pL: exp(-ctrate*t0) -pH - pM
        pL <- dpois(0,lambda*t0)-pH-pM
        
        # PL => pP: two pulses arriving at t<ti
        # 2*exp(-ctrate*ti)  * (1- exp(-ctrate*ti))  : one pulse < ti, the other >ti
        # + (1- exp(-ctrate*ti))*(1- exp(-ctrate*ti)) : tp and tn closer than ti
        pP <- 1-dpois(0,2*lambda*tolsep)  
        # also = ppois(0,2*lambda*tolsep, lower.tail=FALSE) -> given 1 photon, which is the prob of finding >0 ph in the double interval?
        pP2 <- 1-dpois(0,lambda*tolsep)  
        pP1 <- pP - pP2
        # IQ => pI
        pI <- 1-pH-pM-pL-pP
        
        return(list("HQ_ph"=Frac.HQ, "MQ_ph"=Frac.MQ,   "LQ_ph"=Frac.LQ, "IQ_ph"=Frac.IQ,
                    "PL_ph"=Frac.PL, "P1L_ph"=Frac.P1L, "P2L_ph"=Frac.P2L,
                    "HQ_Ps"=pH,      "MQ_Ps"=pM,        "LQ_Ps"=pL,      "IQ_Ps"=pI,     
                    "PL_Ps"=pP,      "P1L_Ps"=pP1,      "P2L_Ps"=pP2))
    }
