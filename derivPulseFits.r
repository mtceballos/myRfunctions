derivPulseFits <- function(fitsFile, fitsExt=1, fitsCol="ADC", npulses, 
                           pulseLength, startSample,headas){
    #
    # Funtion to calculate derivative of pulses in a FITS file
    # Input:
    #      fitsFile: FITS file with pulses in a column 
    #      fitsExt:  FITS file extension
    #      fitsCol:  FITS column name
    #      npulses:  numbre of pulses to be used
    #      pulseLength: pulses length for extraction (<= length of pulses in FITS file)
    #      startSample: initial sample from which pulses will be saved
    #      headas: value of HEADAS variable
    # Output:
    #      derivList: list with differentiated pulses derivList["pulse1"]=c(deriv(1), deriv(2)...)
    #      
    src <- paste("source ",headas,"/headas-init.sh",sep="")
    txtFile <- "fitsNoVar.txt"
    comm <- paste("export HEADAS=",headas," && ", src," && fdump infile=",fitsFile, "+",fitsExt," outfile=", txtFile, 
                  " columns=",fitsCol," rows=1-",npulses," prhead=no showcol=no showrow=no showunit=no clobber=yes",sep="")
    system(comm)

    DataPulse <-read.table(txtFile, header=F)
    nsamplesPerPulse <- nrow(DataPulse)/npulses
    istart <- startSample # first sample to extract the pulse (usually 1000)
    
    derivList <- list()
    for (ip in 1:npulses){
        pulseStart <- 1+(ip-1)*nsamplesPerPulse
        pulseEnd <- pulseStart + nsamplesPerPulse-1
        pulseFull <- DataPulse[pulseStart:pulseEnd,1]
        pulse <- pulseFull[istart:(istart+pulseLength-1)]
        pulseID <- paste("pulse",ip,sep="")
        derivList[[pulseID]] <- diff(pulse)
    }
    file.remove(txtFile)
    return(derivList)
}
