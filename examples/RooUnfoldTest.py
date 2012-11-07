#!/usr/bin/env python

from ROOT import gRandom, TH1, TH1D, gROOT, cout

gROOT.LoadMacro( "/home/skluth/unfold/RooUnfold/libRooUnfold.so" )

from ROOT import RooUnfoldResponse
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import RooUnfoldTUnfold
from ROOT import RooUnfoldInvert
from ROOT import RooUnfoldBasisSplines

from ROOT import TProfile, TCanvas, TF1, TMatrixD, TH2D, TVectorD

TH1.SetDefaultSumw2( False )


#  Gaussian smearing, systematic translation, and variable inefficiency:
def smear( xt, mean=-1.0, sigma=0.05, leff=True ):
    # efficiency
    xeff= 0.3 + (1.0-0.3)/10.0*xt
    x= gRandom.Rndm()
    if leff and x > xeff:
        return None
    # bias and smear
    xsmear= gRandom.Gaus( mean, sigma )
    return xt+xsmear

# Training: create response object:
def train( bininfo, gmean=-1.0, gsigma=0.05, leff=True, loufl=False,
           opttfun="exp" ):
    print "============================ TRAIN ============================="
    txt= "Smear mu, s.d.: " + str(gmean) + ", " + str(gsigma) + ", eff.: " + str(leff) + ", o/u-flow: " + str(loufl) + ", function: " + opttfun
    print txt
    # Create response matrix object:
    response= RooUnfoldResponse( bininfo["mbins"],
                                 bininfo["mlo"],
                                 bininfo["mhi"],
                                 bininfo["tbins"],
                                 bininfo["tlo"],
                                 bininfo["thi"] )
    response.UseOverflow( loufl )
    # Set training truth function:
    if "exp" in opttfun:
        trainfun= TF1( "trainexp", "exp(-x/3.0)",
                       bininfo["tlo"], bininfo["thi"] )
    elif "bw" in opttfun:
        trainfun= TF1( "trainbw", "TMath::BreitWigner(x,4.0,1.0)",
                       bininfo["tlo"], bininfo["thi"] )
    elif "gaus" in opttfun:
        trainfun= TF1( "traingaus", "TMath::Gaus(x,5.0,2.5)",
                       bininfo["tlo"], bininfo["thi"] )
    # Run training simulation:
    for i in xrange( 100000 ):
        xt= trainfun.GetRandom()
        x= smear( xt, gmean, gsigma, leff )
        if x != None:
            response.Fill( x, xt )
        else:
            response.Miss( xt )
    return response

# Generate test distributions:
def generateTest( bininfo, testfun, gmean=-1.0, gsigma=0.05, leff=True ):
    hTrue, hMeas= makeTestHistos( bininfo )
    for i in xrange( 10000 ):
        xt= testfun.GetRandom()
        hTrue.Fill( xt )
        x= smear( xt, gmean, gsigma, leff )
        if x != None:
            hMeas.Fill( x )
    return hTrue, hMeas
def makeTestHistos( bininfo ):
    hTrue= TH1D( "true", "Test Truth",
                 bininfo["tbins"],
                 bininfo["tlo"],
                 bininfo["thi"] )
    hMeas= TH1D( "meas", "Test Measured",
                 bininfo["mbins"],
                 bininfo["mlo"],
                 bininfo["mhi"] )
    return hTrue, hMeas

# Create unfolder object:
def unfolderFactory( optunf, response, hMeas, BasisSpline_m0=0 ):
    if "Bayes" in optunf:
        # Bayes unfoldung until chi^2 cut (max 10):
        unfold= RooUnfoldBayes( response, hMeas, 10, False, True )
    elif "SVD" in optunf:
        # SVD unfolding with free regularisation:
        unfold= RooUnfoldSvd( response, hMeas )
    elif "TUnfold" in optunf:
        # TUnfold with fixed regularisation tau=0.002
        unfold= RooUnfoldTUnfold( response, hMeas, 0.002 )
        # unfold= RooUnfoldTUnfold( response, hMeas )
    elif "Invert" in optunf:
        unfold= RooUnfoldInvert( response, hMeas )
    elif "Reverse" in optunf:
        # Use equivalent Bayes with 1 iteration:
        unfold= RooUnfoldBayes( response, hMeas, 1 )
    elif "BasisSplines" in optunf:
        # unfold= RooUnfoldBasisSplines( response, hMeas, 1.3e-6, 32 )
        # unfold= RooUnfoldBasisSplines( response, hMeas, 0.0, 32, 2 )
        # unfold= RooUnfoldBasisSplines( response, hMeas, 0.0, BasisSpline_m0, 2 )
        unfold= RooUnfoldBasisSplines( response, hMeas, 0.0, BasisSpline_m0, 0 )
    return unfold

# Run a test of the unfolding methods:
def rununfoldtest( bininfo, optunf, response,
                   gmean=-1.0, gsigma=0.05,
                   leff=True, optfun="exp" ):
    optunfs= [ "Bayes", "SVD", "TUnfold", "Invert", "Reverse", "BasisSplines" ]
    if not optunf in optunfs:
        txt= "Unfolding option " + optunf + " not recognised" 
        raise ValueError( txt )
    print "============================= TEST =========================="
    txt= "Method: " + optunf + ", smear mu, s.d.: " + str(gmean) + ", " + str(gsigma) + ", eff.: " +str(leff) + ", function: " + optfun
    print txt
    if "exp" in optfun:
        testfun= TF1( "testexp", "exp(-x/3.0)", 0.0, 10.0 )
    elif "bw" in optfun:
        testfun= TF1( "testbw", "TMath::BreitWigner(x,4.0,1.0)", 0.0, 10.0 )
    hTrue, hMeas= generateTest( bininfo, testfun,
                                gmean=gmean, gsigma=gsigma, leff=leff )
    print "=========================== UNFOLD ========================"
    print "Unfolding method:", optunf
    m0= 0
    if optunf == "BasisSplines":
        if optfun == "exp":
            m0= 9
        else:
            m0= 32
    unfold= unfolderFactory( optunf, response, hMeas, m0 )
    return unfold, hTrue, hMeas

# Create text for histo titles:
def funtxts( opttfun, optfun, leff, loufl ):
    if "exp" in opttfun:
        funttxt= "Exp tau=3"
    elif "bw" in opttfun:
        funttxt= "B-W mu=4, s.d.=1"
    elif "gaus" in opttfun:
        funttxt= "Gaus mu=5, s.d.=2.5"
    if "exp" in optfun:
        funtxt= "Exp tau=3"
    elif "bw" in optfun:
        funtxt= "B-W mu=4, s.d.=1"
    if leff:
        funtxt+= ", eff."
    if loufl:
        funtxt+= ", o/u-flow"
    return funttxt, funtxt

# Provide binning for different test functions:
def createBininfo( optfun ):
    if "exp" in optfun:
        bininfo= { "tbins": 10,
                   "tlo": 0.0,
                   "thi": 10.0,
                   "mbins": 20,
                   "mlo": -2.0,
                   "mhi": 10.0 }
    elif "bw" in optfun:
        bininfo= { "tbins": 40,
                   "tlo": 0.0,
                   "thi": 10.0,
                   "mbins": 80,
                   "mlo": -2.0,
                   "mhi": 10.0 }
    return bininfo

# Plot pull distributions:
def plotPulls( optunf="Bayes", ntest=10, leff=True, loufl=False,
               optfun="exp", opttfun="" ):

    if opttfun == "":
        opttfun= optfun
    bininfo= createBininfo( optfun )
    funttxt, funtxt= funtxts( opttfun, optfun, leff, loufl )

    global histos, canv
    histos= []

    canv= TCanvas( "canv", "thruth vs reco pulls", 600, 800 )
    canv.Divide( 1, 3 )

    gmeantrain= -1.0
    gmeantest= gmeantrain
    for sigma, ipad in [ [ 0.1, 1 ], [ 0.3, 2 ], [ 1.0, 3 ] ]:
        hPulls= TProfile( "pulls",
                          optunf + ", smear mu, s.d.= " +
                          str(gmeantrain) + ", " + str(sigma) +
                          ", train: " + funttxt + ", test: " + funtxt +
                          ", " + str(ntest) + " tests",
                          bininfo["tbins"], bininfo["tlo"], bininfo["thi"] )
        hPulls.SetErrorOption( "s" )
        hPulls.SetYTitle( "Thruth reco pull" )
        histos.append( hPulls )
        response= train( bininfo, gmean=gmeantrain, gsigma=sigma, leff=leff,
                         opttfun=opttfun, loufl=loufl )
        for itest in range( ntest ):
            print "Test", itest
            unfold, hTrue, hMeas= rununfoldtest( bininfo, optunf, response,
                                                 gmean=gmeantest, gsigma=sigma,
                                                 leff=leff, optfun=optfun )
            unfold.PrintTable( cout, hTrue, 2 )
            hReco= unfold.Hreco( 2 )
            nbin= hReco.GetNbinsX()
            for ibin in range( nbin+1 ):
                truevalue= hTrue.GetBinContent( ibin )
                recvalue= hReco.GetBinContent( ibin )
                error= hReco.GetBinError( ibin )
                if error > 0.0:
                    pull= ( recvalue - truevalue )/error
                    hPulls.Fill( hReco.GetBinCenter( ibin ), pull )
        canv.cd( ipad )
        hPulls.SetMinimum( -10.0 )
        hPulls.SetMaximum( 10.0 )
        hPulls.Draw()

    fname= "RooUnfoldTestPulls_" + optunf + "_" + opttfun + "_" + optfun
    if loufl:
        fname+= "_oufl"
    canv.Print( fname + ".pdf" )
        
    return

# Plots with varying smearing resolution:
def featureSizePlots( optunf="Bayes", leff=True, optfun="exp", opttfun="",
                      loufl=False ):

    if opttfun == "":
        opttfun= optfun
    bininfo= createBininfo( optfun )
    funttxt, funtxt= funtxts( opttfun, optfun, leff, loufl )

    global hReco, hMeas, hTrue, hPulls, canv, histos
    histos= []
    canv= TCanvas( "canv", "feature size plots", 600, 800 )
    canv.Divide( 1, 3 )

    gmeantrain= -1.0
    gmeantest= gmeantrain
    for sigma, ipad in [ [ 0.1, 1 ], [ 0.3, 2 ], [ 1.0, 3 ] ]:
        response= train( bininfo, gmean=gmeantrain, gsigma=sigma, leff=leff,
                         opttfun=opttfun, loufl=loufl )
        unfold, hTrue, hMeas= rununfoldtest( bininfo, optunf, response,
                                             gmean=gmeantest, gsigma=sigma,
                                             leff=leff,
                                             optfun=optfun )
        hReco= unfold.Hreco( 2 )
        histos.append( [ hTrue, hMeas, hReco ] )
        canv.cd( ipad )
        hReco.SetTitle( optunf + ", smear mu, s.d.= " +
                        str(gmeantrain) + ", " + str(sigma) +
                        ", train: " + funttxt + ", test: " + funtxt )
        hReco.SetYTitle( "Entries" )
        hReco.Draw()
        hMeas.Draw( "same" )
        hTrue.SetLineColor( 8 )
        hTrue.Draw( "same" )

    fname= "RooUnfoldTest_" + optunf + "_" + opttfun + "_" + optfun
    if loufl:
        fname+= "_oufl"
    canv.Print( fname + ".pdf" )

    return

# Make all plots of unfolding tests:
def doAllPlots( optunfs= "Bayes SVD TUnfold Invert Reverse BasisSplines" ):
    optunftokens= optunfs.split()
    for optunf in optunftokens:
        for optfun in [ "bw", "exp" ]:
            for opttfun in [ "", "gaus" ]:
                for loufl in [ False, True ]:
                    featureSizePlots( optunf=optunf,
                                      optfun=optfun, opttfun=opttfun,
                                      loufl=loufl )
                    plotPulls( optunf=optunf,
                               optfun=optfun, opttfun=opttfun,
                               loufl=loufl )
    return



    

# Run as script:
if __name__ == '__main__':
    doAllPlots()

