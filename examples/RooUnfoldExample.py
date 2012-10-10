#!/usr/bin/env python
# ==============================================================================
#  File and Version Information:
#       $Id: RooUnfoldExample.py 302 2011-09-30 20:39:20Z T.J.Adye $
#
#  Description:
#       Simple example usage of the RooUnfold package using toy MC.
#
#  Author: Tim Adye <T.J.Adye@rl.ac.uk>
#
# ==============================================================================

from ROOT import gRandom, TH1, TH1D, cout, gROOT

gROOT.LoadMacro( "/afs/cern.ch/user/k/kirschen/scratch0/RooUnfold/libRooUnfold.so" )

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import RooUnfoldTUnfold
from ROOT import RooUnfoldInvert


from ROOT import TH1
TH1.SetDefaultSumw2( False )


# ==============================================================================
#  Gaussian smearing, systematic translation, and variable inefficiency
# ==============================================================================

def smear(xt):
    # efficiency
    xeff= 0.3 + (1.0-0.3)/20.0*(xt+10.0)
    x= gRandom.Rndm()
    if x > xeff:
        return None
    # bias and smear
    xsmear= gRandom.Gaus( -2.5, 0.2 )
    return xt+xsmear

# ==============================================================================
#  Example Unfolding
# ==============================================================================

def main( optunf="Bayes" ):

    optunfs= [ "Bayes", "SVD", "TUnfold", "Invert", "Reverse" ]
    if not optunf in optunfs:
        txt= "Unfolding option " + optunf + " not recognised" 
        raise ValueError( txt )  

    global hReco, hMeas, hTrue

    print "==================================== TRAIN ===================================="
    # Create response matrix object for 40 measured and 20
    # unfolded bins:
    response= RooUnfoldResponse( 40, -10.0, 10.0, 20, -10.0, 10.0 )

    #  Train with a Breit-Wigner, mean 0.3 and width 2.5.
    for i in xrange(100000):
        # xt= gRandom.BreitWigner( 0.3, 2.5 )
        xt= gRandom.Gaus( 0.0, 5.0 )
        x= smear (xt)
        if x != None:
            response.Fill( x, xt )
        else:
            response.Miss( xt )
            
    print "==================================== TEST ====================================="
    hTrue= TH1D( "true", "Test Truth", 20, -10.0, 10.0 )
    hMeas= TH1D( "meas", "Test Measured", 40, -10.0, 10.0 )
    #  Test with a Gaussian, mean 0 and width 2.
    for i in xrange(10000):
        # xt= gRandom.Gaus( 0.0, 2.0 )
        xt= gRandom.BreitWigner( 0.3, 2.5 )
        x= smear( xt )
        hTrue.Fill( xt )
        if x != None:
            hMeas.Fill( x )

    print "==================================== UNFOLD ==================================="
    print "Unfolding method:", optunf
    if "Bayes" in optunf:
        # Bayes unfoldung with 4 iterations
        # unfold= RooUnfoldBayes( response, hMeas, 4 )
        unfold= RooUnfoldBayes( response, hMeas, 10, False, True )
    elif "SVD" in optunf:
        # SVD unfoding with free regularisation
        # unfold= RooUnfoldSvd( response, hMeas, 20 )
        unfold= RooUnfoldSvd( response, hMeas )
    elif "TUnfold" in optunf:
        # TUnfold with fixed regularisation tau=0.002
        # unfold= RooUnfoldTUnfold( response, hMeas )
        unfold= RooUnfoldTUnfold( response, hMeas, 0.002 )
    elif "Invert" in optunf:
        unfold= RooUnfoldInvert( response, hMeas )
    elif "Reverse" in optunf:
        unfold= RooUnfoldBayes( response, hMeas, 1 )

    hReco= unfold.Hreco()
    # unfold.PrintTable( cout, hTrue )
    unfold.PrintTable( cout, hTrue, 2 )

    hReco.Draw()
    hMeas.Draw("SAME")
    hTrue.SetLineColor(8)
    hTrue.Draw("SAME")

    return

if __name__ == '__main__':
   main()

