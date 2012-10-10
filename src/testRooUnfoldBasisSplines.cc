
// Unit tests for RooUnfoldBasisSplines

#include "RooUnfoldBasisSplines.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex>

// BOOST test stuff:
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE roounfoldbasissplinestest
#include <boost/test/unit_test.hpp>

// Namespaces:
using std::string;
using std::vector;
using std::complex;

// Test fixture for all tests:
class RooUnfoldBasisSplinesTestFixture {
public:
  RooUnfoldBasisSplinesTestFixture(){
    BOOST_MESSAGE( "Create RooUnfoldBasisSplinesTestFixture" );
  }
  virtual ~RooUnfoldBasisSplinesTestFixture() {
    BOOST_MESSAGE( "Tear down RooUnfoldBasisSplinesTestFixture" );
  }
};

// Declare test suite name and fixture class to BOOST:
BOOST_FIXTURE_TEST_SUITE( roounfoldbasissplinessuite, RooUnfoldBasisSplinesTestFixture )

// Test cases:

// Dummy Test
BOOST_AUTO_TEST_CASE( testRooUnfoldBasisSplinesDummyTest ) {
  BOOST_CHECK_EQUAL( 0, 0 );
}

//test of the empty constructor
BOOST_AUTO_TEST_CASE( testRooUnfoldBasisSplinesEmptyConstructor ){
  RooUnfoldBasisSplines emptyConstr;
  BOOST_CHECK_EQUAL( emptyConstr.GetMinParm(), 0 );
  BOOST_CHECK_EQUAL( emptyConstr.GetMaxParm(), 0 );
  BOOST_CHECK_EQUAL( emptyConstr.GetStepSizeParm(), 0 );
  BOOST_CHECK_EQUAL( emptyConstr.GetDefaultParm(), 0 );
}

//test of the named constructor
BOOST_AUTO_TEST_CASE( testRooUnfoldBasisSplinesConstCharConstr ){
  const char* name = "ciccio";
  const char* title = "panza";
  RooUnfoldBasisSplines charConstr(name, title);
  BOOST_CHECK_EQUAL( charConstr.GetName() , name );
  BOOST_CHECK_EQUAL( charConstr.GetTitle() , title );
}

//test of the named constructor with strigs
BOOST_AUTO_TEST_CASE( testRootUnfoldSplinestStringConstr ){
  TString name = "frengo";
  TString title = "sticastica";
  RooUnfoldBasisSplines strConstr(name, title);
  BOOST_CHECK_EQUAL( strConstr.GetName(), name );
  BOOST_CHECK_EQUAL( strConstr.GetTitle(), title );
}

//test of the named constructor with a constructor
BOOST_AUTO_TEST_CASE( testRootUnfoldSplinesTotalConstr ){
  RooUnfoldResponse res;
  Int_t entries=1000;
  TH1F meas("my name is meas!!","", 100,0,100);
  meas.FillRandom("gaus",entries );
  Double_t tau=1.0e-5;
  Int_t m0=1;
  Int_t iauto=1000;
  const char *name="pirla";
  const char* title="yuyu";
  RooUnfoldBasisSplines Total(&res, &meas, tau, m0, iauto, name, title);
  BOOST_CHECK_EQUAL( Total.response(), &res );
  BOOST_CHECK_EQUAL( Total.Hmeasured(), &meas );
  // to be finished, add the check for the other variables....
}

BOOST_AUTO_TEST_SUITE_END()

