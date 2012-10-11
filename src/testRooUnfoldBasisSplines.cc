
// Unit tests for RooUnfoldBasisSplines

#include "RooUnfoldBasisSplines.h"
#include "testHelperRooUnfoldBasisSplines.h"

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
  RooUnfoldBasisSplinesTestFixture()
    : res(), entries(1000), meas("my name is meas!!","", 100,0,100),
      tau(1.0e-5), m0(1), iauto(1000), name("pirla"), title("yuyu"),
      Total(&res, &meas, tau, m0, iauto, name, title),
      testHelper(&Total)
  {
    BOOST_MESSAGE( "Create RooUnfoldBasisSplinesTestFixture" );

    meas.FillRandom("gaus", entries);
  }
  virtual ~RooUnfoldBasisSplinesTestFixture() {
    BOOST_MESSAGE( "Tear down RooUnfoldBasisSplinesTestFixture" );
  }

  RooUnfoldResponse res;
  Int_t entries;
  TH1F meas;
  Double_t tau;
  Int_t m0;
  Int_t iauto;
  const char* name;
  const char* title;
  RooUnfoldBasisSplines Total;
  testHelperRooUnfoldBasisSplines testHelper;
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
  BOOST_CHECK_EQUAL( Total.response(), &res );
  BOOST_CHECK_EQUAL( Total.Hmeasured(), &meas );
  BOOST_CHECK_EQUAL( testHelper.GetTau(), tau );
  BOOST_CHECK_EQUAL( testHelper.GetM0(), m0 );
  BOOST_CHECK_EQUAL( testHelper.GetIauto(), iauto );
  BOOST_CHECK_EQUAL( Total.GetName(), name );
  BOOST_CHECK_EQUAL( Total.GetTitle(), title );
}

//test the copy constructor
BOOST_AUTO_TEST_CASE( testRooUnfoldSplinesCopyConstructor ){

}


BOOST_AUTO_TEST_SUITE_END()

