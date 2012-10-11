
//
// Unit tests for RooUnfold
//

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex>
#include "TH1.h"
#include "TRandom.h"

// BOOST test stuff:
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE roounfoldtest
#include <boost/test/unit_test.hpp>


#include "RooUnfoldResponse.h"
#include "RooUnfold.h"

// Namespaces:
using std::string;
using std::vector;
using std::complex;
//using namespace INIParser;

const Double_t cutdummy= -99999.0;

Double_t smear (Double_t xt)
{
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff) return cutdummy;
  Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}


// Test fixture for all tests:
class RooUnfoldTestFixture {
public:

  RooUnfold* unfold;
  //  RooUnfoldResponse response;

  RooUnfoldTestFixture(){
    BOOST_MESSAGE( "----------------------------" );
    BOOST_MESSAGE( "Create RooUnfoldTestFixture" );



    TH1::SetDefaultSumw2( false );

  // Turn on error calculation from sums of weights
  // Upsets RooUnfold::PrintTable
  // TH1::SetDefaultSumw2();
  TH1::SetDefaultSumw2( false );


  std::cout << "==================================== TRAIN ====================================" << std::endl;

  RooUnfoldResponse response (40, -10.0, 10.0, 20, -10.0, 10.0 );
  
  // Train with a Breit-Wigner, mean 0.3 and width 2.5.
  for (Int_t i= 0; i<100000; i++) {
    Double_t xt= gRandom->BreitWigner (0.3, 2.5);
    Double_t x= smear (xt);
    
    response.Fill (x, xt);
    
  }
  
  std::cout << "==================================== TEST =====================================" << std::endl;
  // TH1D* hTrue= new TH1D ("true", "Test Truth",    40, -10.0, 10.0);
  TH1D* hTrue= new TH1D( "true", "Test Truth", 20, -10.0, 10.0 );

  TH1D* hMeas= new TH1D ("meas", "Test Measured", 40, -10.0, 10.0);
  // Test with a Gaussian, mean 0 and width 2.
  for (Int_t i=0; i<10000; i++) {
    Double_t xt= gRandom->Gaus (0.0, 2.0), x= smear (xt);
    hTrue->Fill(xt);
    if (x!=cutdummy) hMeas->Fill(x);
  }

  std::cout << "==================================== UNFOLD ===================================" << std::endl;
  // RooUnfoldBayes   unfold (&response, hMeas, 4);    // OR
  RooUnfold* unfold= new RooUnfold( &response, hMeas);
//RooUnfoldSvd     unfold (&response, hMeas, 20);   // OR
//RooUnfoldTUnfold unfold (&response, hMeas);

  TH1D* hReco= (TH1D*) unfold->Hreco();

  unfold->Print();
  unfold->PrintTable (std::cout, hTrue);
  //  hReco->Draw();
  //  hMeas->Draw("SAME");
  //  hTrue->SetLineColor(8);
  //  hTrue->Draw("SAME");


  }
  virtual ~RooUnfoldTestFixture() {
    BOOST_MESSAGE( "Tear down RooUnfoldTestFixture" );
    BOOST_MESSAGE( "----------------------------" );
  }
  //  RooUnfold reader;
  //roounford pointer 
};

// Declare test suite name and fixture class to BOOST:
BOOST_FIXTURE_TEST_SUITE( roounfoldsuite, RooUnfoldTestFixture )

// Test cases:

// Dummy Test
BOOST_AUTO_TEST_CASE( testRooUnfoldDummyTest ) {
  BOOST_CHECK_EQUAL( 0, 0 );
}

BOOST_AUTO_TEST_CASE( testConstructor ) {
  BOOST_MESSAGE( "testConstructor" );
}


BOOST_AUTO_TEST_CASE( Print){
  BOOST_MESSAGE("***Print test***");

}

BOOST_AUTO_TEST_CASE( PrintTable){
  BOOST_MESSAGE("***PrintTable test***");

}

BOOST_AUTO_TEST_CASE( GetErrMat){
  BOOST_MESSAGE("***GetErrMat test***");

}

BOOST_AUTO_TEST_CASE(RunToy){
  BOOST_MESSAGE("RunToy test");
}

BOOST_AUTO_TEST_CASE(GetStepSizeParm){
  BOOST_MESSAGE("GetStepSizeParm test");
  //  unfold->Print();
  //  BOOST_CHECK_MESSAGE(unfold->GetStepSizeParm()==0,"bla");
  //  BOOST_CHECK_EQUAL(GetStepSizeParm(), 0.);
}



//   string value, expectedValue;
//   value= reader.get( "user", "multiline", "" );
//   expectedValue= "this is\na\nmultiline";
//   BOOST_CHECK_EQUAL( value, expectedValue );
// }

// Test error state after parsing:
// BOOST_AUTO_TEST_CASE( testParse ) {
//   BOOST_CHECK_EQUAL( reader.parseError(), 0 );
// }

// Test get string:
// BOOST_AUTO_TEST_CASE( testGetString ) {
//   string value, expectedValue;
//   value= reader.get( "user", "name", "" );
//   expectedValue= "Bob Smith";
//   BOOST_CHECK_EQUAL( value, expectedValue );
//   expectedValue= "default";
//   value= reader.get( "user", "???", "default" );
//   BOOST_CHECK_EQUAL( value, expectedValue );
//   value= reader.get( "???", "name", "default" );
//   BOOST_CHECK_EQUAL( value, expectedValue );
// }

// Test get multi-line string:
// BOOST_AUTO_TEST_CASE( testGetMultilineString ) {
//   string value, expectedValue;
//   value= reader.get( "user", "multiline", "" );
//   expectedValue= "this is\na\nmultiline";
//   BOOST_CHECK_EQUAL( value, expectedValue );
// }

// Test tokens from string:
// BOOST_AUTO_TEST_CASE( testGetStringTokens ) {
//   string values= reader.get( "user", "multiline", "" );
//   vector<string> tokens= getTokens( values );
//   vector<string> expected;
//   expected.push_back( "this" );
//   expected.push_back( "is" );
//   expected.push_back( "a" );
//   expected.push_back( "multiline" );
//   for( size_t i= 0; i < expected.size(); i++ ) {
//     BOOST_CHECK_EQUAL( tokens[i], expected[i] );
//   }
// }

// Test get all names sorted in a section:
// BOOST_AUTO_TEST_CASE( testGetNames ) {
//   vector<string> names= reader.getNames( "protocol" );
//   vector<string> expected;
//   expected.push_back( "complex" );
//   expected.push_back( "number" );
//   expected.push_back( "version" );
//   for( size_t i= 0; i < expected.size(); i++ ) {
//     BOOST_CHECK_EQUAL( expected[i], names[i] );
//   }
// }

// Test bools:
// BOOST_AUTO_TEST_CASE( testGetBools ) {
//   bool value;
//   value= reader.getType( "user", "activetrue", false );
//   BOOST_CHECK_EQUAL( value, true );
//   value= reader.getType( "user", "activeyes", false );
//   BOOST_CHECK_EQUAL( value, true );
//   value= reader.getType( "user", "activeon", false );
//   BOOST_CHECK_EQUAL( value, true );
//   value= reader.getType( "user", "active1", false );
//   BOOST_CHECK_EQUAL( value, true );
//   value= reader.getType( "user", "activefalse", true );
//   BOOST_CHECK_EQUAL( value, false );
//   value= reader.getType( "user", "activeno", true );
//   BOOST_CHECK_EQUAL( value, false );
//   value= reader.getType( "user", "activeoff", true );
//   BOOST_CHECK_EQUAL( value, false );
//   value= reader.getType( "user", "active0", true );
//   BOOST_CHECK_EQUAL( value, false );
//   value= reader.getType( "user", "???", true );
//   BOOST_CHECK_EQUAL( value, true );
// }

// Test long:
// BOOST_AUTO_TEST_CASE( testGetInteger ) {
//   int value;
//   value= reader.getType( "protocol", "version", 0 );
//   BOOST_CHECK_EQUAL( value, 6 );
//   value= reader.getType( "protocol", "???", 0 );
//   BOOST_CHECK_EQUAL( value, 0 );
// }

// // Test double:
// BOOST_AUTO_TEST_CASE( testGetDouble ) {
//   double value;
//   value= reader.getType( "protocol", "number", 0.0 );
//   BOOST_CHECK_EQUAL( value, 1.2345e9 );
//   value= reader.getType( "protocol", "???", 0.0 );
//   BOOST_CHECK_EQUAL( value, 0.0 );
// }

// Test complex:
// BOOST_AUTO_TEST_CASE( testGetComplex ) {
//   complex<double> value;
//   value= reader.getType( "protocol", "complex", complex<double>() );
//   BOOST_CHECK_EQUAL( value, complex<double>(1.0,1.0) );
//   value= reader.getType( "protocol", "???", complex<double>() );
//   BOOST_CHECK_EQUAL( value, complex<double>() );
// }

BOOST_AUTO_TEST_SUITE_END()

