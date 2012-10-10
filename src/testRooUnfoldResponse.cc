// Unit tests for RooUnfoldResponse class
// A. Vanhoefer, E. Schlieckau, 10/2012

#include "RooUnfoldResponse.h"

#include "TRandom.h"
#include "TH2D.h"

// BOOST test stuff:
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RooUnfoldResponseTests
#include <boost/test/unit_test.hpp>

// Namespaces:
using std::string;


//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

// Test fixture for all tests:
class RooUnfoldResponseFixture{
public:
  RooUnfoldResponseFixture()
  {
    BOOST_MESSAGE( "Create RooUnfoldResponseFixture" );
    //Initialize RooUnfoldResponse instance with same entries for testing;
    responseFilledWithSomeEntries=RooUnfoldResponse(40, -10.0, 10.0);
    TRandom FixedSeedRandom(111);
    // Train with a Breit-Wigner, mean 0.3 and width 2.5.
    for (Int_t i=0; i<10000; i++) {
      double xt = FixedSeedRandom.BreitWigner(0.3, 2.5);
      double x = xt + smear(xt); //introduce
						       //bias and smear
    if (x!=cutdummy)
      responseFilledWithSomeEntries.Fill (x, xt);
    else
      responseFilledWithSomeEntries.Miss (xt);
    }
    responseSameBinsMeasuredTruth = RooUnfoldResponse(10,0.,100.);
  }
  virtual ~RooUnfoldResponseFixture(){
    BOOST_MESSAGE( "Tear down RooUnfoldResponseFixture" );
  }
  RooUnfoldResponse response;
  RooUnfoldResponse responseSameBinsMeasuredTruth;
  RooUnfoldResponse responseFilledWithSomeEntries;

private:
  Double_t smear (Double_t xt)
  {
    Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
    Double_t x= gRandom->Rndm();
    if (x>xeff) return cutdummy;
    Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
    return xt+xsmear;
  }
  
};



// Declare test suite name and fixture class to BOOST:
BOOST_FIXTURE_TEST_SUITE( RooUnfoldResponseSuite, RooUnfoldResponseFixture )

// Test cases:

// Failing test
//BOOST_AUTO_TEST_CASE(testFails){
// BOOST_FAIL("This test fails!");
//}

BOOST_AUTO_TEST_CASE(testConstructorNumberOfBins){
  int numberOfBins = 10;
  double low = 1;
  double high = 11;
  RooUnfoldResponse responseWithNumberOfBins(numberOfBins, low, high);

  int resultNumberOfBinsMeasured = responseWithNumberOfBins.GetNbinsMeasured();
  int resultNumberOfBinsTruth    = responseWithNumberOfBins.GetNbinsTruth();
  TH2* responseHistogram = responseWithNumberOfBins.Hresponse();
  double responseHistogramLow = responseHistogram->GetXaxis()->GetBinLowEdge(1);
  double responseHistogramHigh = responseHistogram->GetXaxis()->GetBinLowEdge(numberOfBins)+responseHistogram->GetXaxis()->GetBinWidth(numberOfBins);
  BOOST_CHECK_MESSAGE(numberOfBins == resultNumberOfBinsMeasured, "Number of bins measured not on given value: " << resultNumberOfBinsMeasured << " != " << numberOfBins);
  BOOST_CHECK_MESSAGE(numberOfBins == resultNumberOfBinsTruth, "Number of bins truth not on given value: " << resultNumberOfBinsTruth << " != " << numberOfBins);
  BOOST_CHECK_MESSAGE(resultNumberOfBinsMeasured == resultNumberOfBinsTruth, "Number of bins truth not equal to number of bins measured: " << resultNumberOfBinsTruth << " != " << resultNumberOfBinsMeasured);
  BOOST_CHECK_MESSAGE(low == responseHistogramLow, "First bin low edge not taken correctly: " << responseHistogramLow << " != " << low);
  BOOST_CHECK_MESSAGE(high == responseHistogramHigh, "Last bin high edge not taken correctly: " << responseHistogramHigh << " != " << high);
}


BOOST_AUTO_TEST_CASE(testmethodMiss){
  int numberOfBins = 3;
  double low = 0;
  double high = 3;
  RooUnfoldResponse responseWithNumberOfBins(numberOfBins, low, high);
  const TH1* measured = responseWithNumberOfBins.Hmeasured();
  const TH1* fakes = responseWithNumberOfBins.Hfakes();
  const TH1* truth = responseWithNumberOfBins.Htruth();
  const TH2* response = responseWithNumberOfBins.Hresponse();
  responseWithNumberOfBins.Miss(1.5);
  BOOST_CHECK_MESSAGE(0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0 == measured->GetBinContent(2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(2) << " != 0");
  BOOST_CHECK_MESSAGE(0 == fakes->GetBinContent(2), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(2) << " != 0");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(2), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(2) << " != 1");
  BOOST_CHECK_MESSAGE(0 == measured->GetBinContent(0,2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(0,2) << " != 0");
  responseWithNumberOfBins.Miss(0.5);
  BOOST_CHECK_MESSAGE(0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0 == measured->GetBinContent(2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(2) << " != 0");
  BOOST_CHECK_MESSAGE(0 == fakes->GetBinContent(2), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(2) << " != 0");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(2), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(2) << " != 1");
  BOOST_CHECK_MESSAGE(0 == measured->GetBinContent(0,2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(0,2) << " != 0");
  BOOST_CHECK_MESSAGE(0 == measured->GetBinContent(1), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(1) << " != 0");
  BOOST_CHECK_MESSAGE(0 == fakes->GetBinContent(1), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(1) << " != 0");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(1), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(1) << " != 1");
  BOOST_CHECK_MESSAGE(0 == measured->GetBinContent(0,1), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(0,1) << " != 0");
  responseWithNumberOfBins.Miss(1.5);
  BOOST_CHECK_MESSAGE(0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0 == measured->GetBinContent(2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(2) << " != 0");
  BOOST_CHECK_MESSAGE(0 == fakes->GetBinContent(2), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(2) << " != 0");
  BOOST_CHECK_MESSAGE(2 == truth->GetBinContent(2), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(2) << " != 2");
  BOOST_CHECK_MESSAGE(0 == measured->GetBinContent(0,2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(0,2) << " != 0");
  BOOST_CHECK_MESSAGE(0 == measured->GetBinContent(1), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(1) << " != 0");
  BOOST_CHECK_MESSAGE(0 == fakes->GetBinContent(1), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(1) << " != 0");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(1), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(1) << " != 1");
  BOOST_CHECK_MESSAGE(0 == measured->GetBinContent(0,1), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(0,1) << " != 0");
  }

BOOST_AUTO_TEST_CASE(testFill1D){
  //test with default weight
  double xMeasured = 42;
  double xTruth = 74;
  responseSameBinsMeasuredTruth.Fill(xMeasured,xTruth);
  TH1* measuredHistogram = responseSameBinsMeasuredTruth.Hmeasured();
  int entriesMeasured = measuredHistogram->GetEntries();
  //int 
  //BOOST_CHECK_MESSAGE();


  TH1* truthHistogram = responseSameBinsMeasuredTruth.Htruth();
  TH2* responseHistogram = responseSameBinsMeasuredTruth.Hresponse();
}

//Test of UseOverflowStatus
BOOST_AUTO_TEST_CASE(testUseOverflowStatus){
  //  RooUnfoldResponse testObject;
  BOOST_CHECK_MESSAGE(response.UseOverflowStatus()==false,"default constructor does not initialize with overflow set to false");

}


//Test of add-function
BOOST_AUTO_TEST_CASE(testAddFunction){
  RooUnfoldResponse testObject = responseFilledWithSomeEntries;
  int noOfEntriesInMeasuredHist = testObject.Hmeasured()->GetEntries(); 
  //  std::cout << testObject.Hmeasured()->GetEntries() << std::endl;
  testObject.Add(responseFilledWithSomeEntries);
  BOOST_CHECK_MESSAGE(testObject.Hmeasured()->GetEntries()==2*noOfEntriesInMeasuredHist,"Adding the same RooUnfoldResponse did not result in twice the number of entries in Hmeasured");
  //  std::cout << testObject.Hmeasured()->GetEntries() << std::endl;

  //should be extended :)
}



BOOST_AUTO_TEST_SUITE_END()
