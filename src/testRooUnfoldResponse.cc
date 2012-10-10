// Unit tests for RooUnfoldResponse class
// A. Vanhoefer, E. Schlieckau, 10/2012

#include "RooUnfoldResponse.h"

#include "TRandom.h"
#include "TH2D.h"

// BOOST test stuff:
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RooUnfoldResponseTests
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// Namespaces:
using std::string;


//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

// Test fixture for all tests:
class RooUnfoldResponseFixture{
public:
  RooUnfoldResponseFixture():
    responseSameBinsMeasuredTruth(10,0.,100.)
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
  int measuredDimensions = responseWithNumberOfBins.GetDimensionMeasured();
  int truthDimensions = responseWithNumberOfBins.GetDimensionTruth();
  BOOST_CHECK_MESSAGE(measuredDimensions == 1, "Wrong measured dimension, has to be 1 but is: " << measuredDimensions);
  BOOST_CHECK_MESSAGE(truthDimensions == 1, "Wrong truth dimension, has to be 1 but is: " << truthDimensions);

  TH1* measuredHistogram = responseWithNumberOfBins.Hmeasured();
  int entriesMeasured = measuredHistogram->GetEntries();
  BOOST_CHECK_MESSAGE(entriesMeasured==0,"Wrong number of total entries for the measured histogram. Expected 0, but was: "<<entriesMeasured);
}

BOOST_AUTO_TEST_CASE(testmethodMiss){
  int numberOfBins = 3;
  double low = 0;
  double high = 3;
  RooUnfoldResponse responseWithNumberOfBins(numberOfBins, low, high);
  responseWithNumberOfBins.Miss(1.5);
  const TH1* measured = responseWithNumberOfBins.Hmeasured();
  const TH1* fakes = responseWithNumberOfBins.Hfakes();
  const TH1* truth = responseWithNumberOfBins.Htruth();
  const TH2* response = responseWithNumberOfBins.Hresponse();
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

BOOST_AUTO_TEST_CASE(testmethodFake){
  int numberOfBins = 3;
  double low = 0;
  double high = 3;
  RooUnfoldResponse responseWithNumberOfBins(numberOfBins, low, high);
  responseWithNumberOfBins.Fake(1.5);
  const TH1* measured = responseWithNumberOfBins.Hmeasured();
  const TH1* fakes = responseWithNumberOfBins.Hfakes();
  const TH1* truth = responseWithNumberOfBins.Htruth();
  const TH2* response = responseWithNumberOfBins.Hresponse();
  BOOST_CHECK_MESSAGE(1 == measured->GetBinContent(2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(2) << " != 1");
  BOOST_CHECK_MESSAGE(1 == fakes->GetBinContent(2), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(2) << " != 1");
  BOOST_CHECK_MESSAGE(0 == truth->GetBinContent(2), "truth histogram not filled with zero entry. Number of entries found: " << truth->GetBinContent(2) << " != 0");
}

BOOST_AUTO_TEST_CASE(testFill1D){
  //test with default weight
  double xMeasured = 42;
  double xTruth = 74;
  double  weight =7.3;

  responseSameBinsMeasuredTruth.Fill(xMeasured,xTruth,weight);
  //test if measured histogram was filled right
  TH1* measuredHistogram = responseSameBinsMeasuredTruth.Hmeasured();

  int entriesMeasured = measuredHistogram->GetEntries();
  BOOST_CHECK_MESSAGE(entriesMeasured==1,"Wrong number of total entries for the measured histogram. Expected 1, but was: "<<entriesMeasured);
  double entriesWeightedMeasured = measuredHistogram->Integral();
  BOOST_CHECK_MESSAGE(entriesWeightedMeasured==weight,"Wrong number of weighted entries for the measured histogram. Expected "<< weight <<", but was: "<<entriesWeightedMeasured);

  int binMeasured = measuredHistogram->FindBin(xMeasured);
  double binContentMeasured =measuredHistogram->GetBinContent(binMeasured);
  BOOST_CHECK_MESSAGE(binContentMeasured==weight,"Wrong bin was filled for the measured histogram, expected, bin 4 to be filled, but filled bin: "<<binMeasured);

 //test if truth histogram was filled right
  TH1* truthHistogram = responseSameBinsMeasuredTruth.Htruth();

  int entriesTruth = truthHistogram->GetEntries();
  BOOST_CHECK_MESSAGE(entriesTruth==1,"Wrong number of total entries for the truth histogram. Expected 1, but was: "<<entriesTruth);
  double entriesWeightedTruth = truthHistogram->Integral();
  BOOST_CHECK_MESSAGE(entriesWeightedTruth==weight,"Wrong number of weighted entries for the truth histogram. Expected  "<< weight <<", but was: "<<entriesWeightedTruth);

  int binTruth = truthHistogram->FindBin(xTruth);
  double binContentTruth =truthHistogram->GetBinContent(binTruth);
  BOOST_CHECK_MESSAGE(binContentTruth==weight,"Wrong bin was filled for the truth histogram, expected, bin 7 to be filled, but filled bin: "<<binTruth);

   //test if response histogram was filled right
  TH2* responseHistogram = responseSameBinsMeasuredTruth.Hresponse();

  int entriesResponse = responseHistogram->GetEntries();
  BOOST_CHECK_MESSAGE(entriesResponse==1,"Wrong number of total entries for the truth histogram. Expected 1, but was: "<<entriesResponse);
  double entriesWeightedResponse = truthHistogram->Integral();
  BOOST_CHECK_MESSAGE(entriesWeightedResponse==weight,"Wrong number of weighted entries for the truth histogram. Expected  "<< weight <<", but was: "<<entriesWeightedResponse);

  int binResponse = responseHistogram->FindBin(xMeasured,xTruth);
  double binContentResponse =responseHistogram->GetBinContent(binResponse);
  BOOST_CHECK_MESSAGE(binContentResponse==weight,"Wrong bin was filled for the truth histogram, expected, bin 101 to be filled, but filled bin: "<<binResponse);

  //test default value for weight for filling
  responseSameBinsMeasuredTruth.Fill(xMeasured,xTruth);

  int entries = truthHistogram->GetEntries();
  BOOST_CHECK_MESSAGE(entries==2,"Wrong number of total entries. Expected 2, but was: "<<entries);
  double resultDefaultWeight = truthHistogram->Integral()-weight;
  BOOST_CHECK_CLOSE( 1, resultDefaultWeight, 0.0001 );
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
