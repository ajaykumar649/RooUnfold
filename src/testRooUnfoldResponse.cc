// Unit tests for RooUnfoldResponse class
// A. Vanhoefer, E. Schlieckau, 10/2012

#include "RooUnfoldResponse.h"

#include "TH2.h"

// BOOST test stuff:
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RooUnfoldResponseTests
#include <boost/test/unit_test.hpp>

// Namespaces:
using std::string;

// Test fixture for all tests:
class RooUnfoldResponseFixture{
public:
  RooUnfoldResponseFixture():
    responseSameBinsMeasuredTruth(10,0.,100.)
  {
    BOOST_MESSAGE( "Create RooUnfoldResponseFixture" );
  }
  virtual ~RooUnfoldResponseFixture(){
    BOOST_MESSAGE( "Tear down RooUnfoldResponseFixture" );
  }
  RooUnfoldResponse response;
  RooUnfoldResponse responseSameBinsMeasuredTruth;
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
  BOOST_CHECK_MESSAGE(numberOfBins == resultNumberOfBinsMeasured, "Number of bins measured not on given value: " << resultNumberOfBinsMeasured << " != " << numberOfBins);
  BOOST_CHECK_MESSAGE(numberOfBins == resultNumberOfBinsTruth, "Number of bins truth not on given value: " << resultNumberOfBinsTruth << " != " << numberOfBins);
  BOOST_CHECK_MESSAGE(resultNumberOfBinsMeasured == resultNumberOfBinsTruth, "Number of bins truth not equal to number of bins measured: " << resultNumberOfBinsTruth << " != " << resultNumberOfBinsMeasured);

  TH2* responseHistogram = responseWithNumberOfBins.Hresponse();
  double responseHistogramLow = responseHistogram->GetXaxis()->GetBinLowEdge(1);
  double responseHistogramHigh = responseHistogram->GetXaxis()->GetBinLowEdge(numberOfBins)+responseHistogram->GetXaxis()->GetBinWidth(numberOfBins);
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

BOOST_AUTO_TEST_CASE(testFill1D){
  //test with default weight
  double xMeasured = 42;
  double xTruth = 74;

  responseSameBinsMeasuredTruth.Fill(xMeasured,xTruth);
  //test if measured histogram was filled right
  TH1* measuredHistogram = responseSameBinsMeasuredTruth.Hmeasured();

  int entriesMeasured = measuredHistogram->GetEntries();
  BOOST_CHECK_MESSAGE(entriesMeasured==1,"Wrong number of total entries for the measured histogram. Expected 1, but was: "<<entriesMeasured);
  double entriesWeightedMeasured = measuredHistogram->Integral();
  BOOST_CHECK_MESSAGE(entriesWeightedMeasured==1,"Wrong number of weighted entries for the measured histogram. Expected 1, but was: "<<entriesWeightedMeasured);

  int binMeasured = measuredHistogram->FindBin(xMeasured);
  double binContentMeasured =measuredHistogram->GetBinContent(binMeasured);
  BOOST_CHECK_MESSAGE(binContentMeasured==1,"Wrong bin was filled for the measured histogram, expected, bin 4 to be filled, but filled bin: "<<binMeasured);

 //test if truth histogram was filled right
  TH1* truthHistogram = responseSameBinsMeasuredTruth.Htruth();

  int entriesTruth = truthHistogram->GetEntries();
  BOOST_CHECK_MESSAGE(entriesTruth==1,"Wrong number of total entries for the truth histogram. Expected 1, but was: "<<entriesTruth);
  double entriesWeightedTruth = truthHistogram->Integral();
  BOOST_CHECK_MESSAGE(entriesWeightedTruth==1,"Wrong number of weighted entries for the truth histogram. Expected 1, but was: "<<entriesWeightedTruth);

  int binTruth = truthHistogram->FindBin(xTruth);
  double binContentTruth =truthHistogram->GetBinContent(binTruth);
  BOOST_CHECK_MESSAGE(binContentTruth==1,"Wrong bin was filled for the truth histogram, expected, bin 4 to be filled, but filled bin: "<<binTruth);

   //test if response histogram was filled right
  TH2* responseHistogram = responseSameBinsMeasuredTruth.Hresponse();

  int entriesResponse = responseHistogram->GetEntries();
  BOOST_CHECK_MESSAGE(entriesResponse==1,"Wrong number of total entries for the truth histogram. Expected 1, but was: "<<entriesResponse);
  double entriesWeightedResponse = truthHistogram->Integral();
  BOOST_CHECK_MESSAGE(entriesWeightedResponse==1,"Wrong number of weighted entries for the truth histogram. Expected 1, but was: "<<entriesWeightedResponse);

  int binResponse = responseHistogram->FindBin(xMeasured,xTruth);
  double binContentResponse =responseHistogram->GetBinContent(binResponse);
  BOOST_CHECK_MESSAGE(binContentResponse==1,"Wrong bin was filled for the truth histogram, expected, bin 4 to be filled, but filled bin: "<<binResponse);
}

//Test of UseOverflowStatus
BOOST_AUTO_TEST_CASE(testUseOverflowStatus){
  //  RooUnfoldResponse testObject;
  BOOST_CHECK_MESSAGE(response.UseOverflowStatus()==false,"default constructor does not initialize with overflow set to false");
}

BOOST_AUTO_TEST_SUITE_END()
