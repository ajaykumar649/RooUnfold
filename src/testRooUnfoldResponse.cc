// Unit tests for RooUnfoldResponse class
// A. Vanhoefer, E. Schlieckau, 10/2012

#include "RooUnfoldResponse.h"

#include "TH2.h"
#include "TH2D.h"
#include "TH1.h"

// BOOST test stuff:
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RooUnfoldResponseTests
#include <boost/test/unit_test.hpp>

// Namespaces:
using std::string;

// Test fixture for all tests:
class RooUnfoldResponseFixture{
public:
  RooUnfoldResponseFixture()
  {
    BOOST_MESSAGE( "Create RooUnfoldResponseFixture" );
    responseSameBinsMeasuredTruth = RooUnfoldResponse(10,0.,100.);
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
  const TH2D* responsenooverflow = responseWithNumberOfBins.HresponseNoOverflow();
  BOOST_CHECK_MESSAGE(1 == responseWithNumberOfBins.GetDimensionMeasured() , "Dimension measured is wrong : " << responseWithNumberOfBins.GetDimensionMeasured() << " != 1");
  BOOST_CHECK_MESSAGE(1 == responseWithNumberOfBins.GetDimensionTruth() , "Dimension truth is wrong : " << responseWithNumberOfBins.GetDimensionTruth() << " != 1");
  responseWithNumberOfBins.Miss(1.5);
  BOOST_CHECK_MESSAGE(0.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram is filled! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(0 == fakes->GetBinContent(2), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(2) << " != 0");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(2), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(2) << " != 1");
  responseWithNumberOfBins.Miss(0.5);
  BOOST_CHECK_MESSAGE(0.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram is filled! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(2), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(2) << " != 1");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(1), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(1) << " != 1");
  responseWithNumberOfBins.Miss(1.5);
  BOOST_CHECK_MESSAGE(0.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram is filled! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(2 == truth->GetBinContent(2), "Truth histogram not filled with two entries. Number of entries found: " << truth->GetBinContent(2) << " != 2");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(1), "Truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(1) << " != 1");
  BOOST_CHECK_MESSAGE(0 == truth->GetBinContent(3), "Truth histogram is filled. Number of entries found: " << truth->GetBinContent(3) << " != 0");
  responseWithNumberOfBins.Miss(2.5,5);
  BOOST_CHECK_MESSAGE(0.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram is filled! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(2 == truth->GetBinContent(2), "Truth histogram not filled with two entries. Number of entries found: " << truth->GetBinContent(2) << " != 2");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(1), "Truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(1) << " != 1");
  BOOST_CHECK_MESSAGE(5 == truth->GetBinContent(3), "Truth histogram not filled with five entries. Number of entries found: " << truth->GetBinContent(3) << " != 5");
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

BOOST_AUTO_TEST_SUITE_END()
