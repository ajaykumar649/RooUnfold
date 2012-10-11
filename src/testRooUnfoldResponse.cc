// Unit tests for RooUnfoldResponse class
// A. Vanhoefer, E. Schlieckau, 10/2012

#include "RooUnfoldResponse.h"

#include "TRandom.h"
#include "TH2D.h"
#include "TString.h"

// BOOST test stuff:
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RooUnfoldResponseTests
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// Namespaces:
using std::string;


// Test fixture for all tests:
class RooUnfoldResponseFixture{
public:
  Double_t const returnCutDummy() {return cutdummy_;}
  RooUnfoldResponseFixture():
    responseSameBinsMeasuredTruth(10,0.,100.)
  {
    BOOST_MESSAGE( "Create RooUnfoldResponseFixture" );
    //Initialize RooUnfoldResponse instance with same entries for testing;
    responseFilledWithSomeEntries=RooUnfoldResponse(40, -10.0, 10.0);
    FixedSeedRandom =TRandom(111);
    // Train with a Breit-Wigner, mean 0.3 and width 2.5.
    for (Int_t i=0; i<10; i++) {
      double xt = FixedSeedRandom.BreitWigner(0.3, 2.5);
      double x = xt + smear(xt); //introduce
                                 //bias and smear
      if (x!=returnCutDummy())
      responseFilledWithSomeEntries.Fill (x, xt);
    else
      responseFilledWithSomeEntries.Miss (xt);
    }
  }
  virtual ~RooUnfoldResponseFixture(){
    BOOST_MESSAGE( "Tear down RooUnfoldResponseFixture" );
  }

  /**
   * Gets and returns the requested histograms of the given RooUnfoldResponse object. Known histogram names are: measured, truth and response.
   * rooUnfoldResponse: a RooUnfoldResponse object containing the histograms
   * histograms: required histogram names
   * Return value: vector with the requested histograms
   */
  static std::vector<const TH1*> GetRooUnfoldHistograms(const RooUnfoldResponse& rooUnfoldResponse,TString histograms){
     std::vector<const TH1*> result;

     if(histograms.Contains("measured")){
       const TH1* measuredHistogram = rooUnfoldResponse.Hmeasured();
       result.push_back(measuredHistogram);
     }

     if(histograms.Contains("truth")){
       const TH1* truthHistogram = rooUnfoldResponse.Htruth();
       result.push_back(truthHistogram);
     }

     if(histograms.Contains("response")){
       const TH2* responseHistogram = rooUnfoldResponse.Hresponse();  
       result.push_back(responseHistogram);
     }

    return result;
  }

  /**
   * This method tests the number of total entries for the measured histogram, the truth histogram and the response histogram.
   * rooUnfoldResponse: a RooUnfoldResponse object containing the three histograms
   * entriesExpected: expected number of total entries
   */
  static void testTotalEntriesOfHistograms(const RooUnfoldResponse& rooUnfoldResponse, const int&  entriesExpected){
    if(entriesExpected<0){
      BOOST_MESSAGE("Number of total entries of a histogram has to be positive! (entries: "<< entriesExpected<< ")");
      return;
    }
    //Get three histograms from RooUnfoldResponse object
    std::vector<const TH1*> histograms =  GetRooUnfoldHistograms(rooUnfoldResponse,"measured truth response");

    //test the total entries number for all histograms
    for(unsigned int i=0; i<histograms.size(); i++){
      const TH1* histogram = histograms.at(i);
      int entries = histogram->GetEntries();
      BOOST_CHECK_MESSAGE(entriesExpected==entries,"Wrong number of total entries for histogram "<< histogram->GetName() <<". Expected " << entriesExpected<<", but was: "<<entries);
    }
  }

  /**
   * This method tests the number of weighted entries for the measured histogram, the truth histogram and the response histogram.
   * rooUnfoldResponse: a RooUnfoldResponse object containing the three histograms
   * weightedEntriesExpected: expected number of weighted entries
   */
  static void testWeightedEntriesOfHistograms(const RooUnfoldResponse& rooUnfoldResponse, const double&  weightedEntriesExpected){
    //Get three histograms from RooUnfoldResponse object
    std::vector<const TH1*> histograms =  GetRooUnfoldHistograms(rooUnfoldResponse,"measured truth response");

    //test the total entries number for all histograms
    for(unsigned int i=0; i<histograms.size(); i++){
      const TH1* histogram = histograms.at(i);
      double entries = histogram->Integral();
      BOOST_CHECK_CLOSE(weightedEntriesExpected,entries,0.0001);
    }
  }

  /**
   * Tests if the bin content of the given histogram and bin is correct.
   * histogram: histogram which bin content is to be checked
   * bin: bin of histogram to be checked
   * expectedBinContent: expected content of bin to be checked
   */
  static void testBinContent(const TH1& histogram,int bin, double expectedBinContent){
    double binContent = histogram.GetBinContent(bin);
    BOOST_CHECK_CLOSE(expectedBinContent,binContent,0.0001);
  }


  RooUnfoldResponse response;
  RooUnfoldResponse responseSameBinsMeasuredTruth;
  RooUnfoldResponse responseFilledWithSomeEntries;
  

private:
  TRandom FixedSeedRandom;
  static const Double_t cutdummy_= -99999.0;
  Double_t smear (Double_t xt)
  {
    Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
    Double_t x= FixedSeedRandom.Rndm();
    if (x>xeff) return cutdummy_;
    Double_t xsmear= FixedSeedRandom.Gaus(-2.5,0.2);     // bias and smear
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

BOOST_AUTO_TEST_CASE(testmethodMiss1D){
  int numberOfBins = 3;
  double low = 0;
  double high = 3;
  RooUnfoldResponse responseWithNumberOfBins(numberOfBins, low, high);
  //define the pointers to the histograms
  const TH1* measured = responseWithNumberOfBins.Hmeasured();
  const TH1* fakes = responseWithNumberOfBins.Hfakes();
  const TH1* truth = responseWithNumberOfBins.Htruth();
  const TH2* response = responseWithNumberOfBins.Hresponse();
  const TH2D* responsenooverflow = responseWithNumberOfBins.HresponseNoOverflow();
  //sanity test: check if we are testing the correct dimension
  BOOST_CHECK_MESSAGE(1 == responseWithNumberOfBins.GetDimensionMeasured() , "Dimension measured is wrong : " << responseWithNumberOfBins.GetDimensionMeasured() << " != 1");
  BOOST_CHECK_MESSAGE(1 == responseWithNumberOfBins.GetDimensionTruth() , "Dimension truth is wrong : " << responseWithNumberOfBins.GetDimensionTruth() << " != 1");
  //test if the histograms are empty
  BOOST_CHECK_MESSAGE(0.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram  has no-zero entries before doing anything! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(0.0 == truth->GetEntries(), "truth histogram is filled. Number of entries found: " << truth->GetEntries() << " != 0");
  //check all histograms after the first call of Miss()
  responseWithNumberOfBins.Miss(1.5);
  BOOST_CHECK_MESSAGE(0.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram  has no-zero entries while using Miss() method! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(1.0 == truth->GetBinContent(2), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(2) << " != 1.0");
  BOOST_CHECK_MESSAGE(1.0 == truth->GetEntries(), "truth histogram more than once. Number of entries found: " << truth->GetEntries() << " != 1.0");
  //check all histograms after the second call of Miss()
  responseWithNumberOfBins.Miss(0.5);
  BOOST_CHECK_MESSAGE(0.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram  has no-zero entries while using Miss() method! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(1.0 == truth->GetBinContent(2), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(2) << " != 1");
  BOOST_CHECK_MESSAGE(1.0 == truth->GetBinContent(1), "truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(1) << " != 1");
  BOOST_CHECK_MESSAGE(2.0 == truth->GetEntries(), "truth histogram more than twice. Number of entries found: " << truth->GetEntries() << " != 2.0");
  //check all histograms after the 3rd call of Miss()
  responseWithNumberOfBins.Miss(1.5);
  BOOST_CHECK_MESSAGE(0.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram  has no-zero entries while using Miss() method! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(2 == truth->GetBinContent(2), "Truth histogram not filled with two entries. Number of entries found: " << truth->GetBinContent(2) << " != 2");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(1), "Truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(1) << " != 1");
  BOOST_CHECK_MESSAGE(0 == truth->GetBinContent(3), "Truth histogram is filled. Number of entries found: " << truth->GetBinContent(3) << " != 0");
  BOOST_CHECK_MESSAGE(3.0 == truth->GetEntries(), "truth histogram more than three times. Number of entries found: " << truth->GetEntries() << " != 3.0");
  //check all histograms after the 4th call of Miss()
  responseWithNumberOfBins.Miss(2.5,5);
  BOOST_CHECK_MESSAGE(0.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram  has no-zero entries while using Miss() method! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(2 == truth->GetBinContent(2), "Truth histogram not filled with two entries. Number of entries found: " << truth->GetBinContent(2) << " != 2");
  BOOST_CHECK_MESSAGE(1 == truth->GetBinContent(1), "Truth histogram not filled with one entry. Number of entries found: " << truth->GetBinContent(1) << " != 1");
  BOOST_CHECK_MESSAGE(5 == truth->GetBinContent(3), "Truth histogram not filled with five entries. Number of entries found: " << truth->GetBinContent(3) << " != 5");
  BOOST_CHECK_MESSAGE(4.0 == truth->GetEntries(), "truth histogram more than four times. Number of entries found: " << truth->GetEntries() << " != 4.0");
  }

BOOST_AUTO_TEST_CASE(testmethodFake1D){
  int numberOfBins = 3;
  double low = 0;
  double high = 3;
  RooUnfoldResponse responseWithNumberOfBins(numberOfBins, low, high);
  //define the pointers to the histograms
  const TH1* measured = responseWithNumberOfBins.Hmeasured();
  const TH1* fakes = responseWithNumberOfBins.Hfakes();
  const TH1* truth = responseWithNumberOfBins.Htruth();
  const TH2* response = responseWithNumberOfBins.Hresponse();
  const TH2D* responsenooverflow = responseWithNumberOfBins.HresponseNoOverflow();
  //sanity test: check if we are testing the correct dimension
  BOOST_CHECK_MESSAGE(1 == responseWithNumberOfBins.GetDimensionMeasured() , "Dimension measured is wrong : " << responseWithNumberOfBins.GetDimensionMeasured() << " != 1");
  BOOST_CHECK_MESSAGE(1 == responseWithNumberOfBins.GetDimensionTruth() , "Dimension truth is wrong : " << responseWithNumberOfBins.GetDimensionTruth() << " != 1");
  //test if the histograms are empty
  BOOST_CHECK_MESSAGE(0.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries found: " << responseWithNumberOfBins.FakeEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram has no-zero entries before doing anything! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(0.0 == truth->GetEntries(), "truth histogram is filled. Number of entries found: " << truth->GetEntries() << " != 0");
  //check all histograms after the first call of Fake()
  responseWithNumberOfBins.Fake(1.5);
  BOOST_CHECK_MESSAGE(1.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries is wrong: " << responseWithNumberOfBins.FakeEntries() << " != 1");
  BOOST_CHECK_MESSAGE(1.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 1");
  BOOST_CHECK_MESSAGE(1.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 1");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram has no-zero entries while using Fake() method! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(0.0 == truth->GetEntries(), "truth histogram is filled. Number of entries found: " << truth->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(1 == measured->GetBinContent(2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(2) << " != 1");
  BOOST_CHECK_MESSAGE(1 == fakes->GetBinContent(2), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(2) << " != 1");
  //check all histograms after the second call of Fake()
  responseWithNumberOfBins.Fake(0.5);
  BOOST_CHECK_MESSAGE(2.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries is wrong: " << responseWithNumberOfBins.FakeEntries() << " != 2");
  BOOST_CHECK_MESSAGE(2.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 2");
  BOOST_CHECK_MESSAGE(2.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 2");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram has no-zero entries while using Fake() method! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(0.0 == truth->GetEntries(), "truth histogram is filled. Number of entries found: " << truth->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(1 == measured->GetBinContent(1), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(1) << " != 1");
  BOOST_CHECK_MESSAGE(1.0 == measured->GetBinContent(2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(1) << " != 1.0");
  BOOST_CHECK_MESSAGE(1.0 == fakes->GetBinContent(1), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(1) << " != 1.0");
  BOOST_CHECK_MESSAGE(1.0 == fakes->GetBinContent(2), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(2) << " != 1.0");
  //check all histograms after the 3rd call of Fake()
  responseWithNumberOfBins.Fake(1.5);
  BOOST_CHECK_MESSAGE(3.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries is wrong: " << responseWithNumberOfBins.FakeEntries() << " != 3");
  BOOST_CHECK_MESSAGE(3.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 3");
  BOOST_CHECK_MESSAGE(3.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 3");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram has no-zero entries while using Fake() method! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(0.0 == truth->GetEntries(), "truth histogram is filled. Number of entries found: " << truth->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(1.0 == measured->GetBinContent(1), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(1) << " != 1.0");
  BOOST_CHECK_MESSAGE(2.0 == measured->GetBinContent(2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(1) << " != 2.0");
  BOOST_CHECK_MESSAGE(1.0 == fakes->GetBinContent(1), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(1) << " != 1.0");
  BOOST_CHECK_MESSAGE(2.0 == fakes->GetBinContent(2), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(2) << " != 2.0");
  //check all histograms after the 4th call of Fake()
  responseWithNumberOfBins.Fake(2.5,5);
  BOOST_CHECK_MESSAGE(4.0 == responseWithNumberOfBins.FakeEntries() , "Number of fake entries is wrong: " << responseWithNumberOfBins.FakeEntries() << " != 4");
  BOOST_CHECK_MESSAGE(4.0 == measured->GetEntries(), "Measured histogram is filled. Number of entries found: " << measured->GetEntries() << " != 4");
  BOOST_CHECK_MESSAGE(4.0 == fakes->GetEntries(), "fakes histogram is filled. Number of entries found: " << fakes->GetEntries() << " != 4");
  BOOST_CHECK_MESSAGE(0.0 == response->GetEntries(), "Response histogram is filled. Number of entries found: " << response->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(0.0 == responsenooverflow->GetEntries(), "Responsenooverflow histogram has no-zero entries while using Fake() method! Number of entries found: " << responsenooverflow->GetEntries() << " != 0 [this bug had been reported to the RooUnfold creator]");
  BOOST_CHECK_MESSAGE(0.0 == truth->GetEntries(), "truth histogram is filled. Number of entries found: " << truth->GetEntries() << " != 0");
  BOOST_CHECK_MESSAGE(1.0 == measured->GetBinContent(1), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(1) << " != 1.0");
  BOOST_CHECK_MESSAGE(2.0 == measured->GetBinContent(2), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(1) << " != 2.0");
  BOOST_CHECK_MESSAGE(5.0 == measured->GetBinContent(3), "measured histogram not filled with one entry. Number of entries found: " << measured->GetBinContent(3) << " != 5.0");
  BOOST_CHECK_MESSAGE(1.0 == fakes->GetBinContent(1), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(1) << " != 1.0");
  BOOST_CHECK_MESSAGE(2.0 == fakes->GetBinContent(2), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(2) << " != 2.0");
  BOOST_CHECK_MESSAGE(5.0 == fakes->GetBinContent(3), "fakes histogram not filled with one entry. Number of entries found: " << fakes->GetBinContent(3) << " != 5.0");
}

BOOST_AUTO_TEST_CASE(testFill1D){
  //test with weight
  double xMeasured = 42;
  double xTruth = 74;
  double  weight =7.3;

  responseSameBinsMeasuredTruth.Fill(xMeasured,xTruth,weight);

  //test if number of total and weighted entries are right for all histograms
   RooUnfoldResponseFixture::testTotalEntriesOfHistograms(responseSameBinsMeasuredTruth, 1);
   RooUnfoldResponseFixture::testWeightedEntriesOfHistograms(responseSameBinsMeasuredTruth,weight);

  TH1* measuredHistogram = responseSameBinsMeasuredTruth.Hmeasured();
  TH1* truthHistogram = responseSameBinsMeasuredTruth.Htruth();
  TH2* responseHistogram = responseSameBinsMeasuredTruth.Hresponse();
  


  //test if right bin was filled for measured histogram
  int binMeasured = measuredHistogram->FindBin(xMeasured);
  RooUnfoldResponseFixture::testBinContent(*measuredHistogram,binMeasured, weight);

  //test if right bin was filled for truth histogram
  int binTruth = truthHistogram->FindBin(xTruth);
  RooUnfoldResponseFixture::testBinContent(*truthHistogram,binTruth, weight);

  //test if right bin was filled for response histogram
  int binResponse = responseHistogram->FindBin(xMeasured,xTruth);
  RooUnfoldResponseFixture::testBinContent(*responseHistogram,binResponse, weight);

  //test default value for weight for filling
  responseSameBinsMeasuredTruth.Fill(xMeasured,xTruth);

  int entries = truthHistogram->GetEntries();
  BOOST_CHECK_MESSAGE(entries==2,"Wrong number of total entries. Expected 2, but was: "<<entries);
  double resultDefaultWeight = truthHistogram->Integral()-weight;
  BOOST_CHECK_CLOSE( 1, resultDefaultWeight, 0.0001 );
}

BOOST_AUTO_TEST_CASE(testFill2D){
  //initialising of histograms needed for constructor
  int measuredBinX = 10;
  double measuredLowX = 0;
  double measuredHighX = 10;

  int measuredBinY = 20;
  double measuredLowY = 10;
  double measuredHighY = 70;
  TH2F* measuredHistogram = new TH2F("measured2D","measured histogram 2D",measuredBinX,measuredLowX,measuredHighX,measuredBinY,measuredLowY,measuredHighY);

  int truthBinX = 10;
  double truthLowX = 0;
  double truthHighX = 100;

  int truthBinY = 30;
  double truthLowY = 0;
  double truthHighY = 90;
  TH2F* truthHistogram = new TH2F("truth2D","truth histogram 2D",truthBinX,truthLowX,truthHighX,truthBinY,truthLowY,truthHighY);

  RooUnfoldResponse response2D(measuredHistogram,truthHistogram);

  //testing 2D Fill
  double xMeasured = 3.4;
  double yMeasured = 42;
  double xTruth = 27;
  double yTruth = 67;
  double weight = 2.1;

  int maxBin = measuredBinX* measuredBinY*truthBinX*truthBinY;
  int bin = response2D.Fill(xMeasured,yMeasured,xTruth,yTruth,weight);

  //test if number of total and weighted entries are right for all histograms
  RooUnfoldResponseFixture::testTotalEntriesOfHistograms(response2D, 1);
  RooUnfoldResponseFixture::testWeightedEntriesOfHistograms(response2D,weight);

  TH1* resultMeasuredHistogram =response2D.Hmeasured();
  TH1* resultTruthHistogram =response2D.Htruth();
  TH2* responseHistogram = response2D.Hresponse();

    //test if right bin was filled for measured histogram (error in test)
  int binMeasured = measuredHistogram->FindBin(xMeasured,yMeasured);
  RooUnfoldResponseFixture::testBinContent(*measuredHistogram,binMeasured, weight);

  //test if right bin was filled for truth histogram (error in test)
  int binTruth = truthHistogram->FindBin(xTruth,yTruth);
  RooUnfoldResponseFixture::testBinContent(*truthHistogram,binTruth, weight);

  //test if right bin was filled for response histogram
  int binX = RooUnfoldResponse::FindBin(measuredHistogram, xMeasured, yMeasured);
  int binY = RooUnfoldResponse::FindBin(truthHistogram, xTruth, yTruth);
  int binResponse = responseHistogram->FindBin(binX,binY);
  RooUnfoldResponseFixture::testBinContent(*responseHistogram,binResponse, weight);


  //test dimension also test with wrong dimensions

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
  testObject.Add(responseFilledWithSomeEntries);
  BOOST_CHECK_MESSAGE(testObject.Hmeasured()->GetEntries()==2*noOfEntriesInMeasuredHist,"Adding the same RooUnfoldResponse did not result in twice the number of entries in Hmeasured");

  //should be extended :)
}



BOOST_AUTO_TEST_SUITE_END()
