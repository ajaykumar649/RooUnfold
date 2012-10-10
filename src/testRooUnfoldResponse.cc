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
  RooUnfoldResponseFixture(){
    BOOST_MESSAGE( "Create RooUnfoldResponseFixture" );
  }
  virtual ~RooUnfoldResponseFixture(){
    BOOST_MESSAGE( "Tear down RooUnfoldResponseFixture" );
  }
  RooUnfoldResponse response;
};


// Declare test suite name and fixture class to BOOST:
BOOST_FIXTURE_TEST_SUITE( RooUnfoldResponseSuite, RooUnfoldResponseFixture )
//BOOST_AUTO_TEST_SUITE( RooUnfoldResponseSuite )

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
}

//Test of UseOverflowStatus
BOOST_AUTO_TEST_CASE(testUseOverflowStatus){
  //  RooUnfoldResponse testObject;
  BOOST_CHECK_MESSAGE(response.UseOverflowStatus()==false,"default constructor does not initialize with overflow set to false");
}

BOOST_AUTO_TEST_SUITE_END()
