// Unit tests for RooUnfoldResponse class
// A. Vanhoefer, E. Schlieckau, 10/2012

#include "RooUnfoldResponse.h"

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
BOOST_AUTO_TEST_CASE(testFails){
  BOOST_FAIL("This test fails!");
}

//Test of UseOverflowStatus
BOOST_AUTO_TEST_CASE(testUseOverflowStatus){
  //  RooUnfoldResponse testObject;
  BOOST_CHECK_MESSAGE(response.UseOverflowStatus()==false,"default constructor does not initialize with overflow set to false");
}

BOOST_AUTO_TEST_SUITE_END()
