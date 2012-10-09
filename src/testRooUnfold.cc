
// Unit tests for RooUnfold
//

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex>

// BOOST test stuff:
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE roounfoldtest
#include <boost/test/unit_test.hpp>

// Namespaces:
using std::string;
using std::vector;
using std::complex;
//using namespace INIParser;

// Test fixture for all tests:
class RooUnfoldTestFixture {
public:
  RooUnfoldTestFixture(){
    BOOST_MESSAGE( "Create RooUnfoldTestFixture" );
  }
  virtual ~RooUnfoldTestFixture() {
    BOOST_MESSAGE( "Tear down RooUnfoldTestFixture" );
  }
  //  RooUnfold reader;
};

// Declare test suite name and fixture class to BOOST:
//BOOST_FIXTURE_TEST_SUITE( roounfoldsuite, RooUnfoldTestFixture )

// Test cases:

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

// BOOST_AUTO_TEST_SUITE_END()

