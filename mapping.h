#ifndef ALPHA_MAP_H
#define ALPHA_MAP_H

#include <cstdint>
#include <iostream>
#include <vector>
#include <climits>
#include <iterator> 
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include "linear_code.h"

using namespace std;

uint ALPHABET_SIZE = 26;
uint ALPHABET_BASE = 321;

/**
 * A class to map from a linear code to the alphabet
 * @author Jared Allen
 * @version 9 January 2019
 */
class AlphabetMap
{
public:
  /**
   * Constructor specifying map from original message to encoded
   * message
   * @param encoding the codewords from the linear code
   */
  AlphabetMap( vector< uint > encoding );

  /**
   * a function to convert a message in numbers to a message in
   * letters.
   * @param message the message to be converted
   * @return the message in letters
   */
  vector< char > convert_to_letters( vector< uint > message );

  /**
   * a function to convert a message in letters to a message in
   * numbers.
   * @param message the message to be converted
   * @return the message in numbers
   */
  vector< uint > convert_to_numbers( vector< char > message );

  /**
   * a function to print the letters in the alphabet mapping
   */
  void print_alphabet();

private:

  map< uint, char > encoding_map;
  map< char, uint > rev_encoding_map;

};

AlphabetMap::AlphabetMap( vector< uint > encoding )
{
  //initiate the mappings
  for( int i = 0; i < encoding.size(); i++ )
  {
    char letter = static_cast< char >( i + ALPHABET_BASE );
    encoding_map.insert( { encoding.at( i ), letter } );
    rev_encoding_map.insert( { letter, encoding.at( i ) } );
  }
}

void AlphabetMap::print_alphabet()
{
  map< uint, char >::iterator it;
  
  cout << "the letters in the alphabet: " << endl;
  for( it = encoding_map.begin(); it != encoding_map.end(); it++ )
  {
    cout << it -> second << " ";
  }
  cout << endl;
  cout << endl;
}

vector< char > AlphabetMap::convert_to_letters(
  vector< uint > message )
{
  //use mapping for conversion
  vector< char > converted_message;
  for( uint i = 0; i < message.size(); i++ )
  {
    if( encoding_map.count( message.at( i ) ) == 1 )
    {
      converted_message.push_back(
        encoding_map.find( message.at( i ) ) -> second );
    }
    else
    {
      converted_message.push_back( '|' );
      /*
      converted_message.push_back(
        static_cast< char >( message.at( i ) ) );
      */
    }
  }
  
  return converted_message; 
}

vector< uint > AlphabetMap::convert_to_numbers(
  vector< char > message )
{
  //use reverse mapping for conversion
  vector< uint > converted_message;
  for( uint i = 0; i < message.size(); i++ )
  {
    converted_message.push_back(
      rev_encoding_map.find( message.at( i ) ) -> second );
  }
  
  return converted_message;
}

#endif
