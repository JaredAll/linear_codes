#ifndef NOISY_CHANNEL_H
#define NOISY_CHANNEL_H


#include <climits>
#include <iostream>
#include <vector>
#include <cfloat>
#include <stdlib.h>
#include <time.h>

using namespace std;

/*
 * introduces noise "randomly" into the message.
 * @param message the message to be sent
 * @param code_length the length of the code
 * @param errors_per_word the number of errors 
 * randomly introduced into each word.
 */
void random_noise( vector< uint > &message,
                   uint code_length,
                   uint errors_per_word );

uint find_power( uint base, uint exponent );

/*
 * introduces burst noise into the message, within each word.
 * @param message the message to be sent
 * @param code_length the length of the code
 */
void burst_noise( vector< uint > &message, uint code_length );

/*
uint find_power( uint base, uint exponent )
{
  uint result = 1;
  for( uint i = 0; i < exponent; i++ )
  {
    result *= base;
  }
  return result;
}
*/


void random_noise( vector< uint > &message,
                   uint code_length,
                   uint errors_per_word )
{

  srand( time( NULL ) );

  for( uint j = 0; j < errors_per_word; j++ )
  {
    for( uint i = 0; i < message.size(); i++ )
    {
      uint noise_pv = rand() % code_length;
      uint noise = find_power( 2, noise_pv );
      message.at( i ) ^= noise;
    }
  }
  
}

void burst_noise( vector< uint > &message, uint code_length )
{
  srand( time( NULL ) );

  uint num_bursts = rand() % message.size();
  vector< uint > burst_positions;

  //determine the words in which the bursts will be.
  //allow only 1 burst per word.
  uint placed_bursts = 0;
  while( placed_bursts < num_bursts )
  {
    uint burst_word = rand() % message.size();
    bool word_has_burst = false;
    
    for( uint position = 0; position < burst_positions.size(); position++ )
    {
      if( burst_positions.at( position ) == burst_word )
      {
        word_has_burst = true;
      }
    }

    if( !word_has_burst )
    {
      burst_positions.push_back( burst_word );
      placed_bursts++;
    }
  }    
    

  //place the bursts in the appropriate words.
  for( uint i = 0; i < num_bursts; i++ )
  {
    //determine the burst characteristics with random size <=
    //half the length of the code and determine where in
    //the word it starts
    uint burst_size = rand() % ( code_length / 2 );
    uint burst_start_pv = rand() % ( code_length - burst_size );

    //create the burst
    uint burst = 0;
    for( uint j = 0; j < burst_size; j++ )
    {
      burst += find_power( 2, j );
    }
    burst = burst << burst_start_pv;
    message.at( burst_positions.at( i ) ) ^= burst;
  }
}

#endif
