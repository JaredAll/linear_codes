/* A program to model a linear code over a finite field 
 * F of dimension n and order q, where n is variable and q is
 * restricted to 2. 
 * @author Jared Allen
 * @date November 21, 2018
 */ 


#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <algorithm>
#include "linear_code.h"
#include "noisy_channel.h"
#include "mapping.h"

using namespace std;

/**
 * A function to determine the index of a uint in a vector
 * @param vector the vector
 * @param number the number we are looking for the index of
 * @return the index 
 */
uint index_of( vector< uint > vector, uint number );

/* A function to implement the Gauss-Jordan algorithm, which
 * puts a code matrix in rref
 * @param code matrix the matrix to be put in rref
 * @param code_size the number of columns in the matrix
 * @return the rref matrx
 */
vector< uint > find_rref( vector< uint > code_matrix,
                          uint code_size);

/* A function to find and print the code matrix
 * @param code_matrix the code matrix
 * @param code_length the length of the code words
 */
void print_bitwise( vector< uint > code_matrix,
                       uint code_length );

/* A function to determine the transpose of a code matrix
 * @param code_matrix the code_matrix
 * @param code_length the length of the code words
 * @return the transpose of the matrix
 */
vector< uint > find_transpose( vector< uint > code_matrix,
                               uint code_length );

/* A function to find a simple exponent
 * @param base the base
 * @param exponent the exponent
 * @return the result
 */
uint find_power( uint base, uint exponent );

/* A function to permute the columns in the
 * rref matrix as in algorithm 4.3 of
 * Ling and Xing.
 * @param code_matrix the matrix to be permuted
 * @param code_length the length of the code
 * @return the new matrix 
 */
vector< uint > permute_columns( vector< uint > &code_matrix,
                                uint code_length, vector< uint > permutation );

/* A function to determine the parity check matrix
 * for a code.
 * @param g_permuted G' as in Ling & Xing 4.3
 * @param code_length the length of the code
 * @param permutation the permutation to be reversed
 * @return the parity check matrix
 */
vector< uint > find_pc_matrix( vector< uint > g_permuted,
                               vector< uint > rref_matrix,
                               uint code_length,
                               vector< uint > permutation );

/* A function to determine the necessary permutation
 * to create G', as described in algorithm 4.3.
 * @param code_matrix the code matrix in rref
 * @param code_length the length of the code
 * @return the necessary permutation
 */
vector< uint > find_permutation( vector< uint > code_matrix,
                                 uint code_length );

/* A function to determine if a code matrix is the 
 * identity matrix.
 * @param code_matrix the matrix 
 * @param code_length the length of the code
 * @return whether or not it is the identity
 */
bool is_identity( vector< uint > code_matrix,
                  uint code_length );

/* A function to determine G, the k x n matrix
 * consisting of all nonzero rows of the rref
 * @param code_matrix the rref of the matrix
 * @param code_length the length of the code
 * @return the matrix G
 */
vector< uint > find_g_matrix( vector< uint > code_matrix,
                              uint code_length );



uint index_of( vector< uint > vector, uint number )
{
  ptrdiff_t pos =
    find( vector.begin(), vector.end(), number ) - vector.begin();
  
  if( pos >= vector.size() )
  {
    cout << "no such number in vector." << endl;
    return 0;
  }

  return pos;
}



bool is_identity( vector< uint > code_matrix,
                  uint code_length )
{
  bool the_identity = true;
  //if the matrix isn't square, its not the identity
  if( code_matrix.size() != code_length )
  {
    the_identity = false;
    return the_identity;
  }

  //check to see if each entry is the corresponding power of 2.
  for( uint i = 0; i < code_matrix.size(); i++ )
  {
    uint place_value = find_power( 2, code_length - ( i + 1 ) );
    if( code_matrix.at( i ) != place_value )
    {
      the_identity = false;
    }
  }
  return the_identity;
}


vector< uint > find_g_matrix( vector< uint > code_matrix,
                              uint code_length )
{
  //determine G, the k x n matrix of nonzero rows. 
  vector< uint > g_matrix;
  for( uint i = 0; i < code_matrix.size(); i++ )
  {
    if( code_matrix.at( i ) > 0 )
    {
      //eliminate bits in place values greater than
      //code length
      uint mod_rep;
      mod_rep = code_matrix.at( i ) & ( find_power( 2, code_length ) - 1 );
      if( mod_rep > 0 )
      {
        g_matrix.push_back( mod_rep );
      }
    }
  }
  return g_matrix;
}

vector< uint > find_permutation( vector< uint > g_matrix,
                                 uint code_length )
{
  //find transpose of matrix
  vector< uint > g_transpose = find_transpose( g_matrix,
                                                code_length );

  
  //determine necessary permutation to put into (I|X) form,
  //where I is an identity matrix
  vector< uint > permutation;
  uint place_value = g_matrix.size() - 1;
  while( place_value != UINT_MAX )
  {
    bool first_instance = true;
    for( uint i = 0; i < g_transpose.size(); i++ )
    {
      if( g_transpose.at( i ) == find_power( 2, place_value )
          && first_instance )
      {
        permutation.push_back( i );
        first_instance = false;
      }
    }
    place_value--;
  }

  for( uint i = 0; i < g_transpose.size(); i++ )
  {
    if( count( permutation.begin(), permutation.end(), i ) == 0 )
    {
      permutation.push_back( i );
    }
  }
  
  return permutation;
}


vector< uint > find_pc_matrix( vector< uint > g_permuted,
                               vector< uint > rref_matrix,
                               uint code_length,
                               vector< uint > permutation )
{
  //find transposes of appropriate matrices
  vector< uint > gp_transpose =
    find_transpose( g_permuted, code_length  );
  vector< uint > rref_transpose = find_transpose( rref_matrix,
                                                  code_length );

  vector< uint > x_matrix;
  //determine the size of the x matrix
  
  for( uint i = g_permuted.size(); i < code_length; i++ )
  {
    x_matrix.push_back( gp_transpose.at( i ) );
  }
  

  cout << "the matrix X^T " << endl;
  print_bitwise( x_matrix, g_permuted.size());

  //create H', which is (-X^T | I_( n - k ) )
  uint num_redundant = code_length - g_permuted.size();
  for( uint i = 0; i < x_matrix.size(); i++ )
  {
    x_matrix.at( i ) = x_matrix.at( i ) << num_redundant;
    x_matrix.at( i ) += find_power( 2, num_redundant - ( i + 1 ) );
  }

  vector< uint > hp_matrix = x_matrix;

    cout << "the H' form of the matrix" << endl;
  print_bitwise( hp_matrix, code_length );

  //transform the hp matrix into the h matrix
  //determine reverse permutation from rref_transpose and
  //gp_transpose
  vector< uint > reverse_permutation;
  for( uint i = 0; i < permutation.size(); i++ )
  {
    if( permutation.at( i ) != i )
    {
      reverse_permutation.push_back( index_of( permutation, i ) );
    }
    else
    {
      reverse_permutation.push_back( i );
    }
  }

  /* output permutations */

  cout << "permutation: " << endl;
  for( uint col : permutation )
  {
    cout << col << " ";
  }
  cout << endl;

  for( uint i = 0; i < permutation.size(); i++ )
  {
    cout << i << " ";
  }
  cout << endl;
  cout << endl;

  cout << "reverse permutation: " << endl;
  for( uint col : reverse_permutation )
  {
    cout << col << " ";
  }
  cout << endl;

  for( uint i = 0; i < permutation.size(); i++ )
  {
    cout << i << " ";
  }
  cout << endl;
  cout << endl;

  /* end output permutations */
    
      
  permute_columns( hp_matrix,
                   code_length, reverse_permutation );  

  return hp_matrix;
}
  


vector< uint > permute_columns( vector< uint > &code_matrix,
                                uint code_length,
                                vector< uint > permutation )
{
  //find transpose of code matrix
  vector< uint > cm_transpose = find_transpose(
    code_matrix, code_length );
  
  //perform the permutation
  vector< uint > new_matrix;
  for( uint i = 0; i < permutation.size(); i++ )
  {
    if( permutation.at( i ) != i )
    {
      new_matrix.push_back( cm_transpose.at( permutation.at( i ) ) );
    }
    else
    {
      new_matrix.push_back( cm_transpose.at( i ) );
    }
  }
  
  cm_transpose = new_matrix;
  code_matrix = find_transpose( cm_transpose,
                                code_matrix.size() );
      
  return permutation;
}

uint find_power( uint base, uint exponent )
{
  if( exponent == 0 )
  {
    return 1;
  }
  else
  {
    uint result = 1;
    for( uint i = 0; i < exponent; i++ )
    {
      result *= base;
    }
    return result;
  }
}

vector< uint > find_transpose( vector< uint > code_matrix,
                               uint code_length )
{
  vector< uint > code_matrix_transpose;
  uint new_code_length = 0;

  //initialize the transpose matrix
  for( uint col = code_length; col != 0; col-- )
  {
    code_matrix_transpose.push_back( 0 );
  }

  //establish the new code length
  uint num_code_words = code_matrix.size();

  //add the appropriate values to create the
  //transpose matrix and return it
  for( uint old_row = 0; old_row < code_matrix.size(); old_row++ )
  {
    new_code_length++;
    uint this_old_row = code_matrix.at( old_row );
    for( uint col = code_length - 1; col != UINT_MAX; col-- )
    {
      uint new_row = code_length - col - 1;
      uint place_holder = ( this_old_row >> col ) & 1;
      uint place_value = num_code_words - old_row - 1;
      if( place_holder == 1 )
      {
        code_matrix_transpose.at( new_row ) +=
          find_power( 2, place_value );
      }
    }
  }

  return code_matrix_transpose;
}
      

void print_bitwise( vector< uint > code_matrix, uint code_length )
{
  //find and print the bitwise representation of code_matrix
  for( uint i = 0; i < code_matrix.size(); i++ )
  {
    uint this_code_word = code_matrix.at( i );
    vector< uint > code_word_bitwise;
    for( uint j = 0; j < code_length; j++ )
    {
      uint this_bit = ( this_code_word >> j ) & 1;
      code_word_bitwise.push_back( this_bit );
    }

    for( uint j = code_length - 1; j < UINT_MAX; j-- )
    {
      cout << code_word_bitwise.at( j );
    }
    cout << endl;
  }
  cout << endl;
}


vector< uint > find_rref( vector< uint > code_matrix, uint code_size)
{
  uint row = 0;
  uint col_offset = 1;
  bool in_rref = false;

  //do the Gauss-Jordan algorithm until matrix in rref
  while( !in_rref )
  {
    uint col = code_size - col_offset;
    uint pivot = ( code_matrix.at( row ) >> col ) & 1;
    uint next_row = row;
    while( pivot == 0 and col != UINT_MAX)
    {
      while( pivot == 0 and next_row < code_matrix.size() )
      {
        next_row++;
        if( next_row < code_matrix.size() )
        {
          pivot = ( code_matrix.at( next_row ) >> col ) & 1;
        }
      }
      //if pivot is 0, all entries in col are 0, so increment col
      //and try again.
      if( pivot == 0 or next_row > code_matrix.size() )
      {
        next_row = row;
        col_offset++;
        col = code_size - col_offset;
        pivot = ( code_matrix.at( next_row ) >> col ) & 1;
      }
    }

    //if we havent processed the last row or column, proceed to
    //switch rows. Otherwise, the remaining rows are linear
    //combinations of the first rows, so zero them out.
    
    if( next_row < code_matrix.size() and col != UINT_MAX )
    {
      //switch row with the next row that has a nonzero pivot
    
      uint first_row_switch = code_matrix.at( row );
      uint second_row_switch = code_matrix.at( next_row );
      vector< uint > new_matrix;
      for( uint i = 0; i < code_matrix.size(); i++ )
      {
        if( i != row and i != next_row )
        {
          new_matrix.push_back( code_matrix.at( i ) );
        }
        else if( i == row )
        {
          new_matrix.push_back( second_row_switch );
        }
        else
        {
          new_matrix.push_back( first_row_switch );
        }
      }
      code_matrix = new_matrix;
      
      //eliminate all other entries in the j-th column
      for( uint i = 0; i < code_matrix.size(); i++ )
      {        
        if( i != row )
        {
          uint this_col = ( code_matrix.at( i ) >> col ) & 1;
          if( this_col == 1 )
          {
            code_matrix.at( i ) = code_matrix.at( row ) ^
              code_matrix.at( i );
          }
        }
      }
      row++;
      col_offset++;
    }
    else
    {
      if( row < code_matrix.size() )
      {
        for( uint i = row + 1; i < code_matrix.size(); i++ )
        {
          code_matrix.at( i ) = 0;
        }
      }
      in_rref = true;
    }

    if( row == code_matrix.size() )
    {
      in_rref = true;
    }

    if( col == UINT_MAX )
    {
      in_rref = true;
    }
  }
  return code_matrix;
}


int main()
{
  //get S, a nonempty subset of F
  vector< uint > subset_of_F;
  uint subset_element;
  uint code_length;
  cin >> code_length;
  while( cin >> subset_element )
  {
    subset_of_F.push_back( subset_element );
  }

  cout << "The subset of the field F" << endl;
  print_bitwise( subset_of_F, code_length );

  //determine the rref form of the subset of F.
  //the nonzero rows become the generator matrix for the linear code
  vector< uint > matrix_rref = find_rref( subset_of_F, code_length );
  cout << "the rref of the matrix above" << endl;
  print_bitwise( matrix_rref, code_length );

  
  //determine the nonzero rows of the rref matrix
  vector< uint > g_matrix = find_g_matrix( matrix_rref,
                                           code_length);

  //determine permutation for parity check matrix if applicable
  //and test codes
  if( !is_identity( g_matrix, code_length ) )
  {
    
    vector< uint > permutation = find_permutation(
      g_matrix, code_length );
    vector< uint > g_permuted = g_matrix;
    permute_columns( g_permuted, code_length, permutation );

    cout << "the G' form of the rref matrix" << endl;
    print_bitwise( g_permuted, code_length );
    

    
    //determine the parity check matrix, the basis for C_perp
    vector< uint > parity_check_matrix = find_pc_matrix(
      g_permuted, matrix_rref, code_length, permutation );
  
    cout << "the basis of the dual code, C_perp" << endl;
    print_bitwise( parity_check_matrix, code_length );

    
    //create linear code object
    LinearCode this_code = LinearCode( g_matrix,
                                       parity_check_matrix,
                                       code_length );
    //print linear code information
    this_code.print_generator();
    this_code.print_parity_check();
    this_code.print_words();

    //determine the map between words and encoded words
    vector< uint > encoded_words;
    uint num_words =
      find_power( 2, this_code.get_generator().size() );

    for( uint word = 0; word < num_words; word++ )
    {
      encoded_words.push_back( this_code.encode_word( word ) );
    }
    
    AlphabetMap map = AlphabetMap( encoded_words );
    map.print_alphabet();

    
    
    //start tests
    //-------------------------------------------------------

    //read in message from file
    char letter;
    ifstream message_file;
    vector< char > message;
    message_file.open( "alice_message.txt" );
    if( message_file.is_open() )
    {
      while( message_file.get( letter ) )
      {
        if( letter != '\n' )
        {
          message.push_back( letter );
        }
      }
      message_file.close();
    }

    
    //print original message
    cout << "the original message: " << endl;
    for( char letter : message )
    {
      cout << letter;
    }
    cout << endl;
    cout << endl;

    
    //convert message to encoded words
    vector< uint > encoded_message =
      map.convert_to_numbers( message );

    vector< char > og_message = message;

    
    //introduce noise into message
    uint num_errors = 3;
    random_noise( encoded_message,
                  this_code.get_code_length(), num_errors );

    cout << "number of errors per \"word\": " << num_errors << endl;
    cout << endl;

    
    //output the corrupted message
    message = map.convert_to_letters( encoded_message );

    cout << "the received message: " << endl;
    for( char letter : message )
    {
      cout << letter;
    }
    cout << endl;
    cout << endl;

    
    //extract message from received message
    vector< uint > decoded_message;
    for( uint word : encoded_message )
    {
      decoded_message.push_back( this_code.decode_word( word ) );
    }

    vector< char > char_d_message = map.convert_to_letters( decoded_message );
    
    cout << "the decoded received message: " << endl;
    for( char letter : char_d_message )
    {
      cout << letter;
    }
    cout << endl;
    cout << endl;

    
    //determine accuracy
    float letters_identical = 0;
    for( uint i = 0; i < og_message.size(); i++ )
    {
      if( og_message.at( i ) == char_d_message.at( i ) )
      {
        letters_identical++;
      }
    }
    cout << "percent identity: " <<
      ( letters_identical / og_message.size() ) * 100 << endl;

    //-------------------------------------------------------
    //end tests
  }
  else
  {
    cout << "The matrix G is I" << code_length << endl;
  }


}
