/// \file slim_codec_runlength.cpp
/// Implement classes encoder_runlength and decoder_runlength,
/// used for data that have long stings of repeated values.

//  Copyright (C) 2008, 2009 Joseph Fowler
//
//  This file is part of slim, a compression package for science data.
//
//  Slim is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Slim is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with slim.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <cmath>

#include "slim.h"
#include "slim_single_codec.h"

//#define DEBUG_ENCODING


//----------------------------------------------------------------------
/// \class encoder_runlength
/// Derived class for encoding highly repetitive streams with runlength.
/// That is, as (value, repeats) pairs.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// Construct (optionally by output bitstream).
/// \param dt  The data type code.
/// \param deltas Whether to encode successive differences.
/// \param ob  The output bitstream.
encoder_runlength::encoder_runlength(enum data_t dt, bool deltas, obitstream *ob) 
  : encoder(dt, deltas, ob) {
  set_data_type(SLIM_TYPE_U32);
}



//----------------------------------------------------------------------
/// Destructor is simply the base destructor.
encoder_runlength::~encoder_runlength() {;}



//----------------------------------------------------------------------
/// Encode one datum to the output stream.
/// \param datum  The word to be encoded.
void encoder_runlength::encode(uint32_t datum) const {
#ifdef DEBUG_ENCODING
  cout << "encoder_runlength: encoding single " << datum << endl;
#endif

  mexp_golomb_write(out_bs, datum);  // system uses default order=1
  mexp_golomb_write(out_bs, 1);
}



//----------------------------------------------------------------------
/// Encode an array of data to the output stream.
/// \param data    Array of words to be encoded.
/// \param ndata   Number of words to be encoded.
void encoder_runlength::encode_vector(const uint32_t *data, 
				      int ndata) {
  uint32_t repeated_value;  ///< Value seen before and not yet recorded.
  uint32_t value_counter;   ///< Count of consecutive appearances of value.

  if (use_deltas) {
    repeated_value = data[0]-prev_datum;
    prev_datum = repeated_value;
    value_counter = 1;
    for (int i=1; i<ndata; i++) {
      if (data[i]-prev_datum == repeated_value)
	value_counter ++;
      else {
#ifdef DEBUG_ENCODING
	cout << "encoder_runlength DELTA: encoding " << value_counter
             << " repeats of " << repeated_value << endl;
#endif
	
	mexp_golomb_write(out_bs, repeated_value);
	mexp_golomb_write(out_bs, value_counter);
	repeated_value = data[i]-prev_datum;
	value_counter = 1;
      }
      prev_datum = data[i];
    }
  } else {
    repeated_value = data[0];
    value_counter = 1;
    for (int i=1; i<ndata; i++) {
      if (data[i] == repeated_value)
	value_counter ++;
      else {

#ifdef DEBUG_ENCODING
	cout << "encoder_runlength: encoding " << value_counter
             << " repeats of " << repeated_value << endl;
#endif

	mexp_golomb_write(out_bs, repeated_value);
	mexp_golomb_write(out_bs, value_counter);
	repeated_value = data[i];
	value_counter = 1;
      }
    }
  }

  // There's always an unwritten value when we reach end of loop.
#ifdef DEBUG_ENCODING
	cout << "encoder_runlength: encoding " << value_counter
             << " repeats of " << repeated_value << endl;
#endif
  mexp_golomb_write(out_bs, repeated_value);
  mexp_golomb_write(out_bs, value_counter);
}



//----------------------------------------------------------------------
/// Compute the parameters needed for the encoding algorithm.
/// This algorithm has no parameters: it can be a good or a bad idea,
/// but this routine has nothing to do.
/// \param data  Array of data to be analyzed (can be a statistical sample).
/// \param ndata Length of data array.
/// \return Error code, or 0=no error.
int encoder_runlength::compute_params(const uint32_t *data, 
				      const int ndata) {
  uint32_t prevdata = data[0];
  nchanges_data = 1;
  for (int i=1; i<ndata; i++) {
    if (data[i] != prevdata) {
      nchanges_data ++; 
      prevdata = data[i];
    }
  }
  ndata_checked = ndata;
  return 0;
}



//----------------------------------------------------------------------
/// Compute the parameters needed for the encoding algorithm.
/// This algorithm has no parameters: it can be a good or a bad idea,
/// but this routine has nothing to do.
/// \param data  Array of data to be analyzed (can be a statistical sample).
/// \param ndata Length of data array.
/// \return Error code, or 0=no error.
int encoder_runlength::compute_params(const uint16_t *data, 
				      const int ndata) {
  uint16_t prevdata = data[0];
  nchanges_data = 1;
  for (int i=1; i<ndata; i++) {
    if (data[i] != prevdata) {
      nchanges_data ++; 
      prevdata = data[i];
    }
  }
  ndata_checked = ndata;
  return 0;
}



//----------------------------------------------------------------------
/// Write parameters of the encoder to the output bitstream.
/// \return Error code, or 0=no error.
int encoder_runlength::write_params() const {
  if (out_bs == NULL)
    return -1;

  out_bs->writebits(ALGORITHM_CODE, BITS_SLIM_ALG_CODE);
  out_bs->writebits(data_type, BITS_SLIM_TYPE_CODE);
  return 0;
}



//----------------------------------------------------------------------
/// Do we expect this channel to have no effect?
/// \return true if compression is a bad idea; false if it's good.
bool encoder_runlength::expect_zero_compression() const {
  if (nchanges_data*3 > ndata_checked)
    return true;
  return false;
}



encoder *encoder_runlength::replacement_encoder() {
  return new encoder_reduced_binary(data_type, use_deltas);
}




//----------------------------------------------------------------------
/// \class decoder_runlength
/// Derived class for encoding highly repetitive streams with runlength.
/// That is, as (value, repeats) pairs.
//----------------------------------------------------------------------

/// Constructor is simply the base constructor.
/// \param dt  The data type code.
/// \param deltas Whether to encode successive differences.
/// \param ib  The input bitstream.
decoder_runlength::decoder_runlength(enum data_t dt, bool deltas, ibitstream *ib) :
  decoder(dt, deltas, ib) {
  set_data_type(SLIM_TYPE_U32);

  repeated_value = 0;
  uses_remaining = 0u;
}



//----------------------------------------------------------------------
/// Destructor is simply the base destructor.
decoder_runlength::~decoder_runlength() {}



//----------------------------------------------------------------------
/// Load decoder parameters from the bitstream.
/// This algorithm has no stored parameters.
/// \return Error code, or 0=no error. 
int decoder_runlength::read_params() {
  if (in_bs == NULL)
    return -1;
  return 0;
}



//----------------------------------------------------------------------
/// Decode one word from the input stream.
/// \return  The decoded word.
uint32_t decoder_runlength::decode_u32() 
{
  if (uses_remaining <=0) {
    try {
      repeated_value = mexp_golomb_read_u32(in_bs);
      uses_remaining = mexp_golomb_read_u32(in_bs);
    } catch (const char * s) { // This catches end of input.
      return 0u;
    }
  }
  --uses_remaining;
  return repeated_value;
}



//----------------------------------------------------------------------
/// Write decoder parameters to a stream.
/// \param fout  The writeable output stream.
void decoder_runlength::dump_info(ostream &fout) const {
  fout << "     Runlength decoder\n";
}
