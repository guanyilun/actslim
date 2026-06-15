/// \file slim_codec_reduced_binary.cpp
/// Implement classes encoder_reduced_binary and decoder_reduced_binary,
/// the workhorses of slim.

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
#include <iomanip>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <cmath>

#include "slim.h"
#include "bitstream.h"

//#define DEBUG_ENCODING


//----------------------------------------------------------------------
/// \class encoder_reduced_binary
/// Derived class for encoding simple full range of data.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// Constructor (optionally by output bitstream).
/// Set the encoder algorithm parameters to default values.
/// \param dt  The data type code.
/// \param deltas Whether to encode successive differences.
/// \param ob  The output bitstream.
encoder_reduced_binary::encoder_reduced_binary(enum data_t dt, bool deltas, 
					       obitstream *ob) 
  : encoder(dt, deltas, ob) {
  nbits = 8u * slim_type_size[dt];
  max = UINT_MAX;
  offset = 0u;
  Overflow = UINT_MAX;
}



//----------------------------------------------------------------------
/// Destructor is simply the base destructor.
encoder_reduced_binary::~encoder_reduced_binary() {;}



//----------------------------------------------------------------------
/// Encode one datum to the output stream.
/// \param datum  The word to be encoded.
inline void encoder_reduced_binary::encode(uint32_t  datum) const {
  uint32_t u = datum - offset;
  if (u <= max) {

#ifdef DEBUG_ENCODING
  cout << "Datum " << setw(10) << datum;
  cout << ", offset " << setw(10) << (-offset);
  cout << ", -> encode " << u << " (" << nbits << ")\n";
#endif

    out_bs->writebits(u, nbits);
    return;
  }

#ifdef DEBUG_ENCODING
  cout << "Datum " << setw(10) << datum;
  cout << ", offset " << setw(10) << (-offset);
  cout << ", -> OFLOW  " << u << endl;
#endif

  out_bs->writebits(Overflow, nbits);
  out_bs->writebits(datum, data_size_bits);
}



//----------------------------------------------------------------------
/// Encode one datum to the output stream.
/// \param datum  The word to be encoded.
inline void encoder_reduced_binary::encode(uint16_t  datum) const {
  uint16_t u = datum - offset;
  if (u <= max) {
    out_bs->writebits(u, nbits);
    return;
  }

  out_bs->writebits(Overflow, nbits);
  out_bs->writebits(datum, data_size_bits);
}



//----------------------------------------------------------------------
/// Encode one datum to the output stream.
/// \param datum  The word to be encoded.
inline void encoder_reduced_binary::encode(uint8_t  datum) const {
  uint8_t u = datum - offset;
  if (u <= max) {
    out_bs->writebits(u, nbits);
    return;
  }

  out_bs->writebits(Overflow, nbits);
  out_bs->writebits(datum, data_size_bits);
}



//----------------------------------------------------------------------
/// Count the number of overflows for encoding full sampled data set,
/// assuming a given choice of nbits.  (Idea is to try this for all
/// possible nbits and use the best.)
/// \param histogram The 33-element histogram giving distribution of 
///        sizes of offset-subtracted data.
/// \param n The proposed number of bits to use for symbols.
/// \return  Number of data bits wasted on overflows.
int encoder_reduced_binary::overflow_waste(const int histogram[33], 
					   unsigned int n) {
  int num_oflow=0;
  for (unsigned int l=1+n; l<=data_size_bits; l++)
    num_oflow += histogram[l];
  return num_oflow*data_size_bits;
}



//----------------------------------------------------------------------
/// Compute the number of bits for encoding full sampled data set,
/// assuming a given choice of nbits.  (Idea is to try this for all
/// possible nbits and use the best.)
/// \param histogram The 33-element histogram giving distribution of 
///        sizes of offset-subtracted data.
/// \param ndata  Total number of data in histogram.
int encoder_reduced_binary::best_code_length(const int histogram[33],
					     int ndata) {
  int best_length=INT_MAX;
  int best_nbits=data_size_bits;

  for (int i=data_size_bits; i>0; i--) {
    int cl = overflow_waste(histogram,i) + i*ndata;
    if (cl < best_length) {
      best_length = cl;
      best_nbits = i;
    }
  }
  return best_nbits;
}
 


//----------------------------------------------------------------------
/// Log-base2 of the argument.
inline double log2(double x) {
  return log(x)/M_LN2;
}



//----------------------------------------------------------------------
/// Compute the parameters needed for the encoding algorithm.
/// Here assume that the max/min of the data to be encoded is the
/// range that needs to be encoded without an Overflow.
/// \param data  Array of data to be analyzed (can be a statistical sample).
/// \param ndata Length of data array.
/// \return Error code, or 0=no error.
int encoder_reduced_binary::compute_params(const uint32_t *data, 
					   const int ndata) {
  assert (data_type == SLIM_TYPE_I32 ||
	  data_type == SLIM_TYPE_U32);

  const int32_t *idata = reinterpret_cast<const int32_t *>(data);

  // Estimate mean of the data.  Temporarily call this the "offset"
  double avg;
  if (data_type == SLIM_TYPE_U32) {
    compute_mean(avg, data, ndata);
    offset = (uint32_t)(nearbyint(avg));
  } else {
    compute_mean(avg, idata, ndata);
    offset = (uint32_t)(nearbyint(avg));
  }

  // Make histogram of how many bits we'd need for all data if the mean were removed.
  int histogram[33]={0};
  uint32_t set_bits=0;
  for (int i=0; i<ndata; i++) {
    histogram[bit_size(int32_t(data[i]-offset))]++;
    set_bits |= data[i];
  }
  nbits = best_code_length(histogram, ndata);

  if (nbits > data_size_bits)
    nbits = data_size_bits;
  if (nbits < 1u)
    nbits = 1u;

#ifdef DEBUG_ENCODING
  cout << setw(11) << "debug: data mean = " << offset;
#endif

  // We found assuming integer data and offset = avg in middle of range.
  // But we _want_ unsigned data and offset at the bottom of range.
  if (nbits > 1)
    offset -= 1 << (nbits-1u);

  // set Overflow to the largest nbits-bit number
  Overflow = lowestNset[nbits];
  max = Overflow-1;

#ifdef DEBUG_ENCODING
  cout << setw(11) << " min="<<offset;
  cout << setw(11) << " max="<<offset+max;
  cout << setw(2) << " nbits="<<nbits << endl;
#endif

  return 0;
}



//----------------------------------------------------------------------
/// Compute the parameters needed for the encoding algorithm.
/// Here assume that the max/min of the data to be encoded is the
/// range that needs to be encoded without an Overflow.
/// \param data  Array of data to be analyzed (can be a statistical sample).
/// \param ndata Length of data array.
/// \return Error code, or 0=no error.
int encoder_reduced_binary::compute_params(const uint16_t *data, 
					   const int ndata) {
  assert (data_type == SLIM_TYPE_I16 ||
	  data_type == SLIM_TYPE_U16);

  const int16_t *sdata = 
    reinterpret_cast<const int16_t *>(data);

  // Estimate mean of the data
  double avg;
  if (data_type == SLIM_TYPE_I16) {
    compute_mean(avg, sdata, ndata);
    offset = (uint16_t)(nearbyint(avg));
  } else {
    compute_mean(avg, data, ndata);
    offset = (uint16_t)(nearbyint(avg));
  }

  int histogram[33]={0};
  uint32_t set_bits=0;
  for (int i=0; i<ndata; i++) {
    histogram[bit_size((int16_t)(data[i]-offset))]++;
    set_bits |= data[i];
  }
  nbits = best_code_length(histogram, ndata);

  if (nbits > data_size_bits)
    nbits = data_size_bits;
  if (nbits < 1u)
    nbits = 1u;
  
  // We found assuming integer data and offset = avg in middle of range.
  // But we _want_ unsigned data and offset at the bottom of range.
  if (nbits > 1)
    offset -= 1 << (nbits-1u);

  // set Overflow to the largest nbits-bit number
  Overflow = lowestNset[nbits];
  max = Overflow-1;

  return 0;
}

//----------------------------------------------------------------------
/// Compute the parameters needed for the encoding algorithm.
/// Here assume that the max/min of the data to be encoded is the
/// range that needs to be encoded without an Overflow.
/// \param data  Array of data to be analyzed (can be a statistical sample).
/// \param ndata Length of data array.
/// \return Error code, or 0=no error.
int encoder_reduced_binary::compute_params(const uint8_t *data, 
                                           const int ndata) {
  assert (data_type == SLIM_TYPE_I8 ||
          data_type == SLIM_TYPE_U8);

  const int8_t *sdata = 
    reinterpret_cast<const int8_t *>(data);

  // Estimate mean of the data
  double avg;
  if (data_type == SLIM_TYPE_I8) {
    compute_mean(avg, sdata, ndata);
    offset = (uint8_t)(nearbyint(avg));
  } else {
    compute_mean(avg, data, ndata);
    offset = (uint8_t)(nearbyint(avg));
  }

  int histogram[33]={0};
  uint32_t set_bits=0;
  for (int i=0; i<ndata; i++) {
    histogram[bit_size((int8_t)(data[i]-offset))]++;
    set_bits |= data[i];
  }
  nbits = best_code_length(histogram, ndata);

  if (nbits > data_size_bits)
    nbits = data_size_bits;
  if (nbits < 1u)
    nbits = 1u;
  
  // We found assuming integer data and offset = avg in middle of range.
  // But we _want_ unsigned data and offset at the bottom of range.
  if (nbits > 1)
    offset -= 1 << (nbits-1u);

  // set Overflow to the largest nbits-bit number
  Overflow = lowestNset[nbits];
  max = Overflow-1;

  return 0;
}



//----------------------------------------------------------------------
/// Write parameters of the encoder to the output bitstream.
/// \return Error code, or 0=no error.
int encoder_reduced_binary::write_params() const {
  if (out_bs == NULL)
    return -1;

  out_bs->writebits(ALGORITHM_CODE, BITS_SLIM_ALG_CODE);
  out_bs->writebits(data_type, BITS_SLIM_TYPE_CODE);
  out_bs->writebits(offset, data_size_bits);
  out_bs->writebits(nbits-1u, BITS_SLIM_NBITS);
  return 0;
}



///----------------------------------------------------------------------
// Do we expect this channel to have no effect?
bool encoder_reduced_binary::expect_zero_compression() const {
  return (nbits >= data_size_bits);
}




//----------------------------------------------------------------------
/// \class decoder_reduced_binary
/// Derived class for decoding simple full range of data.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// Constructor.
/// \param dt  The data type code.
/// \param deltas Whether to encode successive differences.
/// \param ib  The input bitstream.
decoder_reduced_binary::decoder_reduced_binary(enum data_t dt, bool deltas,
					       ibitstream *ib)
  : decoder(dt, deltas, ib) {
  switch (dt) {
  case SLIM_TYPE_I8:
  case SLIM_TYPE_U8:
    nbits = 8u;
    break;
  case SLIM_TYPE_I16:
  case SLIM_TYPE_U16:
    nbits = 16u;
    break;
  default:
    nbits = 32u;
  }
  max = UINT_MAX;
  offset = 0u;
  Overflow = 0u;
}



//----------------------------------------------------------------------
/// Destructor is simply the base destructor.
decoder_reduced_binary::~decoder_reduced_binary() {;}



//----------------------------------------------------------------------
/// Load decoder parameters from the bitstream.
/// Only some are stored; reconstruct the others by computing them.
/// \return Error code, or 0=no error. 
int decoder_reduced_binary::read_params() {
  if (in_bs == NULL)
    return -1;

  offset = in_bs->readbits(data_size_bits);
  nbits = in_bs->readbits(BITS_SLIM_NBITS);
  nbits++;
  if (nbits <= 0 || nbits > unsigned(data_size_bits))
    throw "Cannot decode: unexpected number of bits read from file.";

  int pow2n=1;
  for (unsigned int i=0; i<nbits; i++)
    pow2n *= 2;
  Overflow = pow2n-1;
  max = Overflow-1;
  return 0;
}



//----------------------------------------------------------------------
/// Decode one word from the input stream.
/// \return  The decoded word.
uint32_t decoder_reduced_binary::decode_u32() 
{
  uint32_t datum;
  try {
    datum = in_bs->readbits(nbits);
    if (datum == Overflow) 
      return in_bs->readbits(data_size_bits);
    else
      return datum+offset;
  } catch (const char * s) { // This catches end of input.
    return 0;
  }
}



//----------------------------------------------------------------------
/// Decode one word from the input stream.
/// \return  The decoded word.
uint16_t decoder_reduced_binary::decode_u16() 
{
  uint16_t datum;
  try {
    datum = in_bs->readbits(nbits);
    if (datum == Overflow) {
      return in_bs->readbits(data_size_bits);
    } else
      return datum+offset;
  } catch (const char * s) { // This catches end of input.
    return 0;
  }
}



//----------------------------------------------------------------------
/// Decode one word from the input stream.
/// \return  The decoded word.
uint8_t decoder_reduced_binary::decode_u8() 
{
  uint8_t datum;
  try {
    datum = in_bs->readbits(nbits);
    if (datum == Overflow) {
      return in_bs->readbits(data_size_bits);
    } else
      return datum+offset;
  } catch (const char * s) { // This catches end of input.
    return 0;
  }
}



//----------------------------------------------------------------------
/// Write decoder parameters to a stream.
/// \param fout  The writeable output stream.
void decoder_reduced_binary::dump_info(ostream &fout) const {

  switch (data_type) {
  case SLIM_TYPE_I32: case SLIM_TYPE_I16: case SLIM_TYPE_I8:
    fout << "  RedBinary: " << setw(2) << nbits;
    fout << " bit, offset " << setw(11) << int(offset);
    break;

  case SLIM_TYPE_U32: case SLIM_TYPE_U16: case SLIM_TYPE_U8:
    fout << "  RedBinary: " << setw(2) << nbits;
    fout << " bit, offset " << setw(11) << (unsigned int)(offset);
    break;

  case SLIM_TYPE_FLOAT: case SLIM_TYPE_DOUBLE:
  default:
    fout << "  RedBinary: " << setw(2) << nbits;
    fout << " bit, offset " << setw(11) << (offset);
    break;
  }
}
