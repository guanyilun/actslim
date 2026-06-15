/// \file slim_codec_default.cpp
/// Implement classes encoder and decoder that do not actually compress.

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
#include "bitstream.h"

//----------------------------------------------------------------------
/// \class encoder
///
/// Base class for all stream encoders.
//----------------------------------------------------------------------

/// Constructor (optionally by output bitstream).
/// \param dt  The data type code.
/// \param deltas Whether to encode successive differences.
/// \param ob  The output bitstream.
encoder::encoder(enum data_t dt, bool deltas, obitstream *ob) :
  use_deltas(deltas)
{
  set_data_type(dt);
  out_bs = ob;
  prev_datum = 0u;
  prev_sdatum = 0u;
  prev_cdatum = 0u;
}



/// Destructor
encoder::~encoder() {
}



/// Fix the data type from the enumerated list of allowed types.
/// Notice that this can only be run once per object.
/// \param dt_in The data type requested.
/// \return True on success, false on failure.
bool encoder::set_data_type(enum data_t dt_in)
{
  data_type = dt_in;
  switch (data_type) {
  case SLIM_TYPE_I32:
  case SLIM_TYPE_U32:
    data_size_bytes = 4;
    data_size_bits = 32;
    break;
  case SLIM_TYPE_I16:
  case SLIM_TYPE_U16:
    data_size_bytes = 2;
    data_size_bits = 16;
    break;
  case SLIM_TYPE_I8:
  case SLIM_TYPE_U8:
    data_size_bytes = 1;
    data_size_bits = 8;
    break;
  default:
    cout << "Type " << data_type << " size not known.\n";
    throw "Unknown data type in encoder::set_data_type";
    break;
  }
  return true;
}



/// Force the data type to be a signed type.
/// This is important to call if it will be encoding deltas.
void encoder::use_signed_data_type()
{
  if (data_type == SLIM_TYPE_U32)
    data_type = SLIM_TYPE_I32;
  else if (data_type == SLIM_TYPE_U16)
    data_type = SLIM_TYPE_I16;
  else if (data_type == SLIM_TYPE_U8)
    data_type = SLIM_TYPE_I8;
}



/// Is this encoder's data type a signed data type?
bool encoder::is_signed() const {
  if (data_type == SLIM_TYPE_U32 ||
      data_type == SLIM_TYPE_U16 ||
      data_type == SLIM_TYPE_U8)
    return false;
  return true;
}



/// Set up the output file using an open obitstream.  Do not own it
/// \return True on success, false on failure.
bool encoder::set_output(obitstream *ob) {
  out_bs = ob;
  return true;
}



/// Encode to the output stream.
/// This is virtual--subclasses will override this with their own method.
/// \param datum  The word to be encoded.
void encoder::encode(uint32_t datum) const {
  out_bs -> writebits(datum, 32);
}



/// Encode to the output stream.
/// This is virtual--subclasses will override this with their own method.
/// \param datum  The word to be encoded.
void encoder::encode(uint16_t datum) const {
  out_bs -> writebits(datum, 16);
}



/// Encode to the output stream.
/// This is virtual--subclasses will override this with their own method.
/// \param datum  The word to be encoded.
void encoder::encode(uint8_t datum) const {
  out_bs -> writebits(datum, 8);
}



/// Encode a data value to the output stream.
/// \param data  Array of words to be encoded.
void encoder::encode_scalar(const uint32_t *data) {
  if (use_deltas) {
    encode(data[0]-prev_datum);
    prev_datum = data[0];
  } else {
    encode(data[0]);
  }
}



/// Encode a data value to the output stream.
/// \param data  Array of words to be encoded.
void encoder::encode_scalar(const uint16_t *data) {
  if (use_deltas) {
    encode(uint16_t(data[0]-prev_sdatum));
    prev_sdatum = data[0];
  } else {
    encode(data[0]);
  }
}



/// Encode a data value to the output stream.
/// \param data  Array of words to be encoded.
void encoder::encode_scalar(const uint8_t *data) {
  if (use_deltas) {
    encode(uint8_t(data[0]-prev_cdatum));
    prev_cdatum = data[0];
  } else {
    encode(data[0]);
  }
}



/// Encode a data vector to the output stream.
/// \param data  Array of words to be encoded.
/// \param ndata Number of words to be encoded.
void encoder::encode_vector(const uint32_t *data, int ndata) {
  if (use_deltas) {
    // It's important to keep prev_datum out of the loop:
    // This saves 3 instructions per data point.
    encode(data[0]-prev_datum);
    for (int i=1; i<ndata; i++) {
      encode(data[i]-data[i-1]);
    }
    prev_datum = data[ndata-1];
  } else {
    for (int i=0; i<ndata; i++)
      encode(data[i]);
  }
}



/// Encode a data vector to the output stream.
/// \param data  Array of words to be encoded.
/// \param ndata Number of words to be encoded.
void encoder::encode_vector(const uint16_t *data, int ndata) {
  if (use_deltas) {
    // It's important to keep prev_datum out of the loop:
    // This saves 3 instructions per data point.
    encode((uint16_t)(data[0]-prev_sdatum));
    for (int i=1; i<ndata; i++) {
      encode((uint16_t)(data[i]-data[i-1]));
    }
    prev_sdatum = data[ndata-1];
  } else {
    for (int i=0; i<ndata; i++)
      encode(data[i]);
  }
}



/// Encode a data vector to the output stream.
/// \param data  Array of words to be encoded.
/// \param ndata Number of words to be encoded.
void encoder::encode_vector(const uint8_t *data, int ndata) {
  if (use_deltas) {
    // It's important to keep prev_datum out of the loop:
    // This saves 3 instructions per data point.
    encode((uint8_t)(data[0]-prev_cdatum));
    for (int i=1; i<ndata; i++) {
      encode((uint8_t)(data[i]-data[i-1]));
    }
    prev_cdatum = data[ndata-1];
  } else {
    for (int i=0; i<ndata; i++)
      encode(data[i]);
  }
}

/// Reset the prev_ data members to 0 so that a new section can be
/// encoded.
void encoder::reset_previous()
{
  prev_datum = 0u;
  prev_sdatum = 0u;
  prev_cdatum = 0u;
}

/// Compute parameters of the encoder based on a sample of data.
/// This is virtual--subclasses will override this with their own method.
/// \param data   Array of data words.
/// \param ndata  Number of data words in array.
/// \return An error code, 0=no error.
int encoder::compute_params(const uint32_t *data, const int ndata)
{
  data_type = SLIM_TYPE_I32;
  return 0;
}



/// Compute parameters of the encoder based on a sample of data.
/// This is virtual--subclasses will override this with their own method.
/// \param data   Array of data words.
/// \param ndata  Number of data words in array.
/// \return An error code, 0=no error.
int encoder::compute_params(const uint16_t *data, const int ndata)
{
  data_type = SLIM_TYPE_I16;
  return 0;
}



/// Compute parameters of the encoder based on a sample of data.
/// This is virtual--subclasses will override this with their own method.
/// \param data   Array of data words.
/// \param ndata  Number of data words in array.
/// \return An error code, 0=no error.
int encoder::compute_params(const uint8_t *data, const int ndata)
{
  data_type = SLIM_TYPE_I8;
  return 0;
}



/// Write parameters of the encoder to the output bitstream.
/// This is virtual--subclasses will override this with their own method.
/// \return An error code, 0=no error.
int encoder::write_params() const {
  if (out_bs == NULL)
    return -1;
  out_bs->writebits(ALGORITHM_CODE, BITS_SLIM_ALG_CODE);
  out_bs->writebits(data_type, BITS_SLIM_TYPE_CODE);
  return 0;
}



/// Do we expect this channel to have no effect?
bool encoder::expect_zero_compression() const {
  return false;
}



/// A replacement encoder to use when the type in question is no good.
/// This will be the default (dummy) encoder WITHOUT DELTAS, whether or not
/// the current encoder uses them.  (It's the job of slim_channel_encode to
/// notice this change in status.)
/// \return  Encoder of the desired replacement type.
encoder *encoder::replacement_encoder() {
  const bool NEVER_USE_DELTAS = false;
  return new encoder(data_type, NEVER_USE_DELTAS);
}



/// A constant encoder to use when the data is unvarying.
/// This will be the default (dummy) encoder WITHOUT DELTAS, whether or not
/// the current encoder uses them.  (It's the job of slim_channel_encode to
/// notice this change in status.)
/// \return  Encoder of the desired replacement type.
encoder *encoder::constant_encoder(int d0) {
  const bool NEVER_USE_DELTAS = false;
  return new encoder_constant(d0, data_type, NEVER_USE_DELTAS);
}



/// Compute the mean of a vector.
/// This is not used by the default encoder, but it's valuable to
/// some or all of the real encoders.  Let them inherit it from 
/// the base class.
/// \param mean       The mean of the data array.
/// \param data       The data array.
/// \param ndata      The length of the data array.
template <typename T>
void encoder::compute_mean(double& mean,
			   const T *data, int ndata) const {
  // Mean first
  mean = 0.;
  for (int i=0; i<ndata; i++)
    mean += data[i];
  mean /= ndata;
}

/// Instantiation for array of int32_t
template void encoder::compute_mean<int32_t>
(double&, const int32_t *, int) const;

/// Instantiation for array of uint32_t
template void encoder::compute_mean<uint32_t>
(double&, const uint32_t *, int) const;

/// Instantiation for array of int16_t
template void encoder::compute_mean<int16_t>
(double&, const int16_t *, int) const;

/// Instantiation for array of uint16_t
template void encoder::compute_mean<uint16_t>
(double&, const uint16_t *, int) const;

/// Instantiation for array of int8_t
template void encoder::compute_mean<int8_t>
(double&, const int8_t *, int) const;

/// Instantiation for array of uint8_t
template void encoder::compute_mean<uint8_t>
(double&, const uint8_t *, int) const;



//----------------------------------------------------------------------
/// \class decoder
///
/// Base class for all stream decoders.
//----------------------------------------------------------------------

/// Constructor.
/// \param dt  The data type code.
/// \param deltas Whether to encode successive differences.
/// \param ib  The input bitstream.
decoder::decoder(enum data_t dt, bool deltas, ibitstream *ib) :
  use_deltas(deltas) {
  set_data_type(dt);
  in_bs = ib;
  prev_datum = 0u;
  prev_sdatum = 0u;
  prev_cdatum = 0u;
}

/// Destructor
decoder::~decoder() {
}

/// Set up the input file as an open ibitstream.  Do NOT own the stream.
/// \return True on success, false on failure.
bool decoder::set_input(ibitstream *ib) {
  in_bs = ib;
  return true;
}

/// Fix the data type from the enumerated list of allowed types.
/// Notice that this can only be run once per object.
/// \param dt_in The data type requested.
/// \return True on success, false on failure.
bool decoder::set_data_type(enum data_t dt_in)
{
  data_type = dt_in;

  switch (data_type) {
  case SLIM_TYPE_I32:
  case SLIM_TYPE_U32:
    data_size_bytes = 4;
    data_size_bits = 32;
    break;
  case SLIM_TYPE_I16:
  case SLIM_TYPE_U16:
    data_size_bytes = 2;
    data_size_bits = 16;
    break;
  case SLIM_TYPE_I8:
  case SLIM_TYPE_U8:
    data_size_bytes = 1;
    data_size_bits = 8;
    break;
  default:
    cout << "Type " << data_type << " size not known.\n";
    throw "Unknown data type in decoder::set_data_type";
    break;
  }
  return true;
}

/// Decode one word from the input stream, ignoring deltas.
/// \return  The decoded word.
inline uint32_t decoder::decode_u32() 
{
  return in_bs -> readbits(32);
}



/// Decode one word from the input stream, ignoring deltas.
/// \return  The decoded word.
inline uint16_t decoder::decode_u16() 
{
  return in_bs -> readbits(16);
}



/// Decode one word from the input stream, ignoring deltas.
/// \return  The decoded word.
inline uint8_t decoder::decode_u8() 
{
  return in_bs -> readbits(8);
}



/// Decode one word from the input stream.
/// \param data   Array of decoded words.
void decoder::decode_scalar(uint32_t *data)
{
  if (use_deltas) {
    // It's important to keep prev_datum out of the loop:
    // This saves 2 instructions per data point.
    data[0] = decode_u32() + prev_datum;
    prev_datum = data[0];
  } else {
    data[0] = decode_u32();
  }
  return;
}



/// Decode one word from the input stream.
/// \param data   Array of decoded words.
void decoder::decode_scalar(uint16_t *data)
{
  if (use_deltas) {
    // It's important to keep prev_datum out of the loop:
    // This saves 2 instructions per data point.
    data[0] = decode_u16() + prev_sdatum;
    prev_sdatum = data[0];
  } else {
    data[0] = decode_u16();
  }
  return;
}



/// Decode one word from the input stream.
/// \param data   Array of decoded words.
void decoder::decode_scalar(uint8_t *data)
{
  if (use_deltas) {
    // It's important to keep prev_datum out of the loop:
    // This saves 2 instructions per data point.
    data[0] = decode_u8() + prev_cdatum;
    prev_cdatum = data[0];
  } else {
    data[0] = decode_u8();
  }
  return;
}



/// Decode several words from the input stream.
/// \param data   Array of decoded words.
/// \param ndata  Number of words to decode.
void decoder::decode_vector(uint32_t *data, int ndata)
{
  if (use_deltas) {
    // It's important to keep prev_datum out of the loop:
    // This saves 2 instructions per data point.
    data[0] = decode_u32() + prev_datum;
    for (int i=1; i<ndata; i++) {
      data[i] = decode_u32() + data[i-1];
    }
    prev_datum = data[ndata-1];
  } else {
    for (int i=0; i<ndata; i++)
      data[i] = decode_u32();
  }
  return;
}



/// Decode several words from the input stream.
/// \param data   Array of decoded words.
/// \param ndata  Number of words to decode.
void decoder::decode_vector(uint16_t *data, int ndata)
{
  if (use_deltas) {
    // It's important to keep prev_datum out of the loop:
    // This saves 2 instructions per data point.
    data[0] = decode_u16() + prev_sdatum;
    for (int i=1; i<ndata; i++) {
      data[i] = decode_u16() + data[i-1];
    }
    prev_sdatum = data[ndata-1];
  } else {
    for (int i=0; i<ndata; i++)
      data[i] = decode_u16();
  }
  return;
}



/// Decode several words from the input stream.
/// \param data   Array of decoded words.
/// \param ndata  Number of words to decode.
void decoder::decode_vector(uint8_t *data, int ndata)
{
  if (use_deltas) {
    // It's important to keep prev_datum out of the loop:
    // This saves 2 instructions per data point.
    data[0] = decode_u8() + prev_cdatum;
    for (int i=1; i<ndata; i++) {
      data[i] = decode_u8() + data[i-1];
    }
    prev_cdatum = data[ndata-1];
  } else {
    for (int i=0; i<ndata; i++)
      data[i] = decode_u8();
  }
  return;
}



/// Load parameters of the decoder from the input bit stream.
/// \return An error code, 0=no error.
int decoder::read_params()
{
  return 0;
}



/// Write decoder parameters to a stream.
/// \param fout  The writeable output stream.
void decoder::dump_info(ostream &fout) const {
  fout << "  Default (dummy) decoder";
}
