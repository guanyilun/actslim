/// \file slim_codec_constant.cpp
/// Implement classes encoder_constant and decoder_constant.

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
#include "slim_single_codec.h"


//----------------------------------------------------------------------
/// \class encoder_constant
/// Derived class for encoding strictly constant data channels.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// Construct (optionally by output bitstream).
/// \param dt  The data type code.
/// \param deltas Whether to encode successive differences.
/// \param ob  The output bitstream.
encoder_constant::encoder_constant(int value, enum data_t dt, bool deltas, 
                                   obitstream *ob) 
  : encoder(dt, deltas, ob) {
  set_data_type(dt);
  fixed_data = value;
  fixed_sdata = value;
  fixed_cdata = value;
}



//----------------------------------------------------------------------
/// Destructor is simply the base destructor.
encoder_constant::~encoder_constant() {;}


//----------------------------------------------------------------------
/// Encode by writing precisely nothing.
void encoder_constant::encode(uint32_t datum) const {
  if (datum != fixed_data) {
    cerr << "Error: encoder_constant::fixed_data=" << fixed_data 
         <<", this datum=" << datum << endl;
    throw "encoder_constant must be asked to write the same data always.";
  }
}


//----------------------------------------------------------------------
/// Encode by writing precisely nothing.
void encoder_constant::encode(uint16_t datum) const {
  if (datum != fixed_sdata) {
    throw "encoder_constant must be asked to write the same data always.";
  }
}


//----------------------------------------------------------------------
/// Encode by writing precisely nothing.
void encoder_constant::encode(uint8_t datum) const {
  if (datum != fixed_cdata) {
    throw "encoder_constant must be asked to write the same data always.";
  }
}


//----------------------------------------------------------------------
/// Write parameters of the encoder to the output bitstream.
/// \return Error code, or 0=no error.
int encoder_constant::write_params() const {
  if (out_bs == NULL)
    return -1;

  out_bs->writebits(ALGORITHM_CODE, BITS_SLIM_ALG_CODE);
  out_bs->writebits(data_type, BITS_SLIM_TYPE_CODE);
  out_bs->writebits(fixed_data, data_size_bits);
  return 0;
}



//----------------------------------------------------------------------
//----------------------------------------------------------------------



//----------------------------------------------------------------------
/// \class decoder_constant
/// Derived class for encoding strictly constant data.
//----------------------------------------------------------------------

/// Constructor is simply the base constructor.
/// \param dt  The data type code.
/// \param deltas Whether to encode successive differences.
/// \param ib  The input bitstream.
decoder_constant::decoder_constant(enum data_t dt, bool deltas, ibitstream *ib) :
  decoder(dt, deltas, ib) {
  set_data_type(dt);

  fixed_data = 0;
  fixed_sdata = 0;
  fixed_cdata = 0;
}



//----------------------------------------------------------------------
/// Destructor is simply the base destructor.
decoder_constant::~decoder_constant() {}



//----------------------------------------------------------------------
/// Load decoder parameters from the bitstream.
/// This algorithm has one parameter: the fixed data value.
/// \return Error code, or 0=no error. 
int decoder_constant::read_params() {
  if (in_bs == NULL)
    return -1;

  fixed_data = in_bs->readbits(data_size_bits);
  fixed_sdata = fixed_data;
  fixed_cdata = fixed_data;
  return 0;
}



//----------------------------------------------------------------------
/// Decode one word from the input stream.
/// \return  The decoded word.
uint32_t decoder_constant::decode_u32()
{
  return fixed_data;
}



//----------------------------------------------------------------------
/// Decode one word from the input stream.
/// \return  The decoded word.
uint16_t decoder_constant::decode_u16()
{
  return fixed_sdata;
}



//----------------------------------------------------------------------
/// Decode one word from the input stream.
/// \return  The decoded word.
uint8_t decoder_constant::decode_u8()
{
  return fixed_cdata;
}



//----------------------------------------------------------------------
/// Write decoder parameters to a stream.
/// \param fout  The writeable output stream.
void decoder_constant::dump_info(ostream &fout) const {
  switch (data_type) {
  case SLIM_TYPE_I32: case SLIM_TYPE_I16: case SLIM_TYPE_I8:
    fout << "  Constant value signed     " << setw(11) 
         << int(fixed_data);
    break;

  case SLIM_TYPE_U32: case SLIM_TYPE_U16: case SLIM_TYPE_U8:
    fout << "  Constant value unsigned   " << setw(11) 
         << (unsigned int)(fixed_data);
    break;

  case SLIM_TYPE_FLOAT: case SLIM_TYPE_DOUBLE:
    fout << "  Constant value float      " << setw(11) 
         << double(fixed_data);
    break;

  default:
    fout << "  Constant value, type unknown: " << setw(11) 
         << fixed_data;
    break;
  }
}
