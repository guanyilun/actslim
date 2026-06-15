/// \file slim_codec_factory.cpp
/// Contains factory functions encoder_generator() and
/// decoder_generator() to return a derived specific co/dec object as
/// a base-class pointer to encoder or decoder.

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
#include <sysexits.h>

#include "slim.h"



/// Generating function for various encoders.
/// Returned object is a pointer to encoder whose specific derivation
/// depends on code algorithm and data type. 
encoder *encoder_generator(enum code_t code, enum data_t data_type, bool deltas) {
  encoder *s=NULL;

  // For now, pretend float-32 data are int32 and hope for the best.
  if (data_type == SLIM_TYPE_FLOAT)
    data_type = SLIM_TYPE_I32;

  switch (code) {
  case SLIM_ENCODER_DEFAULT:
    s = new encoder(data_type, deltas);
    break;

  case SLIM_ENCODER_REDUCED_BINARY:
  case SLIM_ENCODER_CODE_A:  // Treat these 2 as the same.
    switch (data_type) {
    case SLIM_TYPE_I32:
    case SLIM_TYPE_U32:
    case SLIM_TYPE_I16:
    case SLIM_TYPE_U16:
    case SLIM_TYPE_I8:
    case SLIM_TYPE_U8:
      s = new encoder_reduced_binary(data_type, deltas);
      break;

    default:
//       cout << "Data type " << data_type << " not implemented for code "
//            << code << ".\n";
      s = new encoder(data_type, deltas); // fall back on default encoder.
      break;
    }
    break;

  case SLIM_ENCODER_RUNLENGTH:
    switch (data_type) {
    case SLIM_TYPE_I32:
    case SLIM_TYPE_U32:
      s = new encoder_runlength(data_type, deltas);
      break;

    default:
//       cout << "Data type " << data_type << " not implemented for code "
//            << code << " (RUNLENGTH).\n";
      s = new encoder(data_type, deltas); // fall back on default encoder.
      break;
    }
    break;

  case SLIM_ENCODER_CONSTANT:
    cout << "Encoder for constant data cannot be requested at command line.\n";
    break;

  case SLIM_ENCODER_CODE_B:
    cout << "Encoder for code B is no longer part of slim.\n";
    break;

  case SLIM_ENCODER_HUFFMAN:
    cout << "Encoder for Huffman codes is no longer part of slim.\n";
    break;

    
  default:
    cout << "Encoder code number "<< code<< " is not implemented.\n";
    throw "Unknown encoder type";
  }
  
  assert (s != NULL);
  return s;
}

/// Generating function for various decoders.
/// Returned object depends on code algorithm and data type.
decoder *decoder_generator(code_t code, data_t data_type, bool deltas) {
  decoder *s=NULL;
  switch (code) {
  case SLIM_ENCODER_DEFAULT:
    s = new decoder(data_type, deltas);
    break;

  case SLIM_ENCODER_REDUCED_BINARY: // fall through
  case SLIM_ENCODER_CODE_A:  // Code A has same decoder as Reduced_Binary
    switch (data_type) {
    case SLIM_TYPE_I32:
    case SLIM_TYPE_U32:
    case SLIM_TYPE_I16:
    case SLIM_TYPE_U16:
    case SLIM_TYPE_I8:
    case SLIM_TYPE_U8:
      s = new decoder_reduced_binary(data_type, deltas);
      break;
    default:
      cout << "Data type " << data_type << " not implemented for code "
           << code << " (codeA).\n";
    }
    break;

  case SLIM_ENCODER_CODE_B:
    cout << "Slim can no longer read Code B compressed files.\n";
    exit(EX_USAGE);
    break;

  case SLIM_ENCODER_HUFFMAN:
    cout << "Slim can no longer read Huffman-code compressed files.\n";
    exit(EX_USAGE);
    break;
    
  case SLIM_ENCODER_RUNLENGTH:
    switch (data_type) {
    case SLIM_TYPE_I32:
    case SLIM_TYPE_U32:
      s = new decoder_runlength(data_type, deltas);
      break;
    default:
      cout << "Data type " << data_type << " not implemented for code "
           << code << " (RUNLENGTH).\n";
    }
    break;

  case SLIM_ENCODER_CONSTANT:
    switch (data_type) {
    case SLIM_TYPE_I32:
    case SLIM_TYPE_U32:
    case SLIM_TYPE_I16:
    case SLIM_TYPE_U16:
    case SLIM_TYPE_I8:
    case SLIM_TYPE_U8:
      s = new decoder_constant(data_type, deltas);
      break;
    default:
      cout << "Data type " << data_type << " not implemented for code "
           << code << " (CONSTANT).\n";
    }
    break;
    
  default:
    cout << "Decoder code type " << code << " not implemented.\n";
    throw("Unknown decoder type");
  }
  assert (s != NULL);
  s->set_data_type(data_type);
  
  return s;
}

