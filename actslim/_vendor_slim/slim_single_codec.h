// -*-  mode: c++; -*-

/// \file slim_single_codec.h
/// Inline functions to standardize some multi-step systems for
/// coding and decoding single values.
/// Includes the "modified exponential Golomb" method.

#ifndef SLIM_SINGLE_CODEC_H
#define SLIM_SINGLE_CODEC_H

#include "bitstream.h"



/// Write an unsigned value to a bitstream by method mexp_golomb.
/// \param ob    The bitstream to write on.
/// \param u     The number to encode.
/// \param order The order of the code, i.e. the minimum # of value bits.
inline void mexp_golomb_write(obitstream *ob, uint32_t u, 
			      unsigned int order=1) {
  unsigned int n = bit_size(u);
  if (n > order) {
    ob->write_unary(n-order);
    ob->writebits(u, n-1);
  } else {
    ob->write_unary(0);
    ob->writebits(u, order);
  }
}



/// Read an unsigned 32-bit value from a bitstream by method mexp_golomb.
/// \param ib     The bitstream to read.
/// \param order The order of the code, i.e. the minimum # of value bits.
/// \return       The value read.
inline uint32_t mexp_golomb_read_u32(ibitstream *ib,
                                     unsigned int order=1) {

  uint32_t n_minus_order = ib->read_unary();
  if (n_minus_order > 0) {
    int n_minus_1 = n_minus_order + order - 1;

    uint32_t uval = ib->readbits(n_minus_1);
    return uval | bitNset[n_minus_1];
  } else {
    return ib->readbits(order);
  }
}





#endif  // #ifndef SLIM_SINGLE_CODEC_H
