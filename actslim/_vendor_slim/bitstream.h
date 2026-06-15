// -*- mode: c++; -*-

/// \file bitstream.h
/// Include file for the bitstream and derived classes.

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

#ifndef SLIM_BITSTREAM_H
#define SLIM_BITSTREAM_H

#include <stdint.h>
#include <cstdio>

#ifdef HAVE_LIBZZIP
#include <zzip/zzip.h>
#endif

#ifdef HAVE_LIBLZ4
#include <lz4frame.h>
#endif

#include "bit_constants.h"

using namespace std;

class bitstream {
public:
  bitstream(); 
  bitstream(const char *filename, 
	    int buffersize=DEFAULT_IOBUFFER_SIZE);
  virtual ~bitstream();

  virtual void close() = 0;
  virtual bool is_open() const = 0;
  virtual void setupstream();
  virtual void windup()=0;
  virtual int get_bytes_used(); 
  int get_bitptr(); 
  virtual void print() const = 0;  ///< Print properties of stream (pure virt)

protected:
  static const int Bits_per_word = 8*sizeof(Word_t); ///< Bits per buffer word.

  size_t bufsize;         ///< Size of I/O buffer (bytes)
  size_t buf_used;        ///< # words ever I/O between buffer and disk.
  Byte_t *buffer_base;    ///< Pointer to the buffer
  Byte_t *beyondbuffer;   ///< Pointer just beyond buffer (convenience).
  union {
    Byte_t *Bptr;         ///< Pointer to the current word (as Byte_t *).
    Word_t *Dptr;         ///< Pointer to the current word (as Word_t *). 
  } buffptr;              ///< Pointer to the current word.
  int bitptr;             ///< Pointer to the current bits.
  /// Make *Dptr always point to same place as Bptr.
  /// But where Bptr is a (Byte_t *), *Dptr is a (Word_t *)

public:
  enum {DEFAULT_IOBUFFER_SIZE=1024*1024};    ///< Default buffer size.
  enum {MAX_BITSTREAM_BUFSIZE=16*1024*1024}; ///< Maximum buffer size.
};



class obitstream : public bitstream {
private:
  // No private attributes.

protected:
  FILE *fp;               ///< The I/O stream.

public:
  obitstream(FILE *file, int buffersize=DEFAULT_IOBUFFER_SIZE);
  obitstream(const char *filename,
	     int buffersize=DEFAULT_IOBUFFER_SIZE);
  ~obitstream();

  void writebits(uint32_t data, int nbits);
  void writestring(const char *str, bool write_trailing_null=false);
  template <typename T> void writeword(const T data);
  void write_unary(unsigned int value); 
  virtual void print() const;
  virtual void close();
  virtual bool is_open() const;
  void windup();
  void flush(bool flush_trailing_bits);
};




class ibitstream : public bitstream {
private:
  // No private attributes (unless debugging).
#ifdef DEBUG_READBITS
  int c1,c2, c3, c4, c5, c6;
  void print_debug();
#endif

protected:
#ifdef HAVE_LIBZZIP
  ZZIP_FILE *zfp;               ///< The I/O stream.
#endif
  FILE *fp;
#ifdef HAVE_LIBLZ4
  LZ4F_decompressionContext_t lz4_ctx;  ///< LZ4 decompression context
  Byte_t *lz4_buffer;     ///< Point to LZ4 decompression buffer
  bool using_lz4;         ///< Is the file LZ4 compressed?
#endif

public:
  ibitstream(FILE *file, int buffersize=DEFAULT_IOBUFFER_SIZE);
#ifdef HAVE_LIBZZIP
  ibitstream(ZZIP_FILE *file, int buffersize=DEFAULT_IOBUFFER_SIZE);
#endif
  ibitstream(const char *filename, 
	     int buffersize=DEFAULT_IOBUFFER_SIZE);
  ibitstream(int fd, int buffersize=DEFAULT_IOBUFFER_SIZE);
  ~ibitstream();

  virtual void close();
  virtual bool is_open() const;
  void setupstream();
  void windup();
  virtual void print() const;
  virtual int get_bytes_used(); 
  Word_t readbits(int nbits);
  int32_t readbits_int(int nbits);
  Word_t read_unary(); 
  int readstring(char *s, int count=-1);
  //int get_bits_used() { return bitptr + Bits_per_word*buf_used;}

private:
  void next_word();
  int fill();

  Word_t partial_word;
  int partial_word_bitptr;
};



////////////////////////////////////////////////////////////////////////////////
// Inline functions
////////////////////////////////////////////////////////////////////////////////

/// Find size (on [0,32]) of the smallest # that can hold the integer i.
/// By our convention, [-1,0] require 1 bit, [-2,1] require 2 bits, [-4,-3,
/// 2,3] require 3 bits, etc.
/// This was implemented as loops, but we find we save a lot of compute cycles
/// by laying it out as an explicit tree of tests.
/// \param i The number whose size is being checked.
static inline unsigned int bit_size(int32_t i) {
  if (i < 0) 
    i = (-i)-1; // Convert negative int to non-negatives of same size.

  if (i>lowestNset32bits[15]) {
    if (i>lowestNset32bits[23]) {
      if (i>lowestNset32bits[27]) {
	if (i>lowestNset32bits[29]) {
	  if (i>lowestNset32bits[30]) {
	    return 32;
	  } else {
	    return 31;
	  }
	} else {
	  if (i>lowestNset32bits[28]) {
	    return 30;
	  } else {
	    return 29;
	  }
	}
      } else {	    
	if (i>lowestNset32bits[25]) {
	  if (i>lowestNset32bits[26]) {
	    return 28;
	  } else {
	    return 27;
	  }
	} else {
	  if (i>lowestNset32bits[24]) {
	    return 26;
	  } else {
	    return 25;
	  }
	}
      }
    } else {
      if (i>lowestNset32bits[19]) {
	if (i>lowestNset32bits[21]) {
	  if (i>lowestNset32bits[22]) {
	    return 24;
	  } else {
	    return 23;
	  }
	} else {
	  if (i>lowestNset32bits[20]) {
	    return 22;
	  } else {
	    return 21;
	  }
	}
      } else {	    
	if (i>lowestNset32bits[17]) {
	  if (i>lowestNset32bits[18]) {
	    return 20;
	  } else {
	    return 19;
	  }
	} else {
	  if (i>lowestNset32bits[16]) {
	    return 18;
	  } else {
	    return 17;
	  }
	}
      }
    }
  } else {
    if (i>lowestNset32bits[7]) {
      if (i>lowestNset32bits[11]) {
	if (i>lowestNset32bits[13]) {
	  if (i>lowestNset32bits[14]) {
	    return 16;
	  } else {
	    return 15;
	  }
	} else {
	  if (i>lowestNset32bits[12]) {
	    return 14;
	  } else {
	    return 13;
	  }
	}
      } else {	    
	if (i>lowestNset32bits[9]) {
	  if (i>lowestNset32bits[10]) {
	    return 12;
	  } else {
	    return 11;
	  }
	} else {
	  if (i>lowestNset32bits[8]) {
	    return 10;
	  } else {
	    return 9;
	  }
	}
      }
    } else {
      if (i>lowestNset32bits[3]) {
	if (i>lowestNset32bits[5]) {
	  if (i>lowestNset32bits[6]) {
	    return 8;
	  } else {
	    return 7;
	  }
	} else {
	  if (i>lowestNset32bits[4]) {
	    return 6;
	  } else {
	    return 5;
	  }
	}
      } else {	    
	if (i>lowestNset32bits[1]) {
	  if (i>lowestNset32bits[2]) {
	    return 4;
	  } else {
	    return 3;
	  }
	} else {
	  if (i>lowestNset32bits[0]) {
	    return 2;
	  } else {
	    return 1;
	  }
	}
      }
    }
  }
}



/// Find size (on [1,32]) of the smallest # that can hold the unsigned int u.
/// By our convention, [0,1] require 1 bit, [2,3] require 2 bits, [4,5,
/// 6,7] require 3 bits, etc.
/// \param u The number whose size is being checked.
static inline unsigned int bit_size(unsigned int u) {
  for (int bs=1; bs<=32; bs++)
    if (u == (u&lowestNset32bits[bs]))
      return bs;
  throw "Bit size (unsigned int) fails!";
}


#endif  // #ifdef SLIM_BITSTREAM_H
