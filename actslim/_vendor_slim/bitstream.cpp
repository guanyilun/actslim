/// \file bitstream.cpp
/// Implement the abstract base class bitstream and derived classes
/// ibitstream and obitstream, for i/o of data with access to an
/// arbitrary (integer) number of bits at a time.

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


#include <cstdlib>
#include <cassert>
#include <climits>
#include <cstring>
#include <cerrno>
#include <iostream>
#include "bitstream.h"
#include "slim.h"

//#define DEBUG_BITSTREAM
#ifdef DEBUG_BITSTREAM
#include <cctype>
#endif


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// \class bitstream
/// Bit stream base class.
/// Allows you to R/W data one bit at a time with buffered reading/writing.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/// Dummy default constructor.
bitstream::bitstream() {
  buffer_base = NULL;
}


/*
/// Start bitstream using a FILE ptr to an open file.
bitstream::bitstream(FILE *file, int buffersize) 
{
  fp = file;
  bufsize = buffersize;
  setupstream();
}
*/


/// Start bitstream by filename
bitstream::bitstream(const char *filename, int buffersize) {;}


/*
/// Copy constructor
bitstream::bitstream(const bitstream &b) 
{
  bufsize = b.bufsize;
  fp = b.fp;
  setupstream();
}



/// Assignment operator
bitstream & 
bitstream::operator=(const bitstream &b)
{
  if (&b == this) return *this;
  bufsize = b.bufsize;
  fp = b.fp;
  setupstream();
  return *this;
}
*/


/// Destructor deletes output buffer, closes file.
bitstream::~bitstream()
{
    /*
  close();
  */
  delete [] buffer_base;
  buffer_base= NULL;
}


/// Close the IO file.
void ibitstream::close()
{
#ifdef HAVE_LIBZZIP
  if (zfp)
    zzip_fclose(zfp);
  zfp = NULL;
#endif
#ifdef HAVE_LIBLZ4
  if (lz4_ctx)
    LZ4F_freeDecompressionContext(lz4_ctx);
  lz4_ctx = NULL;
#endif
  if (fp) {
    rewind(fp);
    fclose(fp);
  }
  fp = NULL;
}



/// Is the IO file closed?
bool obitstream::is_open() const 
{
  return (fp != NULL);
}

/// Is the IO file closed?
bool ibitstream::is_open() const 
{
#ifdef HAVE_LIBZZIP
  return (zfp != NULL) || (fp != NULL);
#else
  return (fp != NULL);
#endif
}


/// Allocate a buffer and set up all pointers
void bitstream::setupstream()
{
  // Require buffer to be an integer # of Word_t units in size;
  if (bufsize % sizeof(Word_t))
    bufsize += sizeof(Word_t) - bufsize % sizeof(Word_t);

  if (bufsize > MAX_BITSTREAM_BUFSIZE)
    throw "Buffer size is too big.";

  buffer_base= new Byte_t [bufsize];
  memset(buffer_base, 0, bufsize);
  buf_used = 0;
  
  buffptr.Bptr = buffer_base;
  beyondbuffer = buffer_base + bufsize;
  bitptr = 0;
}



/// Return the number of bytes used so far in this stream.
int bitstream::get_bytes_used() {
  return sizeof(Byte_t)*(buf_used + buffptr.Bptr - buffer_base) +
    (bitptr / 8);
}



/// Get the position of the bitptr
int bitstream::get_bitptr() { 
  return bitptr;
}


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// \class obitstream
/// Output bit stream.
/// Allows you to write data one bit at a time with buffered writing.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/// Start obitstream using a FILE ptr to an open file.
obitstream::obitstream(FILE *file, int buffersize)
{
  fp = file;
  bufsize = buffersize;
  setupstream();
}

/// Start outputbitstream by filename.
obitstream::obitstream(const char *filename, int buffersize) 
{
  fp = fopen(filename, "wb");
  if (fp == NULL) {
    throw bad_output_file(filename, "writing");
  }
  bufsize = buffersize;
  setupstream();
}



// Destructor flushed output buffer and calls base destructor.
obitstream::~obitstream()
{
  close();
}



/// Close the output stream by flushing and closing the FILE *.
void obitstream::close()
{
  windup();
  const bool FLUSH_TRAILING_BITS = true;
  flush(FLUSH_TRAILING_BITS);
  if (fp)
    fclose(fp);
  fp = NULL;
}
  


/// Put the current word to the buffer with upper 0 bits as needed
void obitstream::windup()
{
  int byte_bits_used = bitptr % 8;
  if (byte_bits_used) {
    writebits(0, 8-byte_bits_used);
  }
}



/// Write data to the buffer.
/// \param nbits Number of bits to use
/// \param data Data to write (in the lowest nbits bits).
void obitstream::writebits(uint32_t data32, int nbits)
{
  Word_t data = Word_t(data32); // On a 32-bit machine, this is trivial.

  // Can we get our data into the current Word_t word only?
  if (bitptr+nbits < Bits_per_word) {
    Word_t firstdata = data & lowestNset[nbits];
    *buffptr.Dptr |= (firstdata<<bitptr);
    bitptr += nbits;
    return;
  } 

  //  We need to get our data from the current AND NEXT Word_t word.
  int n_first = Bits_per_word - bitptr; // guaranteed to be <= nbits

  // Do first bits.  It is guaranteed that n_first > 0 (assuming nbits > 0).
  *buffptr.Dptr |= (data<<bitptr);

  // Move to next word and flush buffer, if needed.
  buffptr.Dptr ++;
  bitptr = nbits - n_first; //Note okay to pre-set bitptr to save CPU, if flush(false)
  const bool IGNORE_TRAILING_BITS = false;
  if (buffptr.Bptr  >= beyondbuffer)
    flush(IGNORE_TRAILING_BITS);

  // Do last bits, if any.
  if (bitptr) {
    *buffptr.Dptr = (data>>n_first) & lowestNset[bitptr];
  }
}



/// Write 8-bit character strings to the buffer
/// \param str  The string to write. 0-terminated
/// \param write_trailing_null Whether to put the 0 terminator to the output.
void obitstream::writestring(const char *str, bool write_trailing_null) {
  for (; *str != '\0'; str++)
    writebits(*str-'\0', 8);
  if (write_trailing_null)
    writebits(0, 8);
}



/// Write all the bits of a word.
/// \param data  The word to write.
template <typename T> void obitstream::writeword(const T data)
{
  writebits(data, 8*sizeof(data));
}


/// Instantiation for several types
//template void obitstream::writeword<char>(const char data);
//template void obitstream::writeword<int>(const int data);
template void obitstream::writeword<uint32_t>(const uint32_t data);



/// Write a unary code for the value.
/// Code will be (value) 1s followed by a zero. 
/// Since low-order bits are "first" in our convention, this means
/// that the code for 4 is the 5-bit value 01111, interpreted by the
/// reader as 1,1,1,1,0.
/// \param value  The value to be coded.
void obitstream::write_unary(uint32_t value) {
  const uint32_t MAX_UNARY=1024;

  if (value > MAX_UNARY)
    throw "Attempted to write too large a unary value.";

  for (; value>=(uint32_t)Bits_per_word; value -= Bits_per_word)
    writebits(~0, Bits_per_word);

  /// The length of the code for (value) is 1+value bits.
  Word_t rl_code = lowestNset32bits[value]; 
  writebits(rl_code, value+1);
}



/// Print buffer contents.
void obitstream::print() const
{
  if (buffer_base== NULL) {
    cout << "No buffer allocated\n";
    return;
  }
  Word_t *d;
  for (Byte_t *p=buffer_base; 
       p<buffptr.Bptr || (p==buffptr.Bptr && bitptr>0); p+=4) {
    d = reinterpret_cast<Word_t *>(p);
    cout.width(5);
    cout << p-buffer_base << ": ";
    cout.width(8);
    cout << *d << " = ";
    int bits_to_write = 32;
    if (p == buffptr.Bptr)
      bits_to_write = bitptr;
    for (int i=0; i<bits_to_write; i++) {
      cout.width(1);
      cout << ((*d & lowestNset[i+1])>>i);
      if (i%4 == 3)
        cout << " ";
    }
    cout << endl;
  }
}



/// Flush the write buffer and reset for more data.
/// Do we need to flush even the partial words, as at EOF?
void obitstream::flush(bool flush_trailing_bits)
{
  int thiswrite;
  while (bitptr > 0 && flush_trailing_bits) {
    buffptr.Bptr ++;
    bitptr -= 8*sizeof(Byte_t);
  }

  thiswrite = fwrite(buffer_base, sizeof(Byte_t), 
		     (buffptr.Bptr-buffer_base), fp);
  buf_used += thiswrite;

  if (flush_trailing_bits)
    bitptr = 0;
  buffptr.Bptr = buffer_base;
  if (thiswrite)
    memset(buffer_base, 0, thiswrite); // Clear the output buffer for re-use.
}




// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// \class ibitstream
/// Input bit stream.
/// Allows you to read data N bits at a time with buffered reading.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/// Start ibitstream using a FILE ptr to an open file.
ibitstream::ibitstream(FILE *file, int buffersize)
{
  fp = file;
#ifdef HAVE_LIBZZIP
  zfp = NULL;
#endif
#ifdef HAVE_LIBLZ4
  uint32_t magic_bytes = 0;
  using_lz4 = false;
  lz4_ctx = NULL;
  if (fread(&magic_bytes, sizeof(Byte_t), 4, fp) >= 4) {
    if (magic_bytes == 0x184d2204) {
      size_t c = LZ4F_createDecompressionContext(&lz4_ctx, LZ4F_VERSION);
      if (LZ4F_isError(c))
        throw bad_output_file("unknown",
          ("LZ4 error: " + string(LZ4F_getErrorName(c))).c_str());
      using_lz4 = true;
    }
  }
  fseek(fp, 0, SEEK_SET);
#endif
  bufsize = buffersize;
  setupstream();
}


#ifdef HAVE_LIBZZIP
/// Start ibitstream using a FILE ptr to an open file.
ibitstream::ibitstream(ZZIP_FILE *file, int buffersize)
{
  fp = NULL;
  zfp = file;
  bufsize = buffersize;
  setupstream();
}
#endif


/// Start inputbitstream by filename.
ibitstream::ibitstream(const char *filename, int buffersize) 
{
#ifdef HAVE_LIBZZIP
  fp = NULL;
  zfp = zzip_fopen(filename, "rb");
  if (zfp == NULL) {
#else
  fp = fopen(filename, "rb");
  if (fp == NULL) {
#endif
    throw bad_output_file(filename, "reading");
  }
#ifdef HAVE_LIBLZ4
  uint32_t magic_bytes = 0;
  using_lz4 = false;
  lz4_ctx = NULL;
#ifdef HAVE_LIBZZIP
  if (zzip_fread(&magic_bytes, sizeof(Byte_t), 4, zfp) >= 4) {
#else
  if (fread(&magic_bytes, sizeof(Byte_t), 4, fp) >= 4) {
#endif
    if (magic_bytes == 0x184d2204) {
      size_t c = LZ4F_createDecompressionContext(&lz4_ctx, LZ4F_VERSION);
      if (LZ4F_isError(c))
        throw bad_output_file("unknown",
          ("LZ4 error: " + string(LZ4F_getErrorName(c))).c_str());
      using_lz4 = true;
    }
  }
#ifdef HAVE_LIBZZIP
  zzip_seek(zfp, 0, SEEK_SET);
#else
  fseek(fp, 0, SEEK_SET);
#endif
#endif
  bufsize = buffersize;
  partial_word = 0;
  partial_word_bitptr = -1;
  setupstream();
}

/// Start inputbitstream by file descriptor.
ibitstream::ibitstream(int fd, int buffersize)
{
  fp = fdopen(fd, "rb");
  if (fp == NULL) {
    char *fdname = new char[14];
    snprintf(fdname, 14, "fd#%i", fd);

    throw bad_output_file(fdname, "reading");
  }
  bufsize = buffersize;
  partial_word = 0;
  partial_word_bitptr = -1;
  setupstream();
}



/// Destructor only uses base class destructor.
ibitstream::~ibitstream()
{
  close();
  delete [] buffer_base;
  buffer_base= NULL;
#ifdef HAVE_LIBLZ4
  if (lz4_buffer)
    delete [] lz4_buffer;
  lz4_buffer = NULL;
#endif
#ifdef DEBUG_READBITS
  print_debug();
}



/// Debug message to print when ibitstream is destroyed.
void ibitstream::print_debug() {
  cout.width(8); cout << c1 << " Total calls to readbits\n";
  cout.width(8); cout << c2 << " Calls with last bits\n";
  cout.width(8); cout << c3 << " Calls with first bits\n";
  cout.width(8); cout << c4 << " Buffer fills\n";
  cout.width(8); cout << c5 << " Dangling bytes and recursive calls\n";
  cout.width(8); cout << c6 << " Finish routine\n";
#endif
}



/// Return the number of bytes used so far in this stream.
int ibitstream::get_bytes_used() {
  return      sizeof(Byte_t)*(buf_used + buffptr.Bptr - beyondbuffer) +
    (bitptr / 8);
}



/// Allocate a buffer and set up all pointers, then fill buffer.
void ibitstream::setupstream()
{
  bitstream::setupstream();
#ifdef HAVE_LIBLZ4
  lz4_buffer = new Byte_t [bufsize];
  memset(lz4_buffer, 0, bufsize);
#endif
  fill();
}


/// Move to next word in buffer, refilling it if needed.
inline void ibitstream::next_word() {
  buffptr.Dptr ++;

  // Careful: the bitptr should now be set to zero, but we defer this
  // for optimization purposes.

  if (buffptr.Bptr  >= beyondbuffer) {
    fill();
  }
}



/// Put the current word to the buffer with upper 0 bits as needed
void ibitstream::windup()
{
  int byte_bits_used = bitptr % 8;
  if (byte_bits_used) {
    readbits(8-byte_bits_used);
  }
}



/// Read data from the buffer as unsigned ints.
/// \param n_bits Number of bits to use
/// \return Data read (in the lowest n_bits bits).
Word_t ibitstream::readbits(int n_bits)
{
  Word_t data=0;
  // If there are enough bits in this word, be quick about it.
  if (bitptr + n_bits <= Bits_per_word) {
    data = (*buffptr.Dptr >> bitptr) & lowestNset[n_bits];
    bitptr += n_bits;
    return data;
  }
  // Get n_bits_this from current word, n_bits_next from the next word
  int n_bits_this = Bits_per_word - bitptr;
  int n_bits_next = n_bits - n_bits_this;
  if (n_bits_this > 0)
    data = (*buffptr.Dptr >> bitptr) & lowestNset[n_bits_this];
  bitptr = 0;
  next_word();
  if (bitptr==0) {
      data |= ((*buffptr.Dptr) & lowestNset[n_bits_next]) << n_bits_this;
      bitptr = n_bits_next;
  } else {
      data |= ((*buffptr.Dptr >> bitptr) & lowestNset[n_bits_next]) << n_bits_this;
      bitptr += n_bits_next;
  }
  if (bitptr > Bits_per_word)
    throw "Out of data in ibitstream";
  return data;
}



/// Read a null-terminated string of 8-bit characters from the bit stream
/// \param s     Pointer to the string (must be available memory).
/// \param count Maximum size of the string.
/// \return      Number of chars read (including null).
int ibitstream::readstring(char *s, int count) {
  if (count < 0)
    count = INT_MAX;
  int i;
  for (i=0; i<count; i++, s++) {
    *s = (char)readbits(8);
    if (*s == '\0')
      break;
  }
  return i;
}



/// Read data from the buffer as (signed) ints.
/// \param n_bits Number of bits to read.
/// \return  The data read from the stream.
int ibitstream::readbits_int(int n_bits) {
  int data = int(readbits(n_bits));
  const int size = 8*sizeof(int);
  // If the number was < 32 bits, then fill top with sign bits
  if (n_bits < size) {
    data <<= (size-n_bits);
    data >>= (size-n_bits);
  }
  return data;
}



/// Read a single unary-coded value.
/// \return The unary-coded value from the stream.
Word_t ibitstream::read_unary() {
  Word_t value = 0;
  Word_t datum = 0;
  do {
    datum = readbits(1);
    value++;
  } while (datum);
  value --;  // N times through loop means value = N-1
  return value;
}



/// Print buffer contents.
void ibitstream::print() const
{
  if (buffer_base== NULL) {
    cout << "No buffer allocated\n";
    return;
  }
  for (Byte_t *p=buffer_base; p<buffptr.Bptr || 
	 (p==buffptr.Bptr && bitptr>0); p++) {
    cout.width(5);
    cout << p-buffer_base << ": ";
    cout.width(8);
    cout << *p << " = ";
    //for (int i=31; i>=0; i--) ;
    for (int i=0; i<32; i++) {
      cout.width(1);
      cout << ((*p & lowestNset[i+1])>>i);
      if (i%4 == 3)
        cout << " ";
    }
    cout << endl;
  }
}

/// Fill the read buffer and reset for more data.
/// \return  Bytes read from disk.
int ibitstream::fill()
{
  int thisread;
#ifdef HAVE_LIBZZIP
  if (zfp) {
#ifdef HAVE_LIBLZ4
    if (using_lz4)
      thisread = zzip_fread(lz4_buffer, sizeof(Byte_t), bufsize, zfp);
    else
#endif
    thisread = zzip_fread(buffer_base, sizeof(Byte_t), bufsize, zfp);
  } else if (fp) {
#endif
#ifdef HAVE_LIBLZ4
  if (using_lz4)
    thisread = fread(lz4_buffer, sizeof(Byte_t), bufsize, fp);
  else
#endif
    thisread = fread(buffer_base, sizeof(Byte_t), bufsize, fp);
#ifdef HAVE_LIBZZIP
  } else {
    throw bad_output_file("unknown", "reading");
  }
#endif
#ifdef HAVE_LIBLZ4
  if (using_lz4) {
    if (thisread > 0) {
      size_t to_read = thisread;
      size_t this_read = bufsize;
      size_t c = LZ4F_decompress(lz4_ctx, buffer_base, &this_read, lz4_buffer,
        &to_read, NULL);
      if (LZ4F_isError(c))
        throw bad_output_file("unknown",
          ("LZ4 error: " + string(LZ4F_getErrorName(c))).c_str());
#ifdef HAVE_LIBZZIP
      if (zfp)
        zzip_seek(zfp, (int) to_read - thisread, SEEK_CUR);
      else
#endif
      fseek(fp, (int) to_read - thisread, SEEK_CUR);
      thisread = (int) this_read;
    }
  }
#endif

  // Handle partial words at end by saving them for a future, final call.
  if (thisread % sizeof(Word_t) != 0) {
    int partial_word_offset = thisread - (thisread % sizeof(Word_t));
    partial_word = *(Word_t *)(buffer_base + partial_word_offset);
    
    // Put the partial word in the uppermost bytes of the Word_t word.
    partial_word_bitptr = 8*(sizeof(Word_t) - thisread % sizeof(Word_t));
    partial_word <<= partial_word_bitptr;

    thisread -= (thisread % sizeof(Word_t));
  }

  // If no data read, see whether we saved any partial words last time.
  if (thisread == 0) {
    if (partial_word_bitptr >= 0) {
      buf_used += sizeof(Word_t);
      memcpy(buffer_base, &partial_word, sizeof(Word_t));
      buffptr.Bptr = buffer_base;
      beyondbuffer = buffer_base + sizeof(Word_t);
      bitptr = partial_word_bitptr;
      partial_word_bitptr = -1;
      return sizeof(Word_t)-bitptr/8;
    } else
      throw "Out of data in ibitstream";
  }
   
  buf_used += thisread;
  buffptr.Bptr = buffer_base;
  beyondbuffer = buffer_base + thisread;
  bitptr = 0;

#ifdef DEBUG_READBITS
  cout << "Filled buffer with " << thisread << " bytes (" << buf_used
       << " total)\n";
#endif
  return thisread;
}
