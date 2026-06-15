/// \file raw_section.cpp
/// Implements class raw_section, a buffer with easy random
/// access for parameter computing in the encoding stage.

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


#include "slim.h"
#include "crc.h"   // for CRC32 table 
#include <cstring> // for memcpy


// ---------------------------------------------------------------------------
///\class raw_section
/// Read/write buffer for the *raw* data I/O.
/// Allows random access for use in the sampling = parameter computing
/// part of the encoding.
// ---------------------------------------------------------------------------



/// Constructor
/// \param mode_in Whether buffer is write-only or read-only
raw_section::raw_section(enum section_mode_t mode_in) :
mode(mode_in) {
  private_buf_size = buf_size = 0;
  private_buf = buf = NULL;

  n_chan_alloc = 10;
  n_chan_used = 0;
  chan_reps = new unsigned int[n_chan_alloc];
  chan_offset = new int[n_chan_alloc];
  frame_size = 0;
}



/// Destructor.
raw_section::~raw_section() {
  delete [] private_buf;
  delete [] chan_reps;
  delete [] chan_offset;
}



/// Add a channel to the frame by putting repetitions and size onto the lists.
/// Also update the running total of the frame size.
/// \param reps  Repetitions of this channel.
/// \param chan_size Size of one data value, in bytes.
void raw_section::add_channel(int reps, int chan_size) {

  // If all channels are used, double the allocation.
  if (n_chan_used >= n_chan_alloc) {
    n_chan_alloc *= 2;
    unsigned int *tmp_reps = new unsigned int[n_chan_alloc];
    int *tmp_offset = new int[n_chan_alloc];

    // Copy old data into the new buffers
    int n_chan_copy = n_chan_used;
    if (n_chan_copy > n_chan_alloc)
      n_chan_copy = n_chan_alloc;
    for (int i=0; i<n_chan_copy; i++) {
      tmp_reps[i] = chan_reps[i];
      tmp_offset[i] = chan_offset[i];
    }
    delete [] chan_reps;
    delete [] chan_offset;

    chan_reps = tmp_reps;
    chan_offset = tmp_offset;
    assert(n_chan_used < n_chan_alloc);
  }
  chan_reps[n_chan_used] = reps;
  chan_offset[n_chan_used] = frame_size;
  n_chan_used++;
  frame_size += chan_size*reps;
}



/// Clear the list of channels (repetitions and offsets), without
/// reallocating the raw section data buffer.
void raw_section::reset_channels() {
  n_chan_used = 0;
  frame_size = 0;
}



/// Fill the section buffer from a file.
/// Allowed only in compress mode, not decompress mode.
/// \param fp Open readable FILE pointer.
/// \param size Number of bytes to read.
/// \return     Number of bytes read.
size_t raw_section::fill(FILE *fp, size_t size) {
  
  if (mode != SECTION_COMPRESS_MODE)
    throw "Can only call raw_section::fill() for a COMPRESS section.";

  if (buf_size != size) {
    resize(size);
  }

  if (buf == NULL)
    throw "Must resize with raw_section::resize before calling ::fill";

  size_t n = fread(buf, 1, size, fp);
  return n;
}



/// Resize the I/O buffer to hold a given number of bytes.
/// Preserves the buffer's data up to the minimum of the old/new sizes.
/// Note that we are careful to add MAX_GHOST_BYTES of padding at the
/// end for handling partial words at the end of file.
/// \param size The number of bytes in the raw section.
/// The size need not be a multiple of the frame_size.
/// Return the new size in bytes.
size_t raw_section::resize(size_t size) {

  // Only resize if using private buffer.
  // If we let user resize private while using an external 
  // buffer, it would surely cause confusion.
  if (private_buf != buf)
    throw "Cannot resize a user buffer.";

  // No change.
  if (size == private_buf_size)
    return size;

  // Setting to zero size.
  if (size <= 0) {
    private_buf = buf = NULL;
    num_frames = 0;
    return 0;
  }

  // Reducing the size.
  if (size < private_buf_size) {
    private_buf_size = buf_size = size;
    num_frames = 1+buf_size / frame_size;
    if (num_frames == 0)
      num_frames = 1;
    return size;
  }

  // Increasing the size.
  if (size > MAX_SECTION_LENGTH)
    throw "Cannot resize buffer beyond MAX_SECTION_LENGTH";

  // Make a new buffer and copy the old to it (up to the old size).
  Buffer_t *new_buf = new Buffer_t[size + MAX_GHOST_BYTES];
  if (private_buf)
    memcpy(new_buf, private_buf, buf_size+MAX_GHOST_BYTES);
  delete [] private_buf;

  private_buf = buf = new_buf;
  private_buf_size = buf_size = size;
  num_frames = buf_size / frame_size;
  return buf_size;
}



/// Switch to using a buffer from user space.
/// This only makes sense if the section is for compression, so assert
/// that mode.  Also, a section to be written to disk can safely
/// cast away the const-ness of the buffer, because copying to disk is
/// the only thing that will happen.
/// \param buffer Pointer to the user's buffer.
/// \param length Size of the user's buffer.
void raw_section::use_external_buffer(const unsigned char *buffer, size_t length) {
  if (mode != SECTION_COMPRESS_MODE)
    throw "Can only use a external (user) buffer in COMPRESS mode.";
  buf = const_cast<unsigned char *>(buffer);
  buf_size = length;
  num_frames = buf_size / frame_size;
}



/// Switch to using the section's internal (private) buffer.
void raw_section::use_internal_buffer() {
  buf = private_buf;
  buf_size = private_buf_size;
}



/// Flush the section buffer to a file.
/// Allowed only in expand mode, not compress mode.
/// \param fp Open writeable FILE pointer.
/// \param size Number of bytes to read.
/// \return     Number of bytes read.
size_t raw_section::flush(FILE *fp, size_t size) {
  if (mode != SECTION_EXPAND_MODE)
    throw "Can only call raw_section::fill() for an EXPAND section.";
  if (size == 0)
    size = num_frames * frame_size;
  if (size <= 0 || num_frames <= 0)
    return 0;

  return fwrite(buf, 1, size, fp);
}



/// Compute the CRC-32-IEE 802.3 value.
/// \param  length  Compute CRC on the first _length_ bytes only (0 means 
///                 compute for full section buffer).
/// \return The cyclic redundancy check.
unsigned long raw_section::crc(size_t length) const {
  if (length <= 0)
    length = buf_size;

  const unsigned char *s = buf;
  register unsigned long c = 0xffffffffL;

  if (length) do {
    c = slim_crc_32_table[((int)c ^ (*s++)) & 0xff] ^ (c >> 8);
  } while (--length);

  return c ^ 0xffffffffL;       /* (instead of ~c for 64-bit machines) */
}
  


/// Set num_frames.  
/// This is not meant to be primary: it's only used in ::flush when
/// called with a zero size argument.
/// \param i New number of frames.
/// \return  New number of frames.
int raw_section::set_num_frames(int i) {
  return (num_frames = i);
}



/// Return (as int) the given datum.
/// This version is slower than the 3-argument version, but it's appropriate
/// for the data sampling phase, where the outside world doesn't know or
/// want to know the data point # as a (frame, point) pair.
/// \param ichan  The channel number.
/// \param idat   The data number within the section.
/// \return The data value, as int.  
int32_t& raw_section::ival(int ichan, int idat) {
  int reps = chan_reps[ichan];
  int iframe = idat/reps;
  return reinterpret_cast<int32_t &>
    (buf[iframe*frame_size+chan_offset[ichan]+(idat%reps)*sizeof(int32_t)]);
}



/// Return (as int) a reference to the given datum.
/// This version is faster than the 2-argument version and should be used
/// during the encoding/decoding phase.
/// \param ichan  The channel number.
/// \param iframe The frame number.
/// \param i_inframe The data number within the frame.
/// \return Reference to the data value, as int.  
///         Can be used as an lvalue.
int32_t& raw_section::ival(int ichan, int iframe, int i_inframe) {
  return reinterpret_cast<int32_t &>
    (buf[iframe*frame_size+chan_offset[ichan]+i_inframe*sizeof(int32_t)]);
}



/// Return (as int16_t) the given datum.
/// This version is slower than the 3-argument version, but it's appropriate
/// for the data sampling phase, where the outside world doesn't know or
/// want to know the data point # as a (frame, point) pair.
/// \param ichan  The channel number.
/// \param idat   The data number within the section.
/// \return The data value, as short int.  
int16_t& raw_section::sval(int ichan, int idat) {
  int reps = chan_reps[ichan];
  int iframe = idat/reps;
  return reinterpret_cast<int16_t &>
    (buf[iframe*frame_size+chan_offset[ichan]+(idat%reps)*sizeof(int16_t)]);
}



/// Return (as int16_t) a reference to the given datum.
/// \param ichan  The channel number.
/// \param iframe The frame number.
/// \param i_inframe The data number within the frame.
/// \return Reference to the data value, as int16_t.  
///         Can be used as an lvalue.
int16_t& raw_section::sval(int ichan, int iframe, int i_inframe) {
  return reinterpret_cast<int16_t &>
    (buf[iframe*frame_size+chan_offset[ichan]+i_inframe*sizeof(int16_t)]);
}



/// Return (as int8_t) the given datum.
/// This version is slower than the 3-argument version, but it's appropriate
/// for the data sampling phase, where the outside world doesn't know or
/// want to know the data point # as a (frame, point) pair.
/// \param ichan  The channel number.
/// \param idat   The data number within the section.
/// \return The data value, as short int.  
int8_t& raw_section::cval(int ichan, int idat) {
  int reps = chan_reps[ichan];
  int iframe = idat/reps;
  return reinterpret_cast<int8_t &>
    (buf[iframe*frame_size+chan_offset[ichan]+(idat%reps)*sizeof(int8_t)]);
}



/// Return (as int8_t) a reference to the given datum.
/// \param ichan  The channel number.
/// \param iframe The frame number.
/// \param i_inframe The data number within the frame.
/// \return Reference to the data value, as int8_t.  
///         Can be used as an lvalue.
int8_t& raw_section::cval(int ichan, int iframe, int i_inframe) {
  return reinterpret_cast<int8_t &>
    (buf[iframe*frame_size+chan_offset[ichan]+i_inframe*sizeof(int8_t)]);
}



/// Return (as uint32_t) the given datum.
/// This version is slower than the 3-argument version, but it's appropriate
/// for the data sampling phase, where the outside world doesn't know or
/// want to know the data point # as a (frame, point) pair.
/// \param ichan  The channel number.
/// \param idat   The data number within the section.
/// \return The data value, as uint32_t.  
uint32_t& raw_section::uval(int ichan, int idat) {
  int reps = chan_reps[ichan];
  int iframe = idat/reps;
  return reinterpret_cast<uint32_t &>
    (buf[iframe*frame_size+chan_offset[ichan]+(idat%reps)*sizeof(uint32_t)]);
}



/// Return (as uint32_t) a reference to the given datum.
/// \param ichan  The channel number.
/// \param iframe The frame number.
/// \param i_inframe The data number within the frame.
/// \return Reference to the data value, as uint32_t. 
///         Can be used as an lvalue.
uint32_t& raw_section::uval(int ichan, int iframe, int i_inframe) {
  return reinterpret_cast<uint32_t &>
    (buf[iframe*frame_size+chan_offset[ichan]+i_inframe*sizeof(uint32_t)]);
}



/// Return a pointer to a channel's area of one data frame
/// \param ichan The channel number.
/// \param iframe The frame number.
/// \return      Points the the requested area of the section buffer.
unsigned char * raw_section::ptr(int ichan, int iframe) const {
  return &buf[iframe*frame_size+chan_offset[ichan]];
}
