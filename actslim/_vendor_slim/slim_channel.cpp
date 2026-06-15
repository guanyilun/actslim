/// \file slim_channel.cpp
/// Implements classes slim_channel_array, abstract base slim_channel,
/// and its derived classes slim_channel_encode and slim_channel_decode.

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
#include <cstring>

#include "slim.h"
#include "bitstream.h"

// ----------------------------------------------------------------------
/// \class slim_channel_array
/// Container for all slim_channel objects in a file
// ----------------------------------------------------------------------

/// Constructor.
/// \param n_alloc  How long to make the internal arrays initially.
slim_channel_array::slim_channel_array(int n_alloc) {
  offsets_in_frame = new int[n_alloc];
  chan_array = new slim_channel *[n_alloc];
  assert (offsets_in_frame != NULL && chan_array != NULL);
  num_allocated = n_alloc;

  memset(offsets_in_frame, 0, n_alloc*sizeof(int));
  memset(chan_array, 0, n_alloc*sizeof(slim_channel *));

  num_chan = n_alloc; // this will be corrected when we clear.
  clear();
}



/// Destructor.
slim_channel_array::~slim_channel_array() {
  for (int i=0; i<num_chan; i++)
    delete chan_array[i];
  delete [] chan_array;
  delete [] offsets_in_frame;
}



/// Find the offset within a frame of (the first rep of) a channel.
/// \param i The number of the channel to test.
/// \return  Offset in the frame (in bytes).
int slim_channel_array::offset(int i) const {
  if (i<0 || i >= num_chan)
    return -1;
  return offsets_in_frame[i];
}



/// Add a channel onto the list.
/// We now own this channel, and it will be deleted when this object is.
/// Will grow the existing arrays, if needed.
/// \param c  Pointer to the channel to add.  (Will take ownership of it.)
/// \param frame_offset  The offset of this channel in a frame.
void slim_channel_array::push(slim_channel *c, size_t frame_offset) {

  // Grow the arrays, if needed.
  if (num_chan == num_allocated)
    resize_arrays(num_chan * 2);
  assert (num_chan < num_allocated);

  // Update the circular linked list pointers
  if (num_chan > 0) {
    c->next_chan = chan_array[0];  // new one points to first in list
    chan_array[num_chan-1]->next_chan = c;  // Old last one points to new one
  } else
    c->next_chan = c;

  // Store the channel and its offset in the appropriate arrays
  chan_array[num_chan] = c;
  offsets_in_frame[num_chan] = frame_offset;
  num_chan++;
}



/// Implement a version of realloc for the arrays we keep.
/// Copies all existing data to the minimum of the new size and num_chan.
/// \param n The desired size of the arrays. 
void slim_channel_array::resize_arrays(int n) {
  assert (n > num_chan);
  int *tmp_offsets = new int[n];
  slim_channel **tmp_chan = new slim_channel *[n];
  assert (tmp_offsets != NULL && tmp_chan != NULL);

  // Copy data.  # elements = min(n, num_chan).
  int n_copy = num_chan;
  if (n < num_chan)
    n_copy = n;
  for (int i=0; i<n_copy; i++) {
    tmp_offsets[i] = offsets_in_frame[i];
    tmp_chan[i] = chan_array[i];
  }

  // Free old arrays and point to the new ones.
  delete [] offsets_in_frame;
  delete [] chan_array;
  
  offsets_in_frame = tmp_offsets;
  chan_array = tmp_chan;
  num_allocated = n;
}



/// The array subscript operator gives access to the channels.
/// \param i  Array element desired.
/// \return   The slim_channel pointer for index i, or NULL if bad index.
slim_channel *slim_channel_array::operator[](int i) const
{
  if (i < 0 || i >= num_chan)
    return NULL;
  return chan_array[i];
}



/// Delete all outstanding decoder channels, emptying the list.
void slim_channel_array::clear() {
  for (int i=0; i<num_chan; i++)
    delete chan_array[i];
  for (int i=0; i<num_chan; i++) {
    chan_array[i] = NULL;
    offsets_in_frame[i] = 0;
  }
  num_chan = 0;
}



// ----------------------------------------------------------------------
/// \class slim_channel
/// Base class for the encode/decode derived channel classes.
/// A channel represents one stream of numbers from a single instrument,
/// and its properties are assumed to be consistent over time, but not
/// the same as those of other channels.
// ----------------------------------------------------------------------

//----------------------------------------------------------------------
/// Constructor.
/// \param reps  The number of repetitions of data per frame.
/// \param size  The size (bytes) of a single data word.
/// \param deltas Whether this channel will encode data or its deltas.
/// \param permit Whether this channel will allow bit-rotation.
slim_channel::slim_channel(unsigned int reps, 
			   int size, bool deltas) :
  repetitions(reps),
  raw_size(size),
  frame_size(size*reps),
  encode_deltas(deltas)
{
  bit_rotation = bit_unrotation = 0;
  next_chan = NULL;
}

//----------------------------------------------------------------------
/// Destructor does nothing.
//----------------------------------------------------------------------
slim_channel::~slim_channel() {;}



// ----------------------------------------------------------------------
/// \class slim_channel_encode
/// Derived class for encoding the data from a single channel.
/// A channel represents one stream of numbers from a single instrument,
/// and its properties are assumed to be consistent over time, but not
/// the same as those of other channels.
// ----------------------------------------------------------------------

//----------------------------------------------------------------------
/// Constructor
/// \param reps  The number of repetitions of data per frame.
/// \param size  The size (bytes) of a single data word.
/// \param deltas Whether this channel will encode data or its deltas.
/// \param permit Whether this channel will allow bit-rotation.
//----------------------------------------------------------------------
slim_channel_encode::slim_channel_encode(int reps, int size, 
					 bool deltas, bool permit) :
  slim_channel(reps, size, deltas),
  permit_rotation(permit)
{
  ndata_sampled = 0;
  enc = NULL;
  usual_encoder = NULL;
}



//----------------------------------------------------------------------
/// Destructor deletes the encoder object.
//----------------------------------------------------------------------
slim_channel_encode::~slim_channel_encode() {
  delete enc;
  delete usual_encoder;
}



//----------------------------------------------------------------------
/// Clear the "history" value.
//----------------------------------------------------------------------
void slim_channel_encode::reset_previous() {
  if (enc) 
    enc->reset_previous();
  ndata_sampled = 0;
}



//----------------------------------------------------------------------
/// Find how many low-order bits are constant through the data sample.
/// \param data  The data sample array.
/// \param ndata Number of elements in the array.
/// \return The number of lowest-order bits that do not vary in the sample.
//----------------------------------------------------------------------
template <typename T>
int slim_channel_encode::constant_low_bits(const T *data, 
					   int ndata) const {
  T bor = 0u;
  T band = ~T(0u);
  for (int i=0; i<ndata; i++) {
    bor |= data[i];
    band &= data[i];
  }
  T varying_bits = bor ^ band;
  // When no bits ever change, this is going to be a constant channel.
  if (varying_bits == T(0) )
    return 0;

  for (int i=8*sizeof(T)-1; i>= 0; i--)
    if ((varying_bits & lowestNset[i]) == 0)
      return i;
  return 0;
}

/// Instantiation for an array of uint32_t
template int slim_channel_encode::constant_low_bits<uint32_t>
(const uint32_t *data, int ndata) const;

/// Instantiation for an array of uint16_t
template int slim_channel_encode::constant_low_bits<uint16_t>
(const uint16_t *data, int ndata) const;

/// Instantiation for an array of uint8_t
template int slim_channel_encode::constant_low_bits<uint8_t>
(const uint8_t *data, int ndata) const;



//----------------------------------------------------------------------
/// Rotate data word by bit_rotation bits.
/// \param u The data to rotate.
/// \return  The rotated value.
//----------------------------------------------------------------------
inline uint32_t slim_channel_encode::rotate(uint32_t u) const {
  return (u>>bit_rotation) ^ (u<<bit_unrotation);
}



//----------------------------------------------------------------------
/// Rotate data word by bit_rotation bits.
/// \param i The data to rotate.
/// \return  The rotated value.
//----------------------------------------------------------------------
inline uint32_t slim_channel_encode::rotate(int32_t i) const {
  return (i>>bit_rotation) ^ (i<<bit_unrotation);
}



//----------------------------------------------------------------------
/// Pass an array of channel data to the encoder to compute its parameters.
/// \param data  The array.
/// \param ndata Number of elements in the array.
/// \return Negative error code.  0=no error
//----------------------------------------------------------------------
template <typename T>
int slim_channel_encode::compute_params(T *data, int ndata) {
  if (!enc || ndata < MIN_SAMPLES) {
    return -1;
  }
  ndata_sampled = ndata;
  if (encode_deltas)
    enc->use_signed_data_type();

  if (permit_rotation) {
    bit_rotation = constant_low_bits(data, ndata);
    bit_unrotation = 8*raw_size - bit_rotation;
    if (bit_rotation) {
      if (enc->is_signed()) {
	for (int i=0; i<ndata; i++)
	  data[i] = rotate(int32_t(data[i]));
      } else
	for (int i=0; i<ndata; i++)
	  data[i] = rotate(data[i]);
    }
  } else
    bit_rotation = bit_unrotation = 0;
  return enc->compute_params(data, ndata);
}

/// Instantiation for an array of uint32_t
template int slim_channel_encode::compute_params<uint32_t>(uint32_t *data, int ndata);

/// Instantiation for an array of uint16_t
template int slim_channel_encode::compute_params<uint16_t>(uint16_t *data, int ndata);

/// Instantiation for an array of uint8_t
template int slim_channel_encode::compute_params<uint8_t>(uint8_t *data, int ndata);



//----------------------------------------------------------------------
/// Have encoder object write its parameters.
//----------------------------------------------------------------------
int slim_channel_encode::write_params() const {
  unsigned int delta_code = encode_deltas? 1 : 0;
  ob->writebits(delta_code, 1);
  ob->writebits(bit_rotation, BITS_SLIM_NBITS);

  int retval = enc->write_params();
  return retval;
}



//----------------------------------------------------------------------
/// Encode one frame of data from a buffer to a file.
/// \param buf  The (allocated) buffer containing the raw data.
/// \return Number of raw words actually encoded.
//----------------------------------------------------------------------
size_t slim_channel_encode::encode_frame(void *buf) {

  uint32_t *dptr;
  uint16_t *sptr;
  uint8_t *cptr;

  /// For chans with only one rep per frame, faster to handle specially.
  if (repetitions == 1) {
    return encode_frame_singlevalue(buf);
  }

  unsigned int i;
  unsigned int nwords = repetitions;

  switch (raw_size) {
  case 4:
    dptr = reinterpret_cast<uint32_t *>(buf);
    if (bit_rotation)
      for (i=0; i<nwords; i++)
	dptr[i] = rotate(dptr[i]);
    enc->encode_vector(dptr, nwords);
    return nwords*4;
    break;

  case 2:
    sptr = reinterpret_cast<uint16_t *>(buf);
    if (bit_rotation)
      for (i=0; i<nwords; i++)
	sptr[i] = rotate(sptr[i]);
    enc->encode_vector(sptr, nwords);
    return nwords*2;
    break;

  case 1:
    cptr = reinterpret_cast<uint8_t *>(buf);
    if (bit_rotation)
      for (i=0; i<nwords; i++)
	cptr[i] = rotate(cptr[i]);
    enc->encode_vector(cptr, nwords);
    return nwords;
    break;

  default:
    cerr << "Oops: channel has raw_size=" << raw_size << "\n";
    assert (raw_size == 4 || raw_size == 2 || raw_size == 1);
    return 0;
    break;
  }
  return 0;
}



//----------------------------------------------------------------------
/// Encode one frame of data from a buffer to a file.
/// \param buf  The (allocated) buffer containing the raw data.
/// \param size The maximum number of bytes to read from the buffer 
///             (may use less).
/// \return Number of raw words actually encoded.
//----------------------------------------------------------------------
size_t slim_channel_encode::encode_partial_frame(void *buf, size_t size) {

  uint32_t *dptr;
  uint16_t *sptr;
  uint8_t *cptr;

  /// For chans with only one rep per frame, faster to handle specially.
  if (repetitions == 1 || size == raw_size) {
    return encode_frame_singlevalue(buf);
  }

  unsigned int i;
  unsigned int nwords;
  if (size > frame_size)
    nwords = repetitions;
  else
    nwords = size/raw_size;

  switch (raw_size) {
  case 4:
    dptr = reinterpret_cast<uint32_t *>(buf);
    if (bit_rotation)
      for (i=0; i<nwords; i++)
	dptr[i] = rotate(dptr[i]);
    enc->encode_vector(dptr, nwords);
    return nwords*4;
    break;

  case 2:
    sptr = reinterpret_cast<uint16_t *>(buf);
    if (bit_rotation)
      for (i=0; i<nwords; i++)
	sptr[i] = rotate(sptr[i]);
    enc->encode_vector(sptr, nwords);
    return nwords*2;
    break;

  case 1:
    cptr = reinterpret_cast<uint8_t *>(buf);
    if (bit_rotation)
      for (i=0; i<nwords; i++)
	cptr[i] = rotate(cptr[i]);
    enc->encode_vector(cptr, nwords);
    return nwords;
    break;

  default:
    cerr << "Oops: channel has raw_size=" << raw_size << "\n";
    assert (raw_size == 4 || raw_size == 2 || raw_size == 1);
    return 0;
    break;
  }
  return 0;
}



//----------------------------------------------------------------------
/// Encode one data value from a buffer to a file.
/// Will be called when repetitions==1 to remove useless loops.
/// It's important that this be defined as an inline function.
/// \param buf  The (allocated) buffer containing the raw data.
/// \return Number of raw words actually encoded.
//----------------------------------------------------------------------
inline size_t slim_channel_encode::encode_frame_singlevalue(void *buf) {

  uint32_t *dptr;
  uint16_t *sptr;
  uint8_t *cptr;

  switch (raw_size) {
  case 4:
    dptr = reinterpret_cast<uint32_t *>(buf);
    if (bit_rotation)
      dptr[0] = rotate(dptr[0]);
    enc->encode_scalar(dptr);
    return 4;
    break;

  case 2:
    sptr = reinterpret_cast<uint16_t *>(buf);
    if (bit_rotation)
      sptr[0] = rotate(sptr[0]);
    enc->encode_scalar(sptr);
    return 2;
    break;

  case 1:
    cptr = reinterpret_cast<uint8_t *>(buf);
    if (bit_rotation)
      cptr[0] = rotate(cptr[0]);
    enc->encode_scalar(cptr);
    return 1;
    break;

  default:
    cerr << "Oops: channel has raw_size=" << raw_size << "\n";
    assert (raw_size == 4 || raw_size == 2);
    return 0;
    break;
  }
  return 0;
}



//----------------------------------------------------------------------
/// Set the channel's encoder and take ownership.
//----------------------------------------------------------------------
void slim_channel_encode::set_encoder(encoder *e) {
  if (enc)
    delete enc;
  enc = e;
}



//----------------------------------------------------------------------
/// Set output bitstream by filename.
/// \param out_name The filename to use for bit output
/// \return True on success, false on error.
//----------------------------------------------------------------------
bool slim_channel_encode::set_output(char * const out_name) {
  ob = new obitstream(out_name);
  return enc->set_output(ob);
}



//----------------------------------------------------------------------
/// Set output bitstream to an existing obitstream object
/// \param ob_in  The existing obitstream.
/// \return True on success, false on error.
//----------------------------------------------------------------------
bool slim_channel_encode::set_output(obitstream *ob_in)  {
  ob = ob_in;
  return enc->set_output(ob);
}



//----------------------------------------------------------------------
/// Does the encoder expect a worthless compression?
//----------------------------------------------------------------------
bool slim_channel_encode::expect_zero_compression() const {
  if (ndata_sampled < MIN_SAMPLES)
    return true;
  return enc->expect_zero_compression();
}



//----------------------------------------------------------------------
/// Temporarily replace the usual encoder for this section only.
/// \return  Pointer to the new temporary encoder.
//----------------------------------------------------------------------
encoder *slim_channel_encode::replace_encoder() {
  if (enc == NULL)
    return NULL;

  usual_encoder = enc;
  usual_deltas = encode_deltas;
  enc =  usual_encoder->replacement_encoder();
  enc->set_output(ob);
  encode_deltas = enc->uses_deltas();
  return enc;
}



//----------------------------------------------------------------------
/// Temporarily replace the usual encoder for this section only
/// with the strictly-constant encoder.
/// \param   d0   The value of the constant data.
/// \return  Pointer to the new temporary encoder.
//----------------------------------------------------------------------
encoder *slim_channel_encode::replace_constant(int d0) {
  if (enc == NULL)
    return NULL;

  bit_rotation = bit_unrotation = 0;
  usual_encoder = enc;
  usual_deltas = encode_deltas;
  enc = usual_encoder->constant_encoder(d0);
  enc->set_output(ob);
  encode_deltas = enc->uses_deltas();
  return enc;
}



//----------------------------------------------------------------------
/// Return to using the usual encoder.
/// Undoes the work of replace_encoder().
/// \return  Pointer to the usual encoder.
//----------------------------------------------------------------------
encoder * slim_channel_encode::restore_encoder() {
  if (usual_encoder == NULL) 
    return NULL;

  delete (enc);
  enc = usual_encoder;
  encode_deltas = usual_deltas;
  usual_encoder = NULL;
  return enc;
}



//----------------------------------------------------------------------
/// \class slim_channel_decode
/// Derived class for decoding the data from a single channel.
/// A channel represents one stream of numbers from a single instrument,
/// and its properties are assumed to be consistent over time, but not
/// the same as those of other channels.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// Constructor
//----------------------------------------------------------------------
slim_channel_decode::slim_channel_decode(int reps, int size, bool deltas) :
  slim_channel(reps, size, deltas)
{
  dec = NULL;
}

//----------------------------------------------------------------------
/// Destructor
//----------------------------------------------------------------------
slim_channel_decode::~slim_channel_decode() {
  delete dec;
}

//----------------------------------------------------------------------
/// Set the channel's decoder and take ownership.
//----------------------------------------------------------------------
void slim_channel_decode::set_decoder(decoder *d) {
  if (dec)
    delete dec;
  dec = d;
}

//----------------------------------------------------------------------
/// Set the input ibitstream by name.
/// \param in_name The filename to use for a new ibitstream.
/// \return True on success, false on error.
//----------------------------------------------------------------------
bool slim_channel_decode::set_input(char * const in_name) {
  ib = new ibitstream(in_name);
  return dec->set_input(ib);
}

//----------------------------------------------------------------------
/// Set the input ibitstream to existing object.
/// \param ib_in The ibitstream to use (but not to own).
/// \return True on success, false on error.
//----------------------------------------------------------------------
bool slim_channel_decode::set_input(ibitstream *ib_in)  {
  ib = ib_in;
  return dec->set_input(ib);
}



//----------------------------------------------------------------------
/// Rotate data word by bit_rotation bits.
/// \param u The data to rotate.
/// \return  The rotated value.
//----------------------------------------------------------------------
inline uint32_t slim_channel_decode::rotate(uint32_t u) const {
  return (u>>bit_unrotation) ^ (u<<bit_rotation);
}



//----------------------------------------------------------------------
/// Rotate data word by bit_rotation bits.
/// \param u The data to rotate.
/// \return  The rotated value.
//----------------------------------------------------------------------
inline uint32_t slim_channel_decode::rotate(int u) const {
  return (u>>bit_unrotation) ^ (u<<bit_rotation);
}



//----------------------------------------------------------------------
/// Have the decoder object read in its parameters.
//----------------------------------------------------------------------
int slim_channel_decode::read_params(int bit_rotation_in) {
  bit_rotation = bit_rotation_in;
  bit_unrotation = 8*raw_size - bit_rotation;

  if (!dec)
    return -1;
  dec->read_params();
  
  return 0;
}

//----------------------------------------------------------------------
/// Decode one frame of data from a file and store to a buffer.
/// \param buf  The (allocated) buffer for storing the data.
/// \param size The maximum number of bytes to put in the buffer 
///             (may use less).
/// \return Size (bytes) of the raw data actually decoded.
//----------------------------------------------------------------------
size_t slim_channel_decode::decode_frame(void *buf, size_t size) {

  /// For chans with only one rep per frame, faster to handle specially.
  if (repetitions == 1) {
    return decode_frame_singlevalue(buf);
  }else if (size <= raw_size) {
    decode_frame_singlevalue(buf);
    return size;
  } 

  uint32_t *dptr;
  uint16_t *sptr;
  uint8_t *cptr;

  unsigned int nwords = size/raw_size;
  if (nwords > repetitions)
    nwords = repetitions;

  switch (raw_size) {
  case 4:
    dptr = reinterpret_cast<uint32_t *>(buf);
    
    dec->decode_vector(dptr, nwords);
    if (bit_rotation) {
      for (unsigned int i=0; i<nwords; i++)
	dptr[i] = rotate(dptr[i]);
    }
    return nwords*raw_size;
    break;

  case 2:
    sptr = reinterpret_cast<uint16_t *>(buf);
    dec->decode_vector(sptr, nwords);
    if (bit_rotation) {
      for (unsigned int i=0; i<nwords; i++)
	sptr[i] = rotate(sptr[i]);
    }
    return nwords*raw_size;
    break;

  case 1:
    cptr = reinterpret_cast<uint8_t *>(buf);
    dec->decode_vector(cptr, nwords);
    if (bit_rotation) {
      for (unsigned int i=0; i<nwords; i++)
	cptr[i] = rotate(cptr[i]);
    }
    return nwords*raw_size;
    break;

  default:
    cerr << "Oops: channel has raw_size=" << raw_size << "\n";
    assert (raw_size == 4 || raw_size == 2 || raw_size == 1);
  }
  return 0;
}



//----------------------------------------------------------------------
/// Decode one sample of data from a file and store to a buffer.
/// It is appropriate to call this when repetitions = 1.
/// \param buf  The (allocated) buffer for storing the data.
/// \return Size (bytes) of the raw data actually decoded.
//----------------------------------------------------------------------
size_t slim_channel_decode::decode_frame_singlevalue(void *buf) {

  uint8_t  * cptr = reinterpret_cast<uint8_t *>(buf);
  uint16_t * sptr = reinterpret_cast<uint16_t *>(buf);
  uint32_t * dptr = reinterpret_cast<uint32_t *>(buf);

  switch (raw_size) {
  case 4:
    dec->decode_scalar(dptr);
    if (bit_rotation) {
      dptr[0] = rotate(dptr[0]);
    }
    return raw_size;
    break;

  case 2:
    dec->decode_scalar(sptr);
    if (bit_rotation) {
      sptr[0] = rotate(sptr[0]);
    }
    return raw_size;
    break;

  case 1:
    dec->decode_scalar(cptr);
    if (bit_rotation) {
      cptr[0] = rotate(cptr[0]);
    }
    return raw_size;
    break;

  default:
    cerr << "Oops: channel has raw_size=" << raw_size << "\n";
    assert (raw_size == 4 || raw_size == 2 || raw_size == 1);
  }
  return 0;
}



//----------------------------------------------------------------------
/// Dump channel information by passing request to the decoder.
/// \param fout  The file stream to which to write the dump.
//----------------------------------------------------------------------
void slim_channel_decode::dump_info(ostream &fout) const {
  if (dec)
    dec->dump_info(fout);
}



