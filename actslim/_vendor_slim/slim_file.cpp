/// \file slim_file.cpp
/// Implements the classes slim_compressor_t and slim_expander_t.
/// They are used respectively for compressing raw data (from memory to file)
/// and restoring it from compressed (slim) form (file to memory).

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
#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <utime.h>

#include "slim.h"
#include "bitstream.h"

/// Markers at the end of sections to tell whether more follow.
enum section_foot_markers_t {
  NOT_LAST_SECTION = 0x8,   ///< Marker for all non-final sections in file.
  LAST_SECTION = 0xf,       ///< Marker for the last section in a file.
  BITS_SECTION_FOOT = 4,    ///< Number of bits to write for sect. foot markers.
};

/// Fail an assertion if this is NOT a TWOS-COMPLEMENT machine.
/// The program assumes throughout that it's on a twos-complement
/// machine.  Failing these assertions means that this is NOT.  It
/// might either be a ones-complement machine (recogize when i and ~i
/// are negatives of each other) or a sign-bit machine (recognize when
/// i and -i differ only in the highest bit).
static inline void verify_twos_complement() {

  // This first test is what we really need: UINT_MAX (0xffffffff) 
  // converts to the signed integer -1.
  assert (int(UINT_MAX) == -1);
  assert (INT_MAX + INT_MIN == -1);  
  assert(int(0u-1u)== -1);  

  int i=7;
  assert (i + (~i) == -1);  // i+(~i)==0 on a ones-complement machine.
}



/// Alter the access and modification time on a closed file.
/// \param filename  The file to alter.
/// \param mtime     The new value for filename's atime and mtime.
/// \return          Returns 0 on success or -1 on error (and errno is set).
inline int alter_mtime(const char *filename, time_t mtime) {
  struct utimbuf timbuf;
  timbuf.actime = timbuf.modtime = mtime;
  return utime(filename, &timbuf);
}



// ----------------------------------------------------------------------
/// \class slim_compressor_t
/// A file compression object, with channel encoders for all channels.
/// Capable of compressing from memory.  Use under the fatfile class in order
/// to compress directly from a separate file.
// ----------------------------------------------------------------------

//----------------------------------------------------------------------
/// Constructor
/// \param out_name The output (compressed) data file path.
/// \param flags_in Standard file flags.
/// \param deltas   Should all channels encode deltas?
/// \param samplepct_in What percent to use in sampling data.
//----------------------------------------------------------------------
slim_compressor_t::slim_compressor_t(const char *out_name, 
				     char flags_in, bool deltas,	
				     int samplepct_in) :
  flags(flags_in), 
  sample_pct(samplepct_in),
  encode_deltas(deltas)
{
  
  verify_twos_complement();
  
  frame_size = 0;
  raw_size = 0;
  mtime=0;
  sections_written = 0;
  sec_bytes_stored = 0;
  total_bytes_compressed = 0;
  num_frames = 1;

  ob = new obitstream(out_name);

  // Save output file name
  int len = strlen(out_name);
  out_filename = new char [len+1];
  strncpy(out_filename, out_name, len);
  out_filename[len] = '\0';

  section = new raw_section(SECTION_COMPRESS_MODE);
  curptr = section->ptr(0,0);
  max_frames_per_section = INT_MAX;
}



//----------------------------------------------------------------------
/// Destructor
//----------------------------------------------------------------------
slim_compressor_t::~slim_compressor_t() {

  if (ob->is_open()) {
    close_output();
  }

  delete section;
  delete ob;

  delete [] out_filename;
}



//----------------------------------------------------------------------
/// Close output (compressed) file
//----------------------------------------------------------------------
void slim_compressor_t::close_output() {
  write_last_section_foot();
  ob->close();
}


//----------------------------------------------------------------------
/// Take ownership of an existing encoding channel.
/// \param c A functioning slim_channel_encode object to use.
//----------------------------------------------------------------------
slim_channel_encode * slim_compressor_t::add_channel(slim_channel_encode *c) {
  channels.push(c, frame_size);
  c->set_output(ob);

  // Add this channel's frame_size to total.
  frame_size += c->get_frame_size();

  // Inform the raw section buffer about it.
  section->add_channel(c->get_repetitions(), c->get_raw_size());
  return c;
}



//----------------------------------------------------------------------
/// Create a new encoding channel and keep ownership.
/// \param reps Number of repetitions.
/// \param code ID # of the encoding method to use.
/// \param data_type ID # for the data type.
/// \param deltas Whether to encode successive differences.
/// \param rotate Whether channel can bit-rotate data to move up low bits.
//----------------------------------------------------------------------
slim_channel_encode * slim_compressor_t::add_channel(int reps,
					   enum code_t code, 
					   enum data_t data_type,
					   bool deltas,
					   bool rotate) {
  size_t size=slim_type_size[data_type];

  slim_channel_encode *c = 
    new slim_channel_encode(reps, size, encode_deltas, rotate);
  encoder *e = encoder_generator(code, data_type, deltas);
  c->set_encoder(e);
  add_channel(c);
  return c;
}



//----------------------------------------------------------------------
/// How many channels have been recorded?
/// \return Number of channels owned by this object.
//----------------------------------------------------------------------
inline int slim_compressor_t::num_channels() const {
  return int(channels.size());
}



//----------------------------------------------------------------------
/// Reset the list of channels.
//----------------------------------------------------------------------
void slim_compressor_t::reset_channels() {
  channels.clear();
  section->reset_channels();
  frame_size = 0;
}



//----------------------------------------------------------------------
/// Are all channels repeated only once per frame?
//----------------------------------------------------------------------
bool slim_compressor_t::no_reps() const {
  int nchan = num_channels();
  for (int c=0; c<nchan; c++)
    if (channels[c]->get_repetitions() > 1)
      return false;
  return true;
}



//----------------------------------------------------------------------
/// Confirm that the flags passed in agree with the object.
/// Assert that they are set only if they are correct.
/// Then set ONECHAN and NOREPS if needed.
//----------------------------------------------------------------------
void slim_compressor_t::confirm_flags() {
  if (flags & FLAG_ONECHAN)
    assert (num_channels() == 1);
  // if (num_channels() == 1)
  //   flags |= FLAG_ONECHAN;

  if (flags & FLAG_NOREPS)
    assert (no_reps());
  if (no_reps())
    flags |= FLAG_NOREPS;
}



//----------------------------------------------------------------------
/// Write the file header to the slim (output) file.
/// \param in_filename  The 
/// \return 0 on success, or error code.
//----------------------------------------------------------------------
int slim_compressor_t::write_file_header(const char *in_filename = NULL) {
  const bool with_nullchar = true;

  confirm_flags();
  ob->writestring(FILE_MAGIC);
  ob->writebits(mtime, 32);
  ob->writebits(flags, 8);
  if (flags & FLAG_SIZE)
    ob->writebits(raw_size, 32);
  if (flags & FLAG_NAME) {
    if (in_filename)
      ob->writestring(in_filename, with_nullchar);
    else
      ob->writestring("", with_nullchar);
  }
  if (flags & FLAG_XTRA)
    ob->writebits(0, 16); // XTRA not supported yet.
  assert(! (flags & FLAG_TOC));  // TOC not supported yet.
  
  return 0;
}



//----------------------------------------------------------------------
/// Write the section header to the slim (output) file.
/// \return 0 on success, or error code.
//----------------------------------------------------------------------
int slim_compressor_t::write_section_header() {
  ob->windup();  // Byte-align the section header.

  // Size in bytes (32 bits)
  if (section == NULL)
    throw "Cannot write_section_header for a NULL section.";
  size_t this_sect_size = section->get_size();
  ob->writebits(this_sect_size, BITS_SLIM_SECT_SIZE);

  // Table-of-contents forward pointer.  NOT IMPLEMENTED
  assert(! (flags & FLAG_TOC));  // TOC not supported yet.

  //  # of channels (24 bits)
  int nchan = num_channels();
  if (flags & FLAG_ONECHAN) {
    assert (nchan == 1);
  } else {
    if (nchan != int(nchan & lowestNset[BITS_SLIM_NUM_CHAN]))
      throw "Cannot write number of channels in allowed number of bits.";
    ob->writebits(nchan, BITS_SLIM_NUM_CHAN);
  }

  // A parameter set for each channel.
  bool no_repeats = (flags & FLAG_NOREPS);
  for (int c=0; c<nchan; c++) {
    if ((nchan>1) && !no_repeats) {
      int reps = channels[c]->get_repetitions();
      assert (reps >= 0);
      if (reps != int(reps & lowestNset[BITS_SLIM_REPETITIONS]))
	throw "Cannot write number of reps in allowed number of bits." ;
      ob->writebits(reps, BITS_SLIM_REPETITIONS);
    }
    reinterpret_cast<slim_channel_encode *>(channels[c])->write_params();
  }

  clear_channel_history();  // Reset the "prev" values used for deltas.
  return 0;
}



//----------------------------------------------------------------------
/// Compute the number of data in the input for a given channel.
/// It's okay if this is just a lower limit, as we use this only
/// in saying how many values are available for the sampling phase.
/// \param chan_num Array index of the channel in question.
/// \param section_size Number of bytes available.
/// \return Estimated number of data in the file.
//----------------------------------------------------------------------
int slim_compressor_t::num_data(int chan_num, int section_size) const {
  if (chan_num < 0 || chan_num >= num_channels())
    return 0;

  int frames_available = section_size / frame_size;
  if (frames_available > 0)
    return frames_available * channels[chan_num]->get_repetitions();

  // No complete frames exist.  For now, we know how to solve only if
  // all channels have the same size.
  size_t chan0_size = channels[0]->get_raw_size();
  if(num_channels()>1) {
    for (int i=1; i<num_channels(); i++) {
      assert(chan0_size == channels[i]->get_raw_size());
    }
  }
  return section_size / (num_channels() * chan0_size);
}



//----------------------------------------------------------------------
/// Set the number of frames that can appear in a section.
/// Note that in practice, the number might be reduced to
/// keep the total size of the section below some limit.
/// \param nf   The maximum number of frames.
//----------------------------------------------------------------------
void slim_compressor_t::set_section_frames(unsigned int nf) {

  // Is the requested number of frames too big for a section?
  if (nf*frame_size > MAX_SECTION_LENGTH ) {

    // Is even one frame too big for a section?
    if (frame_size > MAX_SECTION_LENGTH)
      throw "Frame is too long to fit in a single section.";

    if (frame_size > 0)
      nf = MAX_SECTION_LENGTH / frame_size;
    else
      nf = 1;
  }

  max_frames_per_section = nf;
  size_t section_size = nf * frame_size;
  section->resize(section_size);
  curptr = section->ptr(0,0);
}



//----------------------------------------------------------------------
/// Compute where in input a given data element lives for a given channel.
/// \param i_data   Data element number.
/// \param chan_num Array index of the channel in question.
/// \return Location of data in the file (for fseek to use).
//----------------------------------------------------------------------
long slim_compressor_t::data_offset(int i_data, int chan_num) {
  if (i_data < 0)
    return -1;
  if (chan_num < 0 || chan_num >= num_channels())
    return -1;

  slim_channel *c = channels[chan_num];
  long framenum = i_data / c->get_repetitions();
  long offset_in_frame = channels.offset(chan_num) + 
    (i_data % c->get_repetitions())*c->get_raw_size();
  return framenum*frame_size + offset_in_frame;
}



//----------------------------------------------------------------------
/// Compute encoding parameters of each channel from a sample of all data.
/// There are nchan arrays of sampled data, and each channel gets its own.
/// The slim_channel_encode objects decide how to use the data for setting
/// their own parameters.  This just makes the data available.
/// \param length Number of bytes available.
/// \return 0 on success, or error code.
//----------------------------------------------------------------------
int slim_compressor_t::compute_section_params(size_t length) {
  int chan_num;
  int nchan = num_channels();

  // Allocate arrays for sampled data.
  const int Target_group_size = 20;  // Do this many consecutive data points.
  const int Max_number_groups = 1000; // Don't allow more than this many groups.
  const int samplegroups = (sample_pct*Max_number_groups)/100;
  const int Target_data= samplegroups*Target_group_size;

  uint32_t *data = new uint32_t[Target_data+1];
  uint16_t *sdata = (uint16_t *)data; // They can share a buffer
  uint8_t  *cdata = (uint8_t *)data; // They can share a buffer
  assert (data != NULL);

  assert (length > 0);

  // Loop over all channels
  for (chan_num=0; chan_num<nchan; chan_num++) {
    slim_channel_encode *chan =
      reinterpret_cast<slim_channel_encode *>(channels[chan_num]);
    chan->restore_encoder(); // in case were using a temporary one.
    bool delta = chan->get_deltas();

    int raw_size = channels[chan_num]->get_raw_size();
    if (raw_size != 4 && raw_size != 2 && raw_size != 1) {
      cerr << "Channel " << chan_num << " has bad size "
             << raw_size << "\n";
      assert(raw_size == 4 || raw_size == 2 || raw_size == 1);
    }

    int available_data = num_data(chan_num, length);
    int ngroups;
    int skip_between_groups;
    int group_size = Target_group_size;

    // First, check that the channel isn't strictly constant.
    if (available_data >= 1) {
      int32_t d0 = 0;
      bool strictly_const = true;
      if (raw_size == 4) {
	d0 = section->ival(chan_num, 0);
	for (int i=1; i<available_data; i++) {
	  if (section->ival(chan_num, i) != d0) {
	    strictly_const = false;
	    break;
	  }
	}
      } else if (raw_size == 2) {
	d0 = section->sval(chan_num, 0);
	for (int i=1; i<available_data; i++) {
	  if (section->sval(chan_num, i) != d0) {
	    strictly_const = false;
	    break;
	  }
	}
      } else if (raw_size == 1) {
	d0 = section->cval(chan_num, 0);
	for (int i=1; i<available_data; i++) {
	  if (section->cval(chan_num, i) != d0) {
	    strictly_const = false;
	    break;
	  }
	}
      }

      if (strictly_const) {
	chan->replace_constant(d0);
	continue;
      }
    }


    // How much to sample?
    if (available_data > Target_group_size*Max_number_groups) {
      // There are at least Max_number_groups, so use 10%, 20%
      // ... 90% of the Max_number_groups, as requested..
      ngroups = samplegroups;
      group_size = Target_group_size;
      if (ngroups>1)
        skip_between_groups = 
          (available_data-1 - ngroups*group_size)/(ngroups-1) -1;
      else
        skip_between_groups = 0;

    } else if (available_data < Target_group_size) {
      // One full or partial group: use it all.
      ngroups = 1;
      group_size = available_data-1;
      skip_between_groups = 0; //irrelevant

    } else {
      // There is more than 1 full group, less than a full Max_number_groups.
      // So use 10%, 20%, ... 90% of the available data.
      int available_groups = 
        (available_data + Target_group_size/2)/Target_group_size;
      // Treat sample_pct as a percentage of all groups.
      ngroups = (available_groups*sample_pct)/100;
      group_size = Target_group_size;
      if (ngroups>1)
        skip_between_groups = 
          (available_data-1 - ngroups*group_size)/(ngroups-1) -1;
      else
        skip_between_groups = 0;
    }

    int i_in=0, i_out=0;
    try {
      if (raw_size == 4) {
        for (int ig=0; ig<ngroups; ig++) {
          int first = ig*group_size;
          int last = first + group_size;
	  for (i_out = first; i_out<=last; i_out++, i_in++) {
	    data[i_out] = section->uval(chan_num, i_in);
	  }

	  if (delta)
	    for (int j=first; j<last; j++) {
	      data[j] = data[j+1] - data[j];
	    }
          i_in += skip_between_groups;
        }

      } else if (raw_size == 2) {
        for (int ig=0; ig<ngroups; ig++) {
          int first = ig*group_size;
          int last = first + group_size;
	  for (i_out = first; i_out<=last; i_out++, i_in++) {
	    sdata[i_out] = section->sval(chan_num, i_in);
	  }
	  if (delta)
	    for (int j=first; j<last; j++) {
	      sdata[j] = sdata[j+1] - sdata[j];
            }
          i_in += skip_between_groups;
	}
	
      } else if (raw_size == 1) {
        for (int ig=0; ig<ngroups; ig++) {
          int first = ig*group_size;
          int last = first + group_size;
	  for (i_out = first; i_out<=last; i_out++, i_in++) {
	    cdata[i_out] = section->cval(chan_num, i_in);
	  }
	  if (delta)
	    for (int j=first; j<last; j++) {
	      cdata[j] = cdata[j+1] - cdata[j];
            }
          i_in += skip_between_groups;
	}
	
      }
    } catch (const char *s) {}

    if (raw_size == 4) {
      chan->compute_params(data, i_out - (i_out%group_size));
    } else if (raw_size == 2) {
      chan->compute_params(sdata, i_out - (i_out%group_size));
    } else if (raw_size == 1) {
      chan->compute_params(cdata, i_out - (i_out%group_size));
    }

    // If this channel isn't going to get compressed, then replace
    // the compressor class with the trivial ("default") compressor.
    if (chan->expect_zero_compression()) {
      chan->replace_encoder();
      if (raw_size == 4) {
	chan->compute_params(data, i_out - (i_out%group_size));
      } else if (raw_size == 2) {
	chan->compute_params(sdata, i_out - (i_out%group_size));
      } else if (raw_size == 1) {
	chan->compute_params(cdata, i_out - (i_out%group_size));
      }
    }
  }

  delete [] data;
  return 0;
}



//----------------------------------------------------------------------
/// Clear the channel histories (used in computing data deltas).
//----------------------------------------------------------------------
void slim_compressor_t::clear_channel_history() {
  slim_channel_encode *chan;
  for (int i=0; i<num_channels(); i++) {
    chan = reinterpret_cast<slim_channel_encode *>(channels[i]);
    chan->reset_previous();
  }
}


//----------------------------------------------------------------------
/// Write data from a user's buffer, of arbitrary size.
/// This accepts buffers of any length.  It pads out a partially full
/// section (if any), then writes full sections, and then starts a
/// new partial section, as needed to consume the full user's buffer.
/// \param buf    Pointer to the user's buffer.
/// \param length Size of the user's buffer.
/// \return Uncompressed size (bytes) of the written section.  (Should
///  equal the length argument.
//----------------------------------------------------------------------
size_t slim_compressor_t::
write(const unsigned char *buf, size_t length) {

  if (length <= 0 || buf == NULL)
    return 0;

  size_t bytes_remaining = length;

  // If possible, add to existing section.
  const size_t nominal_sect_size = frame_size*num_frames;
  size_t bytes_unused_existing_sect = nominal_sect_size - sec_bytes_stored;
  if (sec_bytes_stored > 0 && bytes_unused_existing_sect > 0) {
    size_t bytes_to_add = bytes_remaining;
    if (bytes_to_add > bytes_unused_existing_sect)
      bytes_to_add = bytes_unused_existing_sect;

    memcpy(curptr, buf, bytes_to_add);
    buf += bytes_to_add;
    curptr += bytes_to_add;
    sec_bytes_stored += bytes_to_add;
    bytes_remaining -= bytes_to_add;
    bytes_unused_existing_sect -= bytes_to_add;
  }

  // If existing is now full, encode_write it.
  if (bytes_unused_existing_sect <= 0) {
    encode_write_section(frame_size * num_frames);
    sec_bytes_stored = 0;
  }

  // If no more to do, return.
  if (bytes_remaining <= 0)
    return length - bytes_remaining;

  // While at least 1 full section still to be written, write full sections
  while (bytes_remaining >= nominal_sect_size) {
    size_t written = write_onesection(buf, nominal_sect_size);
    buf += written;
    bytes_remaining -= written;
  }

  // If no more to do, return.
  if (bytes_remaining <= 0)
    return length - bytes_remaining;

  // Copy remaining bytes to start of a new section.
  memcpy(curptr, buf, bytes_remaining);
  curptr += bytes_remaining;
  sec_bytes_stored += bytes_remaining;
  bytes_remaining = 0;

  return length - bytes_remaining;
}




//----------------------------------------------------------------------
/// Write a complete data section from a user's buffer.
/// \param buf    Pointer to the user's buffer.
/// \param length Size of the user's buffer.
/// \return Uncompressed size (bytes) of the written section.  (Should
///  equal the length argument.
//----------------------------------------------------------------------
size_t slim_compressor_t::
write_onesection(const unsigned char *buf, size_t length) {

  // This call is not allowed if there's a partial section already
  // stored in the private buffer.
  if (curptr != section->ptr(0,0) && section->get_size() > 0)
    throw "Cannot write a full section from user buffer; a partial\n"
      "section is already stored in private buffer.";

  section->use_external_buffer(buf, length);
  size_t written = encode_write_section(length);
  section->use_internal_buffer();
  curptr = section->ptr(0,0);
  return written;
}




//----------------------------------------------------------------------
/// Compress a single section in memory onto disk.
/// If request exceeds section size, no more than the full section size
/// will be written.
/// If request is less than a section, then any partial word at the end
/// will be encoded with padding (ghost bytes).
/// \param length  How much of the raw section buffer (bytes) to write.
/// \return Uncompressed size (bytes) of the written section.
//----------------------------------------------------------------------
size_t slim_compressor_t::encode_write_section(size_t length) {
  uint32_t section_crc=0;

  // Write either the file header, or the previous section's foot.
  // The prev section (if any) is not the last of the file.  Tag it as such.
  if (sections_written == 0)
    write_file_header();
  else {
    ob->writebits(NOT_LAST_SECTION, BITS_SECTION_FOOT);
  }
  sections_written ++;

  curptr = section->ptr(0,0);  // Reset the current buffer ptr.

  size_t this_sect_size = section->get_size();
  if (length < this_sect_size) {
    this_sect_size = length;
    section->resize(length);
  }

  // Compute the encoding parameters and write the section header.
  compute_section_params(this_sect_size);
  write_section_header();

  // If we need the CRC, then it must be computed before chan->encode_frame
  // potentially corrupts the data by applying bit rotations.
  if (flags & FLAG_CRC)
    section_crc = section->crc(this_sect_size);

  // The whole number of frames in this section.
  int frames_this_sect = this_sect_size / frame_size;

  // Loop over all whole frames in this section, encoding and writing.
  unsigned char *buf = section->ptr(0,0);
  slim_channel_encode *chan; // Alias for current channel.
  chan = reinterpret_cast<slim_channel_encode *>(channels[0]);
  int num_chan = num_channels();
  size_t bytes_thiscall;

  for (int frame_num=0; frame_num < frames_this_sect; frame_num++) {
    for (int i=0; i<num_chan; i++) {
      bytes_thiscall = chan->encode_frame(buf);
      buf += bytes_thiscall;
      chan = reinterpret_cast<slim_channel_encode *>(chan->next_chan);
    }
  }
  size_t bytes_written = buf - section->ptr(0,0);

  // If there's a partial frame remaining, this loop writes it, using as
  // many channels as needed.  Here we handle a partial word, too.
  size_t size_request, size_written;
  size_t ghost_bytes = 0;
  for (int i=0; i<num_chan && bytes_written < this_sect_size; i++) {

    // Handle partial words by always rounding up the requested write
    // size to the next complete word.  This rounding will only make a
    // difference when writing the appropriate (final) channel in the file.
    // Clear the ghost memory area to have a unique compressed file (not
    // dependent on uninitialized memory).
    size_request = this_sect_size - bytes_written;
    if (size_request < chan->get_frame_size() &&
        size_request % chan->get_raw_size()) {
      ghost_bytes = chan->get_raw_size() - (size_request % chan->get_raw_size());
      memset(buf + size_request, 0, ghost_bytes);
      assert (ghost_bytes <= MAX_GHOST_BYTES);
    }
    size_written = chan->encode_partial_frame(buf, size_request + ghost_bytes);
    buf += size_written;
    bytes_written += size_written - ghost_bytes;
    chan = reinterpret_cast<slim_channel_encode *>(chan->next_chan);

  }
  
  // Compute and emit the CRC checksum, if requested.
  if (flags & FLAG_CRC) {
    ob->writeword(section_crc);
  }

  sec_bytes_stored = 0; // For next time
  total_bytes_compressed += bytes_written;
  return bytes_written;
}


//----------------------------------------------------------------------
/// Learn the size and modification time of a file to be compressed.
/// \param raw_file_name  Name of the file to be compressed.
//----------------------------------------------------------------------
void slim_compressor_t::get_input_file_stats(const char *raw_file_name) {
  // Save original file Mtime and size.
  struct stat status;
  stat(raw_file_name, &status);
  mtime = status.st_mtime;
  raw_size = status.st_size;
}



//----------------------------------------------------------------------
/// Compress an entire file.
/// \param raw_file_name  Name of the file to be compressed.
/// \return 0 on success, or error code.
//----------------------------------------------------------------------
int slim_compressor_t::compress_from_file(const char *raw_file_name) {

  size_t bytes_written = 0;

  // Find (and save) mtime and size of raw file.
  get_input_file_stats(raw_file_name);
  
  if (raw_size > 0) {

    // Open the input.  Fix the IO buffer size
    FILE *infp = fopen(raw_file_name,"rb");
    if (!infp)
      return -1;
    setvbuf(infp, NULL, _IOFBF, 8192*16);

    // Break up into equal-sized sections (except for last can be short).
    num_frames = max_frames_per_section;
    int breakup_factor = 1;
    if (num_frames * frame_size > MAX_SECTION_LENGTH) {
      breakup_factor = divide_round_up(max_frames_per_section, MAX_SECTION_LENGTH);
      num_frames = divide_round_up(num_frames, breakup_factor);
    }
    int num_sect = divide_round_up(raw_size, frame_size*max_frames_per_section);
    num_sect *= breakup_factor;
    int frames_this_sect;
    int frames_written = 0;

    // Loop over all sections in the output
    for (int i_sect=0; i_sect<num_sect; i_sect++) {

      // Make sure not to miss any frames on the last pass.
      if (i_sect == num_sect-1) {
        int frames_entire_file = divide_round_up(raw_size, frame_size);
        frames_this_sect = frames_entire_file - frames_written;
      } else if (i_sect%breakup_factor == breakup_factor-1) {
        frames_this_sect = max_frames_per_section - (breakup_factor-1)*num_frames;
      } else
        frames_this_sect = num_frames;
    
      size_t this_size = frames_this_sect * frame_size;
      if (this_size + frames_written * frame_size > raw_size) {
        this_size = raw_size - (frames_written * frame_size);
      }
      this_size = section->fill(infp, this_size);
      bytes_written += encode_write_section(this_size);
      frames_written += frames_this_sect;
    }

    fclose(infp);

  } else  // else means file raw_size was zero.
    write_file_header();

  // Handle possible dangling bytes and final section foot.
  close_output();

  // Mark the compressed file with the raw file's modification time.
  alter_mtime(out_filename, mtime);
  
  // Report on the compression ratio
  if (!quiet) {
    struct stat status;
    stat(out_filename, &status);
    size_t final_size = status.st_size;

    double ratio = double(raw_size) / final_size;
    double bitsperword = 32.0/ratio;
    double saved = 100.*(1.-double(final_size)/raw_size);
    cout.width(20);
    cout << raw_file_name << ": ";

    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout.width(8);
    cout.precision(3);
    cout << ratio << ":1, ";
    cout.width(6);
    cout.precision(3);
    cout << bitsperword << " bits/word, ";
    cout.width(6);
    cout.precision(2);
    cout << saved << "% saved.\n";
    cout.precision(6);
  }

  return 0;
}



// ----------------------------------------------------------------------
/// Encode the last end-of-section tag for the file.
/// This used to be complicated before we started using the "ghost
/// bytes" concept to handle partial words at the end of file.
// ----------------------------------------------------------------------
inline void slim_compressor_t::write_last_section_foot() {
  if (sec_bytes_stored) {
    encode_write_section(sec_bytes_stored);
    sec_bytes_stored = 0;
  }
    
  ob->writebits(LAST_SECTION, BITS_SECTION_FOOT);
}



// ----------------------------------------------------------------------
/// \class slim_expander_t
/// A compressed file object, with channel decoders for all channels.
/// Capable of expanding into memory.  Use under slimfile class in order
/// to expand into a separate file.
// ----------------------------------------------------------------------

//----------------------------------------------------------------------
/// Constructor
/// \param in_name   The input (slim) file name.  Will be used read-only.
//----------------------------------------------------------------------
slim_expander_t::slim_expander_t(const char *in_name) {
  verify_twos_complement();

  ib = NULL;
  section = NULL;

  flags = 0;
  mtime = 0;
  raw_size = 0;

  bytes_read = 0;
  sec_bytes_read = 0;
  current_section_size = 0;

  size_t len = strlen(in_name);
  in_filename = new char [len+1];
  strncpy(in_filename, in_name, len);
  in_filename[len] = '\0';
  
  channels.clear();
  num_frames = 0;
  eof_tag_found = false;
  used_read = used_r_onesection = false;
  ignore_crc = false;

  // Find the compressed size
  struct stat status;
  stat(in_name, &status);
  slim_size = status.st_size;

  ib = new ibitstream(in_name);
  if (ib->is_open())
    read_file_header();
}

//----------------------------------------------------------------------
/// Constructor
/// \param in_fd     The input (slim) descriptor.  Will be used read-only.
//----------------------------------------------------------------------
slim_expander_t::slim_expander_t(int in_fd) {
  verify_twos_complement();

  ib = NULL;
  section = NULL;

  flags = 0;
  mtime = 0;
  raw_size = 0;

  bytes_read = 0;
  sec_bytes_read = 0;
  current_section_size = 0;

  in_filename = NULL;

  channels.clear();
  num_frames = 0;
  eof_tag_found = false;
  used_read = used_r_onesection = false;
  ignore_crc = false;

  // Find the compressed size
  struct stat status;
  fstat(in_fd, &status);
  slim_size = status.st_size;

  ib = new ibitstream(in_fd);
  if (ib->is_open())
    read_file_header();
}


//----------------------------------------------------------------------
/// Constructor (actslim addition)
/// \param in_fp        An already-open stdio stream positioned at byte 0.
/// \param in_slim_size Size of the compressed stream (for reporting only).
/// Ownership of in_fp passes to the ibitstream, which fcloses it.
/// This lets callers decode straight from an in-memory buffer via
/// fmemopen(), with no temp file and no extra copy.
//----------------------------------------------------------------------
slim_expander_t::slim_expander_t(FILE *in_fp, size_t in_slim_size) {
  verify_twos_complement();

  ib = NULL;
  section = NULL;

  flags = 0;
  mtime = 0;
  raw_size = 0;

  bytes_read = 0;
  sec_bytes_read = 0;
  current_section_size = 0;

  in_filename = NULL;

  channels.clear();
  num_frames = 0;
  eof_tag_found = false;
  used_read = used_r_onesection = false;
  ignore_crc = false;

  slim_size = in_slim_size;

  ib = new ibitstream(in_fp);
  if (ib->is_open())
    read_file_header();
}


//----------------------------------------------------------------------
/// Destructor
//----------------------------------------------------------------------
slim_expander_t::~slim_expander_t() {
  delete ib;
  delete section;
  if (in_filename)
    delete [] in_filename;
}



//----------------------------------------------------------------------
/// Take ownership of an existing encoding channel.
/// \param c A functioning slim_channel_decode object to use.
/// \param bit_rotation  Number of bits by which raw data are rotated.
//----------------------------------------------------------------------
slim_channel_decode * slim_expander_t::add_channel(slim_channel_decode *c,
					    int bit_rotation) {
  channels.push(c, 0);
  c->set_input(ib);
  c->read_params(bit_rotation);
  return c;
}



//----------------------------------------------------------------------
/// Create a new decoding channel and keep ownership.
/// \param reps Number of repetitions.
/// \param code ID # of the encoding method to use.
/// \param data_type ID # for the data type.
/// \param deltas Whether to encode successive difference data.
/// \param bit_rotation  Number of bits by which raw data are rotated.
//----------------------------------------------------------------------
slim_channel_decode * 
slim_expander_t::add_channel(int reps, enum code_t code, enum data_t data_type, 
		      bool deltas, int bit_rotation) {

  size_t data_size = slim_type_size[data_type];

  slim_channel_decode *c = 
    new slim_channel_decode(reps, data_size, deltas);
  decoder *d = decoder_generator(code, data_type, deltas);
  c->set_decoder(d);
  add_channel(c, bit_rotation);
  return c;
}



//----------------------------------------------------------------------
/// How many channels have been recorded?
/// \return Number of channels owned by this object.
//----------------------------------------------------------------------
inline int slim_expander_t::num_channels() const {
  return int(channels.size());
}



//----------------------------------------------------------------------
/// Is low-level file open?
/// \return Number of channels owned by this object.
//----------------------------------------------------------------------
bool slim_expander_t::is_open() const {
  return ib->is_open();
}



//----------------------------------------------------------------------
/// Read the file header.
/// \return 0 or a negative error code.
//----------------------------------------------------------------------
int slim_expander_t::read_file_header() {
  char file_magic[3]="", 
    orig_filename[256]="";
  ib->readstring(file_magic, 2);
  if (strcmp(file_magic, FILE_MAGIC))
    throw "file is not a slim file.";

  mtime = ib->readbits(32);
  flags = ib->readbits(8);
  if (flags & FLAG_SIZE)
    raw_size = ib->readbits(32);
  if (flags & FLAG_NAME)
    ib->readstring(orig_filename, 256);

  if (flags & FLAG_XTRA) {
    unsigned short xtra_len;
    xtra_len = ib->readbits(16);
    for (int i=0; i<xtra_len; i++)
      ib->readbits(8);
  }
  assert(! (flags & FLAG_TOC));  // TOC not supported yet.

  return 0;
}



//----------------------------------------------------------------------
/// Read the section header.
/// \return 0 or a negative error code.
//----------------------------------------------------------------------
int slim_expander_t::read_section_header() {
  ib->windup(); // Section headers are byte-aligned.

  if (section == NULL)
    section = new raw_section(SECTION_EXPAND_MODE);
  section->reset_channels();
  channels.clear();

  // Size of the raw section
  current_section_size = ib->readbits(BITS_SLIM_SECT_SIZE);
  sec_bytes_read = 0; // we have read none of it yet


  // Table of contents fwd ptr in the compressed file to next section.
  // NOT IMPLEMENTED
  assert(! (flags & FLAG_TOC));  // TOC not supported yet.

  // Code for # of channels
  int nchan;
  if (flags & FLAG_ONECHAN)
    nchan = 1;
  else
    nchan = ib->readbits(BITS_SLIM_NUM_CHAN);

  // Read channel descriptions
  for (int c=0; c<nchan; c++) {
    code_t algo_code;
    data_t type_code;
    unsigned int repetitions; // Bug if >2^19 and not an unsigned int.
    int data_size;

    if (flags & FLAG_NOREPS || nchan <= 1) {
      repetitions = 1;
    } else {
      repetitions = ib->readbits(BITS_SLIM_REPETITIONS);
    }

    unsigned int delta_code = ib->readbits(1);
    bool deltas = (delta_code != 0);
    int bit_rotation = ib->readbits(BITS_SLIM_NBITS);

    algo_code = code_t(ib->readbits(BITS_SLIM_ALG_CODE));
    type_code = data_t(ib->readbits(BITS_SLIM_TYPE_CODE));
    data_size = slim_type_size[type_code];

    // For nchan == 1, repetitions are implicit.
    if (nchan <= 1) {
      repetitions = (current_section_size) / data_size;
      if (repetitions <= 0)
        repetitions = 1;
    }

    slim_channel_decode *cdec = 
      add_channel(repetitions, algo_code, type_code, deltas, bit_rotation);

    // Inform the raw section buffer about it.
    section->add_channel(cdec->get_repetitions(), cdec->get_raw_size());
  }
  
  // Compute # of frames.
  size_t framesize = section->get_framesize();
  num_frames = current_section_size / framesize;
  if (current_section_size % framesize)
    num_frames++;
  section->set_num_frames(num_frames);
  section->resize(current_section_size);
  
  return 0;
}



//----------------------------------------------------------------------
/// Read data from compressed file up to given size.
/// This call will cover multiple sections, if needed.
/// Note that you can't mix calling read() and read_onesection().
/// \param  buf Buffer to fill with expanded (raw) data.
///             It is valid for buf to be NULL, which causes data to be 
///             decoded and discarded.
/// \param  max Size to fill (in bytes).
/// \return Number of bytes read.  The number should always equal the
///  request, unless the file is completely used up.
//----------------------------------------------------------------------
size_t slim_expander_t::read(unsigned char *buf, const size_t max) {
  size_t bytes_decoded=0;
  size_t bytes_thiscall=0, request;

  // Ensure user doesn't mix read() and read_onesection() methods.
  if (used_r_onesection)
    throw "Cannot call slim_expander_t::read() after ::read_onesection.";
  used_read = true;

  // First, see if there are previously-read but unused bytes.
  size_t bytes_unused = current_section_size - sec_bytes_read;
  if (bytes_unused > 0) {
    // This time, get all the unused bytes, up to max
    bytes_thiscall = bytes_unused;
    if (max < bytes_unused)
      bytes_thiscall = max;

    // Copy the wanted bytes.
    if (buf) {
      memcpy(buf, curptr, bytes_thiscall);
      buf += bytes_thiscall;
    }
    sec_bytes_read += bytes_thiscall;

    // Is this enough to meet the request?
    // If so, update the curptr and return.
    if (max == bytes_thiscall) {
      curptr += bytes_thiscall;
      return bytes_thiscall;
    }

    bytes_decoded = bytes_thiscall;
  }

  // The previously-read bytes, if there were any, did not complete the request.

  while (bytes_decoded < max) {
    bytes_thiscall = load_decode_section();
    if (bytes_thiscall == 0)
      break;

    request = max - bytes_decoded;
    if (bytes_thiscall > request) {
      bytes_thiscall = request;
    }

    sec_bytes_read = bytes_thiscall;
    if (buf) {
      memcpy(buf, section->ptr(0,0), bytes_thiscall);
      buf += bytes_thiscall;
    }
    bytes_decoded += bytes_thiscall;
  }
  if (section)
    curptr = section->ptr(0,0) + bytes_thiscall;
  else
    curptr = NULL;
  return bytes_decoded;
}


//----------------------------------------------------------------------
/// Read data from compressed file, precisely one section.
/// This is a private version that will only read from the beginning until 
/// the end of a slim file section.
/// Use ::read() to span section boundaries.
/// \return Number of bytes read.
//----------------------------------------------------------------------
size_t slim_expander_t::load_decode_section() {

  // Did we find the end-of-file tag on the last call?
  if (eof_tag_found) {
    return 0;
  }

  // Recognize EOF by failure to read a section header.
  try {
    read_section_header();
  } catch (const char *s) {
    delete section;
    section = NULL;
    return 0;
  }
  
  // Now we have a section header read, section is prepared.
  // It is possible that some bytes were consumed on previous call.
  // Find out how much to read.
  size_t bytes_thiscall;
  int bytes_remaining=current_section_size;  // signed: can be - if ghost bytes.
  unsigned char *buf = section->ptr(0,0);
  slim_channel_decode *chan = 
    reinterpret_cast<slim_channel_decode *>(channels[0]);

  while (bytes_remaining > 0) {
    // Try to fill the entire request on this channel.
    bytes_thiscall = chan->decode_frame(buf, bytes_remaining);

    // Couldn't fill.  Update and move to next channel in the list.
    buf += bytes_thiscall;
    bytes_remaining -= bytes_thiscall;
    if (bytes_remaining < int(chan->get_raw_size()) && 
        bytes_thiscall < chan->get_frame_size()) {
      if (bytes_remaining > 0) {
        bytes_thiscall = chan->decode_frame(buf, chan->get_raw_size());
        assert (int(bytes_thiscall) > bytes_remaining);
        buf += bytes_remaining;
        bytes_remaining = 0;
      }
      break;
    }

    chan = reinterpret_cast<slim_channel_decode *>
      (chan->next_chan);
  }

  // Check CRC, if they are stored in the current file
  if ((flags & FLAG_CRC) && !ignore_crc) {
    unsigned long crc = section->crc();
    unsigned long expected_crc = ib->readbits(32);
    if (crc != expected_crc) {
      cerr << in_filename << ": CRC-32 error.  compute " << crc <<
        ", file says " << expected_crc << "\n";

      size_t actual_size = current_section_size-bytes_remaining;
      if (actual_size != current_section_size) {
        cerr << in_filename << ": CRC-32 error.  Section size " <<
          actual_size << " (expected " << current_section_size << ")\n";
      } else {
        cerr << in_filename << ": CRC-32 error.  Section size " <<
          actual_size << " (as expected)\n";
      }
      throw "CRC error.";
    }
  }

  // We have read every byte in the section, so read the End-Section block.
  // If we get an error thrown, then we assume it was the end of file.
  // (This will happen if the EOF tag lands at the end of a native word.)
  try {
    unsigned int eof_tag = ib->readbits(BITS_SECTION_FOOT);
    if (eof_tag == LAST_SECTION)
      eof_tag_found = true;
  } catch (const char *s) {
    eof_tag_found = true;
  }
  
  
  // If we read some "ghost bytes", remove them
  size_t bytes_thissect = buf - section->ptr(0,0);
  if (bytes_remaining < 0)
    bytes_thissect += bytes_remaining;
  if (bytes_thissect != current_section_size) {
    cerr <<  "The uncompressed section was "<<bytes_thissect
         <<", not the expected size of "<<current_section_size <<".\n",
    throw "The uncompressed section was not the expected size.";
  }
    
  return bytes_thissect;

}



//----------------------------------------------------------------------
/// Load exactly one section and let user have access to it as a const array.
/// Note that you can't mix calling read() and read_onesection().
/// \param bufptr Called-by-reference allows returning the array location.
/// \return Size of the section read, or 0 if none.
//----------------------------------------------------------------------
size_t slim_expander_t::read_onesection(const unsigned char **bufptr) {

  // Ensure user doesn't mix read() and read_onesection() methods.
  if (used_read)
    throw "Cannot call slim_expander_t::read_onesection() after ::read.";
  used_r_onesection = true;

  size_t sect_size = load_decode_section();
  if (section && sect_size > 0)
    *bufptr = section->ptr(0,0);
  else {
    *bufptr = NULL;
    sect_size = 0;
  }
  return sect_size;
}



//----------------------------------------------------------------------
/// Expand the entire file to another file.
/// \param raw_file_name  A file where we put the raw data.
/// \return 0 or a negative error code.
//----------------------------------------------------------------------
int slim_expander_t::expand_to_file(const char *raw_filename) {
  FILE *fp = fopen(raw_filename, "wb");
  assert (fp != NULL);

  if (!quiet) {
    cout.width(20);
    cout << in_filename << ":\t";
    if (raw_size) {
      cout.width(6);
      cout.precision(1);
      ios_base::fmtflags old = cout.setf(ios_base::fixed, ios_base::floatfield);
      cout << 100.*(1.-double(slim_size)/raw_size);
      cout.setf(old, ios::floatfield);

    } else
      cout << " ???? ";
    cout << "% -- replacing with " << raw_filename << "\n";
  }

  // Loop over all sections, reading and writing them
  while (1) {
    size_t bytes_thissect = load_decode_section();
    if (bytes_thissect == 0)
      break;
    section->flush(fp, bytes_thissect);
  }
  fclose(fp);
  alter_mtime(raw_filename, mtime);
  return 0;
}



//----------------------------------------------------------------------
/// Expand the entire file to standard output.
/// \return 0 or a negative error code.
//----------------------------------------------------------------------
int slim_expander_t::expand_to_stdout() {
  FILE *fp = stdout;
  assert (fp != NULL);

  // Loop over all sections, reading and writing them
  while (1) {
    size_t bytes_thissect = load_decode_section();
    if (bytes_thissect == 0)
      break;
    section->flush(fp, bytes_thissect);
  }
  return 0;
}


//----------------------------------------------------------------------
/// Print slim file infomation to stdout.
/// \return 0 or a negative error code.
//----------------------------------------------------------------------
int slim_expander_t::dump_sliminfo() {
  bool decode_error = false;

  cout << ("-----------------------------------------------------------"
	 "----------------\n");
  cout << "Slim file     "<< in_filename << "\n";
  cout << "Original time "<< ctime(&mtime);  // Note newline built into ctime
  cout << "Slim size     "<< slim_size << "\n";
  if (flags&FLAG_SIZE) {
    cout << "Raw size      "<< raw_size << "\n";
    cout << "Compression   "<< 100.*(1-slim_size*1.0/raw_size) << "%  "<< raw_size<<"\n";
  }

  if (flags&FLAG_NAME)
    cout << "Raw name      present\n";

  if (flags&FLAG_XTRA)
    cout << "XTRA header data present\n";

  if (flags&FLAG_TOC)
    cout << "Table of Contents present.  Yugh!\n";

  if (flags&FLAG_ONECHAN)
    cout << "File contains only 1 channel at a time.\n";
  else
    cout << "Multiple data channels are allowed per section.\n";

  if (flags&FLAG_NOREPS)
    cout << "File channels never repeat before giving way "
	   "to others in a frame.\n";
  else
    cout << "File channels may repeat in a frame.\n";

  if (flags&FLAG_CRC)
    cout << "CRC-32 checksums present.\n";
  else
    cout << "CRC-32 checksums not used.\n";
  
  cout << "-----------------------------------------------------------"
	 "----------------\n";

  int isect=0;
  size_t section_size=0;
  try {
    section_size = load_decode_section();
  } catch (const char *s) {
    decode_error = true;
  }
  while (section_size > 0 || decode_error) {
  
    cout << "SECTION "<<isect<<" (size: "<<section_size<<" bytes):\n";
    cout << "Number of frames:   " << num_frames << "\n";
    cout << "Number of channels: " << num_channels() << "\n";
    slim_channel_decode *chan = NULL; 
    for (int i=0; i<num_channels(); i++) {
      chan = reinterpret_cast<slim_channel_decode *>(channels[i]);
      cout << "Chan " << setw(4) << i << " ";
      cout << "Repeat " << setw(4) << chan->get_repetitions() <<" ";
      cout << "size " << chan->get_raw_size() << " ";
      if (chan->get_deltas())
        cout << "DELT";
      else
        cout << "    ";
      chan->dump_info();
      if (chan->get_bit_rotation())
        cout << " (rot=" << chan->get_bit_rotation() << ")";
      cout << endl;
    }
    cout << endl;
    if (decode_error)
      throw "Decoding error";

    try {
      section_size = load_decode_section();
    } catch (const char *s) {
      decode_error = true;
    }

    isect++;
  }

  cout << "-----------------------------------------------------------"
	 "----------------\n";
  return 0;
}

