/// \file slim_control.cpp
/// Implements to slim_control class to parse command-line arguments.

//  Copyright (C) 2008, 2009, 2013 Joseph Fowler
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
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>

#include "slim.h"
#define RAW_SUFFIX ".raw"   ///< Suffix for all raw files when preserving original.


/// Single-character options
const char short_opt[] = "m:c:r:F:dnpXxSC0bPokquivsyfg987654321V?B:";

/// The long options slim recognizes.
const static struct option long_opt[] = {
  {"method",            required_argument, NULL, 'm'},
  {"num-chan",          required_argument, NULL, 'c'},
  {"repeats",           required_argument, NULL, 'r'},
  {"frames",            required_argument, NULL, 'F'},
  {"deltas",            no_argument,       NULL, 'd'},
  {"filename",          no_argument,       NULL, 'n'},
  {"preserve",          no_argument,       NULL, 'p'},
  {"compress",          no_argument,       NULL, 'X'},
  {"expand",            no_argument,       NULL, 'x'},
  {"rawsize",           no_argument,       NULL, 'S'},
  {"compute-crc32",     no_argument,       NULL, 'C'},
  {"ignore-crc32",      no_argument,       NULL, '0'},
  {"permit-bitrotation",no_argument,       NULL, 'b'},
  {"practice",          no_argument,       NULL, 'P'},
  {"stdout",            no_argument,       NULL, 'o'},
  {"force",             no_argument,       NULL, 'k'},
  {"quiet",             no_argument,       NULL, 'q'},
  {"unsigned",          no_argument,       NULL, 'u'},
  {"int",               no_argument,       NULL, 'i'},
  {"ushort",            no_argument,       NULL, 'v'},
  {"short",             no_argument,       NULL, 's'},
  {"char",              no_argument,       NULL, 'y'},
  {"float",             no_argument,       NULL, 'f'},
  {"double",            no_argument,       NULL, 'g'},
  {"best",              no_argument,       NULL, '9'},
  {"18-pct",            no_argument,       NULL, '9'},
  {"16-pct",            no_argument,       NULL, '8'},
  {"14-pct",            no_argument,       NULL, '7'},
  {"12-pct",            no_argument,       NULL, '6'},
  {"10-pct",            no_argument,       NULL, '5'},
  {"8-pct",             no_argument,       NULL, '4'},
  {"6-pct",             no_argument,       NULL, '3'},
  {"4-pct",             no_argument,       NULL, '2'},
  {"2-pct",             no_argument,       NULL, '1'},
  {"fast",              no_argument,       NULL, '1'},
  {"version",           no_argument,       NULL, 'V'},
  {"help",              no_argument,       NULL, '?'},
  {"debug-buffer",      required_argument, NULL, 'B'},
  {0, 0, 0, 0}
};

// ---------------------------------------------------------------------------
/// \class slim_control
/// Handles command-line arguments and keeps track of control settings.
// ---------------------------------------------------------------------------
///


/// Constructor.
slim_control::slim_control() {
  set_defaults();
}



/// User called this with name "unslim"
void slim_control::unslim() {
  mode = SLIM_DECODE;
}



/// User called this with name "slimcat"
void slim_control::slimcat() {
  use_stdout = true;
  preserve_input = true;
}



/// Parse command-line arguments.
/// \param argc Number of tokens in command line.
/// \param argv Array of ptr to these tokens.
void slim_control::process_args(int argc, char * const argv[])
{
  bool print_version=false, print_help=false;
  const int END_OF_OPTIONS = -1;

  // cout << "Processing "<<argc<<" arguments"<<endl;
  for (opterr = 0;;) {
    int longindex=0;
    int curind = optind;
    int c = getopt_long(argc, argv, short_opt, long_opt, &longindex);
    // cout << "Working with optind="<<curind<<" and c="<<c<<endl;
    if (c == END_OF_OPTIONS)
      break;

    // cout <<"Processed option #"<<longindex<<" which is ascii("<<c<<")='"
    // <<char(c) <<"'"<<endl;
    switch (c) {
    case 0:
      if (strncmp(long_opt[longindex].name, "help",
		   strlen(long_opt[longindex].name))==0)
	usage();
      break;

    case 'm':
      code_method = (enum code_t)(atoi(optarg));
      break;

    case 'c':
      nchan = atoi(optarg);
      onechan = (nchan == 1);
      break;

    case 'F':
      nframes = atoi(optarg);
      break;

    case 'o':
      slimcat();
      break;

    case 'r':
      repeats = atoi(optarg);
      break;

    case 'B':
      debug_buf_size = atoi(optarg);
      break;

    case 'd':
      deltas = true;
      break;

    case 'n':
      save_filename = true;
      break;

    case 'p':
      preserve_input = true;
      break;

    case 'S':
      save_rawsize = true;
      break;

    case 'X':
      mode = SLIM_ENCODE;
      break;

    case 'x':
      //mode = SLIM_DECODE;
      unslim(); // equiv to mode = SLIM_DECODE;
      break;

    case 'V':
      print_version = true;
      break;

    case 'C':
      crc = true;
      break;

    case '0':
      ignore_crc = true;
      break;

    case 'P':
      practice=true;
      break;

    case 'k':
      force_clobber = true;
      break;

    case 'q':
      quiet = true;
      break;

    case 'u':
      data_type = SLIM_TYPE_U32;
      break;

    case 'b':
      permit_bitrotation= true;
      break;

    case 'i':
      data_type = SLIM_TYPE_I32;
      break;

    case 'v':
      data_type = SLIM_TYPE_U16;
      break;

    case 's':
      data_type = SLIM_TYPE_I16;
      break;

    case 'y':
      data_type = SLIM_TYPE_I8;
      break;

    case 'f':
      //data_type = SLIM_TYPE_FLOAT;
      data_type = SLIM_TYPE_I32;
      break;

    case 'g':
      data_type = SLIM_TYPE_DOUBLE;
      break;

    case '9': case '8': case '7':
    case '6': case '5': case '4':
    case '3': case '2': case '1':
      sample_pct = 2*(c-'0');
      assert (sample_pct >= 1 && sample_pct <= 100);
      break;

    case '?':
      print_help = true;
      if (optopt)
        cout << "Bad short opt '"<<optopt << "'"<< endl;
      else
        cout << "Bad long opt '"<<argv[curind] << "'"<< endl;
      break;

    default:
      throw "Error processing options with getopt_long:"
        " unrecognized return value";
      break;
    }
  }

  // Version or help messages are printed, then exit.
  if (print_version)
    version();
  if (print_help)
    usage();
  if (print_version || print_help)
    exit(0);

  if (practice) {
    mode = SLIM_ENCODE;
    preserve_input = true;
  }
}




/// Print a usage message.
void slim_control::usage() const {
  cout << "usage: slim    [options] filename ...\n";
  cout << "usage: unslim  [options] filename ...\n";
  cout << "usage: slimcat [options] filename ...\n";
  cout << "Allowed options are:\n";

  usage_printoptions();
}


/// Print the options portion of a usage message.
void slim_control::usage_printoptions() const {
  const char noarg[]="",
    hasarg[]="<arg required>",
    optional[]="[arg optional]";
  const char *argval[3]={noarg, hasarg, optional};

  ios_base::fmtflags old = cout.setf(ios_base::floatfield);
  cout.setf(ios_base::left);
  for (const struct option *opt = long_opt; opt->name; opt++) {
    cout << "-" << char(opt->val) << ", --";
    cout.width(20);
    cout << opt->name << " " << argval[opt->has_arg] << "\n";
  }
  cout << "Someone really ought to expand this.  JWF 23 July 2007.\n";
  cout.setf(old);
}



/// Print version message.
void slim_control::version() const
{
  string versionstr=SLIM_VERSION;
  if (versionstr[0] == ' ')
    cout << "This is slim, the physics data compressor, untagged version.\n";
  else {
    cout << "This is slim, the physics data compressor, version "
         << versionstr << ".\n";
  }
}



/// Set members to default values.
void slim_control::set_defaults()
{
  deltas = false;
  nchan = 1;
  nframes = 0;
  repeats = 0;
  debug_buf_size = 0;
  sample_pct = 10;
  mode = SLIM_MODE_UNKNOWN;
  code_method = SLIM_ENCODER_REDUCED_BINARY;
  data_type = SLIM_TYPE_I32;
  force_clobber = false;
  preserve_input = false;
  permit_bitrotation = false;
  use_stdout = false;

  save_rawsize = true;
  save_filename = false;
  have_xtra = false;
  have_toc = false;
  onechan = true;
  noreps = false;
  crc = false;
  ignore_crc = false;
  reserved0 = false;
  practice = false;
  quiet = false;
}



/// A flags character suitable for writing to bitstream.
/// Contains bits for the FLAG character in the slim file header.
/// \return Value of FLAG given all relevant booleans.
char slim_control::flags() const
{
  char r = 0;
  if (save_rawsize) r |= FLAG_SIZE;
  if (save_filename) r |= FLAG_NAME;
  if (have_xtra) r |= FLAG_XTRA;
  if (have_toc) r |= FLAG_TOC;
  if (onechan) r |= FLAG_ONECHAN;
  if (noreps) r |= FLAG_NOREPS;
  if (crc) r |= FLAG_CRC;
  return r;
}



/// Set all boolean flags using the FLAG character from the slim file header.
/// \param in  Value of FLAG, which sets all relevant booleans.
/// \return    the input character.
char slim_control::flags(char in)
{
  save_filename = in & FLAG_NAME;
  save_rawsize =  in & FLAG_SIZE;
  have_xtra =     in & FLAG_XTRA;
  have_toc =      in & FLAG_TOC;
  onechan =       in & FLAG_ONECHAN;
  noreps =        in & FLAG_NOREPS;
  crc =           in & FLAG_CRC;

  return in;
}



/// Decide whether to encode/decode based on file suffix.
/// \param fname  The name of the file to be encoded or decoded.
/// \return SLIM_ENCODE or SLIM_DECODE.
enum slim_mode_t slim_control::detect_file_mode(const char *fname) const
{
  if (strstr(fname,".slm") ||
      strstr(fname,".SLM"))
    return SLIM_DECODE;

  return SLIM_ENCODE;
}



/// Handle one file, by encoding or decoding it.
/// \param fname  The name of the file to be encoded or decoded.
void slim_control::handle_one_file(const char *fname)
{
  enum slim_mode_t thismode = mode;

  if (thismode == SLIM_MODE_UNKNOWN)
    thismode = detect_file_mode(fname);
  assert(thismode == SLIM_ENCODE || thismode == SLIM_DECODE);

  // Make sure the file exists and doesn't have multiple hard links
  struct stat st;
  int ret = stat(fname, &st);
  if (ret)
    throw bad_file(fname, ": does not exist.");
  if (S_ISDIR(st.st_mode)) {
    throw bad_file(fname, " is a directory -- ignored.");
  }
  if (!S_ISREG(st.st_mode)) {
    throw bad_file(fname, ": not a regular file.");
  }
  if (st.st_nlink > 1 && !force_clobber && !preserve_input) {
    throw bad_file(fname, ": has more than one hard link (use -k to force).");
  }

  // Make sure that we have read permission on the input file
  if (access(fname, R_OK)) {
    if (errno == EACCES)
      throw bad_file(fname, ": read permission denied.");
    else
      throw bad_file(fname, ": access() call failed; not with EACCES error.");
  }

  if (thismode == SLIM_ENCODE) {
    if (debug_buf_size == 0)
      compress_one_file(fname);
    else {
      debug_compress_from_memory(fname);
    }
  }

  if (thismode == SLIM_DECODE) {
    if (debug_buf_size == 0)
      expand_one_file(fname);
    else
      debug_expand_from_memory(fname);
  }
}



/// Compress one file using slim_compressor_t::compress_from_file()
/// \param rawname  The name of the file to be encoded.
void slim_control::compress_one_file(const char *rawname) {

  // Verify that file exists
  struct stat st;
  int ret = stat(rawname, &st);
  if (ret)
    throw bad_file(rawname, ": does not exist.");
  const size_t raw_file_size = st.st_size;

  // Create the output slim-file name, create compressor object.
  ostringstream outstr;
  outstr << rawname << '.' << SLIM_SUFFIX;
  string outname_string = outstr.str();
  const char *outname = outname_string.c_str();

  if (!force_clobber) {
    ret = stat(outname, &st);
    if (ret == 0)
      throw bad_file(outname, ": slim file exists (use -k to force).");
  }

  slim_compressor_t *compressor =
    new slim_compressor_t(outname, flags(), deltas, sample_pct);

  // If no cmd-line indications, just let there be 1 frame,
  // and all data go in multiple repeats within that section and frame.
  int repeats_thisfile = repeats, nframes_thisfile = nframes;
  if (repeats <= 0 && nframes <= 0) {
    repeats_thisfile = raw_file_size/(sizeof(int)*nchan);
    nframes_thisfile = 1;
  } else if (nframes <= 0) {
    if (repeats == 1)
      nframes_thisfile =  divide_round_up(raw_file_size, sizeof(int)*nchan);
    else
      nframes_thisfile = 1;
  } else if (repeats <= 0) {
    repeats_thisfile = raw_file_size/(sizeof(int)*nchan*nframes);
  }

  // Treat 1-channel sections as a single frame, for efficient compress time.
  // Then break up that frame only if it exceeds the MAX_SECTION_LENGTH
  if (nchan == 1) {
    if (nframes_thisfile > 1)
      repeats_thisfile *= nframes_thisfile;
    nframes_thisfile = divide_round_up(repeats_thisfile,
                                       MAX_SECTION_LENGTH / sizeof(int));
    if (nframes_thisfile > 1)
      repeats_thisfile = MAX_SECTION_LENGTH / sizeof(int);
  }
  noreps = (repeats_thisfile == 1);

  // Build the channel list
  for (int i=0; i<nchan; i++)
    compressor->add_channel(repeats_thisfile, code_method, data_type,
		    deltas, permit_bitrotation);

  compressor->set_section_frames(nframes_thisfile);
  compressor->set_quiet(quiet);
  compressor->compress_from_file(rawname);
  delete compressor;

  // Clean up unwanted files
  struct stat st_raw, st_out;
  int ret1 = stat(rawname, &st_raw);
  int ret2 = stat(outname, &st_out);
  if (ret1 || ret2) return;
  if (!force_clobber && (st_raw.st_size <= st_out.st_size)) {
    if (!quiet)
      cerr << "slim: " << rawname <<
        " expanded when slimmed (use -k to force).\n";
    unlink(outname);
  } else if (practice) {
    unlink(outname);
  } else if (!preserve_input)
    unlink(rawname);

}



/// Expand one compressed file using slim_expander_t::expand_to_file().
/// \param compname  The name of the file to be decoded.
void slim_control::expand_one_file(const char *compname) {

  // Verify that file exists
  struct stat st;
  int ret = stat(compname, &st);
  if (ret)
    throw bad_file(compname, ": does not exist.");

  // Create appropriate raw file name.
  // Normally, this means stripping any trailing dot-suffix.
  // But if we have --preserve, then append a ".raw" suffix instead.
  const char *last_suffix = strrchr(compname,'.');
  size_t baselen;
#ifdef HAVE_LIBLZ4
  if (last_suffix &&
    (!strcmp(last_suffix, ".lz4") || !strcmp(last_suffix, ".LZ4"))) {
    // Strip off extra LZ4 suffix if found
    char *period = (char*) last_suffix;
    ((char *) last_suffix)[0] = '\0';
    last_suffix = strrchr(compname, '.');
    period[0] = '.';
  }
#endif
  if (last_suffix)
    baselen = (last_suffix-compname);
  else
    baselen=strlen(compname);

  const size_t RAWLEN=baselen+strlen(RAW_SUFFIX)+5;
  char *rawname = new char[RAWLEN];
  // Fill with end-string bytes so that the strncpy produce valid strings.
  for (size_t i=0; i<RAWLEN; i++)
    rawname[i] = '\0';
  strncpy(rawname, compname, baselen);
  if (last_suffix == NULL)
    strncpy(rawname+baselen, RAW_SUFFIX, strlen(RAW_SUFFIX));
  else if (preserve_input)
    strncpy(rawname+baselen, RAW_SUFFIX, strlen(RAW_SUFFIX));
  else
    rawname[baselen] = '\0';

  // Make sure not to clobber the output file unless desired.
  if (!force_clobber) {
    ret = stat(rawname, &st);
    if (ret == 0) {
      bad_file bf = bad_file(rawname, ": raw file exists (use -k to force).");
      delete [] rawname;
      throw bf;
    }
  }

  // Do the expansion, deleting the compressed file (unless --preserve).
  slim_expander_t *expander = new slim_expander_t(compname);
  if (ignore_crc)
    expander->set_ignore_crc();
  expander->set_quiet(quiet);

  if (use_stdout)
    expander->expand_to_stdout();
  else
    expander->expand_to_file(rawname);

  delete expander;
  delete [] rawname;
  if (!preserve_input)
    unlink(compname);
}



//////////////////////////////////////////////////////////////////////////////
// Debugging code
// These methods do compression/expansion of a a full file, but they do it
// using the library slim_compressor_t::write() and slim_expander_t::read()
// methods, so we can test them.
//////////////////////////////////////////////////////////////////////////////


/// Compress a single file using slim_compressor_t::write() to test that code.
/// \param rawname  The path of the file to be compressed.
void slim_control::debug_compress_from_memory(const char *rawname) {

  // Verify that file exists
  struct stat st;
  stat(rawname, &st);
  if (!S_ISREG(st.st_mode)) {
    cerr << "slim: " << rawname << ": No such file\n";
    return;
  }
  const size_t raw_file_size = st.st_size;

  // Create the output slim-file name, create compressor object.
  ostringstream outstr;
  outstr << rawname << '.' << SLIM_SUFFIX;
  string outname = outstr.str();

  slim_compressor_t *compressor =
    new slim_compressor_t(outname.c_str(), flags(), deltas, sample_pct);
  compressor->mtime = st.st_mtime;
  compressor->raw_size = st.st_size;

  // If no cmd-line indications, just let there be 1 frame, 1 section,
  // and all data go in multiple repeats within that section and frame.
  int repeats_thisfile = repeats;
  if (repeats <= 0 && nframes <= 0) {
    repeats_thisfile = raw_file_size/sizeof(int);
    nframes = 1;
  } else if (nframes <= 0) {
    if (repeats == 1)
      nframes = INT_MAX;
    else
      nframes = 1;
  } else if (repeats <= 0) {
    repeats_thisfile = 1;
  }

  // Build the channel list
  for (int i=0; i<nchan; i++)
    compressor->add_channel(repeats_thisfile, code_method, data_type,
		    deltas, permit_bitrotation);

  compressor->set_section_frames(nframes);

  FILE *fp = fopen(rawname, "rb");
  if (fp == 0)
    return;

  unsigned char *buf = new unsigned char[debug_buf_size];
  size_t raw_bytes_read=0;
  while (1) {
    raw_bytes_read = fread(buf, 1, debug_buf_size, fp);
    if (raw_bytes_read <= 0)
      break;
    compressor->write(buf, raw_bytes_read);
  }
  fclose(fp);

  delete [] buf;
  delete compressor;

  if (practice)
    unlink(outname.c_str());

  if (!preserve_input)
    unlink(rawname);
}



/// Expand a single file  using slim_expander_t::read() to test that code.
/// \param rawname  The path of the file to be expanded.
void slim_control::debug_expand_from_memory(const char *compname) {

  // Verify that file exists
  struct stat st;
  stat(compname, &st);
  if (!S_ISREG(st.st_mode)) {
    cerr << "slim: " << compname << ": No such file\n";
    return;
  }

  // Create appropriate raw file name.
  // Normally, this means stripping any trailing dot-suffix.
  // But if we have --preserve, then append a ".raw" suffix instead.
  char *rawname;
  size_t len = 1+strrchr(compname,'.')-compname;
  rawname = new char[len+strlen(RAW_SUFFIX)];
  strncpy(rawname, compname, len);
  if (preserve_input) {
    strncpy(strrchr(rawname,'.'), RAW_SUFFIX, strlen(RAW_SUFFIX));
  } else
    rawname[len-1] = '\0';

  // Do the expansion, deleting the compressed file (unless --preserve).
  slim_expander_t *expander = new slim_expander_t(compname);
  if (ignore_crc)
    expander->set_ignore_crc();

  FILE *fp = fopen(rawname, "wb");
  if (fp == NULL) {
    cerr << "slim: " << rawname << ": Unable to open file for writing\n";
    delete [] rawname;
    delete expander;
    return;
  }
  unsigned char *buf = new unsigned char[debug_buf_size];

  size_t written=0, read=0;
  do {
    read = expander->read(buf, debug_buf_size);
    if (read <= 0)
      break;
    written = fwrite(buf, 1, read, fp);
  } while (written > 0);

  fclose(fp);
  delete [] buf;
  delete [] rawname;
  delete expander;

  if (!preserve_input)
    unlink(compname);
}
