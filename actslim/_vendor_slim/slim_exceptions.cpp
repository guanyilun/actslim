/// \file slim_exceptions.cpp
/// Define exceptions thrown in slim.

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
#include <errno.h>


static const size_t MAXLEN=256;

bad_file::bad_file(const char *filename, const char *message) {
  str=new char[MAXLEN];

  snprintf(str, MAXLEN, "Cannot open file %s%s", filename, message);
}


bad_file::bad_file(const bad_file &bf) {
  str = new char[1+strlen(bf.str)];
  strcpy(str, bf.str);
}

bad_file::~bad_file() {
  delete [] str;
  str = NULL;
}

void bad_file::mesg() const {
  cerr << str << endl;
}


bad_output_file::bad_output_file(const char *filename, const char *read_write_gerund) :
  bad_file(filename, "") 
{
  size_t len = strlen(str);
  snprintf(str+len, MAXLEN, " for %s: ", read_write_gerund);
  len = strlen(str);

  if (len+1 < MAXLEN)
    switch (errno) {
    case EACCES:
      snprintf(str+len, MAXLEN-len, "permission denied.");
      break;
    case ENOENT:
      snprintf(str+len, MAXLEN-len, "no such file or directory.");
      break;
    case EIO:
      snprintf(str+len, MAXLEN-len, "I/O error.");
      break;
    case ENOSPC:
      snprintf(str+len, MAXLEN-len, "no space left on device.");
      break;
    case EROFS:
      snprintf(str+len, MAXLEN-len, "read-only file system.");
      break;
    case EBADF:
      snprintf(str+len, MAXLEN-len, "bad file descriptor.");
      break;
    default:
      snprintf(str+len, MAXLEN-len, "(errno=%d).", errno);
      break;
    }
}


bad_output_file::~bad_output_file() 
{
  ;
}

