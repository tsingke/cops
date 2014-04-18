/******************************************************************
  Copyright 2008 by Michael Farrar.  All rights reserved.
  This program may not be sold or incorporated into a commercial product,
  in whole or in part, without written consent of Michael Farrar.  For 
  further information regarding permission for use or reproduction, please 
  contact: Michael Farrar at farrar.michael@gmail.com.
*******************************************************************/

/*
  Written by Michael Farrar, 2008.
  Please send bug reports and/or suggestions to farrar.michael@gmail.com.
*/

#include <stdlib.h>

#ifdef LINUX
#include <sys/time.h>
#endif

#ifdef WIN32
#include <time.h>
#include <sys/timeb.h>
#endif

#include "utils.h"

void *_memalign (size_t size, unsigned int align)
{
  char *ptr;
  void *ret = NULL;

  unsigned long pad;
  unsigned long offset;
  
  /* round the padding to a power of two minimun pad of 16 */
  pad = 1 << 5;
  while (pad < align) {
    pad <<= 1;
  }

  ptr = (char *) malloc(size + sizeof(void *) + (pad-1));
  if (ptr != NULL) {
    offset = (pad - (unsigned long)(ptr + sizeof(void *))) & (pad-1);
    ret = (void *) (ptr + sizeof(void *) + offset);
    *((void **)(ret)-1) = (void *) (ptr);
  }

  return ret;
}

void _memfree (void *ptr)
{
  if (ptr) {
    ptr = *((void **)(ptr)-1);
    free (ptr);
  }
}

#ifdef LINUX
void _stime(double *dtime)
{
  struct timeval tv;

  gettimeofday (&tv, NULL);

  *dtime = (double) tv.tv_sec;
  *dtime = *dtime + (double) (tv.tv_usec) / 1000000.0;
}
#endif
#ifdef WIN32
void _stime (double *dtime)
{
  struct _timeb buf;

  _ftime_s (&buf);

  *dtime = (double) buf.time;
  *dtime = *dtime + (double) buf.millitm / 1000.0;
}
#endif

char *truncElipse (char *str, int len)
{
    int i;
    char *ptr;

    ptr = str + len;
    *ptr-- = '\0';
    for (i = 0; i < 3 && *ptr; ++i, --ptr) {
        if (*ptr) {
            *ptr = '.';
        }
    }

    return str;
}     
