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

#ifndef INCLUDE_UTILS_H
#define INCLUDE_UTILS_H

void *_memalign (size_t size, unsigned int align);
void _memfree (void *ptr);

void _stime (double *dtime);

char *truncElipse (char *str, int len);

#endif /* INCLUDE_UTILS_H */
