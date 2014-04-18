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

#ifndef INCLUDE_STRIPED_H
#define INCLUDE_STRIPED_H

typedef void SW_DATA;

#define ALPHA_SIZE 25
#define MATRIX_SIZE 32

#define MAX_BLK_CNT 32

extern const char AMINO_ACIDS[ALPHA_SIZE];
extern const int AMINO_ACID_VALUE[256];

#define SHORT_BIAS 32768

#ifdef LINUX
#include <stdint.h>
#endif
#ifdef WIN32
typedef          _int16  int16_t;
typedef unsigned _int16 uint16_t;
typedef          _int32  int32_t;
typedef unsigned _int32 uint32_t;

#define __attribute__(x)
#endif

#endif /* INCLUDE_STRIPED_H */
