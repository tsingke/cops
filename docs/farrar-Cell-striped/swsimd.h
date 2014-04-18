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

#ifndef INCLUDE_SWSIMD_H
#define INCLUDE_SWSIMD_H

#ifdef X86
#include <emmintrin.h>
#endif

#include "striped.h"
#include "fastalib.h"

void *simdInit (FASTA_LIB   *query,
                signed char *matrix,
                int          gapOpen,
                int          gapExtend);

void simdSearch (int            count,
                 unsigned char *buffer,
                 uint16_t      *length,
                 uint16_t      *results,
#ifdef X86
                 __m128i       *pvH1,
                 __m128i       *pvH2,
                 __m128i       *pvE,
#endif
                 void          *simdData);

void simdFree (void *simdData);

int simdMaxQueryLen (void);
int simdMaxSeqLen (void);


#endif /* INCLUDE_SWSIMD_H */
