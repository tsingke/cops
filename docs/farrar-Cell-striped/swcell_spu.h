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

#ifndef INCLUDE_SWCELL_SPU_H
#define INCLUDE_SWCELL_SPU_H

#include <stdint.h>

#include "striped.h"

#define MAX_SEQ_SIZE          32*1024

typedef struct {
  char     matrix[MATRIX_SIZE*MATRIX_SIZE]    __attribute__ ((aligned (16)));

  int32_t  esQuerySeq;
  int32_t  queryLen;

  char     gapOpen;
  char     gapExtend;

  char     ceiling;
  char     bias;

} SwSetupInfo;

#endif /* INCLUDE_SWCELL_SPU_H */
