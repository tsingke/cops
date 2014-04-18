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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include <string.h>
#include <errno.h>

#include <libspe2.h>

#include "swsimd.h"
#include "striped.h"
#include "swcell_spu.h"

void *
simdInit(FASTA_LIB   *query,
         signed char *matrix,
         int          gapOpen,
         int          gapExtend)
         
{
  int i, j;

  int segSize;
  int count;

  int bias;
  int ceiling;

  SwSetupInfo *pSwData;

  char *pc;
  char *querySeq;
 
  if (query->seqLength > MAX_SEQ_SIZE) {
    fprintf (stderr, "Query sequence exceeds max %d\n", query->seqLength);
    exit (-1);
  }

  pSwData = (SwSetupInfo *) memalign (16, sizeof (SwSetupInfo));
  if (!pSwData) {
    fprintf (stderr, "Unable to allocate memory for SW data\n");
    exit (-1);
  }

  querySeq = (char *) memalign (16, query->seqLength + 16);
  if (!querySeq) {
    fprintf (stderr, "Unable to allocate query seq %d\n", query->seqLength);
    exit (-1);
  }

  /* reorder the query sequence */
  pc = querySeq;
  segSize = (query->seqLength + 15) / 16;
  count = segSize * 16;
  for (i = 0; i < segSize; ++i) {
    for (j = i; j < count; j += segSize) {
      if (j >= query->seqLength) {
        *pc++ = 0;
      } else {
        *pc++ = query->seqBuffer[j];
      }
    }
  }

  pSwData->esQuerySeq = (uint32_t) querySeq;
  pSwData->queryLen = (uint32_t) query->seqLength;

  /* Find the bias, min and max of the substitution matrix */
  bias = 0;
  ceiling = 0;
  for (i = 0; i < MATRIX_SIZE * MATRIX_SIZE; i++) {
    if (matrix[i] < bias) {
      bias = matrix[i];
    }
    if (matrix[i] > ceiling) {
      ceiling = matrix[i];
    }
  }

  pSwData->bias = (uint32_t) -bias;
  pSwData->ceiling = (uint32_t) ceiling;

  pSwData->gapOpen = (uint32_t) gapOpen;
  pSwData->gapExtend = (uint32_t) gapExtend;

  for (i = 0; i < MATRIX_SIZE * MATRIX_SIZE; i++) {
    pSwData->matrix[i] = matrix[i] - bias;
  }

  return pSwData;
}

void simdFree(void *simdData)
{
    SwSetupInfo *pSwData = (SwSetupInfo *) simdData;

    free ((void *) pSwData->esQuerySeq);
    free (pSwData);
}

void simdSearch (int            count,
                 unsigned char *buffer,
                 uint16_t      *length,
                 uint16_t      *results,
                 void          *context)
{
  int i;
  int status;

  uint32_t mboxData[4];

  mboxData[0] = count;
  mboxData[1] = (uint32_t) length;
  mboxData[2] = (uint32_t) buffer;
  mboxData[3] = (uint32_t) results;

  for (i = 0; i < 4; ++i) {
    status = spe_in_mbox_write (context, 
                                mboxData+i,
                                1, 
                                SPE_MBOX_ALL_BLOCKING);
    if (status == -1) {
      fprintf (stderr, "SPE mailbox write failed: %s\n", strerror (errno));
      exit (1);
    }
  }
}

int simdMaxQueryLen (void)
{
  return MAX_SEQ_SIZE;
}

int simdMaxSeqLen (void)
{
  return MAX_SEQ_SIZE;
}
