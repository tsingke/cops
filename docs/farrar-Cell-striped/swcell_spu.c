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

#include <spu_intrinsics.h>
#include <spu_mfcio.h>
#include <stdint.h>

#include "striped.h"
#include "swcell_spu.h"

#define LIST_SIZE    64
#define MAX_DMA_SIZE 16384

#define GENERIC_DMA_ID        31

unsigned char MATRIX[MATRIX_SIZE*MATRIX_SIZE]    __attribute__ ((aligned (16)));

unsigned char DB_SEQ[MAX_SEQ_SIZE+128]           __attribute__ ((aligned (128)));

unsigned short DATA[MAX_SEQ_SIZE*3]              __attribute__ ((aligned (16)));

#define mfc_ceil16(v) (((v) + 15) & ~0x0f)

static int swWord (const int              queryLen,
                   const unsigned char   *dbSeq,
                   const int              dbLen,
                   vector unsigned short  vBias,
                   vector unsigned short  vMaxInc,
                   vector unsigned short  vGapOpen,
                   vector unsigned short  vGapExtend)
{
  int i, j;
  int iter;

  int score;

  int offset;

  vector unsigned int vJ;
  vector unsigned int vIter;

  vector unsigned short vE, vF, vH;
  vector unsigned short vHEF;
  vector unsigned short vHNext;
  vector unsigned short vETemp, vHTemp;

  vector unsigned short vMaxEF;
  vector unsigned short vMaxScore;

  vector unsigned char  *pvMatrix;
  vector unsigned short *pvData;
  vector unsigned short *pvLast;

  vector unsigned short vMinimums;
  vector unsigned short vTemp;

  vector unsigned short vScore;
  vector unsigned char  vScoreA;
  vector unsigned char  vScoreB;

  vector unsigned short vCeiling;

  vector unsigned short vMask;
  vector unsigned short vMask2;

  vector unsigned int  vCmp;
  vector unsigned char vPattern;

  vector unsigned short vHMin;

  vector unsigned char vInsMask = {  16,  17,   0,   1,   2,   3,   4,   5,
                                      6,   7,   8,   9,  10,  11,  12,  13, };

  vector unsigned char vEven    = { 128,   0, 128,   2, 128,   4, 128,   6,
                                    128,   8, 128,  10, 128,  12, 128,  14, };

  vector unsigned char vOdd     = { 128,   1, 128,   3, 128,   5, 128,   7,
                                    128,   9, 128,  11, 128,  13, 128,  15, };

  iter = ((queryLen + 15) / 16) * 2;

  score = 0;

  /* initialize pointers */
  pvData = (vector unsigned short *) DATA;
  pvMatrix = (vector unsigned char *) MATRIX;

  /* Load the bias to all elements of a constant */
  vBias = (vector unsigned short) spu_extend ((vector signed char) vBias);

  /* Load gap opening penalty to all elements of a constant */
  vGapOpen = (vector unsigned short) spu_extend ((vector signed char) vGapOpen);

  /* Load gap extension penalty to all elements of a constant */
  vGapExtend = (vector unsigned short) spu_extend ((vector signed char) vGapExtend);

  /* initialize the max score */
  vMask = spu_cmpgt (vBias, vGapOpen);
  vMaxScore = spu_sel (vGapOpen, vBias, vMask);
  vMask = spu_cmpgt (vMaxScore, vGapExtend);
  vMaxScore = spu_sel (vGapExtend, vMaxScore, vMask);

  /* initialize the minimum score score */
  vMinimums = spu_shuffle (vMaxScore, vMaxScore, vInsMask);

  /* calculate the ceiling */
  vTemp = spu_cmpeq (vMask, vMask);
  vCeiling = (vector unsigned short) spu_extend ((vector signed char) vMaxInc);
  vCeiling = spu_sub (vTemp, vCeiling);

  /* Convert the data from bytes to shorts */
  j = iter / 2;
  pvLast = pvData + (iter / 2) * 3;
  for (i = 0; i < j; ++ i) {
    vE = *(pvData+2);
    *(pvData+2) = spu_shuffle (vE, vE, vEven);
    *(pvLast+2) = spu_shuffle (vE, vE, vOdd);

    vH = *(pvData+1);
    *(pvData+1) = spu_shuffle (vH, vH, vEven);
    *(pvLast+1) = spu_shuffle (vH, vH, vOdd);

    pvData += 3;
    pvLast += 3;
  }

  pvData = (vector unsigned short *) DATA;
  pvLast = pvData + (iter - 1) * 3;

  /* load the next H value */
  vH = *(pvLast+1);
  vH = spu_shuffle (vH, vMinimums, vInsMask);

  /* load values E and F */
  vE = *(pvData+2);
  vF = vMinimums;

  /* fetch first data asap. */
  offset = dbSeq[0] * 2;
  vScoreA = pvMatrix[offset];
  vScoreB = pvMatrix[offset + 1];

  /* generate the score */
  vPattern = (vector unsigned char) *(pvData+0);
  vScore = (vector unsigned short) spu_shuffle (vScoreA, vScoreB, vPattern);
  vScore = spu_rlmask (vScore, -8);

  vH = spu_add (vH, vScore);
  vH = spu_sub (vH, vBias);

  vMaxEF = vE;

  for (i = 0; i < dbLen; ++i) {

    for (j = 0; j < iter; j++) {

      /* generate the score for the next vector */
      vPattern = (vector unsigned char) *(pvData+3);
      vScore = (vector unsigned short) spu_shuffle (vScoreA, vScoreB, vPattern);
      vScore = spu_rlmask (vScore, -8);

      vHNext = *(pvData+1);

      /* get max from H, E and F */
      vMask = spu_cmpgt (vH, vMaxEF);
      vHEF = spu_sel (vMaxEF, vH, vMask);

      *(pvData+1) = vHEF;

      /* Update highest score encountered this far */
      vMask = spu_cmpgt (vMaxScore, vH);
      vMaxScore = spu_sel (vH, vMaxScore, vMask);

      /* subtract the gap open penalty from H */
      vHTemp = spu_sub (vHEF, vGapOpen);
      vMask = spu_cmpgt (vHTemp, vMinimums);
      vHTemp = spu_sel (vMinimums, vHTemp, vMask);

      /* update E value */
      vETemp = spu_sub (vE, vGapExtend);
      vMask = spu_cmpgt (vETemp, vHTemp);
      vETemp = spu_sel (vHTemp, vETemp, vMask);
      *(pvData+2) = vETemp;

      vE = *(pvData+5);

      /* update F value */
      vF = spu_sub (vF, vGapExtend);
      vMask = spu_cmpgt (vF, vHTemp);
      vF = spu_sel (vHTemp, vF, vMask);

      /* add score to H */
      vHNext = spu_add (vHNext, vScore);
      vH = spu_sub (vHNext, vBias);

      /* get max from E and F */
      vMask = spu_cmpgt (vE, vF);
      vMaxEF = spu_sel (vF, vE, vMask);

      pvData += 3;
    }

    pvData = (vector unsigned short *) DATA;

  lazy_f:
    vJ = spu_splats ((unsigned int) 0);
    vIter = spu_splats ((unsigned int) iter);

    vHMin = spu_add (vMinimums, vGapExtend);
    vTemp = spu_sub (vGapOpen, vGapExtend);

    /* reset pointers to the start of the saved data */
    vH = *(pvData+1);

    /*  the computed F value is for the given column.  since */
    /*  we are at the end, we need to shift the F value over */
    /*  to the next column. */
    vF = spu_shuffle (vF, vMinimums, vInsMask);

    do {
      vHTemp = spu_sub (vH, vTemp);
      vMask = spu_cmpgt (vHTemp, vHMin);
      vHTemp = spu_sel (vHMin, vHTemp, vMask);

      vJ = spu_add (vJ, 1);
      vMask2 = (vector unsigned short) spu_cmpeq (vJ, vIter);

      vMask = spu_cmpgt (vF, vHTemp);
      vMask = spu_andc (vMask, vMask2);
      vCmp = spu_gather (vMask);


      /* update the H value with the new F */
      vMask = spu_cmpgt (vH, vF);
      vH = spu_sel (vF, vH, vMask);
      *(pvData+1) = vH;

      /* update E in case the new H value would change it */
      vE = *(pvData+2);
      vHTemp = spu_sub (vH, vGapOpen);
      vMask = spu_cmpgt (vE, vHTemp);
      vE = spu_sel (vHTemp, vE, vMask);
      *(pvData+2) = vE;

      vH = *(pvData+4);

      /* update F value and saturate */
      vF = spu_sub (vF, vGapExtend);
      vMask = spu_cmpgt (vF, vMinimums);
      vF = spu_sel (vMinimums, vF, vMask);

      pvData += 3;
    } while (spu_extract (vCmp, 0) != 0);

    /* this will be true only in the case where the F value has rolled */
    /* over to the next column.  */
    if (spu_extract (vMask2, 0) != 0) {
      goto lazy_f;
    }

    pvData = (vector unsigned short *) DATA;

    /* fetch first data asap. */
    offset = dbSeq[i + 1] * 2;
    vScoreA = pvMatrix[offset];
    vScoreB = pvMatrix[offset + 1];

    /* generate the score */
    vPattern = (vector unsigned char) *(pvData+0);
    vScore = (vector unsigned short) spu_shuffle (vScoreA, vScoreB, vPattern);
    vScore = spu_rlmask (vScore, -8);

    /* load values E. */
    vE = *(pvData+2);
    vMaxEF = vE;

    /* load the next H value */
    vH = *(pvLast+1);
    vH = spu_shuffle (vH, vMinimums, vInsMask);

    vH = spu_add (vH, (vector unsigned short) vScore);
    vH = spu_sub (vH, vBias);

    /* zero out F value. */
    vF = vMinimums;

    /*  check if we are about to overflow */
    vMask = spu_cmpgt (vMaxScore, vCeiling);
    vCmp = spu_gather (vMask);
    if (spu_extract (vCmp, 0) != 0) {
      break;
    }
  }

  /* find largest score in the vMaxScore vector */
  vTemp = spu_slqwbyte (vMaxScore, 8);
  vMask = spu_cmpgt (vTemp, vMaxScore);
  vMaxScore = spu_sel (vMaxScore, vTemp, vMask);

  vTemp = spu_slqwbyte (vMaxScore, 4);
  vMask = spu_cmpgt (vTemp, vMaxScore);
  vMaxScore = spu_sel (vMaxScore, vTemp, vMask);

  vTemp = spu_slqwbyte (vMaxScore, 2);
  vMask = spu_cmpgt (vTemp, vMaxScore);
  vMaxScore = spu_sel (vMaxScore, vTemp, vMask);

  /* store in temporary variable */
  vMaxScore = spu_sub (vMaxScore, vMinimums);
  score = (int) spu_extract (vMaxScore, 0);

  /* return largest score */
  return score;
}

static int swByte(const int              queryLen,
                  const unsigned char   *dbSeq,
                  const int              dbLen,
                  vector unsigned short  vBias,
                  vector unsigned short  vMaxInc,
                  vector unsigned short  vGapOpen,
                  vector unsigned short  vGapExtend)
{
  int i, j;
  int iter;

  int score;
  int wscore;

  int offset;

  vector unsigned int vJ;
  vector unsigned int vIter;

  vector unsigned short vE, vF, vH;

  vector unsigned short vHEF;
  vector unsigned short vHNext;
  vector unsigned short vETemp, vHTemp;

  vector unsigned short vMaxEF;
  vector unsigned short vMaxScore;

  vector unsigned char  *pvMatrix;
  vector unsigned short *pvData;
  vector unsigned short *pvLast;

  vector unsigned short vMinimums;
  vector unsigned short vTemp;

  vector unsigned char vScore;
  vector unsigned char vScoreA;
  vector unsigned char vScoreB;

  vector unsigned char vCeiling;

  vector unsigned char vMask;
  vector unsigned char vMask2;

  vector unsigned int  vCmp;
  vector unsigned char vPattern;

  vector unsigned short vHMin;

  vector unsigned char vInsMask = { 16,  0,  1,  2,  3,  4,  5,  6,
                                     7,  8,  9, 10, 11, 12, 13, 14, };

  iter = (queryLen + 15) / 16;

  score = 0;
  wscore = 0;

  /* initialize pointers */
  pvMatrix = (vector unsigned char *) MATRIX;

  /* initialize the minimum score */
  vMask = (vector unsigned char) spu_cmpgt (vBias, vGapOpen);
  vMinimums = spu_sel (vGapOpen, vBias, (vector unsigned short) vMask);
  vMask = (vector unsigned char) spu_cmpgt (vMinimums, vGapExtend);
  vMinimums = spu_sel (vGapExtend, vMinimums, (vector unsigned short) vMask);

  /* calculate the ceiling */
  vTemp = spu_cmpeq (vMaxInc, vMaxInc);
  vCeiling = (vector unsigned char) spu_sub (vTemp, vMaxInc);

  /* initialize the maximum score */
  vMaxScore = vMinimums;

  /* Zero out the storage vector */
  pvData = (vector unsigned short *) DATA;
  for (i = 0; i < iter; i++) {
    *(pvData+1) = vMaxScore;
    *(pvData+2) = vMaxScore;
    pvData += 3;
  }

  /* wait for the db sequence transfer to complete */
  mfc_write_tag_mask (1 << GENERIC_DMA_ID);
  mfc_read_tag_status_all ();

  /* fetch first data asap. */
  offset = dbSeq[0] * 2;
  vScoreA = pvMatrix[offset];
  vScoreB = pvMatrix[offset + 1];

  pvData = (vector unsigned short *) DATA;
  pvLast = pvData + (iter - 1) * 3;

  /* generate the score */
  vPattern = (vector unsigned char) *(pvData+0);
  vScore = spu_shuffle (vScoreA, vScoreB, vPattern);

  /* initialize values */
  vH = vMinimums;
  vE = vMinimums;
  vF = vMinimums;

  vH = spu_add (vH, (vector unsigned short) vScore);
  vH = spu_sub (vH, vBias);

  vMaxEF = vMinimums;

  for (i = 0; i < dbLen; ++i) {

    for (j = 0; j < iter; ++j) {

      /* generate the score for the next vector */
      vPattern = (vector unsigned char) *(pvData+3);
      vScore = spu_shuffle (vScoreA, vScoreB, vPattern);

      vHNext = *(pvData+1);

      /* get max from H, E and F */
      vMask = spu_cmpgt ((vector unsigned char) vH,
                         (vector unsigned char) vMaxEF);
      vHEF = spu_sel (vMaxEF, vH, (vector unsigned short) vMask);

      *(pvData+1) = vHEF;

      /* Update highest score encountered this far */
      vMask = spu_cmpgt ((vector unsigned char) vMaxScore,
                         (vector unsigned char) vH);
      vMaxScore = spu_sel (vH, vMaxScore, (vector unsigned short) vMask);

      /* subtract the gap open penalty from H */
      vHTemp = spu_sub (vHEF, vGapOpen);
      vMask = spu_cmpgt ((vector unsigned char) vHTemp,
			 (vector unsigned char) vMinimums);
      vHTemp = spu_sel (vMinimums, vHTemp, (vector unsigned short) vMask);

      /* update E value */
      vETemp = spu_sub (vE, vGapExtend);
      vMask = spu_cmpgt ((vector unsigned char) vETemp,
                          (vector unsigned char) vHTemp);
      vETemp = spu_sel (vHTemp, vETemp, (vector unsigned short) vMask);
      *(pvData+2) = vETemp;

      vE = *(pvData+5);

      /* update F value */
      vF = spu_sub (vF, vGapExtend);
      vMask = spu_cmpgt ((vector unsigned char) vF,
                         (vector unsigned char) vHTemp);
      vF = spu_sel (vHTemp, vF, (vector unsigned short) vMask);

      /* add score to H */
      vHNext = spu_add (vHNext, (vector unsigned short) vScore);
      vH = spu_sub (vHNext, vBias);

      /* get max from E and F */
      vMask = spu_cmpgt ((vector unsigned char) vE,
                         (vector unsigned char) vF);
      vMaxEF = spu_sel (vF, vE, (vector unsigned short) vMask);

      pvData += 3;
    }

  lazy_f:
    vJ = spu_splats ((unsigned int) 0);
    vIter = spu_splats ((unsigned int) iter);

    vHMin = spu_add (vMinimums, vGapExtend);
    vTemp = spu_sub (vGapOpen, vGapExtend);

    pvData = (vector unsigned short *) DATA;

    /* reset pointers to the start of the saved data */
    vH = *(pvData+1);

    /*  the computed F value is for the given column.  since */
    /*  we are at the end, we need to shift the F value over */
    /*  to the next column. */
    vF = spu_shuffle (vF, vMinimums, vInsMask);

    do {
      vHTemp = spu_sub (vH, vTemp);
      vMask  = spu_cmpgt ((vector unsigned char) vHTemp,
                          (vector unsigned char) vHMin);
      vHTemp = spu_sel (vHMin, vHTemp, (vector unsigned short) vMask);

      vJ = spu_add (vJ, 1);
      vMask2 = (vector unsigned char) spu_cmpeq (vJ, vIter);

      vMask = spu_cmpgt ((vector unsigned char) vF,
                         (vector unsigned char) vHTemp);
      vMask = spu_andc (vMask, vMask2);
      vCmp = spu_gather (vMask);


      /* update the H value with the new F */
      vMask = spu_cmpgt ((vector unsigned char) vH,
                         (vector unsigned char) vF);
      vH = spu_sel (vF, vH, (vector unsigned short) vMask);
      *(pvData+1) = vH;

      /* update E in case the new H value would change it */
      vE = *(pvData+2);
      vHTemp = spu_sub (vH, vGapOpen);
      vMask = spu_cmpgt ((vector unsigned char) vE,
                         (vector unsigned char) vHTemp);
      vE = spu_sel (vHTemp, vE, (vector unsigned short) vMask);
      *(pvData+2) = vE;

      vH = *(pvData+4);

      /* update F value and saturate */
      vF = spu_sub (vF, vGapExtend);
      vMask = spu_cmpgt ((vector unsigned char) vF,
                         (vector unsigned char) vMinimums);
      vF = spu_sel (vMinimums, vF, (vector unsigned short) vMask);

      pvData += 3;
    } while (spu_extract (vCmp, 0) != 0);

    /* this will be true only in the case where the F value has rolled */
    /* over to the next column.  */
    if (spu_extract (vMask2, 0) != 0) {
      goto lazy_f;
    }

    pvData = (vector unsigned short *) DATA;

    /* fetch next columns data */
    offset = dbSeq[i + 1] * 2;
    vScoreA = pvMatrix[offset];
    vScoreB = pvMatrix[offset + 1];

    /* generate the score */
    vPattern = (vector unsigned char) *(pvData+0);
    vScore = spu_shuffle (vScoreA, vScoreB, vPattern);

    /* load the next E value */
    vE = *(pvData+2);
    vMaxEF = vE;

    /* load the next H value */
    vH = *(pvLast+1);
    vH = spu_shuffle (vH, vMinimums, vInsMask);

    vH = spu_add (vH, (vector unsigned short) vScore);
    vH = spu_sub (vH, vBias);

    /* zero out F value. */
    vF = vMinimums;

    /*  check if we are about to overflow */
    vMask = spu_cmpgt ((vector unsigned char) vMaxScore, vCeiling);
    vCmp = spu_gather (vMask);
    if (spu_extract (vCmp, 0) != 0) {
      wscore = swWord (queryLen,
                       dbSeq + i + 1,
                       dbLen - i - 1,
                       vBias,
                       vMaxInc,
                       vGapOpen,
                       vGapExtend);
      break;
    }
  }

  /* find largest score in the vMaxScore vector */
  vTemp = spu_slqwbyte (vMaxScore, 8);
  vMask = spu_cmpgt ((vector unsigned char) vTemp,
                      (vector unsigned char) vMaxScore);
  vMaxScore = spu_sel (vMaxScore, vTemp, (vector unsigned short) vMask);

  vTemp = spu_slqwbyte (vMaxScore, 4);
  vMask = spu_cmpgt ((vector unsigned char) vTemp,
                      (vector unsigned char) vMaxScore);
  vMaxScore = spu_sel (vMaxScore, vTemp, (vector unsigned short) vMask);

  vTemp = spu_slqwbyte (vMaxScore, 2);
  vMask = spu_cmpgt ((vector unsigned char) vTemp,
                      (vector unsigned char) vMaxScore);
  vMaxScore = spu_sel (vMaxScore, vTemp, (vector unsigned short) vMask);

  vTemp = spu_slqwbyte (vMaxScore, 1);
  vMask = spu_cmpgt ((vector unsigned char) vTemp,
                      (vector unsigned char) vMaxScore);
  vMaxScore = spu_sel (vMaxScore, vTemp, (vector unsigned short) vMask);

  /* store in temporary variable */
  vMaxScore = spu_sub (vMaxScore, vMinimums);
  score = (int) spu_extract ((vector unsigned char) vMaxScore, 0);

  if (wscore > score) {
    score = wscore;
  }

  return score;
}

int main(uint64_t spuId, uint32_t argp, uint32_t envp)
{
  int i;
  int inx;

  int iter;

  int len;
  int seqLen;

  unsigned char *ptr;

  uint32_t  esWorkPtr;
  uint32_t  esPtr;

  uint32_t  seqCount;
  uint32_t  esSeqBuffer;
  uint32_t  esLengths;
  uint32_t  esResults;

  uint16_t  lengths[LIST_SIZE]                   __attribute__ ((aligned (16)));
  uint16_t  results[LIST_SIZE]                   __attribute__ ((aligned (16)));

  uint32_t  dbPad;

  SwSetupInfo workInfo;

  vector unsigned short vBias;
  vector unsigned short vMaxInc;
  vector unsigned short vGapOpen;
  vector unsigned short vGapExtend;

  vector unsigned short *pvData;
  vector unsigned short *pvLast;

  (void) spuId;
  (void) argp;
  (void) envp;

  /* fetch the work buffer from the mailbox */
  esWorkPtr = spu_read_in_mbox ();

  /* transfer the work information in the local store */
  mfc_get (&workInfo, esWorkPtr, sizeof(workInfo), GENERIC_DMA_ID, 0, 0);
  mfc_write_tag_mask (1 << GENERIC_DMA_ID);
  mfc_read_tag_status_all ();

  /* transfer the query sequence */
  seqLen = workInfo.queryLen;
  esPtr = workInfo.esQuerySeq;
  ptr = (unsigned char *) DATA;
  while (seqLen > 0) {
    len = mfc_ceil16 (seqLen);
    if (len > MAX_DMA_SIZE) {
      len = MAX_DMA_SIZE;
    }

    mfc_get (ptr, esPtr, len, GENERIC_DMA_ID, 0, 0);

    ptr += len;
    esPtr += len;
    seqLen -= len;
  }

  /* copy the matrix data */
  for (i = 0; i < MATRIX_SIZE*MATRIX_SIZE; ++i) {
    MATRIX[i] = (unsigned char) workInfo.matrix[i];
  }

  /* Load gap extension penalty to all elements of a constant */
  vGapExtend = (vector unsigned short) spu_splats (workInfo.gapExtend);

  /* Load gap opening penalty to all elements of a constant */
  vGapOpen = (vector unsigned short) spu_splats (workInfo.gapOpen);
  vGapOpen = spu_add (vGapOpen, vGapExtend);

  /* Load the bias to all elements of a constant */
  vBias = (vector unsigned short) spu_splats (workInfo.bias);

  /* Load the max possible score before overflow */
  vMaxInc = (vector unsigned short) spu_splats (workInfo.ceiling);

  /* wait for the query sequence transfer to complete */
  mfc_write_tag_mask (1 << GENERIC_DMA_ID);
  mfc_read_tag_status_all ();

  /* since the query sequence, H and E values are interleaved, do that now */
  iter = (workInfo.queryLen + 15) / 16;
  pvData = (vector unsigned short *) DATA;
  pvLast = pvData + iter * 3;
  for (i = iter - 1; i >= 0; --i) {
    pvData[i*3] = pvData[i];
    pvLast[i*3] = spu_sl (pvData[i], 8);
  }

  /* notify the host we ready to start */
  spu_write_out_intr_mbox (0);

  /* read info about the block to process */
  seqCount    = spu_read_in_mbox ();

  esLengths   = spu_read_in_mbox ();
  mfc_get (lengths, esLengths, sizeof(lengths), GENERIC_DMA_ID, 0, 0);
  esLengths += sizeof(lengths);

  esSeqBuffer = spu_read_in_mbox ();
  esResults   = spu_read_in_mbox ();

  inx = 0;

  /* wait for the list transfer to complete */
  mfc_write_tag_mask (1 << GENERIC_DMA_ID);
  mfc_read_tag_status_all ();

  /* process all the db sequences */
  while (seqCount--) {

    /* back up the db sequence until it starts on a cache line */
    dbPad = esSeqBuffer & 0x7f;

    ptr = DB_SEQ;
    esPtr = esSeqBuffer - dbPad;
    seqLen = lengths[inx] + dbPad;

    while (seqLen > 0) {
      len = mfc_ceil16 (seqLen);
      if (len > MAX_DMA_SIZE) {
        len = MAX_DMA_SIZE;
      }

      mfc_get (ptr, esPtr, len, GENERIC_DMA_ID, 0, 0);

      ptr += len;
      esPtr += len;
      seqLen -= len;
    }

    results[inx] = swByte(workInfo.queryLen,
                          DB_SEQ + dbPad,
                          lengths[inx],
                          vBias,
                          vMaxInc,
                          vGapOpen,
                          vGapExtend);

    esSeqBuffer += lengths[inx];

    /* check if results and lengths need to be transfered */
    if (++inx >= LIST_SIZE) {
      mfc_put (results, esResults, sizeof(results), GENERIC_DMA_ID, 0, 0);
      esResults += sizeof(results);

      len = sizeof(lengths);
      if (seqCount < LIST_SIZE) {
        len = mfc_ceil16 (seqCount * sizeof (uint16_t));
      }

      mfc_get (lengths, esLengths, len, GENERIC_DMA_ID, 0, 0);
      esLengths += sizeof(lengths);

      inx = 0;

      /* wait for the transfers to complete */
      mfc_write_tag_mask (1 << GENERIC_DMA_ID);
      mfc_read_tag_status_all ();
    }
  }

  /* check if there are any results that have not been written */
  if (inx) {
    len = mfc_ceil16 (inx * sizeof (uint16_t));
    mfc_put (results, esResults, len, GENERIC_DMA_ID, 0, 0);

    /* wait for the results transfer to complete */
    mfc_write_tag_mask (1 << GENERIC_DMA_ID);
    mfc_read_tag_status_all ();
  }

  return 0;
}
