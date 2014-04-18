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

#ifndef INCLUDE_FASTALIB_H
#define INCLUDE_FASTALIB_H

#include <stdio.h>

#define MAX_SEQ_LENGTH (64 * 1024)

typedef struct {
    char *readBuffer;

    char *seqName;
    unsigned char *seqBuffer;
    int seqLength;

    int pos;
    int size;

    FILE *fp;

    int sequences;
    int residues;
} FASTA_LIB;

FASTA_LIB *openFastaLib (char *file);
void closeFastaLib (FASTA_LIB *lib);

unsigned char *nextFastaSeq (FASTA_LIB *lib, int maxLen, int *length);

#define fastaSeqName(LIB) (LIB->seqName)

typedef struct {
  int32_t sequences;
  int32_t residues;

  uint16_t *seqLen;
  uint32_t *hdrPos;

  char *seqBuffer;
  char *hdrBuffer;

  FILE *fpLib;

} FASTA_DB;

FASTA_DB *buildFastaDb (char *file, int maxLen);
void closeFastaDb (FASTA_DB *db);

int  getFastaSeqHdr (FASTA_DB *lib, int index, char *buffer, int bufSize);

#endif /* INCLUDE_FASTALIB_H */
