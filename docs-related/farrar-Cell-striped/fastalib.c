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
#include <ctype.h>
#include <string.h>
#include <malloc.h>

#include "striped.h"
#include "fastalib.h"
#include "utils.h"

#define READ_BUFFER_SIZE (256 * 1024)
#define SEQ_NAME_SIZE    (512)

FASTA_LIB *openFastaLib (char *file)
{
  FILE *fp;

  FASTA_LIB *lib;

  if ((fp = fopen (file, "rb")) == NULL) {
    fprintf (stderr, "Unable to open file %s\n", file);
    exit (-1);
  }

  lib = (FASTA_LIB *) _memalign (sizeof (FASTA_LIB), 16);
  if (!lib) {
    fprintf (stderr, "Unable to allocate memory for library header\n");
    exit (-1);
  }

  lib->readBuffer = (char *) _memalign (READ_BUFFER_SIZE, 16);
  if (!lib->readBuffer) {
    fprintf (stderr, "Unable to allocate memory for read buffer\n");
    exit (-1);
  }

  lib->seqBuffer = (unsigned char *) _memalign (MAX_SEQ_LENGTH, 16);
  if (!lib->seqBuffer) {
    fprintf (stderr, "Unable to allocate memory for sequence\n");
    exit (-1);
  }

  lib->seqName = (char *) _memalign (SEQ_NAME_SIZE, 16);
  if (!lib->seqName) {
    fprintf (stderr, "Unable to allocate memory for sequence name\n");
    exit (-1);
  }

  lib->size = (int) fread (lib->readBuffer, sizeof (char), READ_BUFFER_SIZE, fp);
  if (lib->size == 0 && !feof (fp)) {
    fprintf (stderr, "Error (%d) reading fasta file\n", ferror (fp));
    exit (-1);
  }

  lib->pos = 0;

  lib->fp = fp;

  lib->sequences = 0;
  lib->residues = 0;

  return lib;
}

static int
readNextBlock (FASTA_LIB *lib)
{
  FILE *fp = lib->fp;
  size_t size;

  size = fread (lib->readBuffer, sizeof (char), READ_BUFFER_SIZE, fp);
  if (lib->size == 0 && !feof (fp)) {
    fprintf (stderr, "Error (%d) reading fasta file\n", ferror (fp));
    exit (-1);
  }

  lib->pos = 0;
  lib->size = (int) size;

  return lib->size;
}

unsigned char *
nextFastaSeq (FASTA_LIB *lib, int maxLen, int *length)
{
  int inx;
  int size;
  int done;
  int len = 0;

  char *name;
  unsigned char *seq;

  name = lib->seqName;
  seq = lib->seqBuffer;

  /* check if we are at the end of the library */
  if (lib->size == 0) {
    lib->seqLength = 0;
    if (length) {
      *length = 0;
    }
    return NULL;
  }

  if (lib->pos == lib->size) {
    readNextBlock (lib);
  }

  inx = lib->pos;

  /* check for the start of a sequence */
  if (lib->readBuffer[inx] != '>') {
    fprintf (stderr, "Error parsing fasta file expecting > found %c\n",
             lib->readBuffer[inx]);
    exit (-1);
  }

  ++inx;

  /* read in the sequence name */
  len = 0;
  done = 0;
  do {
    if (inx >= lib->size) {
      size = readNextBlock (lib);
      if (size == 0) {
        if (length) {
          *length = 0;
        }
        return NULL;
      }
      inx = lib->pos;
    } else if (lib->readBuffer[inx] == '\n' || lib->readBuffer[inx] == '\r') {
      *name = '\0';
      done = 1;
    } else if (len < SEQ_NAME_SIZE - 1) {
      *name++ = lib->readBuffer[inx];
      len++;
    }
    ++inx;
  } while (!done);

  lib->pos = inx;

  /* read in the sequence */
  len = 0;
  done = 0;
  do {
    if (inx >= lib->size) {
      size = readNextBlock (lib);
      if (size == 0) {
        *seq = '\0';
        done = 1;
      }
      inx = 0;
    } else if (isspace(lib->readBuffer[inx])) {
      ++inx;
    } else if (lib->readBuffer[inx] == '>') {
      *seq = '\0';
      done = 1;
    } else if (len >= MAX_SEQ_LENGTH) {
      fprintf (stderr, "Sequence %s exceeds maximum length\n", lib->seqName);
      exit (-1);
    } else {
      int value = AMINO_ACID_VALUE[(int)lib->readBuffer[inx]];
      if (value == -1) {
        fprintf (stderr, "Unknown amino acid %c in sequence %s\n",
                 lib->readBuffer[inx], lib->seqName);
        exit (-1);
      }
      *seq++ = (char) value;
      inx++;
      len++;
    }
  } while (!done);

  lib->pos = inx;
  if (len > maxLen) {
    fprintf (stderr, "Truncating sequence %s (%d) exceeds max length\n",
             truncElipse (lib->seqName, 30), len);
    seq = lib->seqBuffer + maxLen;
    len = maxLen;
  }

  lib->sequences++;
  lib->residues += len;

  lib->seqLength = len;
  if (length) {
    *length = len;
  }

  inx = 16 - (len % 16);
  while (inx--) {
    *seq++ = ALPHA_SIZE;
  }
  *seq = '\0';

  return lib->seqBuffer;
}

void closeFastaLib (FASTA_LIB *lib)
{
  fclose (lib->fp);
    
  _memfree (lib->readBuffer);
  _memfree (lib->seqBuffer);
  _memfree (lib->seqName);

  _memfree (lib);
}

#define BEGIN_SEQ      0
#define READ_HEADER    1
#define READ_SEQUENCE  2

FASTA_DB *buildFastaDb (char *file, int maxLen)
{
  FILE *fp;
  int len;

  int state;

  int remaining;
  int offset;
  int count;
  int inx;

  char *seq;
  char *ptr;
  char *block;
  char *start;

  FASTA_LIB *lib;
  FASTA_DB  *db;

  ptr = NULL;

  lib = openFastaLib (file);

  printf ("Scanning library: %s\n", file);

  nextFastaSeq (lib, maxLen, &len);
  while (len > 0) {
    nextFastaSeq (lib, maxLen, &len);
  }

  db = (FASTA_DB *) _memalign (sizeof (FASTA_DB), 16);
  if (!db) {
    fprintf (stderr, "Unable to allocate memory for db header\n");
    exit (-1);
  }

  db->residues = lib->residues;
  db->sequences = lib->sequences;

  db->hdrBuffer = NULL;

  closeFastaLib (lib);

  printf ("Loading library: %s\n", file);

  if ((fp = fopen (file, "rb")) == NULL) {
    fprintf (stderr, "Unable to open file %s\n", file);
    exit (-1);
  }

  db->fpLib = fp;

  /* allocate space for all the sequences */
  db->seqBuffer = (char *) _memalign (db->residues, 16);
  if (!db->seqBuffer) {
    fprintf (stderr, "Unable to allocate sequence buffer %d\n", db->residues);
    exit (-1);
  }

  /* allocate space for the header offsets */
  db->hdrPos = (uint32_t *) _memalign (db->sequences * sizeof(uint32_t), 16);
  if (!db->hdrPos) {
    fprintf (stderr, "Unable to allocate header index %d\n", db->sequences);
    exit (-1);
  }

  /* allocate space for the sequence lengths */
  db->seqLen = (uint16_t *) _memalign (db->sequences * sizeof(uint16_t), 16);
  if (!db->seqLen) {
    fprintf (stderr, "Unable to allocate sequence lengths %d\n", db->sequences);
    exit (-1);
  }

  block = (char *) _memalign (READ_BUFFER_SIZE, 16);
  if (!block) {
    fprintf (stderr, "Unable to allocate read buffer %d\n", READ_BUFFER_SIZE);
    exit (-1);
  }

  state = BEGIN_SEQ;
  remaining = 0;
  offset = 0;
  count = 0;
  inx = 0;

  seq = db->seqBuffer;

  while (count < db->residues) {

    /* check if we need more data to process */
    if (!remaining) {
      if (!(remaining = fread (block, sizeof(char), READ_BUFFER_SIZE, fp))) {
        fprintf (stderr, "Error reading sequence data\n");
        exit (-1);
      }
      ptr = block;
    }

    if (state == BEGIN_SEQ) {
      /* check for the start of a sequence */
      if (*ptr != '>') {
        fprintf (stderr, "Error parsing: expecting '>' found '%c'\n", *ptr);
        exit (-1);
      }
      db->hdrPos[inx] = offset;
      state = READ_HEADER;
    } else if (state == READ_HEADER) {
      if (*ptr == '\n' || *ptr == '\r') {
        db->seqLen[inx] = 0;
        state = READ_SEQUENCE;
        start = seq;
      }
    } else {
      if (isspace(*ptr)) {
        /* do nothing */
      } else if (*ptr == '>') {
        ++inx;
        db->hdrPos[inx] = offset;
        state = READ_HEADER;
      } else {
        if (db->seqLen[inx] < maxLen) {
          db->seqLen[inx]++;
          *seq++ = (char) AMINO_ACID_VALUE[(int)*ptr];
        }
        count++;
      }
    }

    --remaining;
    ++offset;
    ++ptr;
  }

  return db;
}

void closeFastaDb (FASTA_DB *db)
{
  fclose (db->fpLib);

  _memfree (db->seqLen);
  _memfree (db->hdrPos);
  _memfree (db->seqBuffer);

  if (db->hdrBuffer) {
    _memfree (db->hdrBuffer);
  }

  _memfree (db);
}

int  getFastaSeqHdr (FASTA_DB *lib, int index, char *buffer, int bufSize)
{
  int i;
  int size;

  int32_t length = 0;

  char temp[SEQ_NAME_SIZE];
  char *ptr;

  /* seek to the header structure */
  if (fseek (lib->fpLib, lib->hdrPos[index], SEEK_SET) == -1) {
    fprintf (stderr, "Error seeking header index (%d)\n", lib->hdrPos[index]);
    exit (-1);
  }

  /* read a block of the header */
  size = fread (temp, sizeof(char), SEQ_NAME_SIZE, lib->fpLib);
  ptr = temp;

  /* make sure we have a header */
  if (*ptr != '>') {
    fprintf (stderr, "Error parsing: expecting '>' found '%c'\n", *ptr);
    exit (-1);
  }

  ++ptr;
  while (isspace (*ptr)) {
    ptr++;
  }

  i = 0;
  --bufSize;
  while (i < bufSize && size && *ptr != '\n' && *ptr != '\r') {
    *buffer++ = *ptr++;
    length++;
    ++i;
  }
  *buffer = 0;

  return length;
}

