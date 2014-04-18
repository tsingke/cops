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
#include <string.h>
#include <malloc.h>

#include "striped.h"
#include "matrix.h"
#include "fastalib.h"
#include "threads.h"
#include "swsimd.h"
#include "utils.h"

#define GAP_OPEN     10
#define GAP_EXTEND    2

const char AMINO_ACIDS[] = {
    '-', 'A', 'B', 'C', 'D', 'E', 'F', 'G',
    'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 
    'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z',
    '*',
};

const int AMINO_ACID_VALUE[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 24, -1, -1,  0, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, 10, 11, 12, 13, -1,
    14, 15, 16, 17, 18, -1, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1,
    -1,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, 10, 11, 12, 13, -1,
    14, 15, 16, 17, 18, -1, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

void printResults (FASTA_DB *dbLib, int threshold, int rptCount, uint16_t *results);

void printUsage (void)
{
  printf ("Usage: striped [-h] [-i num] [-e num] [-t num] [-c num] "
          "[-T num] [-M num] matrix query db\n");
  printf ("    -h       : this help message\n");
  printf ("    -i num   : gap init penalty (default %d)\n", -GAP_OPEN);
  printf ("    -e num   : gap extension penalty (default %d)\n", -GAP_EXTEND);
  printf ("    -t num   : minimum score threshold (default 20)\n");
  printf ("    -c num   : number of scores to be displayed (default 250)\n");
  printf ("    -T num   : number of threads doing calculations\n");
  printf ("    -M num   : max sequence length to process (default %d)\n", simdMaxSeqLen ());
  printf ("    matrix   : scoring matrix file\n");
  printf ("    query    : query sequence file (fasta format)\n");
  printf ("    db       : sequence database file (fasta format)\n");
}

int main (int argc, char **argv)
{
  int i;

  int penalty;

  int maxLen;

  int threshold = 1;
  int rptCount = 250;

  char *dbFile = NULL;
  char *queryFile = NULL;
  char *matrixFile = NULL;

  signed char *matrix;

  unsigned char *querySeq;
  int queryLen;

  int gapOpen = GAP_OPEN;
  int gapExtend = GAP_EXTEND;

  FASTA_LIB *queryLib;
  FASTA_DB *dbLib;

  int numThreads = 1;

  uint16_t *results;

  void *simdInfo;
  void *threadInfo;

  double startTime;
  double endTime;

  double mcups;

  if (argc < 4) {
    printUsage ();
    exit (-1);
  }

  maxLen = simdMaxSeqLen ();

  i = 1;
  while (i < argc) {
    if (i + 3 == argc) {
      /* should be matrix file name */
      matrixFile = argv[i];

    } else if (i + 2 == argc) {
      /* should be query file name */
      queryFile = argv[i];

    } else if (i + 1 == argc) {
      /* should be matrix file name */
      dbFile = argv[i];

    } else {
      /* process arguements */
      switch (argv[i][1]) {
      case 'h':
        printUsage ();
        break;
      case 'i':
        penalty = atoi (argv[++i]);
        if (penalty > 0 || penalty < -128) {
          fprintf (stderr, "Invalid gap init %d\n", penalty);
          fprintf (stderr, "Valid range is 0 - -128\n");
          exit (-1);
        }
        gapOpen = (unsigned char) -penalty;
        break;
      case 'e':
        penalty = atoi (argv[++i]);
        if (penalty > 0 || penalty < -128) {
          fprintf (stderr, "Invalid gap extension %d\n", penalty);
          fprintf (stderr, "Valid range is 0 - -128\n");
          exit (-1);
        }
        gapExtend = (unsigned char) -penalty;
        break;
      case 't':
        threshold = atoi (argv[++i]);
        if (threshold <= 0) {
          threshold = 1;
        }
        break;
      case 'c':
        rptCount = atoi (argv[++i]);
        if (rptCount < 10) {
          rptCount = 10;
        }
        break;
      case 'T':
        numThreads = atoi (argv[++i]);
        if (numThreads < 1) {
          numThreads = 1;
        }
        break;
      case 'M':
        maxLen = atoi (argv[++i]);
        if (maxLen > simdMaxSeqLen ()) {
          fprintf (stderr, "Max sequence length cannot exceed %d\n", simdMaxSeqLen ());
          exit (-1);
        } else if (maxLen < 1) {
          fprintf (stderr, "Max sequence length must be greater than 0\n");
          exit (-1);
        }
        break;
      default:
        fprintf (stderr, "Invalid option %s\n", argv[i]);
        printUsage ();
        exit (-1);
      }
    }
    ++i;
  }

  if (matrixFile == NULL) {
    fprintf (stderr, "Missing scoring matrix file\n");
    printUsage ();
    exit (-1);
  }

  if (queryFile == NULL) {
    fprintf (stderr, "Missing query sequence file\n");
    printUsage ();
    exit (-1);
  }

  if (dbFile == NULL) {
    fprintf (stderr, "Missing database file\n");
    printUsage ();
    exit (-1);
  }

  matrix = readMatrix (matrixFile);
  if (matrix == NULL) {
    fprintf (stderr, "Error reading matrix\n");
    exit (-1);
  }

  fprintf (stdout, "Striped Smith-Waterman written by Michael Farrar  (c) 2008\n");
  fprintf (stdout, "Please cite: Farrar, M. (2007) Striped Smith-Waterman speeds database\n");
  fprintf (stdout, "             searches six times over other simd implenentations.\n");
  fprintf (stdout, "             Bioinformatics, 23, 156-161\n\n");

  dbLib = buildFastaDb (dbFile, maxLen);

  queryLib = openFastaLib (queryFile);

  /* allocate space for the results list.  add extra space at the end  to */
  /* avoid buffer over runs because the cell implementation will always */
  /* writes the results in blocks with the size a multiple of 16 bytes. */
  results = (uint16_t *) _memalign (dbLib->sequences*sizeof(uint16_t)+16, 128);
  if (!results) {
    fprintf (stderr, "Unable to allocate results list %d\n", dbLib->sequences);
    exit (-1);
  }

  querySeq = nextFastaSeq (queryLib, maxLen, &queryLen);
  if (queryLen == 0) {
    fprintf (stderr, "Empty query sequence\n");
    exit (-1);
  }

  while (queryLen > 0) {

    /* the different architectures require different alignment, structures, */
    /* so call simdInit to do the architectural dependant stuff.  */
    simdInfo = simdInit (queryLib, matrix, gapOpen, gapExtend);

    /* create the smith-waterman threads. */
    threadInfo = initThreads (numThreads, simdInfo, queryLib, dbLib, results);

    printf ("\n\n%s vs %s\n", truncElipse (queryLib->seqName, 30), dbFile);
    printf ("Matrix: %s, Init: %d, Ext: %d\n", matrixFile, -gapOpen, -gapExtend);

    _stime(&startTime);

    runThreads (numThreads, threadInfo);

    _stime(&endTime);

    freeThreads (numThreads, threadInfo);

    printf ("\nScan time: %6.3f   Threads: %d\n", endTime - startTime, numThreads);
  
    printf ("\n");
    printf ("%d residues in query string\n", queryLen);
    printf ("%d residues in %d library sequences\n\n", 
            dbLib->residues, dbLib->sequences);

    mcups = (double) dbLib->residues * (double) queryLen;
    mcups = mcups / (endTime - startTime) / 1000000.0;
    printf ("%d MCUPS\n", (int) mcups);
    printf ("%d MCUPS per thread\n\n",(int) (mcups / (double) numThreads));

    printResults (dbLib, threshold, rptCount, results);

    simdFree (simdInfo);

    querySeq = nextFastaSeq (queryLib, maxLen, &queryLen);
  }

  closeFastaLib (queryLib);
  closeFastaDb (dbLib);

  freeMatrix (matrix);

  _memfree (results);

  return 0;
}

typedef struct _RESULT_STRUCT {
  unsigned int score;
  unsigned int index;
  struct _RESULT_STRUCT *prev;
  struct _RESULT_STRUCT *next;
} RESULT_STRUCT;
  
void printResults (FASTA_DB *dbLib, int threshold, int rptCount, uint16_t *results)
{
  int i;
  int count;

  char seqId[128];

  RESULT_STRUCT *list;
  RESULT_STRUCT *temp, *ptr;

  RESULT_STRUCT top, bottom;

  top.index = (unsigned int) -1;
  top.score = (unsigned int) -1;
  top.prev = NULL;
  top.next = &bottom;

  bottom.index = (unsigned int) -1;
  bottom.score = 0;
  bottom.prev = &top;
  bottom.next = NULL;

  list = (RESULT_STRUCT *) malloc (rptCount * sizeof (RESULT_STRUCT));
  if (!list) {
    fprintf (stderr, "Error allocating results struct %d\n", rptCount);
    exit (-1);
  }

  count = 0;

  for (i = 0; i < dbLib->sequences; ++i) {
    if (results[i] >= threshold) {
      ptr = &top;
      while (ptr && ptr->score >= results[i]) {
        ptr = ptr->next;
      }

      if (ptr != &bottom || count < rptCount) {
        if (count < rptCount) {
          temp = list + count;
          ++count;
        } else {
          temp = bottom.prev;
          bottom.prev = temp->prev;
          temp->prev->next = &bottom;
          if (ptr == temp) {
            ptr = &bottom;
          }
        }
          
        temp->score = results[i];
        temp->index = i;

        temp->prev = ptr->prev;
        temp->next = ptr;
        temp->prev->next = temp;
        ptr->prev = temp;
      }
    }
  }

  printf ("Score  Id\n");

  /* print the score and sequence id */
  ptr = top.next;
  while (ptr != &bottom) {
    getFastaSeqHdr (dbLib, ptr->index, seqId, sizeof (seqId));
    printf ("%5d  %s\n", ptr->score, truncElipse (seqId, 70));
    ptr = ptr->next;
  }

  free (list);
}
