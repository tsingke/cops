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
#ifdef LINUX
#include <stdint.h>
#endif
#include <pthread.h>
#include <string.h>
#include <malloc.h>
#include <errno.h>

#ifdef CELL
#include <libspe2.h>
#endif
#ifdef X86
#include <emmintrin.h>
#endif

#include "striped.h"
#include "fastalib.h"
#include "threads.h"
#include "swsimd.h"
#include "utils.h"

void *workerThread (void *arg);

typedef struct {
  uint32_t   count      __attribute__ ((aligned (16)));
  uint32_t   buffer;
  uint32_t   lengths;
  uint32_t   results;

  uint32_t   simdData;

#ifdef X86
  __m128i   *pvH1;
  __m128i   *pvH2;
  __m128i   *pvE;
#endif

#ifdef CELL
  void      *context;
#endif

  pthread_t  threadId;
} THREAD_INFO;

#ifdef X86
static int startThread = 0;

static pthread_mutex_t startMutex;
static pthread_cond_t startCondVar;
#endif

#ifdef CELL
extern spe_program_handle_t spuSwStriped;
#endif

void *initThreads (int        count, 
                   void      *simdInfo,
                   FASTA_LIB *query,
                   FASTA_DB  *dbLib,
                   uint16_t  *results)
{
  int i;
  int inx;

  int size;
  int block;
  int offset;

#ifdef X86
  size_t bufSize;
#endif

  int status;

  THREAD_INFO *info;
  
  info = (THREAD_INFO *) malloc (count * sizeof (THREAD_INFO));
  if (!info) {
    fprintf (stderr, "Unable to allocate thread info %d\n", count);
    exit (-1);
  }

  /* calculate the block size for each thread */
  block = dbLib->residues / count;
    
  inx = 0;
  size = 0;
  offset = 0;
  for (i = 0; i < count; ++i) {

    info[i].count = 0;
    info[i].buffer = (uint32_t) (dbLib->seqBuffer + offset);
    info[i].lengths = (uint32_t) (dbLib->seqLen + inx);
    info[i].results = (uint32_t) (results + inx);
    info[i].simdData = (uint32_t) simdInfo;

#ifdef X86
    bufSize = (query->seqLength + 16) * sizeof (short);
    info[i].pvH1 = (__m128i *) _memalign (bufSize, 16);
    if (!info[i].pvH1) {
        fprintf (stderr, "Unable to allocate H1 buffer %d\n", bufSize);
        exit (-1);
    }

    info[i].pvH2 = (__m128i *) _memalign (bufSize, 16);
    if (!info[i].pvH2) {
        fprintf (stderr, "Unable to allocate H2 buffer %d\n", bufSize);
        exit (-1);
    }

    info[i].pvE  = (__m128i *) _memalign (bufSize, 16);
    if (!info[i].pvE) {
        fprintf (stderr, "Unable to allocate E buffer %d\n", bufSize);
        exit (-1);
    }
#endif

    /* count the number of sequences and their size */
    size = 0;
    while ((size < block || inx % 8 != 0) && inx < dbLib->sequences) {
      size = size + dbLib->seqLen[inx];
      offset = offset + dbLib->seqLen[inx];
      info[i].count++;
      inx++;
    }
  }

#ifdef X86
  status = pthread_mutex_init (&startMutex, NULL);
  if (status) {
    fprintf (stderr, "Create mutex failed: %s\n", strerror (status));
    exit (1);
  }

  status = pthread_cond_init (&startCondVar, NULL);
  if (status) {
    fprintf (stderr, "Create cond var failed: %s\n", strerror (status));
    exit (1);
  }
#endif

  for (i = 0; i < count; ++i) {

#ifdef CELL
    info[i].context = spe_context_create (0, 0);
    if (!info[i].context) {
      fprintf (stderr, "Create context failed: %s\n", strerror (errno));
      exit (1);
    }

    status = spe_program_load (info[i].context, &spuSwStriped);
    if (status) {
      fprintf (stderr, "SPE program load failed: %s\n", strerror (errno));
      exit (1);
    }
#endif

    status = pthread_create(&info[i].threadId,
                            NULL,
                            workerThread,
                            &info[i]);
    if (status) {
      fprintf (stderr, "Create thread failed: %s\n", strerror (errno));
      exit (1);
    }
  }

#ifdef X86
  status = pthread_mutex_lock (&startMutex);
  if (status) {
    fprintf (stderr, "Lock mutex failed: %s\n", strerror (errno));
    exit (1);
  }

  /* wait for the threads to have started */
  while (startThread < count) {
    status = pthread_cond_wait(&startCondVar, &startMutex);
    if (status) {
      fprintf (stderr, "Wait cond failed: %s\n", strerror (errno));
      exit (1);
    }
  }

  status = pthread_mutex_unlock (&startMutex);
  if (status) {
    fprintf (stderr, "Unlock mutex failed: %s\n", strerror (errno));
    exit (1);
  }
#endif

#if CELL
  /* wait for the SPU to signal they are ready */
  for (i = 0; i < count; ++i) {
    uint32_t buffer;

    buffer = (uint32_t) simdInfo;
    status = spe_in_mbox_write (info[i].context, 
                                &buffer, 
                                1, 
                                SPE_MBOX_ALL_BLOCKING);
    if (status == -1) {
      fprintf (stderr, "SPE mailbox write failed: %s\n", strerror (errno));
      exit (1);
    }

    status = spe_out_intr_mbox_read (info[i].context, 
                                     &buffer, 
                                     1, 
                                     SPE_MBOX_ALL_BLOCKING);
    if (status == -1) {
      fprintf (stderr, "SPE mailbox read failed: %s\n", strerror (errno));
      exit (1);
    }
  }
#endif

  return info;
}

void runThreads (int count, void *data)
{
  int i;
  int status;

  THREAD_INFO *info = (THREAD_INFO *) data;
  
#ifdef X86
  status = pthread_mutex_lock (&startMutex);
  if (status) {
    fprintf (stderr, "Lock mutex failed: %s\n", strerror (errno));
    exit (1);
  }

  /* all the threads have reported.  let them run */
  startThread = 0;

  status = pthread_cond_broadcast (&startCondVar);
  if (status) {
    fprintf (stderr, "Broadcast cond failed: %s\n", strerror (errno));
    exit (1);
  }

  status = pthread_mutex_unlock (&startMutex);
  if (status) {
    fprintf (stderr, "Unlock mutex failed: %s\n", strerror (errno));
    exit (1);
  }
#endif

#ifdef CELL
  /* signal all the SPUs to start */
  for (i = 0; i < count; ++i) {
    simdSearch (info[i].count,
                (unsigned char *) info[i].buffer,
                (uint16_t *) info[i].lengths,
                (uint16_t *) info[i].results,
                (void *) info[i].context);
  }
#endif

  for (i = 0; i < count; ++i) {
    status = pthread_join (info[i].threadId, NULL);
    if (status) {
      fprintf (stderr, "Join failed: %s\n", strerror (errno));
      exit (1);
    }
  } 
}

void freeThreads (int count, void *data)
{
  int i;
  THREAD_INFO *info = (THREAD_INFO *) data;
  
#ifdef CELL
  int status;

  /* destroy all the contexts that were created */
  for (i = 0; i < count; ++i) {
    status = spe_context_destroy ((spe_context_ptr_t) info[i].context);
    if (status == -1) {
      fprintf (stderr, "SPE destroy context failed: %s\n", strerror (errno));
      exit (1);
    }
  }
#endif

#ifdef X86
  pthread_mutex_destroy (&startMutex);
  pthread_cond_destroy (&startCondVar);

  for (i = 0; i < count; ++i) {
      _memfree (info[i].pvH1);
      _memfree (info[i].pvH2);
      _memfree (info[i].pvE);
  }
#endif

  free (info);
}

void *workerThread (void *arg)
{
  int status;

  THREAD_INFO *info = (THREAD_INFO *) arg;
  
#ifdef CELL
  unsigned int entry = SPE_DEFAULT_ENTRY;

  status = spe_context_run (info->context, &entry, 0, 0, 0, 0);
  if (status == -1) {
    fprintf (stderr, "SPE program execution failed: %s\n", strerror (errno));
    exit (1);
  }
#endif

#ifdef X86
  status = pthread_mutex_lock(&startMutex);
  if (status) {
    fprintf (stderr, "Lock mutex failed: %s\n", strerror (errno));
    exit (1);
  }

  /* signal that we have started */
  ++startThread;
  status = pthread_cond_broadcast (&startCondVar);

  /* wait for the signal to start the calculations */
  while (startThread) {
    status = pthread_cond_wait(&startCondVar, &startMutex);
    if (status) {
      fprintf (stderr, "Wait cond failed: %s\n", strerror (errno));
      exit (1);
    }
  }

  status = pthread_mutex_unlock(&startMutex);
  if (status) {
    fprintf (stderr, "Lock mutex failed: %s\n", strerror (errno));
    exit (1);
  }

  simdSearch (info->count,
              (unsigned char *) info->buffer,
              (uint16_t *) info->lengths,
              (uint16_t *) info->results,
              info->pvH1,
              info->pvH2,
              info->pvE,
              (void *) info->simdData);
#endif

  pthread_exit (NULL);
  return NULL;
}

