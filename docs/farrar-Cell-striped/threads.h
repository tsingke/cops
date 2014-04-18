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

#ifndef INCLUDE_THREADS_H
#define INCLUDE_THREADS_H

void *initThreads (int        count, 
                   void      *simdInfo,
                   FASTA_LIB *query,
                   FASTA_DB  *dbLib,
                   uint16_t  *results);

void runThreads (int count, void *data);

void freeThreads (int count, void *data);

#endif /* INCLUDE_THREADS_H */
