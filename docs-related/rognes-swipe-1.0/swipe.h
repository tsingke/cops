/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2011 Torbjørn Rognes, University of Oslo, 
    Oslo University Hospital and Sencel Bioinformatics AS

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjørn Rognes <torognes@ifi.uio.no>, 
    Department of Informatics, University of Oslo, 
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <getopt.h>

//#include <smmintrin.h>
#include <tmmintrin.h>

// Should be 32bits integer
typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

#define WIDTH 32
#define WIDTH_SHIFT 5
#define BLOCKWIDTH 32

#define ext1 ".ssq"
#define ext2 ".ssi"
#define ext3 ".shd"
#define ext4 ".shi"

extern char BIAS;

//#define BIASED

#ifdef BIASED
#define ZERO 0x00
#else
#define ZERO 0x80
#endif

void vector_print(BYTE * vector);
void vector_print_word(WORD * vector);

void * xmalloc(size_t size);

extern int ssse3_present;

void search7(BYTE * * q_start,
	     BYTE gap_open_penalty,
	     BYTE gap_extend_penalty,
	     BYTE * score_matrix,
	     BYTE * dprofile,
	     BYTE * hearray,
	     long sequences,
	     long * seqnos,
	     long * scores,
	     BYTE * adr_psq,
	     UINT32 * seqindex,
	     long qlen);

void search16(WORD * * q_start,
	      WORD gap_open_penalty,
	      WORD gap_extend_penalty,
	      WORD * score_matrix,
	      WORD * dprofile,
	      WORD * hearray,
	      long sequences,
	      long * seqnos,
	      long * scores,
	      BYTE * adr_psq,
	      UINT32 * seqindex,
	      int qlen);

long fullsw(char * dseq,
	    char * dend,
	    char * qseq,
	    char * qend,
	    long * hearray, 
	    long * score_matrix,
	    BYTE gap_open_penalty,
	    BYTE gap_extend_penalty);
