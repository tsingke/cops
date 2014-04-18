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
#include <malloc.h>

#include "striped.h"
#include "matrix.h"

#define BUF_SIZE 512

char *skipSpaces (char *line)
{
  while (isspace (*line)) {
    ++line;
  }

  return line;
}

signed char *readMatrix (char *file)
{
  FILE *fp;
  char line[BUF_SIZE];

  signed char *matrix;

  int mark[ALPHA_SIZE];
  int order[ALPHA_SIZE];

  int done;
  int i;

  int errors = 0;
  int count = 0;

  if ((fp = fopen (file, "r")) == NULL) {
    fprintf (stderr, "Unable to open file %s\n", file);
    exit (-1);
  }

  matrix = (signed char *) malloc (MATRIX_SIZE * MATRIX_SIZE);
  if (!matrix) {
    fprintf (stderr, "Unable to allocate memory for scoring matrix\n");
    exit (-1);
  }

  for (i = 0; i < MATRIX_SIZE * MATRIX_SIZE; i++) {
    matrix[i] = 0;
  }

  /* initialize the order and mark arrays */
  for (i = 0; i < ALPHA_SIZE; ++i) {
    order[i] = -1;
    mark[i] = -1;
  }

  /* read the first line of the matrix giving the amino acid order */
  done = 0;
  while (!done && fgets (line, BUF_SIZE, fp) != NULL) {
    char *ptr = skipSpaces (line);
    if (*ptr && *ptr != '#') {

      while (*ptr && *ptr != '#') {
	int inx = AMINO_ACID_VALUE[*ptr];

	if (inx == -1) {
	  fprintf (stderr, "Unknown amino acid %c in %s\n", *ptr, file);
	  ++errors;
	} else if (mark[inx] != -1) {
	  fprintf (stderr, "Amino acid %c defined twice\n", *ptr);
	  ++errors;
	} else if (count > ALPHA_SIZE) {
	  /* this should not happen, but we will be safe */
	  fprintf (stderr, "Too many amino acids %d\n", count);
	  ++errors;
	} else {
	  order[count++] = inx;
	  mark[inx] = inx;
	}
	ptr = skipSpaces (ptr + 1);
      }

      done = 1;
    }
  }

  /* make sure all amino acids are defined */
  for (i = 0; i < ALPHA_SIZE-2; ++i) {
    if (order[i] < 0) {
      fprintf (stderr, "Missing column for amino acid %c\n", 
	       AMINO_ACIDS[i]);
      ++errors;
    }
  }

  for (i = 0; i < ALPHA_SIZE; ++i) {
    mark[i] = -1;
  }

  if (errors > 0) {
    fprintf (stderr, "Terminating due to errors in matrix file\n");
    exit (-1);
  }

  /* read the scores for the amino acids */
  while (fgets (line, BUF_SIZE, fp) != NULL) {
    signed char *row;
    char *ptr = skipSpaces (line);
    if (*ptr && *ptr != '#') {
      char aminoAcid = *ptr;
      int inx = AMINO_ACID_VALUE[*ptr];
      if (inx == -1) {
	fprintf (stderr, "Unknown amino acid %c in matrix\n", *ptr);
	++errors;
      } else if (mark[inx] != -1) {
	fprintf (stderr, "Row %c defined twice\n", *ptr);
	++errors;
      }

      row = &matrix[inx * MATRIX_SIZE];

      for (i = 0; i < ALPHA_SIZE-2; ++i) {
	int sign = 1;
	int num = 0;

	ptr = skipSpaces (ptr + 1);

	/* check the sign */
	if (*ptr == '-') {
	  sign = -1;
	  ++ptr;
	}

	do {
	  if (*ptr >= '0' && *ptr <= '9') {
	    num = num * 10 + (*ptr - '0');
	    ptr++;
	  } else {
	    char name[16];
	    char *pName;
	    if (isspace (*ptr)) {
	      pName = "space";
	    } else if (*ptr == 0) {
	      pName = "end of line";
	    } else {
	      name[0] = *ptr;
	      name[1] = 0;
	      pName = name;
	    }
	    fprintf (stderr, "Row %c Expecting digit found %s\n", 
		     aminoAcid, pName);
	    exit (-1);
	  }
	} while (*ptr && !isspace (*ptr));

	num = num * sign;

	if (num < -128 || num > 127) {
	  fprintf (stderr, "Weight %d out of range row %c\n", 
		   aminoAcid, num);
	  ++errors;
	  num = 0;
	}

	row[order[i]] = (char) num;
      }

      if (i < ALPHA_SIZE-2) {
	fprintf (stderr, "Amino acid row %c incomplete\n", aminoAcid);
	++errors;
      }

      mark[inx] = 1;
    }
  }

  /* make sure all amino acids are defined */
  for (i = 0; i < ALPHA_SIZE; ++i) {
    if (i == AMINO_ACID_VALUE['-'] || i == AMINO_ACID_VALUE['*']) continue;

    if (mark[i] < 0) {
      fprintf (stderr, "Missing row for amino acid %c\n", 
	       AMINO_ACIDS[i]);
      ++errors;
    }
  }

  if (errors) {
    fprintf (stderr, "Terminating due to errors in matrix %s\n", file);
  }

  fclose (fp);

  return matrix;
}

void freeMatrix (signed char *matrix)
{
  free (matrix);
}

