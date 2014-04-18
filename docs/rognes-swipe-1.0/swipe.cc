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

#include "swipe.h"

/* ARGUMENTS AND THEIR DEFAULTS */

#define DEFAULT_MAXMATCHES 10
#define DEFAULT_MINSCORE 1
#define DEFAULT_QUERYNAME "-"
#define DEFAULT_DATABASENAME ""
#define DEFAULT_GAPOPEN 11
#define DEFAULT_GAPEXTEND 1
#define DEFAULT_MATRIXNAME "BLOSUM62"
#define DEFAULT_THREADS 1
#define DEFAULT_VIEW 0
char * progname;
char * matrixname;
char * databasename;
char * queryname;
long maxmatches;
long minscore;
long gapopen;
long gapextend;
long threads;
long view;

int ssse3_present;

/* Other variables */

#define MAX_THREADS 64

long SCORELIMIT_7;
long SCORELIMIT_8;
long SCORELIMIT_16;
long SCORELIMIT_63;

char BIAS;

int chunk;

#define CHUNKDIVISION 100

pthread_mutex_t hitsmutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t workmutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t progressmutex = PTHREAD_MUTEX_INITIALIZER;

pthread_t pthread_id[MAX_THREADS];

unsigned seqnext = 0;

int fd_psq;
int fd_phr;

UINT32 * adr_pin;
BYTE * adr_psq;
BYTE * adr_phr;
off_t len_psq;
off_t len_pin;
off_t len_phr;

char * dbname;
char * dbdate;
unsigned long totalaa;
unsigned seqcount;
unsigned longest;
unsigned long phroffset;
unsigned long psqoffset;

BYTE gap_open_penalty;
BYTE gap_extend_penalty;

int hits_count;
int compute7 = 0;
int compute8 = 0;
int compute16 = 0;
int compute63 = 0;

long rounds7 = 0;
long rounds8 = 0;
long rounds16 = 0;
long rounds63 = 0;

BYTE * query_sequence = NULL;
int qlen = 0;
char * description = NULL;

char * score_matrix_7 = NULL;
unsigned char * score_matrix_8 = NULL;
short * score_matrix_16 = NULL;
long * score_matrix_63 = NULL;


// BIAS: uses a bias (4) in the score matrix and uses unsigned additions followed by
//         a subtraction of the bias

// Alt.:   used no bias in the score matrix, but uses signed additions
//         Also, the H-, E-, and F-values are based on -128 (0x80)

BYTE map2ncbi[32] = 
  {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  0, 10, 11, 12, 13,  0,
    14, 15, 16, 17, 18,  0, 19, 20, 21, 22, 23,  0,  0,  0,  0,  0 };

struct hits_entry
{
  long score;
  long seqno;
} * hits_list;

BYTE ascii2ncbi(char c)
{
  if ((c & 0xc0) == 0x40)
    return map2ncbi[c & 0x1f];
  else
    return 0;
}

char ncbi2ascii(BYTE b)
{
  char * s = "#ABCDEFGHIKLMNPQRSTVWXYZU#OJ";
  if (b <= 27)
    return s[b];
  else
    return '#';
}

void * xmalloc(size_t size)
{
  const size_t alignment = 16;
  void * t;
  posix_memalign(& t, alignment, size);
  
  // void * t = malloc(size);

  if (t==NULL)
    {
      fprintf(stderr, "Unable to allocate %zu bytes of memory.\n", size);
      exit(1);
    }
  
  //    printf("Memory allocated: %p %zu\n", t, size);

  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  void * t = realloc(ptr, size);
  
  if (t==NULL)
  {
    fprintf(stderr, "Unable to reallocate %zu bytes of memory.\n", size);
    exit(1);
  }

  //  printf("Memory allocated: %p %zu\n", t, size);

  return t;
}

int dots;

inline void show_progress()
{
  if (view==0)
    {  
      int c = 0;
      pthread_mutex_lock(&progressmutex);
      dots += 60;
      while (dots >= (int) seqcount)
	{
	  dots -= seqcount;
	  c++;
	}
      pthread_mutex_unlock(&progressmutex);
      while (c-- > 0)
	printf(".");
      fflush(stdout);
    }
}

void hits_enter(long seqno, long score)
{
  // show_progress();

  //  printf("Entering score %u from sequence no %u.\n", score, seqno);
  
  // find correct place

  pthread_mutex_lock(&hitsmutex);

  if (score < minscore) 
    {
      pthread_mutex_unlock(&hitsmutex);
      return;
    }

  long place = hits_count;

  while ((place > 0) && ((score > hits_list[place-1].score) || 
			 ((score == hits_list[place-1].score) && 
			  (seqno < hits_list[place-1].seqno))))
    place--;
  
  // move entries down
  
  long move = (hits_count < maxmatches ? hits_count : maxmatches - 1) - place;

  //  printf("Inserting at place %d, moving %d.\n", place, move);

  for (long j=move; j>0; j--)
    hits_list[place+j] = hits_list[place+j-1];
  
  // fill new entry

  if (place < maxmatches)
  {
    hits_list[place].score = score;
    hits_list[place].seqno = seqno;
    if (hits_count < maxmatches)
      hits_count++;
  }
  
  if (hits_count == maxmatches)
    minscore = hits_list[maxmatches-1].score;
  
  pthread_mutex_unlock(&hitsmutex);
}

void hits_init()
{
  compute8 = 0;
  compute16 = 0;
  compute63 = 0;
  hits_count = 0;
  hits_list = (struct hits_entry *) xmalloc(maxmatches * sizeof(struct hits_entry));
  dots = 0;
}

void show_d(unsigned seqno, int show)
{
#define HEADERMAX 2048

  unsigned char buffer[HEADERMAX];

  unsigned max = HEADERMAX;
  
  if (show > 0)
    max = show + 10;

  if (max > HEADERMAX)
    max = HEADERMAX;
  
  unsigned offset  = ntohl(adr_pin[phroffset + seqno]);
  unsigned offset2 = ntohl(adr_pin[phroffset + seqno + 1]);

  if (offset + 8 > offset2)
    {
      printf("Problem with header file.\n");
      exit(1);
    }

  unsigned length = offset2 - offset;

  if (length > max)
    length = max;
  
  unsigned got = pread(fd_phr, buffer, length, offset);
  
  if (got < length)
    {
      printf("Unable to read from %s.phr file.\n", databasename);
      exit(1);
    }

  unsigned char expected[7] = { 0x30, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a };
  
  for(int i=0;i<7;i++)
    if (buffer[i] != expected[i])
      {
	printf("Unexpected byte in header at position %d: %02x\n", 
	       i, (unsigned char) buffer[i]);
	exit(1);
      }

  unsigned start = 8;
  unsigned dlen = buffer[7];
  if (buffer[7] & 0x80)
    {
      if (length < 10)
	{
	  printf("Problem with header file.\n");
	  exit(1);
	}
      
      if (buffer[7] == 0x81)
	{
	  dlen = buffer[8];
	  start = 9;
	}
      else if (buffer[7] == 0x82)
	{
	  dlen = (buffer[8] << 8) + buffer[9];
	  start = 10;
	}
      else
	{
	  printf("Unexpected string length in header: %02x.\n", buffer[7]);
	  exit(1);
	}
    }
  
  if (start + dlen > got)
    dlen = got - start;

  if (show > 0)
    {
      if (show < dlen)
	dlen = show;
      
      for(int i=0; i<dlen; i++)
	putchar(buffer[start+i]);
      
      while (dlen < show)
	{
	  putchar(' ');
	  dlen++;
	}
    }
  else
    {
      for(int i=0; i<dlen; i++)
	putchar(buffer[start+i]);
    }
}

void show_descr(unsigned seqno)
{
  show_d(seqno, 74);
}

void show_descr_all(unsigned seqno)
{
  show_d(seqno, 0);
}

void hits_show()
{
  // show hits in list

  if(view == 0)
    {
      printf("Computed (7bit):   %d sequences in %ld rounds\n", compute7, rounds7);
      printf("Computed (16bit):  %d sequences in %ld rounds\n", compute16, rounds16);
      printf("Computed (63bit):  %d sequences in %ld rounds\n", compute63, rounds63);
      printf("\n");

      if (hits_count == 0)
	{
	  printf("\nSorry, no hits.\n");
	}
      else
	{
	  printf("%-74s %5s\n", "High-scoring sequences", "Score");
	  
	  for(int i=0; i<hits_count; i++)
	    {
	      int seqno = hits_list[i].seqno;
	      int score = hits_list[i].score;
	      show_descr(seqno);
	      printf(" %5d\n", score);
	    }
	}
    }
  else if ((view==8)||(view==9))
    {
      for(int i=0; i<hits_count; i++)
	{
	  int seqno = hits_list[i].seqno;
	  int score = hits_list[i].score;
	  printf("%s\t", description);
	  show_descr_all(seqno);
	  printf("\t%d\n", score);
	}
    }
  printf("\n");
}

void query_show()
{
  int linewidth = 60;
  for (unsigned i=0; i<strlen(description); i+=linewidth)
  {
    if (i==0)
      printf("Query description: %-60.60s\n", description+i);
    else
      printf("                   %-60.60s\n", description+i);
  }

#if 0
  for (int j=0; j<qlen; j+=linewidth)
  {
    if (j==0)
      printf("Query sequence:    ");
    else
      printf("                   ");

    for(int k=0; k<linewidth; k++)
    {
      if (j+k < qlen)
	putchar(ncbi2ascii(query_sequence[j+k]));
      else
	break;
    }
    printf("\n");
  }
#endif
}


void query_read(char * queryfilename)
{
  FILE * fp;
  if (strcmp(queryfilename, "-") == 0)
    fp = stdin;
  else
    fp = fopen(queryfilename, "r");
  
  if (!fp)
  {
    fprintf(stderr, "Cannot open query file \"%s\".\n", queryfilename);
    exit(1);
  }

  char line[LINE_MAX];
  
  if((fgets(line, LINE_MAX, fp) == NULL) || (line[0] != '>'))
  {
    fprintf(stderr, "Query file not in proper FASTA format.\n");
    exit(1);
  }
  
  if (line[strlen(line)-1] == '\n')
    line[strlen(line)-1] = 0;

  int dlen = strlen(line+1);
  description = (char*) xmalloc(dlen+1);
  strcpy(description, line+1);
  
  int size = LINE_MAX;
  query_sequence = (BYTE *) xmalloc(size);
  
  while((fgets(line, LINE_MAX, fp) != NULL) && (line[0] != '>'))
  {
    char * p = line;
    while(char c = *p++)
    {
      char m = ascii2ncbi(c);
      if (m)
      {
	if (qlen+1 >= size)
	{
	  size += LINE_MAX;
	  query_sequence = (BYTE*) xrealloc(query_sequence, size);
	}
	query_sequence[qlen++] = m;
      }
    }
  }
  query_sequence[qlen] = 0;

  if (fp != stdin)
    fclose(fp);
}

void args_getstring(int i, int argc, char **argv, char ** result, char * error)
{
  if (i+1 < argc)
  {
    *result = argv[i+1];
  }
  else
  {
    fprintf(stderr, "%s\n", error);
    exit(1);
  }
}

void args_getnum(int i, int argc, char **argv, long * result, char * error)
{
  if (i+1 < argc)
  {
    *result = atol(argv[i+1]);
  }
  else
  {
    fprintf(stderr, "%s\n", error);
    exit(1);
  }
}

void args_usage()
{
  // Free BLAST option letters: Hchjkux

  printf("Usage: swipe [options]\n");
  printf("\n");
  printf("-d s  sequence database file name (base name)  (required)\n");
  printf("-i s  query sequence file name (- means stdin) (default stdin)\n");
  printf("-M s  score matrix file                        (default BLOSUM62)\n");
  printf("-G n  gap opening penalty                      (default 11)\n");
  printf("-E n  gap extension penalty                    (default 1)\n");
  printf("-v n  maximum number of matches to display     (default 10)\n");
  printf("-j n  minimum score of matches to display      (default 1)\n");
  printf("-a n  number of threads (1-64)                 (default 1)\n");
  //  printf("-m    view type (0,8,9)                        (default 0)\n");
  printf("-h    show this help message\n");
}

void args_banner()
{
  char title[] = "SWIPE 1.0 - Smith-Waterman search using Inter-sequence Parallel Execution";
  if (view==0)
    {
      printf("%s\nCompiled %s %s\n\n", title, __DATE__, __TIME__);
    }
  else if ((view==9))
    {
      printf("#%s - Compiled %s %s\n", title, __DATE__, __TIME__);
    }
}

void args_show()
{
  args_banner();
  if (view == 0)
    {

      if (! ssse3_present)
      {
        printf("The performance is reduced because this CPU lacks SSSE3.\n\n");
      }

      printf("Database file:     %s\n", databasename);
      printf("Database title:    %s\n", dbname);
      printf("Database date:     %s\n", dbdate);
      printf("Database size:     %lu residues in %u sequences\n", totalaa, seqcount);
      printf("Longest db seq:    %u residues\n", longest);
      printf("Query file name:   %s\n", queryname);
      printf("Query length:      %d residues\n", qlen);
      query_show();
      printf("Score matrix:      %s\n", matrixname);
      printf("Gap penalty:       %ld+%ldk\n", gapopen, gapextend);
      printf("Min score shown:   %ld\n", minscore);
      printf("Max matches shown: %ld\n", maxmatches);
      printf("Threads:           %ld\n", threads);
    }
}
  
void args_help()
{
  args_banner();
  args_usage();
}

void args_init(int argc, char **argv)
{
  /* Set defaults */
  gapopen = DEFAULT_GAPOPEN;
  gapextend = DEFAULT_GAPEXTEND;
  matrixname = DEFAULT_MATRIXNAME;
  queryname = DEFAULT_QUERYNAME;
  databasename = DEFAULT_DATABASENAME;
  minscore = DEFAULT_MINSCORE;
  maxmatches = DEFAULT_MAXMATCHES;
  threads = DEFAULT_THREADS;
  view = DEFAULT_VIEW;

  progname = argv[0];

  opterr = 1;
  char short_options[] = "d:i:M:G:E:v:j:a:h";

  static struct option long_options[] =
  {
    {"db",               required_argument, NULL, 'd' },
    {"query",            required_argument, NULL, 'i' },
    {"matrix",           required_argument, NULL, 'M' },
    {"gapopen",          required_argument, NULL, 'G' },
    {"gapextend",        required_argument, NULL, 'E' },
    {"num_descriptions", required_argument, NULL, 'v' },
    {"min_score",        required_argument, NULL, 'j' },
    {"num_threads",      required_argument, NULL, 'a' },
    {"help",             no_argument,       NULL, 'h' },
    { 0, 0, 0, 0 }
  };
  
  int option_index = 0;
  int c;
  
  while (1)
    {
      c = getopt_long(argc, argv, short_options, long_options, &option_index);
      if (c == -1)
	break;

      switch(c)
	{
	case 'd':
	  /* database */
	  databasename = optarg;
	  break;
	  
	case 'i':
	  /* query */
	  queryname = optarg;
	  break;
	  
	case 'M':
	  /* matrix */
	  matrixname = optarg;
	  break;
	  
	case 'G':
	  /* gap open */
	  gapopen = atol(optarg);
	  break;
	  
	case 'E':
	  /* gap extend */
	  gapextend = atol(optarg);
	  break;
	  
	case 'v':
	  /* max matches shown */
	  maxmatches = atol(optarg);
	  break;
	  
	case 'j':
	  /* min score threshold */
	  minscore = atol(optarg);
	  break;
	  
	case 'a':
	  /* threads */
	  threads = atol(optarg);
	  break;
	  
	case 'h':
	  args_help();
	  exit(1);
	  break;
	  
	case '?':
	default:
	  printf("Usage: %s [OPTIONS]\n", progname);
	  printf("  -h, --help                 show help\n");
	  printf("  -d, --db=database          database base name (required)\n");
	  printf("  -q, --query=FILE           file name of query sequence (- means stdin)\n");
	  printf("  -M, --matrix=FILE          file name of score matrix\n");
	  printf("  -G, --gapopen=NUM          gap open panalty\n");
	  printf("  -E, --gapextend=NUM        gap extenstion penalty\n");
	  printf("  -v, --num_descriptions=NUM number of sequence descriptions to be shown\n");
	  printf("  -j, --min_score=NUM        minimum score of sequences to be shown\n");
	  printf("  -a, --num_threads=NUM      number of threads to run\n");
	  
	  exit(1);
	  break;
	}
    }
  
  if ((threads < 1) || (threads > MAX_THREADS))
    {
      fprintf(stderr, "Illegal number of threads specified (max %d).\n", MAX_THREADS);
      exit(1);
    }
  else if (strlen(databasename) == 0)
    {
      fprintf(stderr, "No database specified.\n");
      exit(1);
    }
  else if (!((view==0)||(view==8)||(view==9)))
    {
      fprintf(stderr, "Illegal view type.\n");
      exit(1);
    }
  else if ((gapopen < 0) || (gapextend < 0) || ((gapopen + gapextend) < 1))
    {
      fprintf(stderr, "Illegal gap penalties.\n");
      exit(1);
    }

  gap_open_penalty = gapopen + gapextend;
  gap_extend_penalty = gapextend;
}


void score_matrix_dump()
{
  for(int i=0; i<32; i++)
  {
    printf("%2d ", i);
    for(int j=0; j<32; j++)
      {
	printf("%2ld", score_matrix_63[(i<<5) + j]);
      }
    printf("\n");
  }
}

void score_matrix_read(char * matrixfilename)
{
  char line[LINE_MAX];
  char order[LINE_MAX];
  char *p, *q;
  char c;
  int a, b, i;
  int symbols = 0;
  long sc, lo, hi; 
  int read;
  
  score_matrix_7 = (char *) xmalloc(32*32*sizeof(char));
  score_matrix_8 = (unsigned char *) xmalloc(32*32*sizeof(char));
  score_matrix_16 = (short *) xmalloc(32*32*sizeof(short));
  score_matrix_63 = (long *) xmalloc(32*32*sizeof(long));
  memset(score_matrix_63, -1, 32*32*8);
  
  FILE * fp = fopen(matrixfilename, "r");
  
  if (!fp)
  {
    fprintf(stderr, "Cannot open score matrix file \"%s\".\n", matrixfilename);
    exit(1);
  }

  hi = -100;
  lo = 100;

  while(fgets(line, LINE_MAX, fp) != NULL)
  {
    p = line;
    c = *p++;
    
    switch(c)
    {   
    case '\n':
    case '#':
      /* ignore blank lines and comments starting with # */
      break;
      
    case '\t':
    case ' ':
      /* read order of symbols, copy non-whitespace chars */
      q = order;
      while ((c = *p++))
	if (strchr(" \t\n", c) == NULL)
	  *q++ = c;
      *q = 0;
      symbols = strlen(order);
      break;
      
    default:
      /* ordinary lines */
      a = ascii2ncbi(c);
      for (i=0; i<symbols; i++)
	{
	  if (sscanf(p, "%ld%n", & sc, & read) == 0)
	    {
	      fprintf(stderr, "Problem parsing score matrix file.\n");
	      exit(1);
	    }
	  
	  b = ascii2ncbi(order[i]);
	  
	  if ((a==0) || (b==0))
	    sc = -1;
	  
	  if (sc < lo)
	    lo = sc;
	  if (sc > hi)
	    hi = sc;
	  
	  //	  printf("a,b,adr: %d %d %d\n", a, b, (a<<5)+b);

	  score_matrix_63[(a<<5) + b] = sc;
	  
	  p += read;
	}
      break;
    }
  }
  
  BIAS = - lo;
  SCORELIMIT_7  = 128 - hi;
  SCORELIMIT_8  = 256 - hi;
  SCORELIMIT_16 = 65536 - hi;
  
  for(a=0;a<32;a++)
    for(b=0;b<32;b++)
      {
	sc = score_matrix_63[(a<<5) + b];
	
	score_matrix_7[(a<<5) + b] = (char) sc;
	score_matrix_8[(a<<5) + b] = (unsigned char) (BIAS + sc);
	score_matrix_16[(a<<5) + b] = (short) sc;
      }
  
  fclose(fp);

}



void score_matrix_init()
{
  score_matrix_read(matrixname);
  //  score_matrix_dump();
}


void vector_fill(BYTE * vector, BYTE value)
{
  memset(vector, value, 16);
}

void vector_print(BYTE * vector)
{
  for(int i=0; i<16; i++)
    printf(" %02x", vector[i]);
}

void vector_print_word(WORD * vector)
{
  for(int i=0; i<8; i++)
    printf(" %04x", vector[i]);
}

UINT32 getseqlen(UINT32 x)
{
  UINT32 c = psqoffset + x ;
  unsigned a = ntohl(adr_pin[c+1]);
  unsigned b = ntohl(adr_pin[c]);
  return a - b - 1;
}

int getwork(int * first, int * last)
{
  int status = 0;
  pthread_mutex_lock(&workmutex);
  if(seqnext < seqcount)
    {
      * first = seqnext;
      int next = seqnext + chunk;
      if (next > (int) seqcount)
	next = seqcount;
      seqnext = next;
      * last = next - 1;
      status = 1;
    }
  pthread_mutex_unlock(&workmutex);
  return status;
}

void openfiles()
{
  int len = strlen(databasename);
  char * name = (char*)xmalloc(len+5);
  strcpy(name, databasename);

  strcpy(name + len, ".psq");
  fd_psq = open(name, O_RDONLY, 0);
  if (fd_psq < 0)
  {
    fprintf(stderr, "Unable to open file %s.\n", name);
    exit(1);
  }

  strcpy(name + len, ".phr");
  fd_phr = open(name, O_RDONLY, 0);
  if (fd_phr < 0)
  {
    fprintf(stderr, "Unable to open file %s.\n", name);
    exit(1);
  }

  free(name);
}

void closefiles(void)
{
  if (close(fd_psq))
  {
    fprintf(stderr, "Unable to close file %s.psq.\n", databasename);
    exit(1);
  }
  if (close(fd_phr))
  {
    fprintf(stderr, "Unable to close file %s.phr.\n", databasename);
    exit(1);
  }
}

void mmap_seq(int seqfirst, int seqlast,
	      BYTE * * start, BYTE * * base, long * length)
{
  /* mmap area of .psq-file fd including sequences seqfirst to seqlast */
  /* return pointer to mmap in area and base for access in base, plus length */

  UINT32 * seqindex = adr_pin + psqoffset;

  long firstbyte = ntohl(seqindex[seqfirst]);
  long nextbyte = ntohl(seqindex[seqlast+1]);

  long pagesize = getpagesize();
  off_t offset = firstbyte - firstbyte % pagesize;
  size_t len = nextbyte - offset;
  
  BYTE * st = (BYTE*) mmap(NULL, len, PROT_READ, MAP_SHARED, fd_psq, offset);
  
  if (st == MAP_FAILED)
    {
      fprintf(stderr, "Unable to memory map file.\n");
      exit(1);
    }
  
  *start = st;
  *base = st - offset;
  * length = len;

#ifdef DEBUG
  printf("Mapped sequences %d to %d. Memory start: %p. Memory base: %p. Length %d. Offset: %u. Pagesize: %u\n",
	 seqfirst, seqlast, *start, *base, len, offset, pagesize);

#endif

}

void munmap_seq(BYTE * start, long length)
{
  munmap((void*)start,(size_t)length);
}

void * worker(void *)
{
  BYTE* dprofile = (BYTE*) xmalloc(4*16*32);
  BYTE** qtable = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
  BYTE* hearray = (BYTE*) xmalloc(qlen*32);

  for(int i=0; i<qlen; i++)
    {
      qtable[i] = dprofile + 64*query_sequence[i];
    }
    
  int seqfirst;
  int seqlast;

  long * scores = (long*) xmalloc(chunk * sizeof(long));

  long in_count, out_count;

  long * in_list = (long*) xmalloc(chunk * sizeof(long));
  long * out_list = (long*) xmalloc(chunk * sizeof(long));
  long * tmp_list;
      
  while(getwork(&seqfirst, &seqlast))
    {
      BYTE * start;
      BYTE * base;
      long length;
      mmap_seq(seqfirst, seqlast, & start, & base, & length);
      UINT32 * seqindex = adr_pin + psqoffset;
      long sequences = seqlast - seqfirst + 1;
      
      out_count = sequences;
      for (int i=0; i<sequences; i++)
	out_list[i] = seqfirst + i;
      
#if 1

      tmp_list = in_list;
      in_list = out_list;
      out_list = tmp_list;
      in_count = out_count;

      if (in_count > 0)
	{
	  pthread_mutex_lock(&hitsmutex);
	  compute7 += in_count;
	  rounds7++;
	  pthread_mutex_unlock(&hitsmutex);

	  search7(qtable,
		  gap_open_penalty,
		  gap_extend_penalty,
		  (BYTE*) score_matrix_7,
		  dprofile,
		  hearray,
		  in_count,
		  in_list,
		  scores,
		  base,
		  seqindex,
		  qlen);
	  
	  out_count = 0;
	  
	  for (int i=0; i<in_count; i++)
	    {
	      long seqno = in_list[i];
	      long score = scores[i];
	      
	      if (score < SCORELIMIT_7)
		{
		  hits_enter(seqno, score);
		}
	      else
		{
		  out_list[out_count++] = seqno;
		}
	    }
	}
#endif
	  
#if 1

      tmp_list = in_list;
      in_list = out_list;
      out_list = tmp_list;
      in_count = out_count;

      if (in_count > 0)
	{
	  pthread_mutex_lock(&hitsmutex);
	  compute16 += in_count;
	  rounds16++;
	  pthread_mutex_unlock(&hitsmutex);
	  
	  search16((WORD**)qtable,
		   gap_open_penalty,
		   gap_extend_penalty,
		   (WORD*)(score_matrix_16),
		   (WORD*)(dprofile),
		   (WORD*)(hearray),
		   in_count,
		   in_list,
		   scores,
		   base,
		   seqindex,
		   qlen);
	  
	  out_count = 0;

	  for (int i=0; i<in_count; i++)
	    {
	      long seqno = in_list[i];
	      long score = scores[i];
	      
	      if (score < SCORELIMIT_16)
		{
		  hits_enter(seqno, score);
		}
	      else
		{
		  out_list[out_count++] = seqno;
		}
	    }
	}
#endif
      
#if 1

      tmp_list = in_list;
      in_list = out_list;
      out_list = tmp_list;
      in_count = out_count;

      if (in_count > 0)
	{
	  pthread_mutex_lock(&hitsmutex);
	  compute63 += in_count;
	  rounds63++;
	  pthread_mutex_unlock(&hitsmutex);

	  for (int i=0; i<in_count; i++)
	    {
	      long seqno = in_list[i];
	      
	      char * dbegin = (char*) base + ntohl(seqindex[seqno]);
	      char * dend = (char*) base + ntohl(seqindex[seqno+1]) - 1;
	      long score = fullsw(dbegin,
				  dend,
				  (char*) query_sequence, 
				  (char*) query_sequence + qlen, 
				  (long*) hearray,
				  score_matrix_63,
				  gap_open_penalty,
				  gap_extend_penalty);
	      
	      hits_enter(seqno, score);
	    }
	}
      
#endif
      
      munmap_seq(start, length);
    }

  free(dprofile);
  free(qtable);
  free(hearray);
  free(scores);
  free(in_list);
  free(out_list);


  pthread_exit(NULL);
}

void spawn_threads()
{
  long t;
  pthread_attr_t attr;
  void * status;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  chunk = ( seqcount + CHUNKDIVISION*threads) / CHUNKDIVISION / threads;
  //  if (chunk == 0)
  // chunk = 1;

  //  sleep(3);

  for(t=0; t<threads; t++)
    {
      if (pthread_create(pthread_id + t, &attr, worker, (void *)t))
	{
	  fprintf(stderr, "Cannot create thread.\n");
	  exit(1);
	}
    }
  
  
  for(t=0; t<threads; t++) {
    if (pthread_join(pthread_id[t], &status))
      {
	fprintf(stderr, "Cannot join thread.\n");
	exit(1);
      }
  }

  pthread_attr_destroy(&attr);

}

void do_search()
{
  if (view==0)
    {
      printf("Searching...\n");
      //      printf("Searching:         ");
      //fflush(stdout);
      spawn_threads();
    }
  else
    spawn_threads();
}  

void mmap_one(char * name, void * * adr, off_t * len)
{
  int fd = open(name, O_RDONLY, 0);
  if (fd < 0)
  {
    fprintf(stderr, "Unable to open file %s.\n", name);
    exit(1);
  }
  * len = lseek(fd, 0, SEEK_END);
  * adr = mmap(0, *len, PROT_READ, MAP_SHARED, fd, 0);

  //  printf("Mapped file %s to address %p length %ld\n", name, *adr, *len);

  if (*adr == MAP_FAILED)
    {
      fprintf(stderr, "Unable to map file %s of %ld bytes in memory. It may be too large.\n", name, *len);
    exit(1);
  }
}

void mmap_init()
{
  char * name_pin = (char*)xmalloc(strlen(databasename)+5);
  strcpy(name_pin, databasename);
  strcat(name_pin, ".pin");
  mmap_one(name_pin, (void**) & adr_pin, & len_pin);
  free(name_pin);

#if 0
  char * name_phr = (char*)xmalloc(strlen(databasename)+5);
  strcpy(name_phr, databasename);
  strcat(name_phr, ".phr");
  mmap_one(name_phr, (void**) & adr_phr, & len_phr);
  free(name_phr);
#endif

  char * p = (char*) adr_pin;
  // unsigned dummy0 = ntohl(*(UINT32*)p);
  p += 4;
  // unsigned dummy1 = ntohl(*(UINT32*)p);
  p += 4; 
  unsigned namelen = ntohl(*(UINT32*)p);
  p += 4;
  dbname = (char*) xmalloc(namelen+1);
  strncpy(dbname, p, namelen);
  dbname[namelen] = 0;
  p += namelen;
  unsigned datelen = ntohl(*(UINT32*)p);
  p += 4;
  dbdate = (char*) xmalloc(datelen+1);
  strncpy(dbdate, p, datelen);
  dbdate[datelen] = 0;
  p += datelen;
  if ((long)p & 3)
    p++;
  if ((long)p & 3)
    p++;
  if ((long)p & 3)
    p++;
  seqcount = ntohl(*(UINT32*)p);
  p += 4;
  totalaa = *(unsigned long*)p;
  p += 8;
  longest = ntohl(*(UINT32*)p);
  p += 4;
  phroffset = (p - (char*)adr_pin) / 4;
  psqoffset = phroffset + seqcount + 1;
}

void mmap_done()
{
  munmap(adr_pin, len_pin);
  //  munmap(adr_phr, len_phr);
}

#define cpuid(f,a,b,c,d) asm("cpuid": "=a" (a), "=b" (b), "=c" (c), "=d" (d) : "a" (f));


void cpu_features()
{
  unsigned int a __attribute__ ((unused));
  unsigned int b __attribute__ ((unused));
  unsigned int c,d;
  cpuid(1,a,b,c,d);
  //  printf("cpuid: %08x %08x %08x %08x\n", a, b, c, d);

  int has_mmx   = (d >> 23) & 1;
  int has_sse   = (d >> 25) & 1;
  int has_sse2  = (d >> 26) & 1;
  int has_sse3  = (c      ) & 1;
  int has_ssse3 = (c >>  9) & 1;
  int has_sse41 = (c >> 19) & 1;
  int has_sse42 = (c >> 20) & 1;


  if (0)
    {
      printf("MMX:    %d\n", has_mmx);
      printf("SSE:    %d\n", has_sse);
      printf("SSE2:   %d\n", has_sse2);
      printf("SSE3:   %d\n", has_sse3);
      printf("SSSE3:  %d\n", has_ssse3);
      printf("SSE4.1: %d\n", has_sse41);
      printf("SSE4.2: %d\n", has_sse42);
    }
  
  if (! has_sse2)
    {
      printf("Sorry, this program requires a CPU with SSE2.\n");
      exit(1);
    }

  ssse3_present = has_ssse3;
}

int main(int argc, char**argv)
{
  cpu_features();
  args_init(argc,argv);
  query_read(queryname);
  score_matrix_init();
  openfiles();
  mmap_init();
  args_show();
  hits_init();
  do_search();
  hits_show();
  closefiles();
  mmap_done();
}


