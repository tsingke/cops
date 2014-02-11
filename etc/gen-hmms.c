
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#define LINESIZE 70
#define NUMSEQS  10
//#define NLINES	 400
//#define TOTALSIZE  4000

char *chimp = 
"AGACGTGTGGAGTCCCAGCAGAGGCCAACCTGTGTCTCTTCATCTCCGTGGGAAAGGTGCCCCCAAAGTG\
AAAGAGATGGCCTGGTGGAAAGCCTGGATTGAACAGGAGGGTGTCACAGTGAAGAGCAGCTCCCACTTCA\
ACCCAGACCCTGATGCAGAGACCCTCTACAAAGCCATGAAGGGGATCGGGACCAACGAGCAGGCTATCAT\
CGATGTGCTCACCAAGAGAAGCAACACGCAGCGGCAGCAGATCGCCAAGTCCTTCAAGGCTCAGTTCGGC\
AAGGACCTCACTGAGACCTTGAAGTCGGAGCTCAGTGGCAAGTTTGAGAGGCTCATTGTGGCCCTTATGT\
ACCCGCCATACAGATACGAAGCCAAGGAGCTGCATGACGCCATGAAGAGCTTAGGAACCAAGGAGGGTGT\
CATCATTGAGATCCTGGCCTCTCGGACCAAGAACCAGCTGCGGGAGATAATGAAGGCGTATGAGGAAGAC\
TATGGGTCCAGCCTGGAGGAGGACATCCAAGCAGACACAAGTGGCTACCTGGAGAGGATCCTGGTGTGCC\
TCCTGCAGGGCAGCAGGGATGATGTGAGCAGCTTTGTGGACCCGGGACTGGCCCTCCAAGATGCACATGA\
TCTGTATGCGGCAGGCGAGAAGATTCGTGGGACTGATGAGATGAAATTCATCACCATCCTGTGCACGCGC\
AGTGCCACTCACCTGCTGAGAGTGTTTGAAGAGTATGAGAAAATTGCCAACAAGAGCATTGAGGACAGCA\
TCAAGAGTGAGACCCATGGCTCGCTGGAGGAGGCCATGCTCACTGTGGTGAAATGCACCCAAAACCTCCA\
CAGCTACTTTGCAGAGAGACTCTACTATGCCATGAAGGGAGCAGGGACGCGTGATGGGACCCTGATAAGA\
AACATCGTTTCAAGGAGCGAGATTGACTTAAATCTTATCAAGTGTCACTTCAAGAAGATGTATGGCAAGA\
CCCTCAGCAGCATGATCATGGAAGACACCAGCGGCGACTACAAGAACGCCCTGCTGAGCCTGGTGGGCAG\
CGACCCCTGAGGCACAGAAGAACAAGAGCAAAGACCATGAAGCCAGAGTCTCCAGGACTCCTCAGTCAAC\
CTCAGCCATGGACGCAGGTTGGGTGTGAGGGGGGTCCCAGCCTTTCAGTCTTCTATTTCCCTATTTCCAG\
TGCTTTCCAGCTGGGTTTCTGACCCAGAGGGTGGAACCGGCCTGGGCTCCTCTTCCCAACTTCCTCCAGG\
TCATTTCCCAGTGTGAGCACAATGCCAACCTTAGTGTTTCTCCAGCCAGACAGATGCCTCAGCATGAAGG\
GCTTGGGGACTTGTGAGTCATTCCTTCCTCCCTGCAGGAGCTTCCCAAACTGGTCACAGAGTCTCCTGGG\
CACAGGTTATAGAGACCCCAGCCCCATTCCCATCTACTGAAACAGGGTCTCCACAAGAGGGGCCAGGGAA\
TATGGGTTTTTAACAAGCGTCTTACGAAACACTTCTCTATCATGCAGCCAGAGAGCTGGCTGGGAGCCCT\
TTTGTTTTAGAACACACATCCTTCAGCAGCTGAGAAACGAACATGAATCCATCCCAACCAAGATGCCATT\
AACATTCATCTAAAAATATTAGGCTCTAAATGGACGAAAAATTCTCTCGCCATCTTAATAACAAAATAAA\
CTACAAATTCCTGACCCAAGGACACTGTGTTATAAGAGGCGTGGGCTCCCCTGGTGGCTGACCAGGTCAG\
CTGCCCTGGCCTTGCACCCCTCTGCATGCAGCACAGAAGGGTGTGACCATGCCCTCAGCACCACTCTTGT\
CCCCACTGAACGGCAACTGAGACTGGGTACCTGGAGATTCTGAAGTGCCTTTGCTATGATTTTCAAAATA\
ATAAAGATTTGTATTCAACTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";


int main(int argc, char** argv)
{
	char* letters =	"ACDEFGHIKLMNPSTRVQYW....";
					// "ABCDEFGHIJKLMNOPQRSTUVWXYZ......";
					//"ACGT..   ";
			// B J X Z = ambiguidades!
			// U O = special encoded a.a.'s

	int i, j, k;
	int nchars = strlen(letters);
	char data[LINESIZE*2];
	int countempties = 0, countgens = 0;

	if (argc < 2)
		return fprintf(stderr, "ERROR must specify size");

	int TOTALSIZE = atoi(argv[1]); //strlen(chimp); 
	
	if (TOTALSIZE == 0)
		return fprintf(stderr, "ERROR wrong size");

	printf("# STOCKHOLM 1.0\n\n");
	
	srand(time(NULL));

	for (j = 0; j*LINESIZE < TOTALSIZE + countempties; j++)
	{
		k = 0;
//if(0)	for (k = 0; k < NUMSEQS; k++, printf("\n"))

		{
			printf("CONSERV%02d   ", k);
			for (i = 0; i < LINESIZE && i + j*LINESIZE < TOTALSIZE + countempties; i++)
			{	int r = rand() % nchars;
				data[i] = letters[r];
				putc(data[i],stdout);

				if (data[i] == '.')
					countempties++;
			}
		}
		printf("\n");

		for (k++; k < NUMSEQS; k++, printf("\n"))
		{
			printf("CONSERV%02d   ", k);

			for (i = 0; i < LINESIZE && i + j*LINESIZE < TOTALSIZE + countempties; i++)
				if (rand() % 7 == 0)	// mutate
				{
					int r = rand() % nchars;
					char c = letters[r];
					if (c == ' ') continue;
					putc(c,stdout);
//					countempties++;
				}
				else
					 putc(data[i],stdout);
		}
		printf("\n");
	}

	printf("//\n");

	fprintf(stderr, "===== %d empties. Real size: %d =====\n", countempties, TOTALSIZE);
	return 0;
}

