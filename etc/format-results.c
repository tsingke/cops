#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>
#include<ctype.h>
#include<math.h>


#define VITINITIAL	2
#define VITINLINED	3
#define VITPARTITIONS 4
#define VITFILTER	10
#define VITSTREAM	11


typedef struct _aaaa
{
	int prog;
	int nthrs;
	char hmm[30];
	int hmmlen;
	double mcs;
} reg;


int i = 0;

int maxthrs = 0;

int erro(char* msg)
{
	fprintf(stderr, "ERRO: %s na linha %d\n", msg, i);

	return -1;
}

void printreg(reg* r)
{
	printf("%s-%d\t%s:%d\t%.1f\n", 
			(r->prog==VITSTREAM)?"VITSTREAM":"VITFILTER",
			r->nthrs, r->hmm, r->hmmlen, r->mcs);
}


void printexcel(reg** rrs)
{
	int i;
	printf("Model");
	
	int hasPartial = 0, hasThreads = 0;
	
	for (i = 0; rrs[i]; i++)
	{	if (rrs[i]->prog <= VITPARTITIONS)
			hasPartial = 1;
		if (rrs[i]->prog > VITPARTITIONS)
			hasThreads = 1;
	}
	
	if (hasPartial)
	{
		printf("\t%s", "Initial");
		printf("\t%s", "Inlined");
		printf("\t%s", "Partitions");
	}
	
	if (hasThreads)
	{
		for (i = 1; i <= maxthrs; i++)
			if (i != 5 && i != 3 && i != 7)
				printf("\tHMMER-%d", i);		
		for (i = 1; i <= maxthrs; i++)
			if (i != 5 && i != 3 && i != 7)
				printf("\tCOPS-%d", i);
	}
	
	printf("\n");
	
	for (i = 0; rrs[i]; )
	{
		int mlen = rrs[i]->hmmlen;
		printf("%d", mlen);
	
		for( ; rrs[i] && rrs[i]->hmmlen == mlen; i++ )
			if (rrs[i]->nthrs != 3)
			{
				printf("\t%.1f", rrs[i]->mcs);
			}
		
		printf("\n");
	}
}


// -1 para p1 < p2, p1 - p2
// 0  para p1 == p2
// +1 para p1 > p2, p2 - p1

int compregs(const void* a1, const void* a2)
{
	reg* r1 = *((reg**) a1);
	reg* r2 = *((reg**) a2);

	if (r1->hmmlen!= r2->hmmlen)
		return (r1->hmmlen - r2->hmmlen);
	if (r1->prog  != r2->prog)
		return (r1->prog  - r2->prog);
	if (r1->nthrs != r2->nthrs)
		return (r1->nthrs - r2->nthrs);
		
	return -1;
}
		


int main(int argc, char** argv)
{
	int j, k;
	char buf[10010];
	char *ptr;
	reg** regs = calloc(10000, sizeof(reg*));
	
	if (argc < 2)
		return 0;
	
	FILE * fp = fopen(argv[1], "r");
	
	k = 0;
	for (i = 1; fgets(buf, 10000, fp); i++)
	{
		char *tokpat = " \t\n";
	
		//printf("Read line %d\n", i);

		if (buf[0] != 'X')
			continue;
	
		reg* rr = regs[k++] = calloc(1, sizeof(reg));
	
		// Formato:
		// XXXXX	Vitstream	1thrs	hmms-usar1/M0040-1-cysPrx_C.hmm	2864.8
		
		ptr = strtok(buf, tokpat);
		ptr = strtok(NULL,tokpat);
		
		if (!strncmp(ptr, "Vitstream-initial", 15))
			rr->prog = VITINITIAL;
		else
		if (!strncmp(ptr, "Vitstream-inlined", 15))
			rr->prog = VITINLINED;
		else
		if (!strncmp(ptr, "Vitstream-partiti", 15))
			rr->prog = VITPARTITIONS;
		else
		if (!strncmp(ptr, "Vitstream", 9))
			rr->prog = VITSTREAM;
		else
		if (!strncmp(ptr, "Vitfilter", 9))
			rr->prog = VITFILTER;
		else return erro("Programa");
		
		ptr = strtok(NULL, tokpat);
		rr->nthrs = atoi(ptr);
		if (rr->nthrs > maxthrs) 
			maxthrs = rr->nthrs;
		if (rr->nthrs <= 0)
			return erro("Threads");
		
		ptr = strtok(NULL, tokpat);
		
		for (j = 0; ptr[j]; j++)
			if (ptr[j] == '/')
			{	ptr += j;
				j = 0;
			}
		
		strcpy(rr->hmm, ptr);
		if (strlen(rr->hmm) > 25)
			return erro("modelo");
		
		for (j = 0; ptr[j] && !isdigit(ptr[j]) ; j++);
		
		rr->hmmlen = atoi(ptr+j);
		if (rr->hmmlen <= 0)
			return erro("HMM length");
		
		
		ptr = strtok(NULL, tokpat);

		rr->mcs = strtod(ptr, NULL);
		if (!isfinite(rr->mcs) || rr->mcs < 1.0 || rr->mcs > 1e6)
			return erro("MC/s");
	}
	
	
//	printf("READING done\n");
	qsort(regs, k, sizeof(reg*), compregs);
//	printf("SORTING done\n");
	
	//for (j = 0; j < k; j++)	printreg(regs[j]);
	
	printf("\n");
	printexcel(regs);
	
	return 0;
}


