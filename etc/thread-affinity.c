
#define CPUIDgcc(func,ax,bx,cx,dx) \
	__asm__ __volatile__ ("cpuid":	\
		"=a" (ax), "=b" (bx), "=c" (cx), "=d" (dx) :	\
		"a" (func), "c" (0));
//#ifdef _WIN32  __cpuid((int *)regs, (int)i);


void cpuid64bits(unsigned info, unsigned *eax, unsigned *ebx, unsigned *ecx, unsigned *edx)
{
	*eax = info;
	asm volatile
 		("mov %%ebx, %%edi;" /* 32bit PIC: don't clobber ebx */
		"cpuid;"
		"mov %%ebx, %%esi;"
		"mov %%edi, %%ebx;"
		:"+a" (*eax), "=S" (*ebx), "=c" (*ecx), "=d" (*edx)
		: :"edi");
}

/*	Check for Hyper-Threading:
	1. Detect CPU vendor using CPUID function 0
	2. Check for HTT bit 28 in CPU features EDX from CPUID function 1
	3. Get the logical core count from EBX[23:16] from CPUID function 1
	4. Get actual non-threaded CPU core count
		If vendor == 'GenuineIntel' this is 1 plus EAX[31:26] from CPUID function 4
		If vendor == 'AuthenticAMD' this is 1 plus ECX[7:0] from CPUID function 0x80000008
*/

int CPUID_check_hyperthreading()
{
	unsigned eax, ebx, ecx, edx;

	CPUIDgcc(1, eax, ebx, ecx, edx);
	unsigned HTTsupport = ((edx & 0x10000000) != 0);

	if (HTTsupport == 0)
	{	printf("Does not support Hyper-threading\n");
		return 0;
	}
	// Logical cores per CPU
	unsigned ncores = 0, nsiblings = (ebx >> 16) & 0xff; // EBX[23:16]

	// Get vendor
	char vendor[20];
	memset(vendor, 0, sizeof(vendor));
	CPUIDgcc(0, eax, ebx, ecx, edx);
	memcpy(vendor+0, &ebx, 4);
	memcpy(vendor+4, &edx, 4);
	memcpy(vendor+8, &ecx, 4);

	// Physical cores per CPU
	if (!strcmp(vendor, "GenuineIntel")) {
		// Get DCP cache info
		CPUIDgcc(4, eax, ebx, ecx, edx);
		ncores = ((eax >> 26) & 0x3f) + 1; // EAX[31:26] + 1
	}
	else if (!strcmp(vendor, "AuthenticAMD")) {
		// Get NC: Number of CPU cores - 1
		CPUIDgcc(0x80000008, eax, ebx, ecx, edx);
		ncores = ((unsigned)(ecx & 0xff)) + 1; // ECX[7:0] + 1
	}

	printf("CPU %s has %d Physical Cores and %d Logical Cores (siblings)\n", 
			vendor, ncores, nsiblings);
	unsigned HTTenabled = ncores < nsiblings;
	return HTTenabled;
}


	int enabledht = CPUID_check_hyperthreading();
	if (enabledht)
		printf("HAS HYPER THREADING!!\n");


	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);
	if (enabledht)
		CPU_SET(NTHREADS-1, &cpuset);
	else
		CPU_SET(0,&cpuset);

	threads[0] = pthread_self();
	pthread_setaffinity_np(threads[0], sizeof(cpu_set_t), &cpuset);


	
	CPU_ZERO(&cpuset);
	if (enabledht)
	{	int htcore = (i%2==0)*4 + ((i-1)/2);
		printf("SET HTthread %d in %d\n", i, htcore);
		CPU_SET(htcore%2, &cpuset);
	}
	else
		CPU_SET(i%4,&cpuset);

		
	for (i = 1; i < NTHREADS; i++)
	{
		pthr_args *argscopy = calloc(1, sizeof(pthr_args));
		memcpy(argscopy, &args, sizeof(pthr_args));
		argscopy->thrid	= i;
		
		printf("CPUSET %d: %x %d\n", (unsigned) cpuset);
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);

		if (pthread_create(&threads[i], &attr, viterbi_thr_loop, argscopy))
			return fprintf(stderr, "ERROR could not create worker thread\n");
	
	// Alternative method: change affinity after creation
	// pthread_setaffinity_np(threads[i], sizeof(cpu_set_t), &cpuset);
	}

	for (i = 0; i < NTHREADS; i++)
	{
		pthread_getaffinity_np(threads[i], sizeof(cpu_set_t), &cpuset);	
		for (j = 0; j < 8; j++)
			if (CPU_ISSET(j, &cpuset))
				break;
		printf("THR %d running on cpu %d\n", i, j);
	}

	
