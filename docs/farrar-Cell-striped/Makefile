#******************************************************************
#  Copyright 2008 by Michael Farrar.  All rights reserved.
#  This program may not be sold or incorporated into a commercial product,
#  in whole or in part, without written consent of Michael Farrar.  For 
#  further information regarding permission for use or reproduction, please 
#  contact: Michael Farrar at farrar.michael@gmail.com.
#*******************************************************************
#
#  Written by Michael Farrar, 2008.
#  Please send bug reports and/or suggestions to farrar.michael@gmail.com.
#

TARGET ?= x86
COMPILER ?= gcc

OBJS =	striped.o \
	fastalib.o \
	matrix.o \
	threads.o \
	utils.o

CFLAGS  = -O2 -c $(DEFINES) -m32
LFLAGS  = -O2 -lm -lpthread

SPU-CFLAGS = -O2 $(DEFINES) -Wl,-N

ifeq ($(TARGET),cell)
  OBJS    += swcell.o swcell_spu.o

  DEFINES  = -DCELL -DLINUX

ifeq ($(COMPILER),xlc)
  CC          = /opt/ibmcmp/xlc/8.2/bin/ppugxlc
  SPU-CC      = /opt/ibmcmp/xlc/8.2/bin/spugxlc
else
  CC          = /usr/bin/ppu32-gcc
  SPU-CC      = /usr/bin/spu-gcc  
  CFLAGS     += -W -Wall
  SPU-CFLAGS += -W -Wall -Winline -Wno-main
endif

  LFLAGS  += -lspe2 -Wl,-Map,striped.map 

  EMBEDSPU = /usr/bin/ppu-embedspu
endif

ifeq ($(TARGET),sim)
  OBJS    += swcell.o swcell_spu.o

  DEFINES  = -DCELL -DLINUX

ifeq ($(COMPILER),xlc)
  CC          = /opt/ibmcmp/xlc/8.2/bin/ppugxlc
  SPU-CC      = /opt/ibmcmp/xlc/8.2/bin/spugxlc
else
  CC          = /opt/cell/bin/ppu32-gcc
  SPU-CC      = /opt/cell/bin/spu-gcc  
  CFLAGS     += -W -Wall
  SPU-CFLAGS += -W -Wall -Winline -Wno-main
endif

  LFLAGS  += -lspe2 -Wl,-Map,striped.map 

  EMBEDSPU = /opt/cell/bin/ppu-embedspu
endif

ifeq ($(TARGET),x86)
  OBJS    += swsse2.o

  DEFINES  = -DX86 -DLINUX

  CC       = gcc
  CFLAGS  += -msse2
endif

SPU-TIMING = /opt/ibm/cell-sdk/prototype/bin/spu_timing

all: striped

clean:
	rm -f *.o striped *.map swcell_spu *.s *.timing

x86:
	make TARGET=x86 COMPILER=$(COMPILER)

sim:
	make TARGET=sim COMPILER=$(COMPILER)

cell:
	make TARGET=cell COMPILER=$(COMPILER)

striped: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $@

striped.o: striped.c striped.h fastalib.h matrix.h swsimd.h threads.h utils.h
	$(CC) $(CFLAGS) -c $<

fastalib.o: fastalib.c fastalib.h striped.h utils.h
	$(CC) $(CFLAGS) -c $<

matrix.o: matrix.c matrix.h striped.h
	$(CC)  $(CFLAGS) -c $<

swsse2.o: swsse2.c swsimd.h striped.h utils.h
	$(CC)  $(CFLAGS) -c $<

threads.o: threads.c threads.h striped.h
	$(CC)  $(CFLAGS) -c $<

utils.o: utils.c utils.h
	$(CC)  $(CFLAGS) -c $<

swcell.o: swcell.c swsimd.h striped.h swcell_spu.h utils.h
	$(CC)  $(CFLAGS) -c $<

swcell_spu.o: swcell_spu
	$(EMBEDSPU) -m32 spuSwStriped swcell_spu swcell_spu.o

swcell_spu: swcell_spu.c swcell_spu.h
	$(SPU-CC) -o $@ $(SPU-CFLAGS) -Wl,-Map,$@.map $<

timing: swcell_spu.s
	$(SPU-TIMING) -running-count -o swcell_spu.timing $<

swcell_spu.s: swcell_spu.c swcell_spu.h
	$(SPU-CC) -O2 $(DEFINES) -S -o $@ $<
