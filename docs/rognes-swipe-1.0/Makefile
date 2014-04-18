LIBS = -lpthread
CXX=icpc
CC=icpc
CXXFLAGS=-Wall -Wno-missing-declarations -fast -xSSE2

PROG=swipe

all : $(PROG)

OBJS = swipe.o search16.o search7.o search63.o

DEPS = swipe.h Makefile

search63.o : search63.cc $(DEPS)

search16.o : search16.cc $(DEPS)

search7.o : search7.cc $(DEPS)

swipe.o : swipe.cc $(DEPS)

$(PROG) : $(OBJS)
	$(CC) -o $(PROG) $(OBJS) $(LIBS)
	strip swipe

clean :
	rm -f *.o *~ $(PROG)
