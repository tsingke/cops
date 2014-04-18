
VERSION=sse

PROG=decoder

WARN=-Wno-unused-variable -Wno-char-subscripts 

ARG=GALTESYAALGKSSWEEFNANIHTHRFFILVLEIAPAAKDLFSFLKGTSEVVSGPQNNPELQAHAGKKLVYEAIQLLVTGVTVTDATLKNLGSVHVSKGVADGHPPVVKEAILKTIEVVGAKLSEELNSAWTIAYDELAIVQKKEMNLAA

all:
	gcc -Wall $(WARN) -lm -O2 viterbi-$(VERSION).c -o $(PROG)

run: all
	./$(PROG) hmm-globins.txt $(ARG)

grind: all
	valgrind ./$(PROG) hmm-globins.txt $(ARG)

gdb: all
	gdb ./$(PROG) hmm-globins.txt $(ARG)

RESULTS='delete matchx insert'

diff:
	@diff -q -s correct-results-delete.txt results-delete.txt
	@diff -q -s correct-results-matchx.txt results-matchx.txt
	@diff -q -s correct-results-insert.txt results-insert.txt

test: run diff

clean:
	rm -rf  $(PROG)
	rm -rf results-delete.txt results-matchx.txt results-insert.txt

