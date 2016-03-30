.PHONY: all clean

CFLAGS+=
LDFLAGS+=-lm     # link to math library

TARGET=test_hmm

all: $(TARGET)
# type make/make all to compile test_hmm

test_hmm: test_hmm.c hmm.h
	gcc -o test_hmm test_hmm.c 

clean:
	$(RM) $(TARGET)   # type make clean to remove the compiled file