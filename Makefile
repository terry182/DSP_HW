.PHONY: all clean run_train

CFLAGS+= -std=c99
LDFLAGS+=-lm     # link to math library

TARGET=test_hmm train

all: $(TARGET)
# type make/make all to compile test_hmm

test_hmm: test_hmm.c hmm.h
	cc $(CFLAGS) $(LDFLAGS) -o test_hmm test_hmm.c 

train: train.c hmm.h
	cc $(CFLAGS) $(LDFLAGS) -o train train.c

run_train: train
	./train 1 ../model_init.txt ../seq_model_01.txt heq.txt

clean:
	$(RM) $(TARGET)   # type make clean to remove the compiled file
	$(RM) *.o *.gch 
