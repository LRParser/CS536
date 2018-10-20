CC=gcc
CXX=g++ -std=c++11
RM=rm -f

SRCS=CG_hw2.cpp Node.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

default: all

all: CG_hw2

CG_hw2: $(OBJS)
	$(CXX) -o CG_hw2 $(OBJS)

CG_hw2.o: CG_hw2.cpp

clean:
	$(RM) $(OBJS) CG_hw2

debug: all
	./CG_hw2 -f cpts_test2.txt -u 0.0303 -r 0.1

run0: all
	./CG_hw2 -f HW2_test0.txt -u 0.025 -r 0.05 -t 0 -d > out.iv && ivview out.iv

run0a: all
	./CG_hw2 -f HW2_test0.txt -u 0.025 -r 0.05 -t 0.5 -d > out.iv && ivview out.iv

run0b: all
	./CG_hw2 -f HW2_test0.txt -u 0.025 -r 0.05 -t -0.5 -d > out.iv && ivview out.iv


run: all
	./CG_hw2 -f cpts_test3.txt -u 0.0303 -r 0.1 > out.iv && ivview out.iv

run: all
	./CG_hw2 -f cpts_test3.txt -u 0.0303 -r 0.1 > out.iv && ivview out.iv