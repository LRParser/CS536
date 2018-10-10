CC=gcc
CXX=g++ -std=c++11
RM=rm -f

SRCS=CG_hw4.cpp Node.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

default: all

all: CG_hw4

CG_hw4: $(OBJS)
	$(CXX) -o CG_hw4 $(OBJS)

CG_hw4.o: CG_hw4.cpp

clean:
	$(RM) $(OBJS) CG_hw4

debug: all
	./CG_hw4 -f patchPoints.txt -u 2 -v 2 -r 0.1 -F -d

smalloutputonly: all
	./CG_hw4 -f patchPoints.txt -u 2 -v 2 -r 0.1 -F > out.iv

outputonly: all
	./CG_hw4 -f patchPoints.txt -u 9 -v 9 -r 0.1 -F > out.iv


run: all
	./CG_hw4 -S > out.iv && ivview out.iv

run1: all
	./CG_hw4 -A 1.2 -B 1.1 -C 2.3 -r 0.09 -t 1.0 -u 25 -v 20 -S -d > out.iv && ivview out.iv

run3: all
	./CG_hw4 -A 1.9 -B 2.1 -C 0.5 -r 1.9 -t 1.6 -u 23 -v 12 -S > out.iv && ivview out.iv

runsmooth: all
	./CG_hw4 -S -f patchPoints.txt -u 9 -v 9 -r 0.1 > out.iv && ivview out.iv

runold: all
	./CG_hw4 -f patchPoints.txt -u 2 -v 2 -r 0.1 -F > out.iv && ivview out.iv
