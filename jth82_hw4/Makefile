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
	./CG_hw4 > out.iv && ivview out.iv

runsmooth: all
	./CG_hw4 -S -f patchPoints.txt -u 9 -v 9 -r 0.1 > out.iv && ivview out.iv

runold: all
	./CG_hw4 -f patchPoints.txt -u 2 -v 2 -r 0.1 -F > out.iv && ivview out.iv