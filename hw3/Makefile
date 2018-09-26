CC=gcc
CXX=g++ -std=c++11
RM=rm -f

SRCS=CG_hw3.cpp Node.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

default: all

all: CG_hw3

CG_hw3: $(OBJS)
	$(CXX) -o CG_hw3 $(OBJS)

CG_hw3.o: CG_hw3.cpp

clean:
	$(RM) $(OBJS) CG_hw3

run: all
	./CG_hw3 -f patchPoints.txt -u 0.0303 -r 0.1 > out.iv && ivview out.iv