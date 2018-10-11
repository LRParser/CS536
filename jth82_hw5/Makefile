CC=gcc
CXX=g++ -std=c++11
RM=rm -f

SRCS=CG_hw5.cpp Node.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

default: all

all: CG_hw5

CG_hw5: $(OBJS)
	$(CXX) -o CG_hw5 $(OBJS)

CG_hw5.o: CG_hw5.cpp

clean:
	$(RM) $(OBJS) CG_hw5

debug: all
	./CG_hw5 -f patchPoints.txt -u 2 -v 2 -r 0.1 -F -d

run: all
	./CG_hw5 > out.iv && ivview out.iv

