CC=gcc
CXX=g++ -std=c++11
RM=rm -f

SRCS=CG_hw1.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

default: all

all: CG_hw1

CG_hw1: $(OBJS)
	$(CXX) -I /usr/include/Inventor -o CG_hw1 $(OBJS) -lInventor

CG_hw1.o: CG_hw1.cpp

clean:
	$(RM) $(OBJS) CG_hw1

run: all
	./CG_hw1 -f cpts_in.txt -u 0.05 -r 0.1 > out.iv
