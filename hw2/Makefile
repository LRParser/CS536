CC=gcc
CXX=g++ -std=c++11
RM=rm -f

SRCS=CG_hw2.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

default: all

all: CG_hw2

CG_hw2: $(OBJS)
	$(CXX) -I /usr/include/Inventor -o CG_hw2 $(OBJS) -lInventor

CG_hw2.o: CG_hw2.cpp

clean:
	$(RM) $(OBJS) CG_hw2

run: all
	./CG_hw2 -f cpts_in.txt -u 0.05 -r 0.1 > out.iv
