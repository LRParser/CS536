CC=gcc
CXX=g++
RM=rm -f

SRCS=OITest.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: oitest 

oitest: $(OBJS)
	$(CXX) -I /usr/include/Inventor -o oitest $(OBJS) -lInventor

oitest.o: OITest.cpp

support.o: support.hh support.cc

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) oitest
