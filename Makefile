CXX=g++
RM=rm -f
CPPFLAGS=-Wall --std=c++11

SRCS=simulador.cpp

all: CPPFLAGS += -D NDEBUG -O3
all: simulador

debug: CPPFLAGS += -g -Og
debug: simulador

simulador: simulador.cpp
	$(CXX) $(SRCS) $(CPPFLAGS) -o simulador

clean:
	$(RM) simulador