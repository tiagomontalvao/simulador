CXX=g++
RM=rm -f
CPPFLAGS=-O3 -Wall -Isort/include/

SRCS=main.cpp

all: simulador

simulador: main.cpp
	$(CXX) $(SRCS) $(CPPFLAGS) -o simulador

clean:
	$(RM) simulador