CXX=g++
RM=rm -f
CPPFLAGS=-O3 -Wall --std=c++11 -Isort/include/

SRCS=main.cpp

all: get_dependencies simulador

get_dependencies: get_dependencies.sh
	./get_dependencies.sh

simulador: main.cpp
	$(CXX) $(SRCS) $(CPPFLAGS) -o simulador

clean:
	$(RM) simulador