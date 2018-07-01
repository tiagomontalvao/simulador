CXX=g++
RM=rm -f
CPPFLAGS=-O3 -Wall -O3 -Isort/include/ -Itype_traits/include/ -Iconfig/include/ -Istatic_assert/include/ -Icore/include/ -Iserialization/include/ -Impl/include/ -Ipreprocessor/include/ -Irange/include/ -Iiterator/include/

SRCS=main.cpp

all: simulador

simulador: main.cpp
	$(CXX) $(SRCS) $(CPPFLAGS) -o simulador

clean:
	$(RM) main