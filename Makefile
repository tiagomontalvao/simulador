CXX=g++
RM=rm -f
override CPPFLAGS+=-Wall -std=c++11 -Imath/include/

SRCS=simulador.cpp
OUT=simulador

all: CPPFLAGS += -D NDEBUG -O3
all: simulador

simulador: simulador.cpp
	$(CXX) $(SRCS) $(CPPFLAGS) -o $(OUT)

debug: CPPFLAGS += -g -Og
debug: OUT=simulador_dbg
debug: simulador_dbg

simulador_dbg: simulador.cpp
	$(CXX) $(SRCS) $(CPPFLAGS) -o $(OUT)

clean:
	$(RM) simulador simulador_dbg