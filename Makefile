CXX=g++
RM=rm -f
override CPPFLAGS+=-Wall -std=c++11 -Imath/include/ -Iconfig/include/ -Ipredef/include/ -Itype_traits/include/ -Istatic_assert/include/ -Impl/include/ -Ipreprocessor/include/ -Iassert/include/ -Icore/include/ -Ithrow_exception/include/ -Ilexical_cast/include/ -Irange/include/ -Iiterator/include/ -Idetail/include/ -Iconcept_check/include/ -Iutility/include/ -Inumeric_conversion/include/ -Iinteger/include/ -Iarray/include/ -Icontainer/include/ -Imove/include/ -Ismart_ptr/include/

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