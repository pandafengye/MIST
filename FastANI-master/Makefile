CXXFLAGS += -O3 -DNDEBUG -std=c++11 -Isrc -I /usr//include -fopenmp
CPPFLAGS += 

UNAME_S=$(shell uname -s)

ifeq ($(UNAME_S),Darwin)  #macOS clang
	CXXFLAGS += -mmacosx-version-min=10.7 -stdlib=libc++ -Xpreprocessor -fopenmp -lomp
else
	CXXFLAGS += -include src/common/memcpyLink.h -Wl,--wrap=memcpy
	CFLAGS += -include src/common/memcpyLink.h
endif

SOURCES=src/cgi/core_genome_identity.cpp

all : fastANI

fastANI :  
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(SOURCES) -o fastANI -L/usr//lib -lgsl -lgslcblas -lstdc++ -lz -lm  

install : fastANI
	mkdir -p /usr/local/bin/
	cp `pwd`/fastANI /usr/local/bin/

clean :
	-rm -f fastANI
