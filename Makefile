.PHONY: default
default: all
CC=g++
#CXXFLAGS=-std=c++11 -Wall -Wextra -O0 -g # 编译选项
CXXFLAGS=-std=c++11 -Wall -Wextra -O3 # 编译选项
#被依赖的放后面，不能将libhts.a放在后面
LIBFLAGS=-I../third-party/htslib/htslib ../third-party/htslib/libhts.a -lz -lm -lbz2 -llzma -lcurl -lpthread

all: posmix.cpp dataframe.hpp
	$(CC) -o posmix posmix.cpp $(CXXFLAGS) $(LIBFLAGS)

