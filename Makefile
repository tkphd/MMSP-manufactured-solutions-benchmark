# Makefile
# GNU makefile for PFHub Benchmark 7: Method of Manufactured Solutions
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

compiler = /usr/bin/mpic++
flags = -O3 -I $(MMSP_PATH)/include -I /usr/include
deps = manufactured.c
links = -lz -L /usr/lib/x86_64-linux-gnu -lgsl -lgslcblas

all: allen-cahn $(deps)
.PHONY: all clean

allen-cahn: allen-cahn.cpp main.cpp $(deps)
	$(compiler) $(flags) $< $(deps) -o $@ $(links)

manufactured.c: manufactured-sympy.py
	python $<

clean:
	rm -vf allen-cahn manufactured.c manufactured.h
