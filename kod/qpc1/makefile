ifeq ($(OS),Windows_NT)
	PC := python
	RM := del /q
	VR := @echo nie ma valgrinda na windowsie
	RMDIR := rmdir /s /q
	MV := move
	DEVNULL := NUL 2>&1
	SLASH = $(strip \)
else
	ifdef VIRTUAL_ENV
		PC := python
	else
		PC := python3
	endif
	RM := rm -f
	VR = valgrind ./$(EXEC)
	RMDIR := rm -rf
	MV := mv
	DEVNULL := /dev/null
	SLASH := /
endif

.SUFFIXES:

CXX := g++
CFLAGS := -std=c++20 -Wall -Wextra -Werror -Wpedantic -I/usr/include/eigen3 -O3
LIBS := -lm -lgsl -lgslcblas
SRCS := main.cpp def.cpp
HDRS := dek.hpp
OBJS := main.o def.o
EXEC := main.out
DATA := dyspfile.csv kxfile.csv lfile.csv ufile.csv poprzecznefile.csv kx2vfile.csv misc.csv TRfile.csv
MEDIA := dyspfile.png lfile.png uRe.png uIm.png popU.png
SCRIPT := plot.py
TEXNAME := spr
TEX = $(TEXNAME).tex
REPORT := sprawozdanie.pdf
BUILDDIR := build

all: $(REPORT)

%.o: %.cpp $(HDRS)
	$(CXX) $(CFLAGS) -o $@ $< -c

$(EXEC): $(OBJS)
	$(CXX) $^ $(LIBS) -o $@

$(DATA): $(EXEC)
	./$<

$(MEDIA): $(SCRIPT) $(DATA)
	$(PC) $<

$(REPORT): $(TEX) $(MEDIA)
	latexmk -pdf -interaction=nonstopmode -output-directory=build $< > $(DEVNULL)
	$(MV) $(BUILDDIR)$(SLASH)$(TEXNAME).pdf $@ > $(DEVNULL)

gotowe: $(EXEC)
	./$<

safe: $(EXEC)
	$(VR)

clean:
	$(RM) $(OBJS)
	$(RM) $(EXEC)
	$(RM) $(DATA)
	$(RM) $(MEDIA)
	$(RMDIR) $(BUILDDIR)
	$(RM) $(REPORT)

.PHONY: clean gotowe all safe
