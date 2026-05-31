ifeq ($(OS),Windows_NT)
	PC := python
	RM := del /q
	VR := @echo nie ma valgrinda na windowsie
else
	ifdef VIRTUAL_ENV
		PC := python
	else
		PC := python3
	endif
	RM := rm -f
	VR = valgrind ./$(EXEC)
endif

.SUFFIXES:

CXX := g++
CFLAGS := -std=c++20 -Wall -Wextra -Werror -Wpedantic -I/usr/include/eigen3
LIBS := -lm -lgsl -lgslcblas
SRCS := main.cpp def.cpp
HDRS := dek.hpp
OBJS := main.o def.o
EXEC := main.out
DATA := dyspfile.csv kxfile.csv
MEDIA := dyspfile.png
SCRIPT := plot.py

all: $(MEDIA)

%.o: %.cpp $(HDRS)
	$(CXX) $(CFLAGS) -o $@ $< -c

$(EXEC): $(OBJS)
	$(CXX) $^ $(LIBS) -o $@

$(DATA): $(EXEC)
	./$<

$(MEDIA): $(SCRIPT) $(DATA)
	$(PC) $<

gotowe: $(EXEC)
	./$<

safe: $(EXEC)
	$(VR)

clean:
	$(RM) $(OBJS)
	$(RM) $(EXEC)
	$(RM) $(DATA)
	$(RM) $(MEDIA)

.PHONY: clean gotowe all safe
