.SUFFIXES:

CXX := g++
TARGET = main.out
SRCS = main.cpp def.cpp
OBJS = main.o def.o
HDRS = dek.hpp
CFLAGS := -std=c++20 -Wall -Wextra -Werror -Wpedantic
LIBS := -lm

all: gotowe

%.o: %.cpp $(HDRS)
	$(CXX) $(CFLAGS) -o $@ $< -c

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@

gotowe: $(TARGET)
	./$<

safe: $(TARGET)
	valgrind ./$<

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)

.PHONY: clean gotowe all safe
