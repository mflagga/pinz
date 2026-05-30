.SUFFIXES:

CXX := g++
CFLAGS := -std=c++20 -Wall -Wextra -Werror -Wpedantic
LIBS := -lm
SRCS := main.cpp def.cpp
HDRS := dek.hpp
OBJS := main.o def.o
TARGET := main.out

all: gotowe

%.o: %.cpp $(HDRS)
	$(CXX) $(CFLAGS) -o $@ $< -c

$(TARGET): $(OBJS)
	$(CXX) $^ $(LIBS) -o $@

gotowe: $(TARGET)
	./$<

safe: $(TARGET)
	valgrind ./$<

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)

.PHONY: clean gotowe all safe
