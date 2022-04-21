CXX = g++-11
CPPFLAGS = -std=c++20
CXXFLAGS = -O3
TARGET = main

all: $(TARGET)

# compiling main binary 
main: main.o
	@mkdir -p data
	$(CXX) $+ -o $@

# linking object files
.INTERMEDIATE: main.o
main.o: main.cpp
	$(CXX) -c $+ -o $@

.PHONY: clean
clean:
	$(RM) $(TARGET) main.o