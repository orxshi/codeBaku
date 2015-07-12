CC = g++
CFLAGS = -c -std=c++11 -fpic -O3
LDFLAGS = -shared
SOURCES=AFT.cpp\
	Edge.cpp\
	Preset.cpp\
	Trimming.cpp\
	Triangle.cpp\
	Intersection.cpp\
	Point.cpp\
	ReBlanking.cpp\
	FModules.cpp\
	Front.cpp\
	NewGrid.cpp\
	FinalGrid.cpp	
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=libAFT.so

all: $(TARGET)
	
$(TARGET): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(PETSC_LIB) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(TARGET)
