#address = $(CODE_BAKU_SRC)
address1 = ../../../../src/
address2 = ../../../../src/AFT
include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
#include ${PETSC_DIR}/lib/petsc/conf/test
INC1=${PETSC_DIR}/include
INC2=${PETSC_DIR}/$(PETSC_ARCH)/include
CC = g++
CFLAGS = -c -std=c++11 -O3 -I$(INC1) -I$(INC2)
LDFLAGS = -L$(address1) -Wl,-rpath=$(address1)
LDFLAGS += -L$(address2) -Wl,-rpath=$(address2)
SOURCES=main.cpp\
	Init/Init.cpp\
	Init/Read.cpp\
	Init/Output.cpp\
	MovingGrid/OscAirfoil.cpp\
	MovingGrid/StraightMovingAirfoil.cpp\
	Output/Output.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TARGET = out

all: $(SOURCES) $(TARGET)
	
$(TARGET): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(PETSC_LIB) -o $@ -lcb -lAFT

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(TARGET)
