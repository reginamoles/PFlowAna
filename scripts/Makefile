ROOTLIBS = $(shell root-config --libs)
TMVALIBS = -lTreePlayer -lMLP -lMinuit -lThread
LIBS     = $(ROOTLIBS) $(SYSLIBS) -lMinuit -lRooFit

IDIR      = -Isrc -Iinc
ROOTFLAGS = $(shell root-config --cflags)
CXXFLAGS  = -O -Wall -I $(ROOTSYS)/include $(IDIR) $(ROOTFLAGS) -g
LDFLAGS   = $(shell root-config --ldflags)

SRC = ./src/
BIN = ./bin/
EXE = ./bin/

default : obj exe

obj :   $(BIN)PFlowMonitor.o\

exe : 	$(EXE)main_x

$(BIN)%.o : $(SRC)%.cxx
	g++ $(CXXFLAGS) -fPIC -o $@ -c $?

$(LIB)%.so :
	g++ -shared -o $@ $(BIN)*.o

$(EXE)%_x : $(SRC)%.cxx
	g++ -Wno-write-strings $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(BIN)*.o -o $@ $?

clean :   
	rm $(BIN)*.o
	rm $(EXE)*_x
