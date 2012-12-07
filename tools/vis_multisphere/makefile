CXX = g++
CXXFLAGS = $(SVNDEV)

OBJ = vis_multisphere.o molecule.o quaternion.o component.o
OBJD = $(OBJ:%.o=%d.o)

all: vis_multisphere vis_multisphered

vis_multisphere: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o vis_multisphere

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -Wall -O -c $<

vis_multisphered: $(OBJD)
	$(CXX) $(CXXFLAGS) $(OBJD) -o vis_multisphered

%d.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -Wall -g -D DEBUG -c $< -o $*d.o


clean:
	rm *.o
	
clear: clean
	rm vis_multisphere vis_multisphered

tar:
	tar -cvjf 'vis_multisphere.tar.bz2' *.h *.cpp makefile
