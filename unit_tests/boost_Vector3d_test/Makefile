CXX = g++

CXXFLAGS = -g -O0 -I/usr/local/boost_1_55_0

OBJ = boost_Vector3d_test.o ~/my-ls1-branch/MatthiasHeinen/src/utils/Vector3d.o
OBJD = boost_Vector3d_test.do ~/my-ls1-branch/MatthiasHeinen/src/utils/Vector3d.do

HEADER = $(OBJ:%.o=%.h)
CPP = $(OBJ:%.o=%.cpp)
DEPS = $(OBJ:%.o=%.dep)

all: boost_Vector3d_test

boost_Vector3d_test: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o boost_Vector3d_test

boost_Vector3d_testd: $(OBJD)
	$(CXX) $(CXXFLAGS) $(OBJD) -o boost_Vector3d_testd

$(OBJD):
	$(CXX) $(CXXFLAGS) -Wall -g -D DEBUG -c $< -o $@

%.dep: %.cpp %.h
	$(CXX) -MT $*.o -MM $< > $@
	$(CXX) -MT $*.do -MM $< >> $@

clean:
	rm -f $(OBJ) $(OBJD)
        
clear: clean
	rm -f boost_Vector3d_test boost_Vector3d_testd
