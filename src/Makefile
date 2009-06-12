CODE = $(shell find ./ -name "*.cpp")

# only linker and domain decomposition need MPI
all: $(CODE)
	./mkmardyn

parallel/DomainDecomposition.o: parallel/DomainDecomposition.cpp
	$(CXX_PAR) $(CXXFLAGS) -c $< -o $@

tar: clean
	tar -cvjf mardyn_branch.tar.bz2 *.h *.cpp DoxySfile units mkmardyn Makefil* datastructures ensemble integrators md_io molecules parallel plugin tests utils 

clean:
	find -name '*.o' -exec rm -r {} \;
	find -name '*.bz2' -exec rm -r {} \;
	find -name '*~' -exec rm -r {} \;

clear: clean
	rm -f ls1*
