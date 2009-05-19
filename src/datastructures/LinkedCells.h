#ifndef LINKEDCELLS_H_
#define LINKEDCELLS_H_

#include "datastructures/ParticleContainer.h"
#include "utils/Log.h"
//#include "Cell.h"
#include <sstream>

namespace datastructures {
  template<class ParticleType>
  class LinkedCells; 
  
  template<class ParticleType>
  class Cell;
}

//! @brief Linked Cell Data Structure
//!
//! Without any specialized data structure, it needs O(N*N) - where N is
//! the number of particles - time to find all neighbouring pairs of particles.
//! The linked cell data structure is a datastructure which allows to find all
//! neighbouring pairs of particles (neighbouring means particles pairs
//! which have less than a certain distance) in O(N) time.
//! The following picture shows a domain with some particles in it. The blue
//! circle shows the neighbouring area of the red particle
//! \image particles.jpg
//! The problem is that all particles have to be examined to find those within the circle
//! With the linked cell data structure, the domain is divided into cells (using a regular grid).
//! All particles are placed in those cells.
//! For a given cell, neighbouring cells can easily be calculated, so for a given particle,
//! only the particles from neighbouring cells have to be examined. The following
//! picture illustrates this
//! \image particles_with_lc.jpg
//!
//! The spatial domain covered by the linked cells is larger than
//! the bounding box of the (sub)domain. This halo region surrounding
//! the control volume is used for (periodic) boundary conditions
//! and has to be at least as wide as the cutoff radius. \n
//! In total, there are three different cell types:
//! - halo 
//! - boundary
//! - inner
template<class ParticleType>
class datastructures::LinkedCells: public datastructures::ParticleContainer<ParticleType> {
  public:
    //#########################################################################
    //############# methods common for all ParticleContainers #################
    //#########################################################################

    //! @brief initialize the Linked Cell datastructure
    //!  
    //! The constructor sets the following variables:
    //! - cellsPerDimension[3]
    //! - haloWidth[3]
    //! - cellLength[3]
    //! - cellsInCutoffRadius
    //! - lowerCorner and higherCorner
    //!
    //! It resized the cell vector and marks the cells as inner/halo \n
    //! It fills the array innerCellIndices \n
    //! It fills the array with forward and backward neighbour indices \n
    //! The corner parameters for the constructor describe the bounding box
    //! of the phasespace which belongs directly to this process, so they correspond
    //! to a bounding box including inner + boundary cells but excluding halo cells. \n
    //! But the corners of this class have to include the halo cells.
    //! @param bBoxMin lower corner of the bounding box of the domain belonging to this container
    //! @param bBoxMax higher corner of the bounding box of the domain belonging to this container 
    //! @param cutoffRadius distance for which forces have to be calculated
    //! @param cellsInCutoffRadius describes the width of cells relative to the cutoffRadius: \n
    //!        equal (or larger) to the cutoffRadius divided by the length of a cell
    //!        as for the number of cells in each dimension only natural numbers are allowed,
    //!        it can happen that it is not possible to set celllength = cutoffRadius / cellsInCutoffRadius.
    //!        In that case, the celllength is chosen to be the next larger value so that the sum of
    //!        the cell lengths in one dimension equals the length of the phasespace
    //!        Example: phasespacelength=100, cellsInCutoffRadius=2, CutoffRadius=3 \n
    //!        ==> celllength should be: cutoffRadius/cellsInCutoffRadius = 3/2 = 1.5 \n
    //!        ==> cellsPerDimension = phasespacelength/celllength = 100/1.5 = 66.67 cells \n
    //!        ==> cells have to be larger: cellsPerDimension = phasespacelength/celllength = 100/celllength = 66 cells \n
    //!        ==> celllength = 100/66 = 1.5152
    //! @param partPairsHandler specified concrete action to be done for each pair
    LinkedCells(double bBoxMin[3], double bBoxMax[3], double cutoffRadius,
#ifdef COMPLEX_POTENTIAL_SET
                double tersoffCutoffRadius,
#endif
                double cellsInCutoffRadius, datastructures::ParticlePairsHandler<ParticleType>& partPairsHandler);
    
    //! Destructor
    ~LinkedCells();
 
    //! Pointers to the particles are put into cells depending on the spacial position
    //! of the particles. 
    //! Before the call of this method, this distribution might have become invalid.
    //! To ensure, that all Particles (pointers to them) are put into the corresponding cells,
    //! first all cells are cleared and then filled again depending on the spacial position
    //! of the molecules. After the update, exactly one pointer for each particle in this
    //! ParticleContainer is it's corresponding cell.
    void update();    
    
    //! @brief Insert Molecules into the cells.
    //!
    //! Loop over all Molecules and call this->addMolecule(...) for each of them
    //void addMolecules(std::list<ParticleType>& molecules);
    
    //! @brief Insert a single molecule.
    //!
    //! Therefore, first the cell (the index) for the molecule has to be determined,
    //! then the molecule is inserted into that cell.
    // void addMolecule(ParticleType* molecule_ptr);
    // void addParticle(ParticleType* particle);
    void addParticle(ParticleType& particle);
    
    //! @brief calculate the forces between the molecules.
    //!
    //! Only molecules with a distance not larger than the cutoff radius are to be used. \n
    //! Only forces on the Molecules which are in the inner and boundary region have to be calculated
    //! Newton's third law should be used for faster computation:
    //! \li a loop over all inner cells calculates the forces with all forward neighbour cells
    //!     all forward cells have to be used, as none of them can be halo or outside
    //! \li a loop over the boundary cells first calculates forces with all forward cells and all
    //!     backward cells. Here it has to be checked whether the neighbour cell is halo or not.
    //!     If it is Halo, the force is calculated, if it isn't, the force is not calculated, 
    //!     because the same pair of cells has already been processed in one of the other loops. 
    void traversePairs();

    //#########################################################################
    //############# special methods for LinkedCells ###########################
    //#########################################################################
    
    //! @brief Initialze index vectors and cells. 
    //!
    //! Fill the vector with the indices of the inner and boundary cells.
    //! Assign each cell it's region (halo, boundary, inner).
    void initializeCells();
        
    //! @brief Calculate neighbour indices.
    //!
    //! This method is executed once for the molecule container and not for
    //! each cell. E.g. the index (in the cell vector) of the right neighbour of a cell
    //! always equals the index of the cell minus one. This method calculates two vectors
    //! of index offsets, one for positive offsets (forward neighbours) and one for negative 
    //! offsets (backward neighbours). So given a specific cell, the neighbours can be retrieved
    //! by adding to the index of the cell the offsets in the two vectors. 
    //! 
    //! The method works as follows: \n
    //! The loop runs over all potential neighbour cells (bounding box which contains
    //! the cell itself, and in each dimension on the lower and on the higher side as
    //! many cells as the width of the halo strip. E.g. if the haloWidth is 2, a box 
    //! of 5x5x5 cell is considered as potential neighbours
    //! for each of those cells, the minimal possible distance between that cell
    //! and the central cell is calculated (sqrt(x^2+y^2+z^2)). If that distance
    //! is larger than the cutoff radius, the cell can be neglected.
    //! The distance in one dimension is the width of a cell multiplied with the number 
    //! of cells between the two cells (this is received by substracting one of the difference). 
    void calculateNeighbourIndices(); 

    //! @brief Get the index in the cell vector to which this Molecule belong
    //!
    //! each spacial position within the bounding box of the linked cells
    //! belongs unambiguously to one cell. \n
    //! This method determines for a given Molecule the corresponding cell
    //! and returns the index of that cell in the cell vector. \n
    //! If the molecule is not inside the bounding box, an error is printed
    unsigned long getCellIndexOfMolecule(ParticleType* molecule);
    
    //! @brief given the 3D index of a cell, return the index in the cell vector.
    //!
    //! A cell can be identified by a 3D index. \n
    //! This method determines for a given 3D index the corresponding cell
    //! and returns the index of that cell in the cell vector. \n
    //! The method can also be used to get the offset between two cells in the cell
    //! vector when called with the 3D cell index offets (e.g. x: one cell to the left,
    //! y: two cells back, z: one cell up,...)
    unsigned long cellIndexOf3DIndex(int xIndex, int yIndex, int zIndex);

    //! @brief gets the width of the halo region in dimension index
    double get_halo_L(int index);

    //! @return the number of particles stored in the Linked Cells
    unsigned long getNumberOfParticles();
    
    //! @brief returns a pointer to the first particle in the Linked Cells
    //!
    //! Internally, the particles are store in a std::list. To traverse this
    //! list, a iterator (_particleIter) for the list is used.
    //! This method sets this iterator to point to the begin of the list
    //! and return a pointer to the value pointed to by the iterator
    ParticleType* begin();

    //! @brief returns a pointer to the next particle in the Linked Cells
    //!
    //! The iterator _particleIter is first incremented. Then a pointer
    //! to the value pointed to by the iterator is returned. If the
    //! iterator points to the end of the list (which is one element after the last
    //! element), NULL is returned 
    ParticleType* next();

    //! @brief returns NULL
    ParticleType* end();
        
    //! @brief delete all Particles which are not within the bounding box
    void deleteOuterParticles();

    double getCutoff() { return this->_cutoffRadius; }
#ifdef COMPLEX_POTENTIAL_SET
    double getTersoffCutoff() { return this->_tersoffCutoffRadius; }
#endif
    void countParticles(Domain* d);
    //! @brief counts all particles inside the bounding box
    unsigned countParticles(int cid);
    //! @brief counts particles in the intersection of bounding box and control volume
    unsigned countParticles(int cid, double* cbottom, double* ctop);

#ifdef GRANDCANONICAL
    void deleteMolecule(unsigned long molid, double x, double y, double z);
    double getEnergy(ParticleType* m1);
    int localGrandcanonicalBalance() { return this->_localInsertionsMinusDeletions; }
    int grandcanonicalBalance(parallel::DomainDecompBase* comm);
    void grandcanonicalStep(ensemble::ChemicalPotential* mu, double T);
#endif
    
  private:
    //! Logging interface
    static utils::Log _log;

    // datastructures::ParticlePairsHandler<ParticleType>& _particlePairsHandler;

    //! the list contains all molecules from the system
    list<ParticleType> _particles; // CHECKED
    
    //! Iterator to traverse the list of particles (_particles)
    typename std::list<ParticleType>::iterator _particleIter;
    
    //! Vector containing all cells (including halo)
    std::vector<Cell<ParticleType> > _cells;
    
    //! Vector containing the indices (for the cells vector) of all inner cells (without boundary)
    std::vector<unsigned long> _innerCellIndices;
    
    //! Vector containing the indices (for the cells vector) of all boundary cells
    std::vector<unsigned long> _boundaryCellIndices;
    
    //! Neighbours that come in the total ordering after a cell
    std::vector<unsigned long> _forwardNeighbourOffsets;
    
    //! Neighbours that come in the total ordering before a cell
    std::vector<unsigned long> _backwardNeighbourOffsets;

    //! low corner of the bounding box around the linked cells (including halo)
    double _haloBoundingBoxMin[3];
    
    //! high corner of the bounding box around the linked cells (including halo)
    double _haloBoundingBoxMax[3];  
    
    //! Number of Cells in each spatial dimension (including halo)
    int _cellsPerDimension[3];
    //! Halo width (in cells) in each dimension
    int _haloWidthInNumCells[3];
    //! width of the halo strip (in size units)
    double _haloLength[3];
    //! length of the cell (for each dimension)
    double _cellLength[3];
    //! cutoff radius
    double _cutoffRadius;
#ifdef COMPLEX_POTENTIAL_SET
    //! Tersoff cutoff radius
    double _tersoffCutoffRadius;
#endif
    
    //! Number of neighbour cells in one direction.
    double _cellsInCutoffRadius;

#ifdef GRANDCANONICAL
    int _localInsertionsMinusDeletions;
#endif
};

#include "datastructures/LinkedCells.cpph"

#endif /*LINKEDCELLS_H_*/
