#ifndef CELL_H_
#define CELL_H_


#include <list>
#include <vector>
using namespace std;

namespace datastructures {
  template<class ParticleType>
  class Cell; 
}

//! @brief Cell data structure.
//!
//! A Cell stores a list of pointers to the molecules in that cell. \n
//! It also knows to which region it belongs. The following regions are possible: \n
//! \li completele outside (more than the cutoffradius away from all cells that
//!                       belong directly to the MoleculeContainer)
//!                      Such a cell shouldn't exist!!!
//! \li halo region (not belonging directly to the MoleculeContainer, but within
//!                       the cutoffradius of at least one cell of the MoleculeContainer)
//! \li boundary region (belonging directly to the MoleculeContainer, but not more than
//!                           the cutoffradius away from the boundary)
//! \li inner region (more than the cutoffradius away from the boundary)
//! 
//! There are three boolean member variables for the last three regions. \n
//! If more than one of them is true, there must be an error in the code \n
//! If none of them is true, the cell wasn't assigned to any region yet.
//! A cell which is completely outside shouldn't exist, as it completely useless.
//! 
//! For cuboid domains, a cell doesn't need much information
//! (you could even avoid a seperate class completely)
//! Most information (position, neighbours,...) can be calculated
//! from the cell's index in the linked cell vector.
//! e.g. for a domain of 10x10x10 cells, the cell with index (counting from zero)
//! 561 has spacial index (counting from zero)  (1,6,5) and the cell indices of
//! the neighbours are e.g. 450, 451, 452, 460, 461, ... 
//! But for either unregular domains or unregular spacial decompositions
//! in the parallelisation (space filling curves), it's not easy to store that
//! information in the linked cell vector (if it is a vector at all)
//! That's why this class was introduced, so additional information which is
//! specific to a cell should be stored here.
//! @todo Improve class comment
template<class ParticleType>
class datastructures::Cell{
  public:
    //! The constructor sets all cell-states to false 
    Cell();
  
    //! removes all elements from the list molecules
    void removeAllParticles();
        
    //! insert a single molecule into this cell
    void addParticle(ParticleType* particle_ptr);
      
    //! return a reference to the list of molecules (molecule pointers) in this cell
    list<ParticleType*>& getParticlePointers();

#ifdef GRANDCANONICAL
    bool deleteMolecule(unsigned long molid);
#endif
    
    //! Set the flag for a Halo Cell
    void assignCellToHaloRegion();
    
    //! Set the flag for a Boundary Cell
    void assignCellToBoundaryRegion();
    
    //! Set the flag for a Inner Cell
    void assignCellToInnerRegion();

    //! returns true, if the cell is a Halo Cell, otherwise false
    bool isHaloCell();
    
    //! returns true, if the cell is a Boundary Cell, otherwise false
    bool isBoundaryCell();
    
    //! returns true, if the cell is a Inner Cell, otherwise false
    bool isInnerCell();    

    
  private:
    //! each cell contains a list of pointers to the molecules in the cell
    list<ParticleType*> _particlePointers; 
    
    //! true when the cell is in the halo region 
    bool _haloCellState;
    //! true when the cell is in the boundary region
    bool _boundaryCellState;
    //! true when the cell is in the inner region
    bool _innerCellState;
};

#include "datastructures/Cell.cpph"

#endif /*CELL_H_*/
