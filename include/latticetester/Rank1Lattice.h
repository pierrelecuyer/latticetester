// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the occasional supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LATTICETESTER_RANK1LATTICE_H
#define LATTICETESTER_RANK1LATTICE_H

#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/IntLatticeExt.h"

namespace LatticeTester {

  /**
   * This subclass of `IntLatticeExt` defines a general rank 1 lattice rule in \f$d\f$ dimensions,
   * whose points \f$\mathbb{u}_i\f$ are defined by
   * \f{equation}{
   *    \mathbf{u}_i = (i \mathbf{a} \mod m)/m,
   * \f}
   * where \f$\mathbf{a} \in \mathbb{Z}_m^d\f$ is called the generating vector.
   * The lattice is rescaled simply by removing the division by \f$m\f$.
   * Given \f$\mathbf{a}\f$, a basis of the rescaled (integer) lattice is given by
   * \f{align*}{
   *    \mathbf{v}_1 & = a \\
   *    \mathbf{v}_2 & = m \mathbf{e}_2 \\
   *    \vdots & \\
   *    \mathbf{v}_d & = m \mathbf{e}_v \\
   * \f}
   * where \f$\mathbf{e}_i\f$ is the \f$i^\text{th}\f$ unit vector.
   * 
   * A condition that is often required when building a rank 1 lattice is that
   * \f$\gcd(a_i, m) = 1,\ \forall 1\leq i \leq d\f$. When this condition is
   * verified, each lower-dimensional projection of the lattice contains the same number
   * of points (has the same density \f$m\f$) as the full lattice.
   * When searching for lattices that satisfy this condition, one may assume
   * without loss of generality generality that \f$a_1 = 1\f$.
   *
   * ***  CHANGES: added parameter withDual to the constructor.
   *   Removed  init();
   *   I changed buildBasis for the dual to use the direct construction.
   *   incDimBasis needs to be made more efficient, by just updating the basis!
   *   It seems that setLac could be in IntLattice already.
   *
   */
template<typename Int, typename Real>
class Rank1Lattice: public IntLatticeExt<Int, Real> {

      private:
        typedef NTL::vector<Int>  IntVec;
        typedef NTL::matrix<Int>  IntMat;
        typedef NTL::vector<Real> RealVec;

      public:

        /**
         * This constructor takes as input the modulus `m`, the generating vector `aa`,
         * the (maximal) dimension `maxDim`, and the norm used to measure the vector lengths.
         * The length of the vector `aa` should be `maxDim`.
         * This constructor does not build the basis, to leave
         * more flexibility in the dimension when doing so.
         */
        Rank1Lattice (const Int & m, const IntVec & aa, int64_t maxDim,
        	bool withDual=false, NormType norm = L2NORM);

        /**
         * Constructor for the special case of a Korobov lattice.
         * Here the generating vector has the form aa = (1, a, a^2 mod m, a^3 mod m, ...)
         * where a is an integer such that 1 < a < m.
         */
        Rank1Lattice (const Int & m, const Int & a, int64_t maxDim,
        	bool withDual=false, NormType norm = L2NORM);

        /**
         * Copy constructor.
         */
        Rank1Lattice (const Rank1Lattice<Int, Real> & Lat);

        /**
         * Assigns `Lat` to this object.
        */
        Rank1Lattice & operator= (const Rank1Lattice<Int, Real> & Lat);

        /**
         * Destructor.
         */
        ~Rank1Lattice();

        /**
         * Returns the first components of the generating vector \f$\ba\f$ as a string.
         * The number of components in the string will be the current dimension of the lattice.
         */
        std::string toStringCoef() const;

        /**
         * Builds a basis in `d` dimensions. This `d` must not exceed `this->maxDim()`.
         * This initial basis will be upper triangular.
         */
        void buildBasis (int64_t d);
        
        /**
         * Builds only the basis in `d` dimensions. This `d` must not exceed `this->maxDim()`.
         * This dual basis will be lower triangular.
         */
        void buildDualBasis (int64_t d);

        /**
         * Increases the current dimension by 1 and updates the basis.
         * The dimension must be smaller than `maxDim` when calling this function.
         */
        void incDimBasis ();
        
        /**
         * Increases the current dimension of only the dual lattice by 1 and updates only the dual basis.
         * The dimension must be smaller than `maxDim` when calling this function.
         */
        void incDimDualBasis ();

      protected:

        /**
         * Initializes the rank 1 lattice. This just invokes `IntLatticeExt::initProj()`.
         *        So why do we need it ???  check this ....   ?????         *****
         */
        // void initProj();

        /**
         * Vector of multipliers (generating vector) of the rank 1 lattice rule.
         * They are stored for up to `maxDim()` dimensions.
         * The first coordinate has index 0.
         */
        IntVec m_a;
        
        /**
         * This number stores the length of the input vector aa. This is necessary to make
         * an efficient implementation of 'incDimBasis'
         */
        int m_k = 1;
    };


//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice (
         const Int & m, const IntVec & aa, int64_t maxDim, bool withDual, NormType norm):
         IntLatticeExt<Int, Real> (m, maxDim, withDual, norm) {
    this->m_a = aa;
    //this->initProj();
  }

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice (
        const Int & m, const Int & a, int64_t maxDim, bool withDual, NormType norm):
        IntLatticeExt<Int, Real> (m, maxDim, withDual, norm) {
    m_a.SetLength(maxDim);
	Int powa(1);  m_a[0] = powa;
    for (int64_t i=1; i < maxDim; i++) {
    	//powa = (a * powa) % m;
    	m_a[i] = a * m_a[i-1] % m;
    }
    //this->initProj();    // Does not initialize a basis....
  }

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::~Rank1Lattice() {
      this->m_a.kill ();
      // ~(this->m_a);
    }


  //============================================================================

/*
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::initProj() {
    IntLatticeExt<Int, Real>::initProj();
      // for (int64_t r = 1; r < this->getDim(); r++)
      //   this->m_lgVolDual2[r] = this->m_lgVolDual2[r - 1];
    }
*/

  //============================================================================

// Essenti
template<typename Int, typename Real>
Rank1Lattice<Int, Real> & Rank1Lattice<Int, Real>::operator= (
         const Rank1Lattice<Int, Real> & lat) {
      if (this == &lat)
         return *this;
      this->copy (lat);
      //this->initProj();
      this->m_a = lat.m_a;
      return *this;
    }

  //============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice (
        const Rank1Lattice<Int, Real> & lat):
      IntLatticeExt<Int, Real> (
          lat.m_modulo, lat.getDim (), lat.m_withDual, lat.getNormType ()) {
    // MyExit (1, "Rank1Lattice:: constructor is incomplete" );
    //this->initProj();
    this->m_a = lat.m_a;
  }

  //============================================================================

template<typename Int, typename Real>
std::string Rank1Lattice<Int, Real>::toStringCoef ()const {
      return toString (this->m_a, 0, this->getDim ());
    }

  //============================================================================

//template<typename Int, typename Real>
//void Rank1Lattice<Int, Real>::incDimBasis () {
//    	assert(1 + this->getDim() <= this->m_maxDim);
//    	buildBasis (1 + this->getDim ());
//    	this->setNegativeNorm ();
//    	this->setDualNegativeNorm ();
//
//	}

//template<typename Int, typename Real>
//void Rank1Lattice<Int, Real>::incDimBasis () {
//		int64_t d = 1 + this->getDim();
//    	assert(d <= this->m_maxDim);
//    	IntMat temp;
//    	temp.SetDims(d, d);
//    	//Use old basis for first d - 1 dimension
//    	for (int i = 0; i < d-1; i++) {
//            for (int j = 0; j < d-1; j++) {
//               temp[i][j] = this->m_basis[i][j];	
//            }
//        }
//    	// Code for general k
////        //Add d-th component for the first d - 1 vectors 
////        //Recall that indices start with 0!
////    	for (int i = 0; i < d-1; i++) {
////    		temp[i][d-1] = 0;
////    		for (int k = 1; k < this->m_k + 1; k++) {
////    	       temp[i][d-1] = temp[i][d-1] + m_a[k] * temp[i][d-1-k];
////            }
//// 	       temp[i][d-1] = temp[i][d-1] % this->m_modulo;
////    	}
////    	//Set last vector equal m times unit vector
////    	for (int j = 0; j < d-1; j++)
////    		temp[d-1][j] = 0;
////    	temp[d-1][d-1] = this->m_modulo;
//    	
//    	//Fill in the new component
//    	for (int j = 0; j < d-1; j++)
//    		temp[d-1][j] = 0;
//        for (int i = 0; i < d; i++) 
//        	temp[i][d-1] = 0;     
//        temp[0][d-1] = m_a[d-1];
//    	temp[d-1][d-1] = this->m_modulo;
//        
//        this->setDim (d);
//        this->m_basis.SetDims(d, d);        
//        this->m_basis = temp;              
//        this->setNegativeNorm ();        
//
//        if (!this->m_withDual) return;
//
//    	              
//        //Use old basis for first d - 1 dimension
//        for (int i = 0; i < d-1; i++) {
//           for (int j = 0; j < d-1; j++) {
//              temp[i][j] = this->m_dualbasis[i][j];	
//           }
//        }
//        //Add extra coordinate to each vector
//       	for (int j = 0; j < d-1; j++)
//        	temp[d-1][j] = 0;
//        for (int i = 0; i < d; i++) 
//           	temp[i][d-1] = 0;      
//        temp[d-1][0] = m_a[d-1];
//        temp[d-1][d-1] = 1;      
////        //Code for general k
////        Int add;  
////        //Calculate the new basis vector
////        //Proceed columwise
////        for (int j = 0; j < d; j++) {
////           if (d-1 != j) temp[d-1][j] = 0;
////           else temp[d-1][j] = 1;
////           //Sum over the rows
////           for (int i = 0; i < d-1; i++) {
////              add = this->m_basis[i][d-1] * temp[i][j];
////              add = add / this->m_modulo;
////              temp[d-1][j] = temp[d-1][j] - add; 
////           }
////           if (temp[d-1][j] != 0)
////              temp[d-1][j] = temp[d-1][j] % this->m_modulo;
////        }
////        this->m_dualbasis.SetDims(d, d);
//        this->m_dualbasis = temp;
//        this->setDualNegativeNorm ();
//    }

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::incDimBasis () {
		int64_t d = 1 + this->getDim();
    	assert(d <= this->m_maxDim);
    	IntMat temp;
    	temp.SetDims(d, d);
    	//Use old basis for first d - 1 dimension
    	for (int i = 0; i < d-1; i++) {
            for (int j = 0; j < d-1; j++) {
               temp[i][j] = this->m_basis[i][j];	
            }
        }
     	//Fill in the new component
    	for (int j = 0; j < d-1; j++)
    		temp[d-1][j] = 0;
    	temp[d-1][d-1] = this->m_modulo; 
        for (int i = 0; i < d-1; i++) {
        	temp[i][d-1] = this->m_a[d-1] * this->m_basis[i][0];
        	temp[i][d-1] = temp[i][d-1] % this->m_modulo;
        }
        
        this->setDim (d);
        this->m_basis.SetDims(d, d);        
        this->m_basis = temp;              
        this->setNegativeNorm ();        

        if (!this->m_withDual) return;
        
        Int add;
    	              
        //Use old basis for first d - 1 dimension
        for (int i = 0; i < d-1; i++) {
           for (int j = 0; j < d-1; j++) {
              temp[i][j] = this->m_dualbasis[i][j];	
           }
        }
        //Add extra coordinate to each vector
        for (int i = 0; i < d; i++) {
           	temp[i][d-1] = 0;
           	temp[d-1][i] = 0;
        }
        temp[d-1][d-1] = 1;
        
       	for (int j = 0; j < d-1; j++) {
       		add = 0;
       		for (int i = 0; i < d-1; i++) {
       			add = add - this->m_basis[i][d-1] * this->m_dualbasis[i][j];
       		}
   			add = add / this->m_modulo;
   			temp[d-1][j] = temp[d-1][j] + add;
       	}
        this->m_dualbasis = temp;
        this->setDualNegativeNorm ();
    }


//============================================================================

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::incDimDualBasis () {
	int64_t d = 1 + this->getDim();
	assert(d <= this->m_maxDim);
    this->setDim (d);
	IntMat temp;
	temp.SetDims(d, d);
	
    //Use old basis for first d - 1 dimension
    for (int i = 0; i < d-1; i++) {
       for (int j = 0; j < d-1; j++) {
          temp[i][j] = this->m_dualbasis[i][j];	
       }
    }
    //Add extra coordinate to each vector
    for (int i = 0; i < d; i++) {
    	temp[i][d-1] = 0;
    }         
    temp[d-1][0] = -m_a[d-1];
    //temp[d-1][0] = temp[d-1][0] % this->m_modulo;
    temp[d-1][d-1] = 1;
    this->m_dualbasis.SetDims(d, d);
    this->m_dualbasis = temp;
    this->setDualNegativeNorm ();	
    
    //The dimension of the primal lattice has to be increased, as well
    this->m_basis.SetDims(d, d);      
}

//============================================================================

// The basis is built directly, as explained in the guide of LatMRG.
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildBasis (int64_t d) {
      assert(d <= this->m_maxDim);
      this->setDim (d);
      this->m_basis.SetDims(d,d);
      this->m_dualbasis.SetDims(d,d);
      int64_t i, j;

      // This builds an upper triangular basis in a standard way.
      for (j = 0; j < d; j++) {
        this->m_basis[0][j] = this->m_a[j];
      }
      for (i = 1; i < d; i++) {
        for (j = 0; j < d; j++) {
          if (i == j) this->m_basis[i][i] = this->m_modulo;
          else this->m_basis[i][j] = 0;
        }
      }
      this->setNegativeNorm ();

      if (!this->m_withDual) return;
      // If `withDual`, we construct the dual basis also in a direct way.
      
      this->m_dualbasis[0][0] = this->m_modulo;
      for (j = 1; j < d; j++)
         this->m_dualbasis[0][j] = 0;
      for (i = 1; i < d; i++) {
         this->m_dualbasis[i][0] = -this->m_basis[0][i];
         for (j = 1; j < d; j++) {
            if (i == j) this->m_dualbasis[i][i] = 1;
            else this->m_dualbasis[i][j] = 0;
         }
      }
      this->setDualNegativeNorm ();
      
//      //New implementation for general k according to L'Ecuyer and Couture (1997) 
//      for (int i = 0; i < this->m_k; i++) {
//         for (int j = 0; j < d; j++) {
//        	 if (i==j) this->m_dualbasis[i][j] = this->m_modulo; 
//        	 else this->m_dualbasis[i][j] = 0;
//         }
//      }
//      for (int i = this->m_k; i < d; i++) {
//    	  for (int j = 0; j < d; j++) {
//              if (i==j) this->m_dualbasis[i][j] = 1;
//              else this->m_dualbasis[i][j] = 0;
//              if (j < this->m_k) {
//                  for (int k = 0; k < this->m_k; k++) {
//        	    	 this->m_dualbasis[i][j] = this->m_dualbasis[i][j] - m_a[k] * this->m_basis[j][i - k];
//                  }
//              }
//          }
//      }
    }

//============================================================================

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildDualBasis (int64_t d) {
    assert(d <= this->m_maxDim);
    this->setDim (d);
    this->m_basis.SetDims(d,d);
    this->m_dualbasis.SetDims(d,d);
    for (int i = 0; i < d; i++)
    	this->m_dualbasis[i][0] = -m_a[i];
    for (int i = 0; i < d; i++) {
    	if (i<1) this->m_dualbasis[i][i] = this->m_modulo;
    	else this->m_dualbasis[i][i] = 1;
    }
   
//    for (int i = 0; i < this->m_k; i++) {
//       for (int j = 0; j < d; j++) {
//      	 if (i==j) this->m_dualbasis[i][j] = 1; 
//      	 else this->m_dualbasis[i][j] = 0;
//       }
//    }
//    for (int i = this->m_k; i < d; i++) {
//  	  for (int j = 1; j < d+1; j++) {
//  		  if(d-j >= this->m_k) {
//  			  if(i==d-j) this->m_dualbasis[i][d-j] = 1;
//  			  else this->m_dualbasis[i][d-j] = 0;
//  		  } else {
//  			  for (int k = 0; k < this->m_k; k++) {
//  				  this->m_dualbasis[i][d-j] = this->m_dualbasis[i][d-j] + m_a[k+1] * this->m_dualbasis[i-(k+1)][d-j]; 
//  			  }
//  			  this->m_dualbasis[i][d-j] = this->m_dualbasis[i][d-j] % this->m_modulo;
//  		  }
//      }
//   }
//   for (int i = this->m_k; i < d; i++)
//	   for (int j = 0; j < this->m_k; j++)
//		   this->m_dualbasis[i][j] = -this->m_dualbasis[i][j]; 
//   for (int i = 0; i < this->m_k; i++) 
//	   this->m_dualbasis[i][i] = this->m_modulo;
}

//============================================================================

/*
 *  This is already in IntLattice!
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::dualize () {
      IntMat tmps(this->m_basis);
      this->m_basis = this->m_dualbasis;
      this->m_dualbasis = tmps;
    }
*/

//============================================================================

/**
 * Selects and stores a vector of indices with lacunary values.
 */
template<typename Int>
void setLac(const Lacunary<Int>& lac) {           // ??????????
	return;
}


//============================================================================

template class Rank1Lattice<std::int64_t, double>;
template class Rank1Lattice<NTL::ZZ, double>;
template class Rank1Lattice<NTL::ZZ, NTL::RR>;

} // End namespace LatticeTester

#endif
