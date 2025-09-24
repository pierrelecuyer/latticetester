// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
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

#ifndef LATTICETESTER__ENUMTYPES_H
#define LATTICETESTER__ENUMTYPES_H
#include <string>
#include <array>

using namespace std;

/**
 * \file EnumTypes.h
 *
 * Enumeration types and global constants used in LatticeTester.
 * Some of them may be no longer used or deprecated.
 * This file should be cleaned up at some point, when we are done.
 */

namespace LatticeTester {

  /**
   * The available norm types to measure the length of vectors.
   * For \f$X = (x_1,…,x_t)\f$:<br>
   * `SUPNORM` corresponds to \f$\Vert X\Vert= \max(|x_1|,…,|x_t|)\f$.<br>
   * `L1NORM` corresponds to \f$\Vert X\Vert= |x_1|+\cdots+|x_t|\f$.<br>
   * `L2NORM` corresponds to \f$\Vert X\Vert= (x_1^2+\cdots+x_t^2)^{1/2}\f$.<br>
   * `ZAREMBANORM` corresponds to \f$\Vert X\Vert= \max(1, |x_1|)\cdots\max(1, |x_t|)\f$.
   */
  enum NormType { SUPNORM, L1NORM, L2NORM, ZAREMBANORM };

  /**
   * Different choices of output formats.
   *
   * `TERM`: the results will appear only on the terminal screen.<br>
   * `RES`: the results will be in plain text and sent to a `.res` file.<br>
   * `TEX`: the results will be in a LaTeX file with extension `.tex`.<br>
   * `GEN`: a list of retained generators will be sent to a file with extension `.gen`,
   *        in a specific format so this list can be read again for further analysis.
   */
  // enum OutputType { TERM, RES, TEX, GEN };

  /**
   * Types of problems that LatticeTester can handle.
   * Not sure if we still need this.    *******
   */
  // enum ProblemType { BASIS, DUAL, REDUCTION, SHORTEST, MERIT };

  /**
   * This can be supersed by the Real type.
   *
   * Types of precision that the NTL can use for real numbers:
   * `DOUBLE` -- double
   * `QUADRUPLE` -- quad_float (quasi quadruple precision)
   *         this is useful when roundoff errors can cause problems
   * `XDOUBLE` -- xdouble (extended exponent doubles)
   *         this is useful when numbers get too big
   * `RR` -- RR (arbitrary precision floating point).
   * The choice `DOUBLE` is usually the fastest,
   * but may be prone to roundoff errors and/or overflow.
   * See `https://github.com/u-u-h/NTL/blob/master/doc/LLL.txt`.
   */
 // enum PrecisionType { DOUBLE, QUADRUPLE, EXPONENT, ARBITRARY, EXACT };
 // enum PrecisionType { DOUBLE, QUADRUPLE, XDOUBLE, RR };

  /**
   * Indicates whether an integer is prime, probably prime, composite or its
   * status is unknown (or we do not care).
   */
  // enum PrimeType { PRIME, PROB_PRIME, COMPOSITE, UNKNOWN };

  /**
   * Merit criteria to measure the quality of generators or lattices.
   * TO DO: this list is not very clear. Maybe outdated.  ****************
   *
   * `LENGTH`: Only using the length of the shortest vector as a criterion.
   * `SPECTRAL`: figure of merit \f$S_T\f$ based on the spectral test.<br>
   * `BEYER`: figure of merit is the Beyer quotient \f$Q_T\f$.<br>
   * `PALPHA`: figure of merit based on \f$P_{\alpha}\f$.<br>
   * <tt>BOUND_JS</tt>: figure of merit based on
   *     the Joe-Sinescu bound \cite rSIN08a.<br>   ???
   */
  // enum CriterionType { LENGTH, SPECTRAL, BEYER, PALPHA, BOUND_JS };

  /**
   * Different types of normalizations that can be used for shortest-vector lengths.
   * Corresponds to different ways of approximating the Hermite constants `gamma_t`.
   *
   * `BESTLAT`: the value used for \f$d_t^*\f$ corresponds to the best
   * lattice.<br>
   * `BESTUPBOUND`: the value used for \f$d_t^*\f$ corresponds to the best
   * bound known to us.<br>
   * `LAMINATED`: the value used for \f$d_t^*\f$ corresponds to the best
   * *laminated* lattice.<br>
   * `ROGERS`: the value for \f$d_t^*\f$ is obtained from *Rogers’* bound on the
   * density of sphere packing.<br>
   * `MINKL1`: the value for \f$d_t^*\f$ is obtained from the theoretical bounds
   * on the length of the shortest nonzero vector in the lattice using the
   * \f${\mathcal{L}}_1\f$ norm.<br>
   * `MINKHLAW`: the value for \f$d_t^*\f$ is obtained from *Minkowski’*
   * theoretical bounds on the length of the shortest nonzero vector in the
   * lattice using the \f${\mathcal{L}}_2\f$ norm.<br>
   * `NONE`: no normalization will be used.<br>
   */
 // enum NormaType { BESTLAT, BESTUPBOUND, LAMINATED, ROGERS, MINKOWSKI, MINKL1, MINK, L1, L2, NONE };
  enum NormaType { BESTLAT, BESTUPBOUND, LAMINATED, ROGERS, MINKL1, MINKHLAW, NONE };

  /**
   * Indicates which type of calculation is considered for the
   * \f$P_{\alpha}\f$ test. \anchor REF__Const_CalcType_def
   * Is this used anywhere?         ************
   *
   * `PAL` is for the \f$P_{\alpha}\f$ test. <br>
   * `BAL` is for the bound on the \f$P_{\alpha}\f$ test. <br>
   * `NORMPAL` is for the \f$P_{\alpha}\f$ test `PAL`, with the result normalized
   *      over the `BAL` bound. <br>
   * `SEEKPAL` is for the \f$P_{\alpha}\f$ seek, which searches
   *      for good values of the multiplier.
   */
  enum CalcType { PAL, NORMPAL, BAL, SEEKPAL };

  /**
   * A list of all the possible lattice reductions implemented in `LatticeTester`.
   *
   * `PAIR`: Pairwise reductions only.
   * `LLL`: LLL reduction only.
   * `BKZ`: block Korkine-Zolotarev reduction only.
   * `BB`: direct shortest vector search with BB (no pre-red.).
   * `PAIRBB`: Pairwise reduction followed by BB.
   * `LLLBB`: LLL followed by BB.
   * `BKZBB`: BKZ followed by BB.
   */
   // enum PreReductionType { NOPRERED, DIETER, LLL, BKZ, FULL };
  enum ReductionType { PAIR, LLL, BKZ, BB, PAIRBB, LLLBB, BKZBB };

  /**
   * Two possible ways of obtaining a triangular matrix to compute the bounds
   * in the BB algorithm.
   *
   * `CHOLESKY`: use a lower-triangular matrix obtained as the Cholesky decomposition
   *             of the matrix of scalar products.
   * `TRIANGULAR`: use a lower-triangular basis
   */
  enum DecompTypeBB { CHOLESKY, TRIANGULAR };
  
  /**
   * Two possible ways of computing the basis for a projection.
   *
   * `LLLPROJ`: uses LLL reduction.
   * `UPPERTRIPROJ`: use an upper-triangular basis construction.
   */
  enum ProjConstructType { LLLPROJ, UPPERTRIPROJ };
  
  /**
   * Two different types of figures of merit
   *
   * `MERITM`: based on shortest vector.
   * `MERITQ`: based on Beyer quotient.
   */
  enum MeritType { MERITM, MERITQ };


  //============================================================================

  /**
   * The following are functions for printing the `enum` constants in this module.
   * Each function returns the value of the `enum` variable given as input as a string.
   */

  static std::string toStringNorm (NormType norm) {
    switch (norm) {
      case SUPNORM:
        return "SUPNORM";
      case L1NORM:
        return "L1NORM";
      case L2NORM:
        return "L2NORM";
      case ZAREMBANORM:
        return "ZAREMBANORM";
      default:
        return "***** NormType: UNDEFINED CASE ";
     }
  }
  /*
  static std::string toStringOutput (OutputType out) {
    switch (out) {
      case TERM:
        return "TERM";
      case RES:
        return "RES";
      case TEX:
        return "TEX";
      case GEN:
        return "GEN";
      default:
        return "***** OutputType: UNDEFINED CASE ";
    } 
  }
  
  static std::string toStringProblem (ProblemType prob) {
	switch (prob) {
	   case BASIS:
	     return "BASIS";
	   case DUAL:
	     return "DUAL";
	   case REDUCTION:
	     return "REDUCTION";
	   case SHORTEST:
	     return "SHORTEST";
	   case MERIT:
	     return "MERIT";
	   default:
	     return "***** ProblemType: UNDEFINED CASE ";
	}
  }
  
  static std::string toStringPrecision (PrecisionType precision) {
    switch (precision) {
      case DOUBLE:
        return "DOUBLE";
      case QUADRUPLE:
        return "QUADRUPLE";
      case XDOUBLE:
        return "XDOUBLE";
      case RR:
        return "RR";
      default:
        return "***** PrecisionType: UNDEFINED CASE ";
    }
  }

  static std::string toStringPrime (PrimeType prim) {
	switch (prim) {
	  case PRIME:
	     return "PRIME";
	  case PROB_PRIME:
	     return "PROB_PRIME";
	  case COMPOSITE:
	     return "COMPOSITE";
	  default:
	    return "UNKNOWN";
	}
  }

  static std::string toStringCriterion (CriterionType criter) {
    switch (criter) {
      case LENGTH:
        return "LENGTH";
      case SPECTRAL:
        return "SPECTRAL";
      case BEYER:
        return "BEYER";
      case PALPHA:
        return "PALPHA";
      default:
        return "***** CriterionType: UNDEFINED CASE ";
    }
  }
*/
  static std::string toStringNorma (NormaType norma) {
    switch (norma) {
      case BESTLAT:
        return "BESTLAT";
      case BESTUPBOUND:
        return "BESTUPBOUND";
      case LAMINATED:
        return "LAMINATED";
      case ROGERS:
        return "ROGERS";
      case MINKL1:
        return "MINKL1";
      case MINKHLAW:
        return "MINKHLAW";
      case NONE:
        return "NONE";
      default:
        return "***** NormaType: UNDEFINED CASE ";
    }
  }
  /*
  static std::string toStringCalc (CalcType calc) {
    switch (calc) {
      case PAL:
        return "PAL";
      case NORMPAL:
        return "NORMPAL";
      case BAL:
        return "BAL";
      case SEEKPAL:
        return "SEEKPAL";
      default:
        return "***** CalcType: UNDEFINED CASE ";
    }
  }
*/
  static std::string toStringReduction (ReductionType reduct) {
    switch (reduct) {
      case PAIR:
        return "PAIR";
      case LLL:
        return "LLL";
      case BKZ:
        return "BKZ";
      case BB:
        return "BB";
      case PAIRBB:
        return "PAIRBB";
      case LLLBB:
        return "LLLBB";
      case BKZBB:
        return "BKZBB";
      default:
        return "***** ReductionType: UNDEFINED CASE ";
    }
  }
  
  static std::string toStringDecomp(DecompTypeBB decomp) {
    switch (decomp) {
      case CHOLESKY:
        return "CHOLESKY";
      case TRIANGULAR:
        return "TRIANGULAR";
      default:
        return "***** DecompType: UNDEFINED CASE ";
    }
  }
  
  static std::string toStringProjConstruct(ProjConstructType proj) {
    switch (proj) {
      case LLLPROJ:
        return "LLLPROJ";
      case UPPERTRIPROJ:
        return "UPPERTRIPROJ";
      default:
        return "***** ProjConstructType: UNDEFINED CASE ";
    }
  }
  
  static std::string toStringMeritType(MeritType merit) {
    switch (merit) {
      case MERITM:
        return "MERITM";
      case MERITQ:
        return "MERITQ";
      default:
        return "***** MeritType: UNDEFINED CASE ";
    }
  }

};

#endif
