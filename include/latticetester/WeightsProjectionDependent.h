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

#ifndef LATTICETESTER__PROJECTION_DEPENDENT_WEIGHTS_H
#define LATTICETESTER__PROJECTION_DEPENDENT_WEIGHTS_H

#include "latticetester/Weights.h"

#include <map>
#include <vector>

namespace LatticeTester {

/**
 * \class WeightsProjectionDependent
 *
 * Defines projection-dependent weights, for which the weight for any given projection
 * can be set individually by #setWeight().
 * Internally, the weights are regrouped by largest coordinate index in
 * different `std::map` objects.
 */
class WeightsProjectionDependent: public Weights {
protected:

   typedef std::map<Coordinates, Weight> WeightsMap;

   /// Per-projection weights, regrouped by largest coordinate index.
   std::vector<WeightsMap> m_weights;

   /// Used only to return an empty map.
   static const WeightsMap m_emptyWeights;

public:

   /**
    * Constructs projection-dependent weights.
    */
   WeightsProjectionDependent();

   /**
    * Copy constructor.
    */
   WeightsProjectionDependent(const WeightsProjectionDependent&);

   /**
    * Destructor.
    */
   virtual ~WeightsProjectionDependent() {
   }

   /**
    * Returns the weight of the projection specified by `projection`.
    */
   virtual Weight getWeight(const Coordinates &projection) const;

   /**
    * Returns the number of weights.
    */
   virtual uint64_t getSize() const {
      return (uint64_t) m_weights.size();
   }

   /**
    * Returns a map of weights for all projections whose largest index is
    * `largestIndex`.
    */
   virtual const WeightsMap& getWeightsForLargestIndex(Coordinates::value_type largestIndex) const;

   /**
    * Sets the weight of the projection specified by `projection`.
    */
   virtual void setWeight(const Coordinates &projection, Weight weight);

protected:
   /// \copydoc LatticeTester::Weights::format()
   virtual void format(std::ostream &os) const;

   friend std::istream& operator>>(std::istream&, WeightsProjectionDependent&);
};

/**
 * \relates WeightsProjectionDependent
 * Reads formatted projection-dependent weights into the object \c weights.
 *
 * The input should be a sequence of projection-to-weight mappings, of the format:
 * \code
 * <match1>: <weight1>, <match2>: <weight2>, ...
 * \endcode
 * where <tt>\<weight<i>n</i>\></tt> is the weight (a floating point number)
 * associated to the projection-match <tt>\<match<i>n</i>\></tt>, and
 * <tt>\<match<i>n</i>\></tt> is one of:
 * - a set of coordinates, as specified in
 *   #operator>>(std::istream&, LatticeTester::Coordinates&)
 *   to explicitly set the weight for the projection that
 *   correspond to these coordinates;
 * - the string <tt>order <i>m</i></tt> to implicitly set the weights of
 *   projections of order <tt><i>m</i></tt>;
 * - the string \c default to set the default weight for other projections.
 *
 * \remark Comments can be appended to any line after a \c # character.
 *
 * \remark The match-weight pairs can be separated by commas and/or whitespace
 * (including newlines).
 *
 * \remark The colons (\c :) can be replaced with <tt>=\></tt> or <tt>-\></tt>.
 *
 * \sa  #operator>>(std::istream&, LatticeTester::Coordinates&)
 */
std::istream& operator>>(std::istream &is, WeightsProjectionDependent &weights);

}
#endif
