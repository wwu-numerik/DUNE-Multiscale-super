#ifndef __BONDARYCOND__
#define __BONDARYCOND__

#include"../baseTypes.hh"
#include"../Utils/floatCompare.hh"

namespace RandomDiff {
using namespace Dune;
class BoundaryCondition
{
public:
    bool isDirchleBoundry(const BaseTypes::VectorType& point) const;
    bool isNeumannBoundry(const BaseTypes::VectorType& point) const;
    double evalDirichleBoundary(const BaseTypes::VectorType& point) const;
};

} // END RANDOMDIFF
#endif // __BONDARYCOND__
