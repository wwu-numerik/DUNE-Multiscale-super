#include "boundaryCond.hh"

using namespace RandomDiff;

bool BoundaryCondition::isDirchleBoundry(const BaseTypes::VectorType &point) const
{
    bool leftB = floatCompare<BaseTypes::floatType>(point[0], 0, 5e-15);
    bool rightB = floatCompare<BaseTypes::floatType>(point[0], 1, 5e-15);
    if (leftB || rightB) {
        return true;
    }
    return false;
}

bool BoundaryCondition::isNeumannBoundry(const BaseTypes::VectorType &point) const
{
    return (!isDirchleBoundry(point));
}

double BoundaryCondition::evalDirichleBoundary(const BaseTypes::VectorType &point) const
{
    double boundaryVal = 0;
    floatCompare<BaseTypes::floatType> (point[0], 0, 5e-15) ? boundaryVal = 1 : boundaryVal = 0;
    return boundaryVal;
}
