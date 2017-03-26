#ifndef __BASETYPES__
#define __BASETYPES__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<dune/grid/yaspgrid.hh>

namespace RandomDiff {
using namespace Dune;

struct BaseTypes
{
    // Grid, Dim and Integral types
    static constexpr int dimProblem = 2;
    typedef double floatType;
    typedef double rangeType;

    // Vectors, Ranges etc....
    typedef FieldVector<floatType, dimProblem> VectorType;

    // Grid
    typedef YaspGrid<dimProblem> GridType;
    typedef typename YaspGrid<dimProblem>::LeafGridView GridViewType;
    typedef typename GridViewType::Grid::ctype Coord;


    //Iterators
    typedef typename GridViewType::IntersectionIterator IntersectionIterator;                   // Intersection iterator
    typedef typename IntersectionIterator::Intersection::Geometry IntersectionGeometry;         // Intersection geometry
};
}
#endif // __BASETYPES__
