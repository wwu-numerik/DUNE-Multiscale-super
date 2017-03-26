// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __DUNE_GRID_HOWTO_ELEMENT_DATA_HH
#define __DUNE_GRID_HOWTO_ELEMENT_DATA_HH

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "../Permeability/SRUpscale.hh"

//! Parameter for mapper class
/** This class is only here to show what such a class looks like -- it does
    exactly the same as Dune::MCMGElementLayout. */
template<int dimgrid>
struct P0Layout
{
    bool contains (Dune::GeometryType gt)
    {
        if (gt.dim()==dimgrid) return true;
        return false;
    }
};

// demonstrate attaching data to elements
template<class G, class F>
void elementdata (const G& grid, const F& f, const char * name, bool useRenorm, int times)
{
    // the usual stuff
    //const int dim = G::dimension;
    const int dimworld = G::dimensionworld;
    typedef typename G::ctype ct;
    typedef typename G::LeafGridView GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementLeafIterator;
    typedef typename ElementLeafIterator::Entity::Geometry LeafGeometry;

    // get grid view on leaf part
    GridView gridView = grid.leafGridView();

    // make a mapper for codim 0 entities in the leaf grid
    Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout>
            mapper(grid);

    // allocate a vector for the data
    std::vector<double> c(mapper.size());

    // iterate through all entities of codim 0 at the leaves
    for (ElementLeafIterator it = gridView.template begin<0>();
         it!=gridView.template end<0>(); ++it)
    {
        // cell geometry
        const LeafGeometry geo = it->geometry();

        Dune::FieldVector<ct,dimworld> global = geo.center();
        // evaluate functor and store value
        c[mapper.index(*it)] = f.getSchiftPerm(global, useRenorm, times);
        //c[mapper.index(*it)] = f(global);
    }
    // generate a VTK file
    Dune::VTKWriter<typename G::LeafGridView> vtkwriter(gridView);
    vtkwriter.addCellData(c,"permDist");
    vtkwriter.write( name, Dune::VTK::appendedraw );
}

#endif //__DUNE_GRID_HOWTO_ELEMENT_DATA_HH
