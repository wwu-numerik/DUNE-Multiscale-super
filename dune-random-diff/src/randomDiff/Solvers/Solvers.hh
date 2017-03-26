#ifndef __FVSOLVER__
#define __FVSOLVER__

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>

#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>

#include<dune/pdelab/stationary/linearproblem.hh>

//Local Operatos
//---------------------------------------------
#include"../LocalOperators/CellCenteredFV.hh"
//---------------------------------------------


namespace RandomDiff {
using namespace Dune;

//------------------------------------------------------------------------------------------------------------------------------------
//FiniteVolume Solver
//------------------------------------------------------------------------------------------------------------------------------------

template<class PermeabilityGenType, class BoundaryCond>
class FVSolver
{
    //typedef PDELab::istl::VectorBackend<> VBE;                                                                                           // Vector Backend
    //typedef PDELab::istl::BCRSMatrixBackend<> MBE;                                                                                       // Matrix Backend
    typedef PDELab::ISTLVectorBackend<> VBE;
    typedef PDELab::istl::BCRSMatrixBackend<> MBE;
    typedef FVLOperator<PermeabilityGenType,BoundaryCond> LOP;                                                                           // Local operator

    typedef PDELab::P0LocalFiniteElementMap<BaseTypes::Coord, BaseTypes::floatType, BaseTypes::dimProblem> FEM;                          // Finite element map
    typedef PDELab::NoConstraints CON;                                                                                                   // No constraints!
    typedef PDELab::GridFunctionSpace<BaseTypes::GridViewType, FEM, CON, VBE> GFS;                                                       // Grid Function Space
    typedef PDELab::EmptyTransformation CC;                                                                                              // Empty Constrainst Container

    typedef PDELab::GridOperator<GFS, GFS, LOP, MBE, BaseTypes::floatType, BaseTypes::floatType ,BaseTypes::floatType, CC, CC> GO;       // Grid Operator
    typedef PDELab::ISTLBackend_CG_AMG_SSOR<GO> LS;                                                                                      // Lenear Solver
    typedef typename GO::Traits::Domain U;                                                                                               // Deg of freedom
    typedef PDELab::StationaryLinearProblemSolver<GO,LS,U> StationaryLinearProblemSolver;                                                // Solver

public:
    //********************************
    // gfs cotr, copy or reff ???!
    //********************************
    FVSolver(const PermeabilityGenType& reff_perm, BoundaryCond& reff_boundary, const BaseTypes::GridType& reff_Grid, bool useRenormalization, int times, int sparcityPattern = 5)
        : reff_Grid(reff_Grid), reff_perm(reff_perm), reff_boundary(reff_boundary),
          fem(GeometryType(GeometryType::cube, BaseTypes::dimProblem)),
          gfs(reff_Grid.leafGridView(), fem),
          mbe(sparcityPattern),
          lop(reff_perm, reff_boundary, useRenormalization, times),
          go(gfs, gfs, lop, mbe),
          deg(gfs, 0.0)
    {
        iterSolver = 0;
        timeSolver = 0;
        gfs.name("pressure");
    }

    ~FVSolver() {}

    void apply(double tol = 1e-7, int verboseLevel = 0, int maxIter = 5000, double defectRate = 1e-99)
    {
        LS ls(gfs, maxIter, verboseLevel, false, false);
        StationaryLinearProblemSolver slp(go, ls, deg, tol, defectRate, verboseLevel);
        slp.apply();

        // Update Statistics for the solution
        iterSolver = slp.result().linear_solver_iterations;
        timeSolver = slp.result().linear_solver_time + slp.result().assembler_time;
    }

    void writeToVTKFormat(const char * name) const
    {
        VTKWriter<BaseTypes::GridViewType> vtkwriter(reff_Grid.leafGridView(), VTK::conforming);
        PDELab::addSolutionToVTKWriter(vtkwriter, gfs, deg);
        vtkwriter.write(name, VTK::appendedraw);
    }

    const GFS* solProxy() const
    {
        return &gfs;
    }

    const GFS getSol() const
    {
        return gfs;
    }


    const int getNumOfIter()
    {
        return iterSolver;
    }

    const double getTimeToSolve()
    {
        return timeSolver;
    }

    double surfaceFlow(bool upscale = false, int times = 0)
    {
        const int dim = BaseTypes::dimProblem;                                                // Dim world
        typedef typename BaseTypes::GridViewType::ctype ct;                                   // Coordinates type
        typedef typename BaseTypes::GridViewType::IntersectionIterator IntersectionIterator;  // Intersection iterator
        typedef typename IntersectionIterator::Intersection::Geometry IntersectionGeometry;   // Intersection geometry
        typedef typename PDELab::DiscreteGridFunction<GFS,U> DGF;                             // Discrate grid function space
        typedef FieldVector<ct, 1> Scalar;
        typedef FieldVector<ct, dim> Vector;


        DGF dgf(gfs, deg); // Descrate grid function space (Solution)

        ct localFlux = 0;
        ct k_avg;
        ct aproxGrad;


        Scalar val;
        Vector dist;
        Vector atPos;

        const BaseTypes::GridViewType& reff_gv = reff_Grid.leafGridView();
        auto endit = reff_gv.end<0, Dune::Interior_Partition>();
        for ( auto it = reff_gv.begin<0, Dune::Interior_Partition>(); it != endit; ++it) {
            IntersectionIterator isend = reff_gv.iend(*it);
            for (IntersectionIterator isit = reff_gv.ibegin(*it); isit != isend; ++isit) {
                if(isit->boundary() && floatCompare<ct>(isit->geometry().center()[0], 0, 1e-15)) {
                    const IntersectionGeometry igeo = isit->geometry();
                    // dist in global coordinates
                    dist = it->geometry().center() - igeo.center();
                    dist.two_norm();

                    //Position in local coordinates
                    atPos = it->geometry().local(it->geometry().center());
                    dgf.evaluate(*it, atPos, val);
                    aproxGrad = (val - 1) / dist[0];
                    atPos = isit->geometry().center();
                    //**********************************************
                    k_avg = reff_perm.getSchiftPerm(it->geometry().center(), upscale , times);
                    //**********************************************
                    //reff_perm(atPos);
                    //**********************************************

                    localFlux -= aproxGrad * k_avg * igeo.volume();
                }
            }
        }
        return localFlux;
    }
private:

    const BaseTypes::GridType& reff_Grid;
    const PermeabilityGenType& reff_perm;
    const BoundaryCond& reff_boundary;

    int iterSolver;
    double timeSolver;

    FEM fem;
    GFS gfs;
    MBE mbe;
    LOP lop;
    GO go;
    U deg;
};

} // END NAMEPACE RANDOMDIFF
#endif // __FVSOLVER__
