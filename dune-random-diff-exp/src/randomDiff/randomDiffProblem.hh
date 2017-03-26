#ifndef __RANDOM_DIFF__
#define __RANDOM_DIFF__

#include<array>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>


#include "../baseTypes.hh"
#include "../Permeability/permeability.hh"
#include "../Permeability/CEGenerator.hh"
#include "../Permeability/correlation.hh"
#include "../Utils/elementdata.hh"
#include "../Utils/vertexdata.hh"

#include "Solvers/Solvers.hh"

/*
 * To Do:
 * Replace const char* with std::string
 * Duplicate MPI_Comm
 * Input params checks
 */

namespace RandomDiff {
using namespace Dune;

struct RandomDiffTypes
{
    typedef CorrelationTwoNorm<BaseTypes::dimProblem, BaseTypes::VectorType, BaseTypes::floatType> CorrelationType;
    typedef CirculantEmbeding<BaseTypes::dimProblem, BaseTypes::VectorType, BaseTypes::floatType, CorrelationType > PermeabilityType;
    //typedef Permeability<BaseTypes::dimProblem, BaseTypes::VectorType, BaseTypes::floatType, CorrelationType > PermeabilityType;
    typedef BoundaryCondition BoundaryConditionType;
    typedef FVSolver<PermeabilityType, BoundaryCondition> SolverType;

};

class RandomDiffProblem
{
public:

    typedef RandomDiffTypes::PermeabilityType PermeabilityGenType;

    RandomDiffProblem(int log2seg, int overlap, const PermeabilityGenType& reff_bindPerm, int renormTimes);
    ~RandomDiffProblem();

    void init(MPI_Comm comm);

    int getLog2seg() const;
    int getOverlap() const;
    MPI_Comm getComm() const;

    bool useRenormalization() const;
    double getRenormalizationTimes() const;
    double getSolTime() const;
    int getIterTime() const;

    void applaySolver(double tol = 1e-7, int verboseLevel = 0, int maxIter = 5000, double defectRate = 1e-99);

    double surfaceFlow();

    void writeSolToVtKFormat(const char * name);
    void writePermToVtkFormat(const char * name);
    void writePermToVtkFormatVartex(const char * name);
private:
    int log2seg;
    int overlap;

    bool useRenorm;
    double renormTimes;

    MPI_Comm comm;

    const RandomDiffTypes::PermeabilityType& perm;
    const BaseTypes::GridType* unitCube;
    RandomDiffTypes::SolverType* p_solver;
    RandomDiffTypes::BoundaryConditionType boundary;
};

} // EDN RANDOMDIFF
#endif // __RANDOM_DIFF__
