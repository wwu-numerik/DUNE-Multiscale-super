#include "randomDiffProblem.hh"

using namespace RandomDiff;

RandomDiffProblem::RandomDiffProblem(int log2seg, int overlap, const PermeabilityGenType& reff_bindPerm, int renormTimes)
: log2seg(log2seg), overlap(overlap),
  useRenorm(static_cast<bool>(renormTimes)), renormTimes(renormTimes),
  perm(reff_bindPerm), unitCube(nullptr), p_solver(nullptr)
{}

void RandomDiffProblem::init(MPI_Comm comm)
{
    this->comm = comm;

    // Create grid
    std::array<int, BaseTypes::dimProblem> cellsPerDirection;
    int cellNum = (1<<log2seg);
    std::fill(cellsPerDirection.begin(), cellsPerDirection.end(), cellNum);
    BaseTypes::VectorType upperRightCorner(1.0);
    bool periodic = false;

    if (unitCube != nullptr || p_solver != nullptr) {
        throw("Already Inited");
    }
    unitCube = new BaseTypes::GridType(upperRightCorner, cellsPerDirection, periodic, overlap, comm);
    p_solver = new RandomDiffTypes::SolverType(perm, boundary, *unitCube, useRenorm, renormTimes);

}

RandomDiffProblem::~RandomDiffProblem()
{
    if (unitCube != nullptr || p_solver != nullptr) {
        delete this->p_solver;
        delete this->unitCube;
    }
    unitCube = nullptr;
    p_solver = nullptr;
}


int RandomDiffProblem::getLog2seg() const
{
    return log2seg;
}

int RandomDiffProblem::getOverlap() const
{
    return overlap;
}

MPI_Comm RandomDiffProblem::getComm() const
{
    return comm;
}

bool RandomDiffProblem::useRenormalization() const
{
    return useRenorm;
}

double RandomDiffProblem::getRenormalizationTimes() const
{
    return renormTimes;
}

double RandomDiffProblem::getSolTime() const
{
    return p_solver->getTimeToSolve();
}

int RandomDiffProblem::getIterTime() const
{
    return p_solver->getNumOfIter();
}

void RandomDiffProblem::applaySolver(double tol, int verboseLevel, int maxIter, double defectRate)
{
    p_solver->apply(tol, verboseLevel, maxIter, defectRate);
}

double RandomDiffProblem::surfaceFlow()
{

    double flux = 0;
    double localFlux = 0;
    assert(sizeof(MPI_DOUBLE) == sizeof(double));
    localFlux = p_solver->surfaceFlow(useRenorm, renormTimes);
    MPI_Allreduce(&localFlux, &flux, 1, MPI_DOUBLE, MPI_SUM, unitCube->comm());
    return flux;
}

void RandomDiffProblem::writeSolToVtKFormat(const char* name)
{
    p_solver->writeToVTKFormat(name);
}

void RandomDiffProblem::writePermToVtkFormat(const char* name)
{
    elementdata<BaseTypes::GridType, RandomDiffTypes::PermeabilityType>(*unitCube, perm, name, static_cast<bool>(renormTimes), renormTimes);
}

void RandomDiffProblem::writePermToVtkFormatVartex(const char *name)
{
    vertexdata<BaseTypes::GridType, RandomDiffTypes::PermeabilityType>(*unitCube, perm, name);
}
