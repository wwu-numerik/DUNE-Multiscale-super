#include "mlmc.hh"
#include <unistd.h>

namespace MultiLevelMonteCarlo {


/// Terminates program if condition is not satisfied.
/// \param condition   condition to check
/// \param message     message to show when condition is false
/// TODO embed into exception handling
void check(bool condition, const char* message) {
    if(!condition) {
        std::cerr << message << "\n";
        exit(1);
    }
};


bool isPowerOf2(int x) { return (x != 0) && ((x & (x-1)) == 0); }

LevelDiffAccumulator::LevelDiffAccumulator()
    :
      performedRealizations(0),
      quantitySum(0),
      quantitySumSqueard(0),
      computationalTime(0),
      avgSolvingTimeAcc(0),
      p_diffLevels(nullptr),
      procsPerRealization(0),
      levelComm(MPI_COMM_NULL),
      baseSeed(1)
{
    diffGroupsNumber = 0;
    localDiffComm = MPI_COMM_NULL;
    diffCommColor = 0;
    diffQuantitySum = 0;
    difRealizations = 0;
    diffTime = 0;
}

LevelDiffAccumulator::LevelDiffAccumulator(MPI_Comm levelComm)
:
  performedRealizations(0),
  quantitySum(0),
  quantitySumSqueard(0),
  computationalTime(0),
  avgSolvingTimeAcc(0),
  p_diffLevels(nullptr),
  procsPerRealization(1),
  levelComm(levelComm),
  baseSeed(1)
{
    diffGroupsNumber = 0;
    localDiffComm = MPI_COMM_NULL;
    diffCommColor = 0;
    diffQuantitySum = 0;
    diffQuantitySumSqueard = 0;
    difRealizations = 0;
    diffTime = 0;
    assignProcessors();
}

LevelDiffAccumulator::LevelDiffAccumulator(MPI_Comm levelComm, int procsPerRealization)
:
  performedRealizations(0),
  quantitySum(0),
  quantitySumSqueard(0),
  computationalTime(0),
  avgSolvingTimeAcc(0),
  p_diffLevels(nullptr),
  procsPerRealization(procsPerRealization),
  levelComm(levelComm),
  baseSeed(1)
{
    diffGroupsNumber = 0;
    localDiffComm = MPI_COMM_NULL;
    diffCommColor = 0;
    diffQuantitySum = 0;
    diffQuantitySumSqueard = 0;
    difRealizations = 0;
    diffTime = 0;
    assignProcessors();
}

LevelDiffAccumulator::LevelDiffAccumulator(MPI_Comm levelComm, Difference * p_diffLevels, int initialRep)
:
  performedRealizations(0),
  quantitySum(0),
  quantitySumSqueard(0),
  computationalTime(0),
  avgSolvingTimeAcc(0),
  p_diffLevels(p_diffLevels),
  procsPerRealization(1),
  levelComm(levelComm),
  baseSeed(1)
{
    realizationsToDo = initialRep;
    diffGroupsNumber = 0;
    localDiffComm = MPI_COMM_NULL;
    diffCommColor = 0;
    diffQuantitySum = 0;
    diffQuantitySumSqueard = 0;
    difRealizations = 0;
    diffTime = 0;
    assignProcessors();
}


LevelDiffAccumulator::LevelDiffAccumulator(MPI_Comm levelComm, int procsPerRealization, Difference *p_diffLevels, int initialRep)
    :
      performedRealizations(0),
      quantitySum(0),
      quantitySumSqueard(0),
      computationalTime(0),
      avgSolvingTimeAcc(0),
      p_diffLevels(p_diffLevels),
      procsPerRealization(procsPerRealization),
      levelComm(levelComm),
      baseSeed(1)
{
    realizationsToDo = initialRep;
    diffGroupsNumber = 0;
    localDiffComm = MPI_COMM_NULL;
    diffCommColor = 0;
    diffQuantitySum = 0;
    diffQuantitySumSqueard = 0;
    difRealizations = 0;
    diffTime = 0;
    assignProcessors();
}

LevelDiffAccumulator::~LevelDiffAccumulator()
{
    //MPI_Comm_free(&localDiffComm);
}

// Returns number of performed realizations.
int LevelDiffAccumulator::getPerformed() const
{
    return performedRealizations;
}

// Returns average time for realization on present level.
double LevelDiffAccumulator::getAvgComputeTime() const
{
    bool isZero = static_cast<bool> (!performedRealizations);
    if (isZero) {
        return 0; // Check for div by zero
    }
    return avgSolvingTimeAcc / performedRealizations;
}

// Returns average time for realization on present level.
double LevelDiffAccumulator::getTotalTime() const
{
    return computationalTime;
}

// Gets base seed
int LevelDiffAccumulator::getBaseSeed() const
{
    return baseSeed;
}

// Gets defined procs per realization
int LevelDiffAccumulator::getProcsPerRealization() const
{
    return procsPerRealization;
}

// Gets how many realizations is to be done
int LevelDiffAccumulator::getRealizationsToDo() const
{
    return realizationsToDo;
}

// Sets base seed
void LevelDiffAccumulator::setBaseSeed(int seed)
{
    baseSeed = seed;
}

// Sets realizations to do
void LevelDiffAccumulator::setRealizationsToDo(int realizations)
{
    realizationsToDo = realizations;
}

// Returns mean value of realizations so far
double LevelDiffAccumulator::mean() const
{
    if (!static_cast<bool> (performedRealizations)) {
        return 0; // Check div by zero
    }
    return quantitySum / performedRealizations;
}

// Returns empirical variance of realizations so far
double LevelDiffAccumulator::var() const
{
    if (!static_cast<bool> (performedRealizations)) {
        return 0;
    }
    return quantitySumSqueard / performedRealizations - (quantitySum / performedRealizations) * (quantitySum / performedRealizations);
}

// Evaluate difference.
void LevelDiffAccumulator::eval()
{
    assert(!(p_diffLevels == nullptr));
    //assert(realizationsToDo);

    loadBalancer();

    diffTime = MPI_Wtime();
    localEval(difRealizations);
    diffTime = MPI_Wtime() - diffTime;  // Local Computational time
    realizationsToDo = 0;               // Sets Statistical calculator to done state
}

// Eval and synhronize
void LevelDiffAccumulator::evalAndUpdate()
{
    eval();
    updateStatistics();
}

// Returns number of groups on present level.
int LevelDiffAccumulator::groupsCount() const
{
    return diffGroupsNumber;
}

// Return number of processors per group.
int LevelDiffAccumulator::procs() const
{
    return procsPerRealization;
}

// Clear all results so far
void LevelDiffAccumulator::clearResults()
{
    performedRealizations = 0;
    quantitySum = 0;
    computationalTime = 0;

    diffGroupsNumber = 0;
    diffCommColor = 0;
    diffQuantitySum = 0;
    difRealizations = 0;
    diffTime = 0;
}

void LevelDiffAccumulator::initDiff(int offset)
{
    p_diffLevels->init(levelComm, localDiffComm, offset * (diffGroupsNumber + 1));
}

// AssignProcessors, based on procs per problem
void LevelDiffAccumulator::assignProcessors()
{
    assert(levelComm != MPI_COMM_NULL);
    int commRank, commSz;
    MPI_Comm_rank(levelComm, &commRank);
    MPI_Comm_size(levelComm, &commSz);

    assert(commSz % procsPerRealization == 0);

    diffGroupsNumber = commSz / procsPerRealization;
    diffCommColor = commRank / procsPerRealization;
    MPI_Comm_split(levelComm, diffCommColor,commRank, &localDiffComm);
}


void LevelDiffAccumulator::synhronize()
{
    double sendBuf[4] = {0};
    double resvBuf[4] = {0};
    sendBuf[0] = diffTime;
    sendBuf[1] = diffQuantitySum;
    sendBuf[2] = diffQuantitySumSqueard;
    sendBuf[3] = static_cast<double>(difRealizations);
    MPI_Allreduce(&sendBuf, &resvBuf, 4, MPI_DOUBLE, MPI_SUM, levelComm);
    avgSolvingTimeAcc += (resvBuf[0] / procsPerRealization);
    quantitySum += (resvBuf[1] / procsPerRealization);
    quantitySumSqueard += (resvBuf[2] / procsPerRealization);
    performedRealizations += static_cast<int> (resvBuf[3] / procsPerRealization);
    computationalTime = performedRealizations * getAvgComputeTime();  //(performedRealizations / diffGroupsNumber) * getAvgComputeTime()  + (performedRealizations % diffGroupsNumber) * getAvgComputeTime();
}

void LevelDiffAccumulator::loadBalancer()
{
     // Common for all processors, it mai need to be splited
    int additionalRep = realizationsToDo;
    int localRep = additionalRep / diffGroupsNumber;
    int remAdditionalRep = additionalRep % diffGroupsNumber;

    // split it to groups
    difRealizations = localRep;
    // Put remaining repettions to the first rem processors
    if ((diffCommColor + 1) <= remAdditionalRep) {
        ++difRealizations;
    }
}

void LevelDiffAccumulator::localEval(int realizations)
{
    if (MPI_COMM_NULL == localDiffComm) {
        return;
    }

    double realizationEval = 0;
    while (realizations > 0) {
        //diffQuantitySum += test;
        realizationEval = p_diffLevels->eval();
        diffQuantitySum += realizationEval;
        diffQuantitySumSqueard += realizationEval * realizationEval;
        --realizations;
    }
}

void LevelDiffAccumulator::updateStatistics()
{
    synhronize();

    //Zero local variables
    diffQuantitySum = 0;
    diffQuantitySumSqueard = 0;
    difRealizations = 0;
    diffTime = 0;
}


MultiLevelMonteCarloEst::MultiLevelMonteCarloEst()
:
  mean(0),
  totalVar(0),
  totalComputationalTime(0),
  tolerance(0),
  proportionalityConst(0),
  monteCarloComm(MPI_COMM_NULL),
  levelComm(MPI_COMM_NULL)
{}

MultiLevelMonteCarloEst::MultiLevelMonteCarloEst(double tolerance)
:
    mean(0),
    totalVar(0),
    totalComputationalTime(0),
    tolerance(tolerance),
    proportionalityConst(0),
    monteCarloComm(MPI_COMM_NULL),
    levelComm(MPI_COMM_NULL)
{}

MultiLevelMonteCarloEst::MultiLevelMonteCarloEst(MPI_Comm MonteCarloComm, double tolerance)
:
    mean(0),
    totalVar(0),
    totalComputationalTime(0),
    tolerance(tolerance),
    proportionalityConst(0),
    monteCarloComm(MonteCarloComm),
    levelComm(MPI_COMM_NULL)
{
    assignProcessors();
}

void MultiLevelMonteCarloEst::eval()
{
    //print info
    int rank = 0;
    MPI_Comm_rank(monteCarloComm, &rank);
    int parallelCompTime = MPI_Wtime();
    bool converged = false;
    do {
        for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
            StatisticalAccumulators[i].eval();
        }

        for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
            StatisticalAccumulators[i].updateStatistics();
        }
        MPI_Barrier(monteCarloComm);
        computeProportionallityConstant();
        for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
            computeRealizationsOnDiffLevel(i, 0.5);
        }

        for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
            mean += StatisticalAccumulators[i].mean();
            totalComputationalTime += StatisticalAccumulators[i].getTotalTime();
            totalVar += StatisticalAccumulators[i].var();
        }

        // To DO: Add ofstream
        if (rank == 0) {
            for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
                std::cout << "===Level(" << i << ")\n";
                std::cout << "Requierd Add. Real:" << StatisticalAccumulators[i].getRealizationsToDo() << "\n";
                std::cout << "Performed Realizations:" << StatisticalAccumulators[i].getPerformed() << "\n";
                std::cout << "Expected Value:" << StatisticalAccumulators[i].mean() << "\n";
                std::cout << "Var:" << StatisticalAccumulators[i].var() << "\n";
                std::cout << "Time (Total):" << StatisticalAccumulators[i].getTotalTime() << "\n";
                std::cout << "Time (AVG):" << StatisticalAccumulators[i].getAvgComputeTime() << "\n";
                double eff = StatisticalAccumulators[i].getTotalTime();
                double constVal = StatisticalAccumulators[i].getPerformed() % StatisticalAccumulators[i].groupsCount() == 0 ? StatisticalAccumulators[i].getPerformed() / StatisticalAccumulators[i].groupsCount()
                                                                                                                            : StatisticalAccumulators[i].getPerformed() / StatisticalAccumulators[i].groupsCount() + 1;
                eff /=  constVal * (StatisticalAccumulators[i].getAvgComputeTime() * StatisticalAccumulators[i].groupsCount());
                std::cout << "Balancer Eff:" << eff << "\n";
                double avgOnProcs = (1 - eff) * StatisticalAccumulators[i].getAvgComputeTime() * StatisticalAccumulators[i].groupsCount();
                std::cout << "Compute Time lost:" << avgOnProcs  << "\n";
                std::cout << "Compute Time lost(AVG):" << avgOnProcs / StatisticalAccumulators[i].groupsCount()  << "\n";
                std::cout << std::flush;
            }
        }
        // check if converged
        converged = true;
        for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
            converged &= (!static_cast<bool> (StatisticalAccumulators[i].getRealizationsToDo()));
        }
        // Zero computed values
        if (!converged) {
            mean = 0;
            totalComputationalTime = 0;
            totalVar = 0;
        }
    } while(!converged);

    if (rank == 0) {
        parallelCompTime = MPI_Wtime() - parallelCompTime;
        std::cout << "MC computed: E[Q] = " << mean << " Time[s] (Total work time): " << totalComputationalTime  << " Time[s] (mlmc): "  <<  parallelCompTime << "\n";
    }
}

double MultiLevelMonteCarloEst::computeProportionallityConstant()
{
    double alpha = 0;
    int n = static_cast<int>(StatisticalAccumulators.size());
    for(int i=0; i<n; ++i) {
        alpha +=  sqrt(StatisticalAccumulators[i].var() * StatisticalAccumulators[i].getAvgComputeTime());
    }
    alpha /= (tolerance * tolerance);
    proportionalityConst = alpha;
    return alpha;
}

double MultiLevelMonteCarloEst::computeRealizationsOnDiffLevel(int diffLevel, double factor = 1)
{
    assert(diffLevel < StatisticalAccumulators.size());
    double doneSoFar = StatisticalAccumulators[diffLevel].getPerformed();
    double estimatedSoFar =  proportionalityConst * sqrt(StatisticalAccumulators[diffLevel].var() /  StatisticalAccumulators[diffLevel].getAvgComputeTime());
    int requiredAdd = std::max(0.0, (estimatedSoFar + 0.5) - doneSoFar);
//    if (requiredAdd > 100) {
//        requiredAdd = static_cast<int> (factor * requiredAdd);
//        // Fill gapping processors, for better performance
//        requiredAdd += requiredAdd % StatisticalAccumulators[diffLevel].groupsCount();
//    }
    StatisticalAccumulators[diffLevel].setRealizationsToDo(requiredAdd);
    return StatisticalAccumulators[diffLevel].getRealizationsToDo();
}

void MultiLevelMonteCarloEst::addDifference(Difference* diff, int procsPerRealization, int initialRep)
{
    StatisticalAccumulators.push_back(LevelDiffAccumulator(levelComm, procsPerRealization, diff, initialRep));
    int offset =  StatisticalAccumulators.size() - 1;
    StatisticalAccumulators[offset].initDiff(offset);
}

void MultiLevelMonteCarloEst::addDifference(Difference *diff, int initialRep)
{
    addDifference(diff, 1, initialRep);
}

void MultiLevelMonteCarloEst::assignProcessors()
{
    assert(monteCarloComm != MPI_COMM_NULL);
    MPI_Comm_dup(monteCarloComm, &levelComm);
}
} // END namespace MultiLevelMonteCarlo
