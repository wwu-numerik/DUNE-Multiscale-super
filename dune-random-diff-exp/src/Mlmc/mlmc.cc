#include "mlmc.hh"
#include <unistd.h>

namespace MultiLevelMonteCarlo {

// To Do: Remove
//--------------------------------------------------
void check(bool condition, const char* message)
{
    if(!condition) {
        std::cerr << message << "\n";
        exit(1);
    }
}
// To Do: Move to Utils
//--------------------------------------------------
bool isPowerOf2(int x)
{
    return (x != 0) && ((x & (x-1)) == 0);
}


StaticAccRawData::StaticAccRawData()
{
    quantitySum = 0;
    quantitySumSqueard = 0;
    computationalTime = 0;
    avgSolvingTimeAcc = 0;
    performedRealizations = 0;
    realizationsToDo = 0;
    tag = -1;
}

ostream& operator<<(ostream& os, const StaticAccRawData& sD)
{
    os << "quantitySum:           " << sD.quantitySum << "\n";
    os << "quantitySumSqueard:    " << sD.quantitySumSqueard << "\n";
    os << "computationalTime:     " << sD.computationalTime << "\n";
    os << "avgSolvingTimeAcc:     " << sD.avgSolvingTimeAcc << "\n";
    os << "performedRealizations: " << sD.performedRealizations << "\n";
    os << "realizationsToDo:      " << sD.realizationsToDo << "\n";
    return os;
}

SolutionStatistics::SolutionStatistics()
{
    clear();
}

void SolutionStatistics::clear()
{
    iterationsCoarse = 0;
    iterationsFine  = 0;
    timeCoarse  = 0;
    timeFine = 0;

    totalTime = 0;
    permInitTime = 0;
    permGenTime = 0;
}

ostream& operator<<(ostream& os, const SolutionStatistics& sS)
{
    os << "iterationsCoarse:" << sS.iterationsCoarse;
    os << "iterationsFine:  " << sS.iterationsFine;
    os << "timeCoarse:      " << sS.timeCoarse;
    os << "timeFine:        " << sS.timeFine;

    os << "totalTime:       " << sS.totalTime;
    os << "permInitTime:    " << sS.permInitTime;
    os << "permGenTime:     " << sS.permGenTime;

    return os;
}


LevelDiffAccumulator::LevelDiffAccumulator()
{
    cotrEmpty();
}

LevelDiffAccumulator::LevelDiffAccumulator(MPI_Comm levelComm) throw(invalid_argument)
{
    cotrEmpty();
    if (levelComm == MPI_COMM_NULL) {
        invalid_argument(string("Trying to Construct LevelDiffAccumulator, with Invalid Comm"));
    }
    this->levelComm = levelComm;
}

LevelDiffAccumulator::LevelDiffAccumulator(MPI_Comm levelComm, int procsPerRealization) throw (invalid_argument)
{
    cotrEmpty();
    if (levelComm == MPI_COMM_NULL) {
        invalid_argument(string("Trying to Construct LevelDiffAccumulator, with Invalid Comm"));
    }
    if (procsPerRealization < 1) {
        invalid_argument(string("procsPerRealization must be possitive"));
    }

    // To Do Comm Dup
    this->levelComm = levelComm;
    this->procsPerRealization = procsPerRealization;

    assignProcessors();
}

LevelDiffAccumulator::LevelDiffAccumulator(MPI_Comm levelComm, Difference * p_diffLevels, int initialRep) throw (invalid_argument)
{
    cotrEmpty();
    assert(dynamic_cast<Difference*> (p_diffLevels));
    if (levelComm == MPI_COMM_NULL) {
        invalid_argument(string("Trying to Construct LevelDiffAccumulator, with Invalid Comm"));
    }

    if (initialRep < 0) {
        invalid_argument(string("initialRep must be positive"));
    }

    this->levelComm = levelComm;
    this->p_diffLevels = p_diffLevels;
    this->realizationsToDo = initialRep;

    assignProcessors();
}


LevelDiffAccumulator::LevelDiffAccumulator(MPI_Comm levelComm, int procsPerRealization, Difference *p_diffLevels, int initialRep) throw (invalid_argument)
{
    cotrEmpty();
    assert(dynamic_cast<Difference*> (p_diffLevels));

    if (levelComm == MPI_COMM_NULL) {
        invalid_argument(string("Trying to Construct LevelDiffAccumulator, with Invalid Comm"));
    }
    if (procsPerRealization < 1) {
        invalid_argument(string("procsPerRealization must be possitive"));
    }

    if (initialRep < 0) {
        invalid_argument(string("initialRep must be positive"));
    }
    this->levelComm = levelComm;
    this->procsPerRealization = procsPerRealization;
    this->p_diffLevels = p_diffLevels;
    this->realizationsToDo = initialRep;

    assignProcessors();
}

LevelDiffAccumulator::LevelDiffAccumulator(int procsPerRealization, Difference *p_diffLevels, int initialRep) throw(invalid_argument)
{
    cotrEmpty();
    assert(dynamic_cast<Difference*> (p_diffLevels));

    if (procsPerRealization < 1) {
        invalid_argument(string("procsPerRealization must be possitive"));
    }

    if (initialRep < 0) {
        invalid_argument(string("initialRep must be positive"));
    }

    this->procsPerRealization = procsPerRealization;
    this->p_diffLevels = p_diffLevels;
    this->realizationsToDo = initialRep;
}

LevelDiffAccumulator::~LevelDiffAccumulator()
{
    if (localDiffComm != MPI_COMM_NULL) {
        MPI_Comm_free(&localDiffComm);
    }
}

// Returns number of performed realizations.
int LevelDiffAccumulator::getPerformed() const
{
    return performedRealizations;
}

// Returns average time for realization on present level.
double LevelDiffAccumulator::getAvgComputeTime() const
{
    if (!static_cast<bool> (performedRealizations)) {
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

// Gets how many realizations is has to be done
int LevelDiffAccumulator::getRealizationsToDo() const
{
    return realizationsToDo;
}

// Sets Comm and assigns procs
void LevelDiffAccumulator::setComm(MPI_Comm newComm)
{

    if (newComm == MPI_COMM_NULL) {
        invalid_argument(string("Trying to Construct LevelDiffAccumulator, with Invalid Comm"));
    }

    if (localDiffComm != MPI_COMM_NULL) {
        MPI_Comm_free(&localDiffComm);
    }

    levelComm = newComm;
    assignProcessors();
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
    assert((p_diffLevels != nullptr) && (localDiffComm != MPI_COMM_NULL));

    loadBalancer();

    diffTime = MPI_Wtime();
    difRealizations = difRealizations - localEval(difRealizations); // Get the actual
    diffTime = MPI_Wtime() - diffTime;  // Local Computational time
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

void LevelDiffAccumulator::setAccParams(const StaticAccRawData& params)
{
    avgSolvingTimeAcc = params.avgSolvingTimeAcc;
    computationalTime = params.computationalTime;
    performedRealizations = params.performedRealizations;
    quantitySum = params.quantitySum;
    quantitySumSqueard =params.quantitySumSqueard;
    realizationsToDo = params.realizationsToDo;
}
void LevelDiffAccumulator::interupComputation()
{
    int myRank, sz;
    MPI_Comm_rank(levelComm, &myRank);
    MPI_Comm_size(levelComm, &sz);

    // Case where master is the only processor
    if (sz == 1) {
        interupSignal = true;
        return;
    }

    // Send signal
    if (myRank == localMasterRank) {
        MPI_Request req = MPI_REQUEST_NULL;
        for (int i = 0; i < sz; ++i) {
            if (myRank != i) {
                MPI_Isend(nullptr, 0, MPI_INT, i, localMasterRank, levelComm, &req);

            }
        }
        // Complete ISend;
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        interupSignal = true;
    }
    else {
        MPI_Recv(nullptr, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, levelComm, MPI_STATUS_IGNORE);
        interupSignal = true;
    }
}

void LevelDiffAccumulator::bindMasterComm(MPI_Comm mastercomm, int masterRank)
{
    // To Do:Add checks
    masterComm = mastercomm;
    localMasterRank = masterRank;
}

void LevelDiffAccumulator::resetInteruped()
{
    interupSignal = false;
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
    interupSignal = false;
}

ostream& operator<<(ostream& os, const LevelDiffAccumulator& lDA)
{
    os << "E[Q] =        " << lDA.mean() << "\n";
    os << "V[Q] =        " << lDA.var() << "\n";
    os << "T[Avg] =      " << lDA.getAvgComputeTime() << "\n";
    os << "T[Total] =    " << lDA.getTotalTime() << "\n";
    os << "Samples done: " << lDA.getPerformed() << "\n";
    os << "Samples rem.: " << lDA.getRealizationsToDo() << "\n";
    return os;
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

    // To Do: Add base seed
    p_diffLevels->init(levelComm, localDiffComm, 0);
}

void LevelDiffAccumulator::synhronize()
{
    constexpr int bufSz = 4;
    double sendBuf[bufSz] = {0};
    double resvBuf[bufSz] = {0};

    double perforemd = 0;

    assert(sizeof(double) == sizeof(MPI_DOUBLE));
    assert(sizeof(sendBuf) / sizeof(double) == sizeof(resvBuf) / sizeof(double));

    sendBuf[0] = diffTime;
    sendBuf[1] = diffQuantitySum;
    sendBuf[2] = diffQuantitySumSqueard;
    sendBuf[3] = static_cast<double>(difRealizations);

    if (localMasterRank != -1) {

        MPI_Reduce(&sendBuf, &resvBuf, sizeof(sendBuf) / sizeof(double), MPI_DOUBLE, MPI_SUM, localMasterRank, levelComm);
    }
    else {
        MPI_Allreduce(&sendBuf, &resvBuf, sizeof(sendBuf) / sizeof(double), MPI_DOUBLE, MPI_SUM, levelComm);
    }

    avgSolvingTimeAcc += (resvBuf[0] / procsPerRealization);
    quantitySum += (resvBuf[1] / procsPerRealization);
    quantitySumSqueard += (resvBuf[2] / procsPerRealization);

    perforemd = static_cast<int> (resvBuf[3] / procsPerRealization);
    performedRealizations += perforemd;
    realizationsToDo -= perforemd;

    computationalTime = (performedRealizations / diffGroupsNumber) * getAvgComputeTime()  + (performedRealizations % diffGroupsNumber) * getAvgComputeTime();
    //computationalTime /= diffGroupsNumber;
}

void LevelDiffAccumulator::loadBalancer()
{
     // Common for all processors, it may need to be splited
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

int LevelDiffAccumulator::localEval(int realizations)
{
    int myRank = 0;
    MPI_Comm_rank(levelComm, &myRank);
    double realizationEval = 0;

    while (realizations > 0 && !interupSignal) {
        if (masterComm != MPI_COMM_NULL && myRank == localMasterRank) {
            int hasMsg = 0;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, masterComm, &hasMsg, MPI_STATUS_IGNORE);
            if (hasMsg) {
                // Propagate Singal
                interupComputation();
                break;
            }
        }
        else {
            int hasMsg = 0;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, levelComm, &hasMsg, MPI_STATUS_IGNORE);
            if (hasMsg) {
                interupComputation();
                break;
            }
        }
        realizationEval = p_diffLevels->eval();
        diffQuantitySum += realizationEval;
        diffQuantitySumSqueard += realizationEval * realizationEval;
        --realizations;
    }
    return realizations;
}

void LevelDiffAccumulator::updateStatistics()
{
    synhronize();

    //Zero local variables
    diffQuantitySum = 0;
    diffQuantitySumSqueard = 0;
    difRealizations = 0;
    diffTime = 0;
    //interupSignal = false;
}

StaticAccRawData LevelDiffAccumulator::getRawData() const
{
    StaticAccRawData toReturn;
    toReturn.avgSolvingTimeAcc = avgSolvingTimeAcc;
    toReturn.computationalTime = computationalTime;
    toReturn.quantitySum = quantitySum;
    toReturn.quantitySumSqueard = quantitySumSqueard;
    toReturn.performedRealizations = performedRealizations;
    toReturn.realizationsToDo = realizationsToDo;
    return toReturn;
}


void LevelDiffAccumulator::cotrEmpty()
{
    // Global Data
    performedRealizations = 0;
    quantitySum = 0;
    quantitySumSqueard = 0;
    computationalTime = 0;
    avgSolvingTimeAcc = 0;
    procsPerRealization = 0;

    p_diffLevels = nullptr;
    levelComm = MPI_COMM_NULL;
    baseSeed = 1;
    localMasterRank = -1;

    // Local Data
    localDiffComm = MPI_COMM_NULL;
    diffGroupsNumber = 0;
    diffCommColor = 0;
    diffQuantitySum = 0;
    difRealizations = 0;
    diffTime = 0;
    interupSignal = false;
}

MultiLevelMonteCarloEst::MultiLevelMonteCarloEst()
:
  mean(0),
  totalVar(0),
  proportionalityConst(0),
  totalComputationalTime(0),
  tolerance(0),
  monteCarloComm(MPI_COMM_NULL),
  levelComm(MPI_COMM_NULL),
  isLeader(false),
  tag(-1)
{}

MultiLevelMonteCarloEst::MultiLevelMonteCarloEst(double tolerance)
:
    mean(0),
    totalVar(0),
    proportionalityConst(0),
    totalComputationalTime(0),
    tolerance(tolerance),
    monteCarloComm(MPI_COMM_NULL),
    levelComm(MPI_COMM_NULL),
    isLeader(false),
    tag(-1)
{}

MultiLevelMonteCarloEst::MultiLevelMonteCarloEst(MPI_Comm MonteCarloComm, double tolerance)
:
    mean(0),
    totalVar(0),
    proportionalityConst(0),
    totalComputationalTime(0),
    tolerance(tolerance),
    monteCarloComm(MonteCarloComm),
    levelComm(MPI_COMM_NULL),
    isLeader(false),
    tag(-1)
{
    MPI_Comm_size(MonteCarloComm, &MCCommSz);
    masterRank = 0;
    //assignProcessors();
}

void MultiLevelMonteCarloEst::eval()
{
    //Temp!!!
    mean = 0;
    totalVar = 0;
    totalComputationalTime = 0;

    bool converged = false;
    bool needSamples = true;

    int myRank = 0;
    MPI_Comm_rank(monteCarloComm, &myRank);
    do{
        do {
            //Redistribute work
            loadBalancer();
            createGroups();
            createComs();

            // Do Parallel Comp. on level 0 paralell tree
            if (levelComm != MPI_COMM_NULL) {
                StatisticalAccumulators[tag].setComm(levelComm);
                // Do parallel Comp. on level 1 parallel tree
                StatisticalAccumulators[tag].bindMasterComm(groupLeadersComm, levelLeaderRank); //Should be split
                StatisticalAccumulators[tag].evalAndUpdate();
            }
            stopSignal();
            gatherComputedData();
            MPI_Bcast(&mlmcResults[0], static_cast<int>(mlmcResults.size()) * sizeof(StaticAccRawData) / sizeof(double), MPI_DOUBLE, masterRank, monteCarloComm);

            freeComms();
            freeGroups();

            // Populate results
            for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
                StatisticalAccumulators[i].setAccParams(mlmcResults[i]);
                StatisticalAccumulators[i].resetInteruped();
            }

            // Check if there are some samples left to do interup
            needSamples = false;
            //double totalRemReal = 0;
            for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
                //totalRemReal = StatisticalAccumulators[i].getRealizationsToDo());
                needSamples |= static_cast<bool> (StatisticalAccumulators[i].getRealizationsToDo());
            }
        }
        while(needSamples);

        // Compute add realizations
        computeProportionallityConstant();
        for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
            computeRealizationsOnDiffLevel(i, 0.8);
        }

        converged = true;
        for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
            converged &= (!static_cast<bool> (StatisticalAccumulators[i].getRealizationsToDo()));
        }

        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

        if (myRank == 0) {
            for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
                std::cout << "===Level (" << i << ")" << "\n";
                std::cout << StatisticalAccumulators[i] << "\n" << std::flush;
                std::cout << std::flush;
            }

            std::cout << "MSE:" << errorMSE() << "\n";
            std::cout << "RMS:" << errorRMS() << "\n";
            std::cout << std::flush;

        }
    } while(!converged);
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
    int requiredAdd = std::max(0.0, estimatedSoFar - doneSoFar);
    StatisticalAccumulators[diffLevel].setRealizationsToDo(requiredAdd);
    return StatisticalAccumulators[diffLevel].getRealizationsToDo();
}

double MultiLevelMonteCarloEst::errorMSE() const
{
    double mse = 0;
    for (int i = 0; i < StatisticalAccumulators.size(); ++i) {
        mse += StatisticalAccumulators[i].var() / StatisticalAccumulators[i].getPerformed();
    }
    return mse;
}

double MultiLevelMonteCarloEst::errorRMS() const
{
    return errorMSE() / tolerance;
}

void MultiLevelMonteCarloEst::addDifference(Difference* diff, int procsPerRealization, int initialRep)
{
    StatisticalAccumulators.push_back(LevelDiffAccumulator(procsPerRealization, diff, initialRep));
    mlmcResults.push_back(StaticAccRawData());
    procsPerLevel.push_back(0);
}

void MultiLevelMonteCarloEst::addDifference(Difference *diff, int initialRep)
{
    addDifference(diff, 1, initialRep);
}

void MultiLevelMonteCarloEst::assignProcessors()
{
    assert(monteCarloComm != MPI_COMM_NULL);
    //MPI_Comm_dup(monteCarloComm, &levelComm);

}

void MultiLevelMonteCarloEst::loadBalancer()
{
    // Branch & Bound Should be used
    int levels = static_cast<int> (StatisticalAccumulators.size());
    int remainProcs = MCCommSz;
    int realizationsOnLevel = 0;
    std::fill(procsPerLevel.begin(), procsPerLevel.end(), 0);

    for (int i = levels - 1; i > -1;  --i) {
        realizationsOnLevel = StatisticalAccumulators[i].getRealizationsToDo();
        if (realizationsOnLevel > 0) {
            if (remainProcs <= realizationsOnLevel) {
                //Put all procs on this level
                procsPerLevel[i] = remainProcs;
                remainProcs = 0;
            }
            else {
                procsPerLevel[i] = realizationsOnLevel;
                remainProcs = remainProcs - realizationsOnLevel;
            }
        }
    }

//    int myRank = 0;
//    MPI_Comm_rank(monteCarloComm, &myRank);
//    if (myRank == 0) {
//        std::cout << "===Procs remain:" << remainProcs << "\n" << std::flush;
//    }

    if (remainProcs > 0) {
        procsPerLevel[0] += remainProcs;
        remainProcs = 0;
    }
}

void MultiLevelMonteCarloEst::stopSignal()
{
    int myRank = 0;
    int leaderCommSz = 0;
    if (isLeader) {
        //Distribute State to leaders
        MPI_Comm_rank(groupLeadersComm, &myRank);
        MPI_Comm_size(groupLeadersComm, &leaderCommSz);
        MPI_Request req = MPI_REQUEST_NULL;

        // Notify master
        if (myRank != masterRank) {
            MPI_Send(nullptr, 0, MPI_INT, masterRank, myRank, groupLeadersComm);
            MPI_Recv(nullptr, 0, MPI_INT, masterRank, MPI_ANY_TAG, groupLeadersComm, MPI_STATUS_IGNORE);
        }
        else {
            // Master
            for (int i =0; i < leaderCommSz; ++i) {
                if (i != masterRank) {
                   MPI_Isend(nullptr, 0, MPI_INT, i, masterRank, groupLeadersComm, &req);
                }
            }

            // Recive acc. from the slaves
            for (int i =0; i < leaderCommSz -1; ++i) {
                MPI_Recv(nullptr, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, groupLeadersComm, MPI_STATUS_IGNORE);
            }

            MPI_Wait(&req, MPI_STATUS_IGNORE);
        }
    }
}

void MultiLevelMonteCarloEst::gatherComputedData()
{
    if (isLeader) {
        int myRank = 0;
        MPI_Comm_rank(groupLeadersComm, &myRank);

        StaticAccRawData computedSoFar = StatisticalAccumulators[tag].getRawData();
        computedSoFar.tag = tag; //Set the tag

        vector<StaticAccRawData> resvBuf(mlmcResults.size(), StaticAccRawData());

        int count = sizeof(StaticAccRawData) / sizeof(double);
        MPI_Gather(&computedSoFar, count, MPI_DOUBLE, &resvBuf[0], count, MPI_DOUBLE, masterRank, groupLeadersComm);
        if (masterRank == myRank) {
            //Sort out the data (Possibale when there is not enough processors)
            int sz = static_cast<int> (resvBuf.size());
            for (int i = 0; i < sz; ++i) {
                if (static_cast<int> (resvBuf[i].tag) != -1) {
                    int level = resvBuf[i].tag;
                    mlmcResults[level] = resvBuf[i];
                    mlmcResults[level].tag = -2;
                }
            }
            for (int i = 0; i < sz; ++i) {
                if (static_cast<int> (mlmcResults[i].tag) != -2) {
                    mlmcResults[i] = StatisticalAccumulators[i].getRawData();
                }
                // return flag to -1
                mlmcResults[i].tag = -1;
            }
        }
    }
}

void MultiLevelMonteCarloEst::createGroups()
{
    // Extract Group that wiil be worked on
    MPI_Group monteCarloCommGrp;
    MPI_Comm_group(monteCarloComm, &monteCarloCommGrp);

    //Extract the sz of Statistical estimators
    int sz = static_cast<int> (StatisticalAccumulators.size());

    // Temp container holders for ranks
    vector<vector<int>> ranksInGroup(sz, vector<int>());
    vector<int> rankIngroupLeaders;
    vector<int> mapRanks;

    // Extract current rank of proc in monteCralo comm
    int myRank = 0;
    MPI_Group_rank(monteCarloCommGrp, &myRank);

    // Put ranks to Groups
    int procsCounter = 0;
    for (int i = 0; i < sz; ++ i) {
        for (int j = 0; j < procsPerLevel[i]; ++j) {
            if(procsCounter == myRank) {
                tag = i;
            }
            ranksInGroup[i].push_back(procsCounter);
            ++procsCounter;
        }
    }

    assert(procsCounter == MCCommSz);

    // Set leaders
    for (int i = 0; i < ranksInGroup.size(); ++ i) {
        if (ranksInGroup[i].size() > 0) {
            rankIngroupLeaders.push_back(ranksInGroup[i][0]);
            mapRanks.push_back(MPI_UNDEFINED);
        }
    }

    myRank == ranksInGroup[tag][0] ? isLeader = true : isLeader = false;
    // Create Groups
    MPI_Group_incl(monteCarloCommGrp, procsPerLevel[tag], &ranksInGroup[tag][0], &levelGroup);
    MPI_Group_incl(monteCarloCommGrp, static_cast<int> (rankIngroupLeaders.size()), &rankIngroupLeaders[0], &groupLeadersGroup);

    // Translate Ranks form MonteCarloComm to levelComm and set levelCommLeaderRank
    MPI_Group_translate_ranks(monteCarloCommGrp, static_cast<int>(rankIngroupLeaders.size()),  &rankIngroupLeaders[0], levelGroup, &mapRanks[0]);
    levelLeaderRank = MPI_UNDEFINED;
    for (int i = 0; i < mapRanks.size(); ++i) {
        if (mapRanks[i] != MPI_UNDEFINED) {
            levelLeaderRank = mapRanks[i];
            break;
        }
    }
}

void MultiLevelMonteCarloEst::createComs()
{
    MPI_Comm_create(monteCarloComm, levelGroup, &levelComm);
    MPI_Comm_create(monteCarloComm, groupLeadersGroup, &groupLeadersComm);
}

void MultiLevelMonteCarloEst::freeComms()
{
    MPI_Comm_free(&levelComm);
    if (isLeader) {
        MPI_Comm_free(&groupLeadersComm);
    }

    levelComm = MPI_COMM_NULL;
    groupLeadersComm = MPI_COMM_NULL;

}

void MultiLevelMonteCarloEst::freeGroups()
{
    MPI_Group_free(&levelGroup);
    if (isLeader) {
        MPI_Group_free(&groupLeadersGroup);
    }
    levelGroup = MPI_GROUP_NULL;
    groupLeadersGroup = MPI_GROUP_NULL;
}
} // END namespace MultiLevelMonteCarlo
