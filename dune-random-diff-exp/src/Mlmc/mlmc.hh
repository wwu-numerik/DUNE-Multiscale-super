#ifndef MLMC_H
#define MLMC_H

#include <mpi.h>
#include <fftw3-mpi.h>

#include<assert.h>
#include<exception>
#include<iostream>
#include<string>

#include<vector>
#include<array>
#include<math.h>
#include<algorithm>


namespace MultiLevelMonteCarlo {
//-------------------------------------
using std::vector;
using std::array;

using std::exception;
using std::string;

using std::invalid_argument;

using std::ofstream;
using std::ostream;
//-------------------------------------


void check(bool condition, const char* message);
bool isPowerOf2(int x);

// Raw Data Structure for StatAcc
struct StaticAccRawData
{
    StaticAccRawData();

    double quantitySum;
    double quantitySumSqueard;
    double computationalTime;
    double avgSolvingTimeAcc;
    double performedRealizations;
    double realizationsToDo;
    double tag;

    friend ostream& operator<<(ostream& os, const StaticAccRawData& sD);
};


// Raw Data Structure to hold Statistics for the Slover
struct SolutionStatistics
{
    SolutionStatistics();

    void clear();

    double iterationsCoarse;
    double iterationsFine;
    double timeCoarse;
    double timeFine;

    double totalTime;
    double permInitTime;
    double permGenTime;

    friend ostream& operator<<(ostream& os, const SolutionStatistics& sS);
};



class Difference
{
public:
    Difference() {}

    virtual void init(MPI_Comm global, MPI_Comm local, int seed) = 0;
    virtual double eval() = 0;
    virtual const int getLog2Seg() = 0;
    virtual const int getOverlap() = 0;
    virtual const double getCorrLen() = 0;
    virtual const double getSgima() = 0;
    virtual const double getRenormTimes() = 0;
    virtual const int getPermLog2Seg() = 0;

    struct SolutionStatistics solStat;

    virtual ~Difference(){}
};


class LevelDiffAccumulator
{
public:
    LevelDiffAccumulator();
    LevelDiffAccumulator(MPI_Comm levelComm) throw(invalid_argument);
    LevelDiffAccumulator(MPI_Comm levelComm, int procsPerRealization) throw(invalid_argument);
    LevelDiffAccumulator(MPI_Comm levelComm, Difference *p_diffLevels, int initialRep) throw(invalid_argument);
    LevelDiffAccumulator(MPI_Comm levelComm, int procsPerRealization, Difference *p_diffLevels, int initialRep) throw(invalid_argument);

    LevelDiffAccumulator(int procsPerRealization, Difference *p_diffLevels, int initialRep) throw(invalid_argument);
// TO DO *************************
//    Shollow Copy
//    LevelDiffAccumulator(const LevelDiffAccumulator& rhs);
//    LevelDiffAccumulator(const LevelDiffAccumulator&& rhs);
//    LevelDiffAccumulator operator= (const LevelDiffAccumulator& rhs);
//    LevelDiffAccumulator operator= (const LevelDiffAccumulator&& rhs);
//**********************************

    virtual ~LevelDiffAccumulator();            // Free Communicatos

    int getPerformed() const;                   // Returns number of performed realizations.
    double getAvgComputeTime() const;           // Returns average time for realization on present level.
    double getTotalTime() const;                // Returns total time for realization on present level.
    int getBaseSeed() const;                    // Gets base seed;
    int getProcsPerRealization() const;         // Gets defined procs per realization
    int getRealizationsToDo() const;            // Gets how many realizations is to be done

    void setBaseSeed(int seed);                 // Set base seed
    void setRealizationsToDo(int realizations); // Sets realizations to do over level comm
    void setComm(MPI_Comm newComm);             // Resets the Comm

    double mean() const;                        // Returns mean value of realizations on present level.
    double var() const;                         // Returns empirical variance of realizations on present level.
    void eval();                                // Evaluate difference.
    void evalAndUpdate();                       // Eval and update statistics

    void updateStatistics();                    // Updates values on demand, synhronize processes
    StaticAccRawData getRawData() const;      // Extract raw data from acc


    int groupsCount() const;                    // Return the number of the groups
    int procs() const;                          // Return number of processors per group.

    void interupComputation();                  // Interup computation
    void setAccParams(const StaticAccRawData &params); // Set Params explicitly
    void bindMasterComm(MPI_Comm masterComm, int masterRank);   // Bind a Master comm
    void resetInteruped();                      // Set to false again
    void clearResults();                        // Clear all results so far

    friend ostream& operator<<(ostream& os, const LevelDiffAccumulator& lDA);

protected:
    virtual void assignProcessors();            // Assigns processors for realizations
    virtual void synhronize();                  // Synhronize statistics
    virtual void loadBalancer();                // Default load balancer
    virtual int localEval(int realizations);   // Local evaluation

private:

    void cotrEmpty();                           // Construct empty obj

    int performedRealizations;                  // Number of performed realizations
    int realizationsToDo;                       // Number of realizations to do

    double quantitySum;                         // Sum of results so far
    double quantitySumSqueard;                  // Sum of results^2 so far
    double computationalTime;                   // Total computation time so far
    double avgSolvingTimeAcc;                   // Accumulator for solving time of a realization

    Difference * p_diffLevels;                  // Pointer to Difference object

    int procsPerRealization;                    // Number of processors per realization
    int diffGroupsNumber;                       // Number of diff groups
    MPI_Comm levelComm;                         // Communicator to work on
    MPI_Comm localDiffComm;                     // Local diff communicator
    MPI_Comm masterComm;                        // Local diff communicator

    int diffCommColor;                          // Local color (id) of diff comm
    int interupSignal;                          // Signal to stop computations
    int interupMsg;                             // Signal to stop computations
    int localMasterRank;

    double diffQuantitySum;                     // Local QuantitySum
    double diffQuantitySumSqueard;              // Local QuantitySum^2
    int difRealizations;                        // Local number of repettiosions needed to be done
    double diffTime;                            // Local time

    int baseSeed;                               // Base seed to pass to diffLevels
};

class MultiLevelMonteCarloEst
{
public:
    MultiLevelMonteCarloEst();
    MultiLevelMonteCarloEst(double tolerance);
    MultiLevelMonteCarloEst(MPI_Comm MonteCarloComm, double tolerance);


    void eval();
    void addDifference(Difference *diff, int procsPerRealization, int initialRep);
    void addDifference(Difference *diff, int initialRep);

    vector<StaticAccRawData> mlmcResults;

protected:
    void assignProcessors();
    void loadBalancer();
    void stopSignal();
    void gatherComputedData();

    void createGroups();
    void createComs();
    void freeGroups();
    void freeComms();

private:
    double computeProportionallityConstant();
    double computeRealizationsOnDiffLevel(int diffLevel, double factor);

    double errorRMS() const;
    double errorMSE() const;

    double proportionalityConst;
    vector<LevelDiffAccumulator> StatisticalAccumulators;

    double mean;
    double totalVar;
    double totalComputationalTime;
    double tolerance;

    MPI_Comm monteCarloComm;
    MPI_Comm levelComm;
    MPI_Comm groupLeadersComm;

    MPI_Group levelGroup;
    MPI_Group groupLeadersGroup;

    bool isLeader;
    int masterRank;
    int tag;
    int MCCommSz;
    int levelLeaderRank;
    vector<int> procsPerLevel;
    vector<int> groupLeaders;

};
//************************************************************************************************************

//****************************************************************************************************
/// Class for doing statistics on a given level.
class Level {
public:

    /// Constructs empty level.
    Level() : _n(0), _N(0), _p(0), _g(0), _sumX(0), _sumX2(0), _T(0),
        _diff(NULL), _masters(MPI_COMM_NULL) {}

    /// Constructs level from Difference object.
    /// \param diff     Difference object
    /// \param minProc  minimal number of processors needed to compute solution
    Level(Difference& diff, int minProc=1) : _n(0), _N(0), _p(minProc), _g(0),
        _sumX(0), _sumX2(0), _T(0), _diff(&diff), _masters(MPI_COMM_NULL){
        check(minProc>0,"minProc must be positive.");
    }


    /// Assigns communicator and group
    /// \param world   global communicator
    void assignProcessors(MPI_Comm world=MPI_COMM_WORLD, int offset = 0) {
        check(_diff!=NULL,"No Difference object set.");
        int size, rank;
        int range[1][3];
        MPI_Group all, group, masters;
        MPI_Comm_size(world, &size);
        MPI_Comm_rank(world, &rank);
        MPI_Comm_group(world, &all);
        range[0][0] = rank - rank%_p;
        range[0][1] = range[0][0]+_p-1;
        range[0][2] = 1;
        MPI_Group_range_incl(all, 1, range, &group);
        MPI_Comm_create(world, group, &_group);
        _g = size/_p;
        _diff->init(world,_group, offset * (_g + 1));
        range[0][0] = 0;
        range[0][1] = size-1;
        range[0][2] = _p;
        MPI_Group_range_incl(all, 1, range, &masters);
        MPI_Comm_create(world, masters, &_masters);
        _isMaster = (rank%_p)==0;
    }

    /// Sets number of repetitions.
    /// \param N   number of repetitions
    void setRepetitions(double N) {
        _N = ceil(N);
    }

    /// Returns time of next break.
    /// \param iBreak   index of next break;
    /// \param nBreak   total number of breaks
    double nextBreakTime(int iBreak, int nBreak) const {

        double dT = (_T>0)? (_N-_n)*(_T/(_n*_g)) : 0.1*nBreak;  // TODO
        if(iBreak<nBreak)
            dT = dT * iBreak / nBreak;
        return MPI_Wtime() + dT;
    }

    /// Returns communicator of masters
    MPI_Comm getMasters() const { return _masters; }

    /// Indicates master
    bool isMaster() const { return _isMaster; }

    /// Returns communicator of masters
    MPI_Comm getGroup() const { return _group; }

    /// Clears statistics.
    void clear() { _n=0; _sumX=0; _sumX2=0; _T=0; }

    /// Updates statistics.
    void update(int n, double sumX, double sumX2, double T) {
        _n+=n; _sumX+=sumX; _sumX2+=sumX2; _T+=T;
    }

    /// Returns number of performed realizations.
    int performed() const { return _n; }

    /// Returns number of required realizations.
    int required() const { return _N; }

    /// Returns mean value of realizations on present level.
    double mean() const { return _sumX/_n; }

    /// Returns empirical variance of realizations on present level.
    double var()  const { return _sumX2/_n - mean()*mean(); }

    /// Returns average time for realization on present level.
    double time() const { return _T/_n; }

    /// Returns number of groups on present level.
    int groups() const { return _g; }

    /// Return number of processors per group.
    int procs() const { return _p; }

    /// Evaluate difference.
    double eval() { return _diff->eval(); }

    const SolutionStatistics getSolStats() { return _diff->solStat;}

private:
    int  _n;        ///< number of performed realizations
    int  _N;        ///< number of required realizations
    int  _p;        ///< number of processors per realization
    int  _g;        ///< number of groups acting in parallel
    double _sumX;   ///< sum of results so far
    double _sumX2;  ///< sum of squared results so far
    double _T;      ///< total computation time so far
    Difference *const _diff;  ///< pointer to Difference object
    MPI_Comm _group;          ///< communicator of group
    MPI_Comm _masters;        ///< communicator of masters
    bool _isMaster; ///< flag indicating master
};
//****************************************************************************************************


/// Environment for computing expected values by the
/// Multi Level Monte Carlo Method.
class MLMC {
public:

    /// Default constructor
    MLMC() { _output = &(std::cout);}

    /// Set output stream
    void setOutput(std::ostream &output) {
        _output = &output;
    }

    /// Adds difference object.
    /// \param diff     Difference object
    /// \param minProc  minimal number of processors needed to compute solution
    void addDifference(Difference& diff, int minProc) {
        check(isPowerOf2(minProc),"minProc must be power of 2.");
        _level.push_back(Level(diff, minProc));
    }

    /// Computes optimal number of repetitions per level .
    /// \param tol  absolute tolerance of expected value
    void setRepetitions(double tol) {
        double alpha = 0;
        int n = static_cast<int>(_level.size());
        for(int i=0; i<n; ++i)
            alpha += sqrt( _level[i].time() * _level[i].var() );
        alpha /= tol*tol;
        for(int i=0; i<n; ++i)
            _level[i].setRepetitions( alpha * sqrt(_level[i].var() / _level[i].time() ) );
    }

    /// Computes the expected value of a scalar random variable up to a given
    /// tolerance by the Multi Level Monte Carlo approach.
    /// \param tol    absolute tolerance of expected value
    /// \param nBreak number of breaks for recomputing required realizations
    /// \param world  world communicator
    double expectation(double tol, int nBreak, MPI_Comm world=MPI_COMM_WORLD) {
        double duration = MPI_Wtime();  // TODO remove


        int nLevel = static_cast<int>(_level.size());
        int rank,size;
        MPI_Comm_rank(world,&rank);
        MPI_Comm_size(world,&size);
        check(isPowerOf2(size),"Number of processors must be power of 2.");

        // Check input.
        check(nLevel>0,"No levels set.");
        check(tol>0, "tol must be positive.");
        check(nBreak>1, "nBreak must be greater than 1.");

        // Create groups and communicators for different levels.
        for(int i=0; i<nLevel; ++i) {
            _level[i].setRepetitions(8); // TODO
            _level[i].assignProcessors(world, i);
        }

        // Loop over breaks and levels
        bool ready = false;
        for(int iBreak=1; !ready; ++iBreak) {
            ready = true;
            for(int iLevel=0; iLevel<nLevel; ++iLevel) {

                // Moments of groups.
                Level* l = &(_level[iLevel]);
                if(l->required() <= l->performed())
                    continue;

                ready = false;
                double grpSumX = 0;
                double grpSumX2 = 0;
                double grpT = MPI_Wtime();
                int grpN = 0;
                double tBreak = l->nextBreakTime(iBreak,nBreak);
                while(MPI_Wtime()<tBreak || grpN*l->groups()<6) { //TODO 6
                    double x = l->eval();
                    grpSumX += x;
                    grpSumX2 += x*x;
                    ++grpN;
                }
                grpT = MPI_Wtime()-grpT;

                // Accumulate group results.
                double sumX = 0;
                double sumX2 = 0;
                double T = 0;
                int n = 0;
                MPI_Comm masters = l->getMasters();
                MPI_Comm group = l->getGroup();
                // Accumulate over group masters
                if(l->isMaster()) {
                    MPI_Allreduce(&grpSumX,&sumX,1,MPI_DOUBLE,MPI_SUM,masters);
                    MPI_Allreduce(&grpSumX2,&sumX2,1,MPI_DOUBLE,MPI_SUM,masters);
                    MPI_Allreduce(&grpT,&T,1,MPI_DOUBLE,MPI_SUM,masters);
                    MPI_Allreduce(&grpN,&n,1,MPI_INT,MPI_SUM,masters);
                }
                // Distribute master results over groups
                MPI_Bcast(&sumX,1,MPI_DOUBLE,0,group);
                MPI_Bcast(&sumX2,1,MPI_DOUBLE,0,group);
                MPI_Bcast(&T,1,MPI_DOUBLE,0,group);
                MPI_Bcast(&n,1,MPI_INT,0,group);
                // Update moments.
                l->update(n,sumX,sumX2,T);
            }
            if(ready) break;
            // Compute optimal number of repetitions per algorithm.
            setRepetitions(tol);

            if(rank==0) {
                (*_output) << "Break " << iBreak << "\n";
                for(int i=0; i<nLevel; ++i)
                    (*_output) << "  lev" << i << ":"
                               << " N " << _level[i].required()
                               << " n " << _level[i].performed()
                               << " t " << _level[i].time()
                               << " v " << _level[i].var()
                               << "\n";
                (*_output).flush();

                double e = 0;
                for(int i=0; i<nLevel; ++i)
                    e += _level[i].mean();

            }
            //******************
            //break;
            //******************
        }

        // Compute expected value from mean values of differences.
        double e = 0;
        for(int i=0; i<nLevel; ++i) {
            e += _level[i].mean();
            if (MPI_COMM_NULL != _level[i].getMasters()) {
                _levelStat.push_back(updateStatistics(_level[i].getMasters(), i));
            }
        }

        if(rank==0) {
            duration = MPI_Wtime()-duration;
            (*_output) << "Expected value: " << e << " Duration: " << duration << " s\n";

            (*_output) << "===Statistical infromation\n";
            for (int i = 0; i < _levelStat.size(); ++i) {
                (*_output) << "==Level: " << i << "\n";
                (*_output) << "iterationsCoarse:\t" << _levelStat[i].iterationsCoarse << "\n";
                (*_output) << "iterationsFine:  \t" << _levelStat[i].iterationsFine << "\n";
                (*_output) << "permGenTime:     \t" << _levelStat[i].permGenTime << "\n";
                (*_output) << "permInit:        \t" << _levelStat[i].permInitTime << "\n";
                (*_output) << "timeCoarse:      \t" << _levelStat[i].timeCoarse << "\n";
                (*_output) << "timeFine:        \t" << _levelStat[i].timeFine << "\n";
                (*_output) << "totalTime:       \t" << _levelStat[i].totalTime << "\n";
            }
        }
        return e;
    }

private:
    std::vector<Level> _level;           ///< vector of levels
    std::ostream *_output;               ///< output stream

    std::vector<SolutionStatistics> _levelStat;

    SolutionStatistics updateStatistics(MPI_Comm comm, int level)
    {
        SolutionStatistics avgStat;
        avgStat.clear();
        SolutionStatistics levelStat = _level[level].getSolStats();
        MPI_Allreduce(&levelStat.iterationsCoarse, &avgStat.iterationsCoarse, 1, MPI_DOUBLE,MPI_SUM, comm);
        MPI_Allreduce(&levelStat.iterationsFine, &avgStat.iterationsFine, 1, MPI_DOUBLE,MPI_SUM, comm);
        MPI_Allreduce(&levelStat.permGenTime, &avgStat.permGenTime, 1, MPI_DOUBLE,MPI_SUM, comm);
        MPI_Allreduce(&levelStat.permInitTime, &avgStat.permInitTime, 1, MPI_DOUBLE,MPI_SUM, comm);
        MPI_Allreduce(&levelStat.timeCoarse, &avgStat.timeCoarse, 1, MPI_DOUBLE,MPI_SUM, comm);
        MPI_Allreduce(&levelStat.timeFine, &avgStat.timeFine, 1, MPI_DOUBLE,MPI_SUM, comm);
        MPI_Allreduce(&levelStat.totalTime, &avgStat.totalTime, 1, MPI_DOUBLE,MPI_SUM, comm);

        avgStat.iterationsCoarse /= _level[level].performed();
        avgStat.iterationsFine /= _level[level].performed();
        avgStat.permGenTime /= _level[level].performed();
        avgStat.permInitTime /= _level[level].groups();
        avgStat.timeCoarse /= _level[level].performed();
        avgStat.timeFine /= _level[level].performed();
        avgStat.totalTime /= _level[level].performed();

        if (_level[level].performed() == _level[level].groups()) {
            avgStat.permGenTime /= 2;
        }
        return avgStat;
    }
};
}


#endif

