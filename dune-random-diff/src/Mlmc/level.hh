#ifndef __LEVEL__
#define __LEVEL__

#include "mlmc.hh"
#include "../Permeability/SRUpscale.hh"
#include "../randomDiff/randomDiffProblem.hh"

using namespace MultiLevelMonteCarlo;
using namespace RandomDiff;

class CroaseDiff: public Difference
{
public:
    CroaseDiff(int log2seg, int overlap, double corrLen, double sigma, int renormTimes, int permLog2Seg)
        : log2seg(log2seg), overlap(overlap), corrLen(corrLen), sigma(sigma), renormTimes(renormTimes), permLog2Seg(permLog2Seg),
          corr(corrLen, sigma),
          perm(),
          estimatorCroase(log2seg, overlap, perm, renormTimes)
    {}

    void init(MPI_Comm world, MPI_Comm group, int seed)
    {
        this->group = group;

        int rank;
        MPI_Comm_rank(world,&rank);
        perm.init(group, corr, permLog2Seg, seed + rank + 1 , overlap);

        Difference::solStat.permInitTime += perm.getInitTime();

        estimatorCroase.init(group);
    }

    double eval()
    {
        perm.create();
        Difference::solStat.permGenTime += perm.getGeneratonTime();

        estimatorCroase.applaySolver();
        double flow = estimatorCroase.surfaceFlow();

        Difference::solStat.iterationsCoarse += estimatorCroase.getIterTime();
        Difference::solStat.timeCoarse += estimatorCroase.getSolTime();
        return flow;
    }

    inline const int getLog2Seg() {return log2seg;}
    inline const int getOverlap() {return overlap;}
    inline const double getCorrLen() {return corrLen;}
    inline const double getSgima() {return sigma;}
    inline const double getRenormTimes() {return renormTimes;}
    inline const int getPermLog2Seg() {return permLog2Seg;}

private:
    int log2seg;
    int overlap;

    double corrLen;
    double sigma;
    int renormTimes;
    int permLog2Seg;
    RandomDiffTypes::CorrelationType corr;
    RandomDiffTypes::PermeabilityType perm;
    RandomDiffProblem estimatorCroase;

    MPI_Comm group;

};


class TelescopicDiff: public Difference
{
public:
    TelescopicDiff(int log2seg, int overlap, double corrLen, double sigma, int renormTimes, int permLog2Seg)
        : log2seg(log2seg), overlap(overlap), corrLen(corrLen), sigma(sigma), renormTimes(renormTimes), permLog2Seg(permLog2Seg),
          corr(corrLen, sigma),
          perm(),
          estimatorFine(log2seg, 2 * overlap, perm, renormTimes),
          estimatorCroase(log2seg - 1, 2 * overlap, perm, renormTimes + 1)
    {}

    void init(MPI_Comm world, MPI_Comm group, int seed)
    {
        this->group = group;
        int rank;
        MPI_Comm_rank(world,&rank);
        perm.init(group, corr, permLog2Seg, seed + rank + 1 , 2 * overlap);

        Difference::solStat.permInitTime += perm.getInitTime();

        estimatorFine.init(group);
        estimatorCroase.init(group);

    }

    double eval()
    {
        perm.create();
        Difference::solStat.permGenTime += perm.getGeneratonTime();

        estimatorFine.applaySolver();
        estimatorCroase.applaySolver();
        double flow = estimatorFine.surfaceFlow();
        flow -= estimatorCroase.surfaceFlow();

        Difference::solStat.iterationsCoarse += estimatorCroase.getIterTime();
        Difference::solStat.timeCoarse += estimatorCroase.getSolTime();

        Difference::solStat.iterationsFine += estimatorFine.getIterTime();
        Difference::solStat.timeFine += estimatorFine.getSolTime();

        Difference::solStat.totalTime = Difference::solStat.timeFine + Difference::solStat.timeCoarse;

        return flow;
    }

    inline const int getLog2Seg() {return log2seg;}
    inline const int getOverlap() {return overlap;}
    inline const double getCorrLen() {return corrLen;}
    inline const double getSgima() {return sigma;}
    inline const double getRenormTimes() {return renormTimes;}
    inline const int getPermLog2Seg() {return permLog2Seg;}


private:
    int log2seg;
    int overlap;

    double corrLen;
    double sigma;
    int renormTimes;
    int permLog2Seg;

    RandomDiffTypes::CorrelationType corr;
    RandomDiffTypes::PermeabilityType perm;
    RandomDiffProblem estimatorFine;
    RandomDiffProblem estimatorCroase;

    MPI_Comm group;
};

class RenromDump: public Difference
{
public:
    RenromDump(int log2seg, int overlap, double corrLen, double sigma, int renormTimes, int permLog2Seg)
        : log2seg(log2seg), overlap(overlap), corrLen(corrLen), sigma(sigma), renormTimes(renormTimes), permLog2Seg(permLog2Seg),
          corr(corrLen, sigma),
          perm(),
          estimatorFine(this->log2seg, 2 * overlap, perm, renormTimes),
          estimatorCroase(log2seg - 1, 2 * overlap, perm, renormTimes + 1),
          estimatorCroasest(log2seg - 2, 2 * overlap, perm, renormTimes + 2)
    {}

    void init(MPI_Comm world, MPI_Comm group)
    {
        this->group = group;
        int rank;
        MPI_Comm_rank(world,&rank);
        perm.init(group, corr, permLog2Seg, 1, 2 * overlap);
        estimatorFine.init(group);
        estimatorCroase.init(group);
        estimatorCroasest.init(group);
    }

    double eval()
    {

        perm.create();
        estimatorFine.writePermToVtkFormat("PermField_r0");
        estimatorCroase.writePermToVtkFormat("PermField_r1");
        estimatorCroasest.writePermToVtkFormat("PermField_r2");
        estimatorFine.applaySolver();
        estimatorCroase.applaySolver();
        estimatorCroasest.applaySolver();

        estimatorFine.writeSolToVtKFormat("Fine");
        estimatorCroase.writeSolToVtKFormat("Coarse");
        estimatorCroasest.writeSolToVtKFormat("Coarsest");
        return -1;
    }

    inline const int getLog2Seg() {return log2seg;}
    inline const int getOverlap() {return overlap;}
    inline const double getCorrLen() {return corrLen;}
    inline const double getSgima() {return sigma;}
    inline const double getRenormTimes() {return renormTimes;}
    inline const int getPermLog2Seg() {return permLog2Seg;}


private:
    int log2seg;
    int overlap;

    double corrLen;
    double sigma;
    int renormTimes;
    int permLog2Seg;

    RandomDiffTypes::CorrelationType corr;
    RandomDiffTypes::PermeabilityType perm;
    RandomDiffProblem estimatorFine;
    RandomDiffProblem estimatorCroase;
    RandomDiffProblem estimatorCroasest;
    MPI_Comm group;
};


#endif
