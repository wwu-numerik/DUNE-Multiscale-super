#ifndef __DUNE_RANDOM_DIFF__
#define __DUNE_RANDOM_DIFF__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "randomDiff/randomDiffProblem.hh"
#include "Mlmc/level.hh"
#include "Mlmc/mlmc.hh"


template< class CroaseType, class TelescopicType>
class MlmcComposer
{
public:
    static const int MAX_POSSIBALE_LEVELS = 4;
    MlmcComposer(int levels, int fineLog2seg, int overlap, double corrLen, double sigma, double tol, int breaks, const int* procsPerLevel)
        :levels(levels), fineLog2seg(fineLog2seg), overlap(overlap), corrLen(corrLen), sigma(sigma), tol(tol), breaks(breaks)
    {
        if (MAX_POSSIBALE_LEVELS < levels) {
            throw("exceed max Possibale Levels threshold");
        }

        for (int i = 0; i < MAX_POSSIBALE_LEVELS; ++i) {
            diffObj[i] = nullptr;
        }

        diffObj[0] = new CroaseType (fineLog2seg - levels, overlap, corrLen, sigma, levels, fineLog2seg);
        this->procsPerLevel[0] = procsPerLevel[0];
        mlmc.addDifference(*diffObj[0], this->procsPerLevel[0]);

        for (int i = 1; i < levels + 1; ++i) {
            //RTTI check...
            diffObj[i] = new TelescopicType(fineLog2seg - (levels - i), overlap, corrLen, sigma, (levels - i), fineLog2seg);
            this->procsPerLevel[i] = procsPerLevel[i];
            mlmc.addDifference(*diffObj[i], procsPerLevel[i]);
        }

    }

    void expectation(MPI_Comm comm = MPI_COMM_WORLD, std::string filename = "")
    {
        std::ofstream results;
        //Not the best sol
        if (filename.size()) {
            mlmc.setOutput(results);
            int rank;
            MPI_Comm_rank(comm, &rank);
            results.open(filename);
            if (!results.is_open()) {
                throw ("Fatal error, results file cannot be opend, terminating the program");
            }

            if (rank == 0) {
                results << "<levels> <fineGridSz> <overlap> <corrLen> <sigma> <tol> <braks>\n"
                        << levels + 1 << "\t" << fineLog2seg << "\t" << overlap << "\t"
                        << corrLen << "\t" << sigma << "\t" << tol << "\t" << breaks << "\n";

                for (int i = 0; i < levels + 1; ++i) {
                    results << "Data for <level>: " << i + 1 << "\n";
                    results << "<permSz> " << diffObj[i]->getPermLog2Seg() << "\t"
                            << "<renorm times fine>" << diffObj[i]->getRenormTimes() << "\t"
                            << "<gridSz>" << diffObj[i]->getLog2Seg() << "\n";
                }
                results << "\n";
            }
        }

        mlmc.expectation(tol,breaks, comm);

        if (filename.size()) {
            results.close();
        }

    }

    ~MlmcComposer() {
        for (int i = 0; i < levels + 1; ++i) {
            delete diffObj[i];
        }
    }

 private:
    int procsPerLevel[MAX_POSSIBALE_LEVELS];
    int levels;
    int fineLog2seg;
    int overlap;
    double corrLen;
    double sigma;
    double tol;
    int breaks;

    Difference* diffObj[MAX_POSSIBALE_LEVELS];
    MLMC mlmc;
};

typedef MlmcComposer<CroaseDiff, TelescopicDiff> SREstimator;

#endif
