#include "dune-random-diff.hh"

int main(int argc, char** argv)
{

    // Parameters
     int    pFast;       // processors per fast solve
     int    pMedium;     // processors per medium solve
     int    pSlow;       // processors per accurate solve
     int    log2seg;     // number of segments per dimension
     int    overlap;     // overlap in domain decomposition
     double corrLen;     // correlation length
     double sigma;       // sigma of correlation distribution
     double tol;         // tolerance computing mean flux
     int    breaks;      // number of breaks for reestimating variances
     int levels;         // number of levels

     try{
         Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
         if (!(argc!=11 || argc != 12)) {
             if(helper.rank()==0) {
                 std::cout
                         << "usage: dune_mlmc <pFast> <pMedium> <pSlow> <log2seg> <overlap> "
                         <<"<corrLen> <sigma> <tol> <breaks> <levels> [results file]\n";
             }
             return 1;
         }


         // Read parameters
         pFast   = atoi(argv[1]);

         pMedium = atoi(argv[2]);
         pSlow   = atoi(argv[3]);
         log2seg = atoi(argv[4]);
         overlap = atoi(argv[5]);
         corrLen = atof(argv[6]);
         sigma   = atof(argv[7]);
         tol     = atof(argv[8]);
         breaks  = atoi(argv[9]);
         levels  = atoi(argv[10]);

         int procs[3] = {1, 1, 1}; // Temp !!!
         procs[0] = pFast;
         procs[1] = pMedium;
         procs[2] = pSlow;
         std::string resultsFile;
         //Some checks should be done
         if (argc == 12 ) {
             resultsFile.assign(argv[11]);
         }

        CroaseDiff coarsediff(log2seg - 2, 1, corrLen,sigma, 2, log2seg);
        TelescopicDiff middle(log2seg - 1, 1, corrLen, sigma, 1, log2seg);
        TelescopicDiff slow(log2seg,1, corrLen,sigma,0,log2seg);
        MultiLevelMonteCarloEst estimator(MPI_COMM_WORLD, tol);
        estimator.addDifference(&coarsediff,300);
        estimator.addDifference(&middle,12);
        estimator.addDifference(&slow,6);
        int myRank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        MPI_Barrier(MPI_COMM_WORLD);
        double time = MPI_Wtime();
        estimator.eval();
        time = MPI_Wtime() - time;
        MPI_Barrier(MPI_COMM_WORLD);
        if (myRank == 0) {
            std::cout << "Compute time[s]:" << time;
        }

        //RenromDump test(log2seg, overlap, corrLen, sigma, 0, log2seg);
        //test.init(MPI_COMM_WORLD, MPI_COMM_WORLD);
        //test.eval();

//        SREstimator mlmcEst(levels, log2seg, overlap, corrLen, sigma, tol, breaks, procs);
//        mlmcEst.expectation(MPI_COMM_WORLD, resultsFile);

//        RandomDiff::RandomDiffTypes::CorrelationType corr(corrLen, sigma);
//        RandomDiff::RandomDiffTypes::PermeabilityType perm;
//        MPI_Comm local;
//        int rank;
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//        MPI_Comm_split(MPI_COMM_WORLD,rank, rank, &local);

    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::flush << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;

    }
    catch (const char * msg) {
        std::cerr << msg << std::flush << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl << std::flush;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
