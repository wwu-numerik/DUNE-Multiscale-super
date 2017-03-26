#ifndef __CEGENERATOR__
#define __CEGENERATOR__

#include<array>
#include<vector>
#include<random>
#include<assert.h>
#include<math.h>
#include<mpi.h>
#include<fftw3-mpi.h>


#include"SRUpscale.hh"
#include"../Utils/floatCompare.hh"

template<int DIM, typename X, typename R, typename COR>
class CirculantEmbeding
{
public:
    enum TypeOfGen {Vartex = 1, CellCenter = 0};

    CirculantEmbeding()
        :comm(MPI_COMM_NULL), localBuffsz(nullptr), dispArr(nullptr), fft_plan(nullptr), ifft_plan(nullptr), p_base(nullptr), p_layer(nullptr)
    {
        minimal = 0;
        local_size = local_start = alloc_local = 0;
        seed = 0;
        base_size = 0;
        segment_size = 0;
        step = 0;
        typeOfGen = static_cast<TypeOfGen> (-1);
        leftPermFields = 0;
        permFieldIndex = 0;
    }

    CirculantEmbeding(const CirculantEmbeding&) = delete;
    CirculantEmbeding(const CirculantEmbeding&&) = delete;

    //CirculantEmbeding(MPI_Comm comm, const COR& corrFunction, int directionSz, int seed, TypeOfGen genType = TypeOfGen::CellCenter , double minimal=1e-8)


    void init(MPI_Comm comm, const COR& corrFunction, int directionSz, int seed, int overlap=1, TypeOfGen genType = TypeOfGen::CellCenter, double minimal=1e-8)
    {
        timeForInit = MPI_Wtime();
        // Set up FFTW lib
        fftw_mpi_init();

        int commSz,rank;
        MPI_Comm_size(comm, &commSz);
        MPI_Comm_rank(comm, &rank);
        localBuffsz = new int[commSz];
        dispArr = new int[commSz];

        this->comm = comm;
        corr_functor = corrFunction;
        this->seed = seed;
        this->minimal = minimal;
        typeOfGen = genType;

        //Set seed and random engine
        rand = std::default_random_engine(seed);
        normal = std::normal_distribution<double>(0,1);

        base_size = (1 << (directionSz + genType));
        segment_size = base_size << 1;
        step = 1.0 / base_size;

        // Alloc memory
        alloc_local = fftw_mpi_local_size_2d(segment_size, segment_size, comm, &local_size,  &local_start);
        p_base = fftw_alloc_complex(alloc_local);
        p_layer = fftw_alloc_complex(alloc_local);
        permeabilityField.resize(segment_size * segment_size);

        // Populate boundries of permability, should be one call, or compute it by local information
        MPI_Allgather(&local_size, 1, MPI_INT, localBuffsz, 1, MPI_INT, this->comm);
        MPI_Allgather(&local_start, 1, MPI_INT, dispArr, 1, MPI_INT, this->comm);

        //Update buffers
        for (int i = 0; i < commSz; ++i) {
            localBuffsz[i] *= 2 * segment_size;
            dispArr[i] *= 2 * segment_size;
        }

        //Create plans
        fft_plan = fftw_mpi_plan_dft_2d(segment_size, segment_size, p_base, p_base, this->comm, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT);
        ifft_plan = fftw_mpi_plan_dft_2d(segment_size, segment_size, p_layer, p_layer, this->comm, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);

        setUpBasis();
        fftw_execute(fft_plan);
        setFactor();
        computeLayers();
        timeForInit = MPI_Wtime() - timeForInit;

        timeForGeneration = MPI_Wtime();
        fftw_execute(ifft_plan);
        generateLocalPerm();
        populate();
        timeForGeneration = MPI_Wtime() - timeForGeneration;

        leftPermFields = 1;
        permFieldIndex = 0;
    }

    ~CirculantEmbeding()
    {
        delete []localBuffsz;
        delete []dispArr;
        fftw_free(p_base);
        fftw_free(p_layer);
        fftw_destroy_plan(fft_plan);
        fftw_destroy_plan(ifft_plan);
    }

    void create()
    {
        perpNewEval();
    }

    double getGeneratonTime()
    {
        return this->timeForGeneration;
    }

    double getInitTime()
    {
        return this->timeForInit;
    }

    R operator() (const X& x) const
    {
        if (x[0] < 0 || x[0] > 1 || x[1] < 0 || x[1] > 1) {
            throw("Bad point recived");
        }
        X rem = x - calcLowerLeftPoint(x);
        double dataPoints[4] = {0, 0, 0, 0};
        double yArr[2] = {0, 0};
        // On Vartex
        if (floatCompare<double>(rem[0], 0.0, 1e-15) && floatCompare<double>(rem[1], 0.0, 1e-15)) {
            return permeabilityField[mapPointToIndex(calcLowerLeftPoint(x))][permFieldIndex];
        }

        X tmp = calcLowerLeftPoint(x);
        if (tmp[0] == 1) {
            tmp[0] -= step;
        }
        if (tmp[1] == 1) {
            tmp[1] -= step;
        }
        dataPoints[0] = permeabilityField[mapPointToIndex(tmp)][permFieldIndex]; //Q11
        tmp[0] += step;
        dataPoints[1] = permeabilityField[mapPointToIndex(tmp)][permFieldIndex]; //Q21
        tmp[1] += step;
        dataPoints[3] = permeabilityField[mapPointToIndex(tmp)][permFieldIndex]; //Q22
        tmp[0] -= step;
        dataPoints[2] = permeabilityField[mapPointToIndex(tmp)][permFieldIndex]; //Q12

        double k;
        X lLPoint = calcLowerLeftPoint(x);
        auto aprox1D = [] (double x1, double x2, double f1, double f2, double x) -> double { return (x2 - x) * f1 + (x - x1) * f2;};
        yArr[0] = aprox1D (lLPoint[0], lLPoint[0] + step, dataPoints[0], dataPoints[1], x[0]);
        yArr[1] = aprox1D (lLPoint[0], lLPoint[0] + step, dataPoints[2], dataPoints[3], x[0]);
        k = aprox1D (lLPoint[1], lLPoint[1] + step, yArr[1], yArr[0], x[1]) / (step * step);
        return std::max(k, minimal);
    }

    // Only for centers of cells
    R getSchiftPerm(const X& pointCoordinates, const bool useRenormalization,const int times, const SimplifiedRenorm::FlowDirection direction  = SimplifiedRenorm::FlowDirection::X) const
    {
        return getRenromPerm(times, pointCoordinates);
    }
    const std::vector<std::array<double, 2>> * getPermData()
    {
        return &permeabilityField;
    }

    void perpNewEval(int newSeed = -1)
    {

        if (newSeed != -1) {
            seed = newSeed;
        }
        // Switch to next ready perm
        if (leftPermFields > permFieldIndex) {
            ++permFieldIndex;
        }
        else {
            timeForGeneration += MPI_Wtime();
            regenerate();
            timeForGeneration = MPI_Wtime() - timeForGeneration;
        }
    }

private:
    // Maps 2D to 1D representation of the array
    int mapArrInd (int i, int j)
    {
        if (i * segment_size + j > segment_size * segment_size) {
            throw("Bad mapping between point and perm index");
        }
        return i * segment_size + j;
    }

    X calcLowerLeftPoint (const X& x) const
    {
        X LLPoint;
        LLPoint[0] = static_cast<int>(x[0] / step) * step;
        LLPoint[1] = static_cast<int>(x[1] / step) * step;
        return LLPoint;
    }

    int mapPointToIndex (const X& x) const
    {
        int indX = x[0] / step;
        int indY = x[1] / step;
        return indX * segment_size + indY;
    }

    void setUpBasis()
    {
        auto calcPoint = [this] (int i, int j) -> X
        {
            int rank = 0;
            MPI_Comm_rank(comm, &rank);
            X point;
            point[0] = (i > (base_size)) ? ((segment_size) - i) * step : i * step;
            point[1] = (j > (base_size)) ? ((segment_size) - j) * step : j * step;
            return point;
        };

        for (int i = 0; i < local_size; ++i) {
            for (int j = 0; j < segment_size; ++j) {
                //Setup coorelation
                p_base[mapArrInd(i, j)][0] = corr_functor(calcPoint(i + local_start, j));
                p_base[mapArrInd(i, j)][1] = 0;
            }
        }
    }

    void setFactor()
    {
        double factor = 1.0 / (segment_size * segment_size);
        for(int i=0; i< local_size; ++i) {
            for (int j = 0; j < segment_size; ++j) {
                p_base[mapArrInd(i, j)][0] = std::sqrt(factor * std::abs(p_base[mapArrInd(i, j)][0]));
                p_base[mapArrInd(i, j)][1] = 0;
            }
        }
    }

    void computeLayers()
    {
        // Multiply coefficients with N(0,1)-random numbers
        for(int i=0; i< local_size; ++i) {
            for (int j = 0; j < segment_size; ++j) {
                p_layer[mapArrInd(i, j)][0] = p_base[mapArrInd(i, j)][0] * normal(rand);
                p_layer[mapArrInd(i, j)][1] = p_base[mapArrInd(i, j)][0] * normal(rand);
            }
        }
    }

    void generateLocalPerm()
    {
        for(int i=0; i< local_size; ++i) {
            for (int j = 0; j < segment_size; ++j) {
                p_layer[mapArrInd(i, j)][0] = std::exp(p_layer[mapArrInd(i, j)][0]);
                p_layer[mapArrInd(i, j)][1] = std::exp(p_layer[mapArrInd(i, j)][1]);
            }
        }
    }

    void populate()
    {

        int workerID, sz;
        MPI_Comm_size(comm,&sz);
        MPI_Comm_rank(comm, &workerID);

        auto chekPerm = [this](double val) -> double { return val < minimal ? minimal : val; };
        for(int i= 0; i< local_size; ++i) {
            for (int j = 0; j < segment_size; ++j) {
                permeabilityField[mapArrInd(i, j)][0] = chekPerm(p_layer[mapArrInd(i, j)][0]);
                permeabilityField[mapArrInd(i, j)][1] = chekPerm(p_layer[mapArrInd(i, j)][1]);
            }
        }
        //To Do: Set up localdistribution
        MPI_Allgatherv(&permeabilityField[0], localBuffsz[workerID], MPI_DOUBLE, &permeabilityField[0], localBuffsz, dispArr, MPI_DOUBLE, comm);
    }

    void regenerate()
    {
        permFieldIndex = 0;

        computeLayers();
        fftw_execute(ifft_plan);
        generateLocalPerm();
        populate();
    }

    double getRenromPerm(int times, X cellCenter, const SimplifiedRenorm::FlowDirection direction = SimplifiedRenorm::FlowDirection::X) const
    {
        double uPperms[4] = {0,0,0,0};
        double translate = 0;
        if (times == 1) {
            translate = step / 2.0;
        }
        else if (times >1) {
            translate = step;
        }
        X schiftToCenter (times * translate);
        X topR = cellCenter + schiftToCenter;
        X bottomL = cellCenter - schiftToCenter;

        schiftToCenter[0] *= -1;
        X topL = cellCenter + schiftToCenter;
        X bottomR = cellCenter - schiftToCenter;

        if (times == 0) {
            double upPerm = permeabilityField[mapPointToIndex(cellCenter)][permFieldIndex];
            if (upPerm < 0) {
                throw("SR faild");
            }
            return upPerm;
        }
        else {
            // Shift to center of old cells
            uPperms[0] = getRenromPerm(times - 1, bottomL, direction);
            uPperms[1] = getRenromPerm(times - 1, bottomR, direction);
            uPperms[2] = getRenromPerm(times - 1, topL, direction);
            uPperms[3] = getRenromPerm(times - 1, topR, direction);
        }
        return SimplifiedRenorm::upscale2D(uPperms,4, direction);
    }


    MPI_Comm comm;                           // MPI-communicator
    int* localBuffsz;                        // array of local portions sizes
    int* dispArr;                            // array of displasments
    double minimal;                          // minimal permeability
    ptrdiff_t local_size;                    // local num. of segments along 1st dim
    ptrdiff_t local_start;                   // 1st local index along 1st dim.
    ptrdiff_t alloc_local;                   // local allocation
    fftw_plan fft_plan;                      // plan for fast fourier transform
    fftw_plan ifft_plan;                     // plan for inverse fft
    COR corr_functor;                        // correlation as function of x-y

    int seed;
    std::default_random_engine rand;          // random number generator
    std::normal_distribution<double>  normal; // normal distribution

    fftw_complex * p_base;                    // basis functions
    fftw_complex * p_layer;                   // local layers of perm
    int base_size;                            // base sz
    int segment_size;
    double step;                              // Step
    TypeOfGen typeOfGen;                      // Vartex or cell centered
    int leftPermFields;                       // number of generated perm fileds
    int permFieldIndex;                       // current perm field
    typedef std::array<double, 2 > permPair;
    std::vector<permPair> permeabilityField; // full perm field

    double timeForGeneration;
    double timeForInit;
};

#endif
