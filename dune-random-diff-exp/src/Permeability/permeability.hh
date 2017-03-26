#ifndef PERMEABILITY_HH
#define PERMEABILITY_HH

#include<exception>
#include<array>
#include<vector>
#include<random>
#include<iostream>
#include<assert.h>
#include<cmath>
#include<mpi.h>

#include<fftw3-mpi.h>
#include"SRUpscale.hh"
#include"../Utils/floatCompare.hh"

#include<fstream>

template< int DIM, typename X, typename R, typename COR >
class Permeability
{
    typedef std::array<double,2> complex;
    typedef std::vector<complex> cvec;
    typedef std::array<int,DIM> iarr;

public:
    Permeability()
    {
        _fft = nullptr;
        _ifft = nullptr;
    }

    Permeability( const Permeability&) = delete;

    ~Permeability() {
        fftw_destroy_plan(_fft);
        fftw_destroy_plan(_ifft);
    }

    Permeability(MPI_Comm comm, const COR& corr, int log2Seg, int seed, int overlap=1, double minimal=1e-8)
    {
        _fft = nullptr;
        _ifft = nullptr;
        init(comm, corr, log2Seg, seed, overlap, minimal);
    }

    void init(MPI_Comm comm, const COR& corr, int log2Seg, int seed, int overlap=1, double minimal=1e-8)
    {

        permInitTime = MPI_Wtime();

        // Initialize
        _comm = comm;
        _corr = corr;
        _N    = 1<<log2Seg;
        _rand = std::default_random_engine(seed);
        _normal = std::normal_distribution<double>(0,1);
        _overlap = overlap;
        _minimal = minimal;
        _part = 1;

        if (DIM == 3) {
            calcDirectionArray3D(translationVecArray3D, 8);
        }

        if(_fft != nullptr) {
            fftw_destroy_plan(_fft);
            fftw_destroy_plan(_ifft);
        }

        ptrdiff_t local_size;
        double h = 1.0 / _N;
        int _2N = 2 * _N;
        X x;

        MPI_Comm_size(_comm,&_nProc);
        MPI_Comm_rank(_comm,&_iProc);
        fftw_mpi_init();
        setRange();


        // Create basis functions in 2D
        if(DIM == 2) {
            local_size = fftw_mpi_local_size_2d(_2N, _2N, _comm, &_n0, &_start);

            assert(_n0==_2N / _nProc);
            _base  = cvec(local_size);

            fftw_complex* base = (fftw_complex*) _base.data();
            _layer = cvec(local_size);

            fftw_complex *layer = (fftw_complex*) _layer.data();

            _fft   = fftw_mpi_plan_dft_2d(_2N,_2N,base,base,_comm,FFTW_FORWARD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT);
            _ifft  = fftw_mpi_plan_dft_2d(_2N,_2N,layer,layer,_comm,FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);

            int p=0;
            for(int i=0; i<_n0; ++i) {
                x[0] = (i +_start > _N) ? (_2N -_start - i) * h : (_start + i) * h;
                for(int j=0; j<_2N; ++j, ++p) {
                    x[1] = (j > _N) ? ( _2N - j) * h : j * h;
                    base[p][0] = _corr(x);
                    base[p][1] = 0;
                }
            }
            fftw_execute(_fft);

            double factor = 1.0 / (_2N*_2N);
            for(int i=0; i<_n0 * _2N; ++i) {
                //assert(base[i]>=0); //TODO: clarify
                base[i][0] = sqrt(factor * fabs(base[i][0]) );
                base[i][1] = 0;
            }
        }


        // Create basis functions in 3D
        else if(DIM == 3) {
            local_size = fftw_mpi_local_size_3d(_2N,_2N,_2N,_comm,&_n0,&_start);
            assert(_n0 == _2N / _nProc);
            _base  = cvec(local_size);
            fftw_complex* base = (fftw_complex*) _base.data();
            _layer = cvec(local_size);
            fftw_complex* layer = (fftw_complex*) _layer.data();
            _fft   = fftw_mpi_plan_dft_3d(_2N,_2N,_2N,base,base,_comm,FFTW_FORWARD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT);
            _ifft  = fftw_mpi_plan_dft_3d(_2N,_2N,_2N,layer,layer,_comm,FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);

            int p=0;
            for(int i=0; i <_n0; ++i) {
                x[0] = (i+_start>_N)? (_2N-_start-i) * h : (_start+i) * h;
                for(int j=0; j < _2N; ++j) {
                    x[1] = (j>_N) ? (_2N-j) * h : j * h;
                    for(int k=0; k< _2N; ++k, ++p) {
                        x[2] = ( k > _N) ? ( _2N - k) * h : k * h;
                        base[p][0] = _corr(x);
                        base[p][1] = 0;
                    }
                }
            }

            fftw_execute(_fft);

            double factor = 1.0 / (_2N*_2N*_2N);
            for(int i=0; i <_n0 * _2N * _2N; ++i) {
                //assert(base[i]>=0); //TODO: clarify
                base[i][0] = sqrt(factor * fabs(base[i][0]));
                base[i][1] = 0;
            }
        }
        else {
            throw("Only dimensions 2 and 3 are implemented.");
        }
        permInitTime = MPI_Wtime() - permInitTime;
    }

    /// Compute number of processors per dimension
    /// \param   nProc       total number of processors
    /// \param   procPerDim  number of processors per dimension
    /// \returns false, if number of processors is not a power of 2
    static bool partition(int nProc, int* procPerDim)
    {
        double h = std::log2(nProc);
        int    log2P = int(h);

        //        if(log2P!=h) return false;
        auto isPowerOf2 = [](int x) { return (x != 0) && ((x & (x-1)) == 0); };
        if(isPowerOf2(nProc)) return false;

        int nAll = log2P/DIM;
        int nMore = log2P - nAll*DIM;
        for(int i=0; i<DIM; ++i)
            procPerDim[i] = 1 << (nAll + (i<nMore));
        return true;
    }

    /// Create random permeability field from basis.
    void create() {

        // toggle between real and imaginary part of permeability field.
        // Recompute only if next real part is requested.
        _part = 1-_part;
        if(_part==1) {
            permGenTime = 0;
            return;
        }
        permGenTime = MPI_Wtime();
        ptrdiff_t i,n;
        n = (DIM==2)? _n0*2*_N : _n0*4*_N*_N;
//********************************************
        //dummy_1.clear();
        //dummy_2.clear();
//********************************************
        // Multiply coefficients with N(0,1)-random numbers and apply IFFT
        for(i=0; i<n; ++i) {
            double random_p_1 = _normal(_rand);
            double random_p_2 = _normal(_rand);
            _layer[i][0] =  random_p_1 * _base[i][0];
            _layer[i][1] = random_p_2 * _base[i][0];
//********************************************
            dummy_1.push_back(random_p_1);
            dummy_2.push_back(random_p_2);
//********************************************
        }
        fftw_execute(_ifft);

        // k = exp(Y)
        int rank = 0;
        MPI_Comm_rank(this->_comm,&rank);
        for(i=0; i<n; ++i) {
            _layer[i][0] = exp(_layer[i][0]);
            _layer[i][1] = exp(_layer[i][1]);
        }

        // Redistribute permeability field
        redistribute();
        permGenTime = MPI_Wtime() - permGenTime;
    }

    /// Evaluate permeability field.
    /// \param x position
    /// \return permeability at x
    R operator() (const X& x) const {

        for (int i = 0; i < DIM; ++i) {
            if (x[i] < 0 || x[i] > 1) {
                throw("Bad point recived");
            }
        }

        int cell = 0;
        int c00,c01,c10,c11;
        double t[DIM];
        double k;
        for(int i=0; i<DIM; ++i) {
            double p = x[i]*_N;
            if(p<_iMin[i] || p>_iMax[i])
                throw "outside";
            p    -= _iMin[i];
            int j = floor(p);
            t[i]  = p-j;
            cell  = (_size[i])*cell + j;
        }
        if(DIM==2) {
             c00 = cell;
             c01 = c00 + 1;
             c10 = c00 + (_size[1]);
             c11 = c10 + 1;

            k = ( (1-t[1]) * _perm[c00][_part] + t[1] * _perm[c01][_part]) * (1-t[0])
                    +((1-t[1]) * _perm[c10][_part] + t[1] * _perm[c11][_part])*t[0];
        }
        else if (DIM==3) {

            int c000 = cell;
            int c001 = c000+1;
            int c010 = c000+_size[2];
            int c011 = c010+1;
            int c100 = c000+_size[1]*_size[2];
            int c101 = c100+1;
            int c110 = c100+_size[2];
            int c111 = c110+1;
            k = ( ( (1-t[2])*_perm[c000][_part]
                  + t[2]*_perm[c001][_part])*(1-t[1])
                    + ((1-t[2])*_perm[c010][_part]
                    + t[2]*_perm[c011][_part])*t[1] ) * (1-t[0])
                    + ( ( (1-t[2])*_perm[c100][_part]
                        + t[2]*_perm[c101][_part])*(1-t[1])
                    + ((1-t[2])*_perm[c110][_part]
                    + t[2]*_perm[c111][_part])*t[1] ) * t[0];
        }
        else {
            std::cerr << "Not implemented, yet.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

//********************************************
        if (k != k) {
            std::cerr << "NaN in Perm.\n" << "point:" << "(" << x[0] << "," << x[1] << ") "
                      << c00 << " : "
                      << c01 << " : "
                      << c10 << " : "
                      << c11 << " Perm_mem_sz: " << _perm.size()  << " size: " << _size[1];
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
//********************************************
        return std::max(k,_minimal);
    }

// **********************************************************************************************
    double getGeneratonTime()
    {
        return permGenTime;
    }

    double getInitTime()
    {
        return permInitTime;
    }

    R getMaxPermVal() const
    {
        double max = 0;
        for (int i = 0; i < this->_perm.size(); ++ i) {
            if (max < _perm[i][_part]) {
                max = _perm[i][_part];
            }
        }
        return max;
    }

    R getMinPermVal() const
    {
        double min = 1e+100;
        for (int i = 0; i < this->_perm.size(); ++ i) {
            if (min > _perm[i][_part]) {
                min = _perm[i][_part];
            }
        }
        return min;
    }

    // Only for centers of cells
    R getSchiftPerm(const X& pointCoordinates, const bool useRenormalization ,const int timesSr = 0, const SimplifiedRenorm::FlowDirection direction = SimplifiedRenorm::FlowDirection::X ) const
    {
        if (timesSr < 0) {
            throw ("Bad nernorm params");
        }
        double upPerm = 0;
        X schift ((1.0 / _N) / 2);
        if (timesSr == 0) {
            // Go to the shifted point
            return _perm[mapPointToInx(pointCoordinates - schift, _N)][_part];
        }
        else {
           if (DIM == 2) {
               upPerm = getRenromPerm(timesSr - 1, pointCoordinates, direction);
               if (upPerm < 0) {
                   throw("SR faild");
               }
           }
           if (DIM == 3) {
               upPerm = getRVECell3d(timesSr - 1, pointCoordinates, direction);
               if (upPerm < 0) {
                   throw("SR faild");
               }
           }
           return upPerm;
        }
    }
// **********************************************************************************************
private: 
// **********************************************************************************************
    // Maps point to ind of perm
    int mapPointToInx (const X x, int N) const
    {
        double step = 1.0 / N;
        int ind[DIM];
        for (int i = 0; i < DIM; ++i) {
            ind[i] = x[i] / step;
        }
        int permInd = 0;
        if (DIM == 3) {
            permInd = ind[2] * (N + 1) * (N + 1) + ind[1] * (N + 1) + ind[0];
        }
        else if (DIM == 2) {
            permInd = ind[1] * (N + 1) + ind[0];
        }

        if (_perm.size() < permInd) {
            throw ("Ind out of bounds");
        }
        return permInd;
    }

    double getRenromPerm(int times, X cellCenter, const SimplifiedRenorm::FlowDirection direction = SimplifiedRenorm::FlowDirection::X) const
    {
        double upPerm = 0;
        double uPperms[4] = {0,0,0,0};
        if (times == 0) {
           X schiftToCenter ((1.0 / (2 * _N)) );
           X topR = cellCenter + schiftToCenter;
           X bottomL = cellCenter - schiftToCenter;

           schiftToCenter[0] *= -1;
           X topL = cellCenter + schiftToCenter;
           X bottomR = cellCenter - schiftToCenter;
           for (int i = 0; i < DIM; ++i) {
              if (topL[i] < 0 || topR[i] < 0 || bottomL[i] < 0 || bottomR[i] < 0 ||
                      topL[i] > 1 || topR[i] > 1 || bottomL[i] > 1 || bottomR[i] > 1) {
                  throw("Bad schift");
              }
          }
           X schift ((1.0 / (2 * _N)) );
           upPerm = SimplifiedRenorm::upscale(
                       _perm[mapPointToInx(topL - schift, _N)][_part],
                       _perm[mapPointToInx(topR - schift, _N)][_part],
                       _perm[mapPointToInx(bottomL - schift, _N)][_part],
                       _perm[mapPointToInx(bottomR - schift, _N)][_part],
                       direction
                       );

           if (upPerm < 0) {
               throw("SR faild");
           }
           return upPerm;
        }
        else {
            // Shift to center of old cells
            X schiftToCenter ( times * (1.0 / _N));
            X topR = cellCenter + schiftToCenter;
            X bottomL = cellCenter - schiftToCenter;

            schiftToCenter[0] *= -1;
            X topL = cellCenter + schiftToCenter;
            X bottomR = cellCenter - schiftToCenter;

            uPperms[0] = getRenromPerm(times - 1, bottomL, direction);
            uPperms[1] = getRenromPerm(times - 1, bottomR, direction);
            uPperms[2] = getRenromPerm(times - 1, topL, direction);
            uPperms[3] = getRenromPerm(times - 1, topR, direction);
        }
        return SimplifiedRenorm::upscale2D(uPperms,4, direction);
    }

// To Do: Has to be checked
// ******************************
    double getRVECell3d (int times, X cellCenter, const SimplifiedRenorm::FlowDirection direction = SimplifiedRenorm::FlowDirection::X) const
    {
        double shiftCoord[8] = {0};
        if (times == 0) {
            for (int i = 0; i < 8; ++i) {
                shiftCoord[i] = _perm[mapPointToInx(cellCenter + translationVecArray3D[i],_N)][_part];
            }
        }
        else {
            for (int i = 0; i < 8; ++i) {
                X tmp;
                for (int j = 0; j < tmp.size(); ++j) {
                    tmp[j] =  cellCenter[j] + translationVecArray3D[i][j] * times * 2;
                }
                shiftCoord[i] = getRVECell3d(times - 1, tmp  ,direction);
            }
        }
        return SimplifiedRenorm::upscale3D(shiftCoord, 8);
    }

    void calcDirectionArray3D(X* transformArray, int sz)
    {
        int vecInx = 0;
        for (int i = -1; i < 2; i += 2) {
            for (int j = -1; j < 2; j += 2) {
                for (int k = -1; k < 2; k += 2) {
                    transformArray[vecInx][0] = static_cast<double> (k) / (2 * _N); // X coord
                    transformArray[vecInx][1] = static_cast<double> (j) / (2 * _N); // Y coord
                    transformArray[vecInx][2] = static_cast<double> (i) / (2 * _N); // Z coord
                    ++vecInx;
                }
            }
        }
    }
// **********************************************************************************************

    /// Compute coordinates of processor on cartesian processor grid
    /// \param iProc       processor index
    /// \param procPerDim  number of processors per dimension
    /// \param pos         position on processor grid
    static void gridPosition(int iProc, const int* procPerDim, int* pos) {
        int off[DIM];
        off[0] = 1;
        for(int i=1; i<DIM; ++i) {
            off[i] = off[i-1]*procPerDim[i-1];
        }
        for(int i=DIM-1; i>=0; --i) {
            pos[i] = iProc/off[i];
            iProc -= pos[i]*off[i];
        }
    }

    /// Compute index range of local part of permeability field
    void setRange() {
        int size = 1;
        int procPerDim[DIM];
        int pos[DIM];
        partition(_nProc,procPerDim);
        gridPosition(_iProc,procPerDim,pos);
        for(int i=0; i<DIM; ++i) {
            int len = _N / procPerDim[i];
            _iMin[i] = std::max(0,len*pos[i]-_overlap);
            _iMax[i] = std::min(_N,len*(pos[i]+1)+_overlap);
            _size[i] = _iMax[i]-_iMin[i]+1;
            size    *= _size[i];
        }
        _perm = cvec(size);
    }

    /// Redistribute permeability from layers to blocks.
    void redistribute() {
        MPI_Win win;
        int nCplx = sizeof(complex);
        MPI_Win_create(_layer.data(), nCplx*_layer.size(), nCplx, MPI_INFO_NULL, _comm, &win);
        MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE, win);
        complex *dest = _perm.data();
        if(DIM==2) {
            int stepDest = _size[1];
            int stepSource = 2*_N;
            for(int i=_iMin[0]; i<=_iMax[0]; ++i) {
                int rankSource = i/_n0;
                int source = (i%_n0)*stepSource + _iMin[1];
                MPI_Get(dest, stepDest, MPI_DOUBLE_COMPLEX, rankSource, (MPI_Aint) source, stepDest, MPI_DOUBLE_COMPLEX, win);
                dest+= stepDest;
            }
        }
        else if(DIM==3) {
            int stepDest = _size[2];
            int stepSource = 4*_N*_N;
            for(int i=_iMin[0]; i<=_iMax[0]; ++i) {
                int rankSource = i/_n0;
                int source = (i%_n0)*stepSource + 2*_N*_iMin[1] + _iMin[2];
                for(int j=_iMin[1]; j<=_iMax[1]; j++) {
                    MPI_Get(dest, stepDest, MPI_DOUBLE_COMPLEX, rankSource, (MPI_Aint) source, stepDest, MPI_DOUBLE_COMPLEX, win);
                    dest+= stepDest;
                    source += 2*_N;
                }
            }
        }
        else {
            throw("Only dimensions 2 and 3 are implemented.");
        }
        MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED, win);
        MPI_Win_free(&win);
    }




    //--- Members ------------------------------------------------------------
    int _iProc;                               ///< index of processor
    int _nProc;                               ///< total number of processors
    int _overlap;                             ///< overlap of subgrids
    int _N;                                   ///< number of segments on [0,1]
    double _minimal;                          ///< minimal permeability
    ptrdiff_t _n0;                            ///< local num. of segments along 1st dim
    ptrdiff_t _start;                         ///< 1st local index along 1st dim.
    COR _corr;                                ///< correlation as function of x-y
    MPI_Comm _comm;                           ///< MPI-communicator
    fftw_plan _fft;                           ///< plan for fast fourier transform
    fftw_plan _ifft;                          ///< plan for inverse fft
    cvec _base;                               ///< base functions
    cvec _layer;                              ///< local layer of permeability field
    cvec _perm;                               ///< permeability field
    int  _part;                               ///< 0: real, 1: imaginary part of _perm
    iarr _iMin;                               ///< minimal global indices of subgrid
    iarr _iMax;                               ///< maximal global indices of subgrid
    iarr _size;                               ///< number of nodes per dim in subgrid
    std::default_random_engine _rand;         ///< random number generator
    std::normal_distribution<double> _normal; ///< normal distribution
    std::vector<double> dummy_1, dummy_2;


    //**************************
    X translationVecArray3D[8];
    double permGenTime;
    double permInitTime;
};

#endif
