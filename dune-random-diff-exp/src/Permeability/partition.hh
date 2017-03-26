#ifndef MLMC_PARTITION
#define MLMC_PARTITION
/// Class for partitioning domain. The number of processors
/// must be an integer power of 2.
/// \tparam DIM   dimension
#include<dune/grid/yaspgrid.hh>
#include"../Utils/floatCompare.hh"
template <int DIM>
class Partition : public Dune::YLoadBalance<DIM>
{
public:

    //typedef Dune::FieldVector<int, DIM>  iTupel;
    typedef std::array<int, DIM> iTupel;
    virtual ~Partition() {}

    /// Find partitioning.
    /// \param size   numbers of segments per dimension of global grid
    /// \param P      number of processors
    /// \param dims   numbers of processors per dimension

    virtual void loadbalance(const iTupel& size, int P, iTupel& dims) const
    {

        double h = std::log2(P);

        int    log2P = int(h);
        if (floatCompare<double>(log2P, h, 1e-15))
//        if(log2P!=h)
            DUNE_THROW(Dune::ParallelError,
                       "Number of processors must be power of 2.");
        int nAll = log2P/DIM;
        int nMore = log2P - nAll*DIM;
        for(int i=0; i<DIM; ++i)
            dims[i] = 1 << (nAll + (i<nMore));
    }
};
#endif
