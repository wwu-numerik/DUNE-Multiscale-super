#ifndef MLMC_CORRELATION
#define MLMC_CORRELATION
/// Class defining exponential correlation function.
/// \tparam DIM   dimension
/// \tparam X     point type
/// \tparam R     scalar type
///
/// \author jan.mohring@itwm.fraunhofer.de
/// \date 2014
template<int DIM, typename X, typename R>
class CorrelationOneNorm {
public:
  /// Constructor
  /// \param corrLen   correlation length
  /// \param sigma     standard deviation     
  CorrelationOneNorm(R corrLen=0.1, R sigma=1.0)
    : _corrLen(corrLen), _sigma2(sigma*sigma) {}

  CorrelationOneNorm(const CorrelationOneNorm& old) {
    _corrLen = old._corrLen;
    _sigma2  = old._sigma2;
  }

  /// Evaluation
  /// \param d   difference of points to take corretation of
  R operator() (X d) const {
      R corr = 0;
      for (int i = 0; i < DIM; ++i) {
          corr += std::abs(d[i]);
      }
      return _sigma2 * exp(-corr / _corrLen);
  }

private:
  R _corrLen;   //< correlation length
  R _sigma2;    //< standard deviation
};


template<int DIM, typename X, typename R>
class CorrelationTwoNorm {
public:
  /// Constructor
  /// \param corrLen   correlation length
  /// \param sigma     standard deviation
  CorrelationTwoNorm(R corrLen=0.1, R sigma=1.0)
    : _corrLen(corrLen), _sigma2(sigma*sigma) {}

  CorrelationTwoNorm(const CorrelationTwoNorm& old) {
    _corrLen = old._corrLen;
    _sigma2  = old._sigma2;
  }

  R operator() (X d) const {
      R sumX2 = 0;
      for(int i=0; i<DIM; ++i) {
         sumX2 += d[i]*d[i];
      }
      return _sigma2 * exp(-sqrt(sumX2)/_corrLen);
    }


private:
  R _corrLen;   //< correlation length
  R _sigma2;    //< standard deviation
};

#endif
