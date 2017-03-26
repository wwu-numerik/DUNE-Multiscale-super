#ifndef __CELLCENTERFV__
#define __CELLCENTERFV__

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

#include"../boundaryCond.hh"

namespace RandomDiff {
using namespace Dune;

template<class PermeabilityGenerator, class BoundaryCond>
class FVLOperator:
        public PDELab::NumericalJacobianApplyVolume<FVLOperator<PermeabilityGenerator, BoundaryCond> >,
        public PDELab::NumericalJacobianVolume<FVLOperator<PermeabilityGenerator, BoundaryCond> >,
        public PDELab::NumericalJacobianApplySkeleton<FVLOperator<PermeabilityGenerator, BoundaryCond> >,
        public PDELab::NumericalJacobianSkeleton<FVLOperator<PermeabilityGenerator, BoundaryCond> >,
        public PDELab::NumericalJacobianApplyBoundary<FVLOperator<PermeabilityGenerator, BoundaryCond> >,
        public PDELab::NumericalJacobianBoundary<FVLOperator<PermeabilityGenerator, BoundaryCond> >,
        public PDELab::FullSkeletonPattern,
        public PDELab::FullVolumePattern,
        public PDELab::LocalOperatorDefaultFlags
{
public:
    // pattern assembly flags
    enum { doPatternVolume = true };
    enum { doPatternSkeleton = true };

    // residual assembly flags
    enum { doAlphaVolume  = false };
    enum { doAlphaSkeleton  = true };
    enum { doAlphaBoundary  = true };

    FVLOperator(const PermeabilityGenerator& reff_permGen, const BoundaryCond& reff_bCond, bool useRenrom, int renormTimes)
        :reff_perm(reff_permGen), reff_bondary(reff_bCond), useRenormalization(useRenrom), renormTimes(renormTimes) {}

    FVLOperator(const FVLOperator&) = delete;
    FVLOperator(const FVLOperator&&) = delete;
    FVLOperator& operator=(const FVLOperator& other) = delete;
    FVLOperator& operator=(const FVLOperator&& other) = delete;

    // No reaction term, emtpy here
    //template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    //void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const {}

    // skeleton integral depending on test and ansatz functions, each face is only visited ONCE!
    template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_skeleton (const IG& ig, const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s, const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, R& r_s, R& r_n) const
    {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef FieldVector<DF, IG::dimension> Point; // Point on the grid
        typedef typename IG::Geometry FaceGeometry;

        Point insideGlobal, outsideGlobal;
        DF distance, faceVolume, kAvg, k1, k2;

        // Distance between cell centers in global coordinates, between self (cell) and neighbor (cell)
        insideGlobal = ig.inside().geometry().center();
        outsideGlobal = ig.outside().geometry().center();
        insideGlobal -= outsideGlobal;
        distance = insideGlobal.two_norm();

        faceVolume = ig.geometry().volume();


        //kAvg = reff_perm.getSchiftPerm(ig.inside().geometry().center(), useRenormalization, renormTimes);
        //kAvg += reff_perm.getSchiftPerm(ig.outside().geometry().center(), useRenormalization, renormTimes);
        //kAvg /= 2;

        k1 = reff_perm.getSchiftPerm(ig.inside().geometry().center(), useRenormalization, renormTimes);
        k2 = reff_perm.getSchiftPerm(ig.outside().geometry().center(), useRenormalization, renormTimes);

        kAvg = 2 * (k1 * k2) / (k1 + k2);


        r_s.accumulate(lfsu_s, 0, -(x_n(lfsu_n,0) - x_s(lfsu_s,0)) * kAvg *  faceVolume / distance);
        r_n.accumulate(lfsu_n, 0, (x_n(lfsu_n,0) - x_s(lfsu_s,0)) * kAvg * faceVolume / distance);
    }

    // skeleton integral depending on test and ansatz functions
    template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s, R& r_s) const
    {
        // Domain, Range, Point types, etc...
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef FieldVector<DF, IG::dimension> Point; // Point on the grid
        typedef typename IG::Geometry FaceGeometry;

        DF bondaryVal, distance, faceVolume, kAvg;
        Point insideGlobal;
        Point faceCenter = ig.geometry().center();

        if (reff_bondary.isDirchleBoundry (faceCenter)) {
            bondaryVal = reff_bondary.evalDirichleBoundary(faceCenter);
            kAvg = reff_perm.getSchiftPerm(ig.inside().geometry().center(), useRenormalization, renormTimes);

            insideGlobal = ig.inside().geometry().center();
            insideGlobal -= faceCenter;
            distance = insideGlobal.two_norm();
            faceVolume = ig.geometry().volume();

            r_s.accumulate(lfsu_s,0, -(bondaryVal - x_s(lfsu_s,0)) * kAvg * faceVolume / distance);
        }
        else
        {
            r_s.accumulate(lfsu_s,0, 0);
        }
    }

private:
    const PermeabilityGenerator& reff_perm;
    const BoundaryCond& reff_bondary;
    bool useRenormalization;
    int renormTimes;
};

} // END NAMESPACE RANDOMDIFF
#endif // __CELLCENTERFV__
