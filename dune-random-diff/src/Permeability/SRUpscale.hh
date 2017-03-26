#ifndef __SRUpscale__
#define __SRUpscale__

#include<math.h>
#include<memory>
// If error occur returns -1
template<typename PermValType = double>
class RVERenormalization
{
public:
    RVERenormalization() = delete;
    RVERenormalization(const RVERenormalization& ) = delete;
    RVERenormalization(const RVERenormalization&& ) = delete;
    RVERenormalization operator=(const RVERenormalization& ) = delete;
    RVERenormalization operator=(const RVERenormalization&& ) = delete;

    enum FlowDirection{X, Y, Z};


    static PermValType upscale(PermValType topL, PermValType topR, PermValType bottomL, PermValType bottomR, FlowDirection direction = X)
    {
        PermValType toArray[4] = {0};
        toArray[0] = bottomL;
        toArray[1] = bottomR;
        toArray[2] = topL;
        toArray[3] = topR;
        return upscale2D(toArray, 4, direction);
    }

    static PermValType upscale(const PermValType* zOLayer, FlowDirection direction = X)
    {
        PermValType toArray[4] = {0};
        for (int i = 0; i < 4; ++i) {
            toArray[i] = zOLayer[i];
        }
        return upscale2D(toArray, 4, direction);
    }

    static PermValType upscale(PermValType z0_topL, PermValType z0_topR, PermValType z0_bottomL, PermValType z0_bottomR,
                               PermValType z1_topL, PermValType z1_topR, PermValType z1_bottomL, PermValType z1_bottomR,
                               FlowDirection direction = X)
    {
        PermValType toArray[8] = {0};
        toArray[0] = z0_bottomL;
        toArray[1] = z0_bottomR;
        toArray[2] = z0_topL;
        toArray[3] = z0_topR;

        toArray[4] = z1_bottomL;
        toArray[5] = z1_bottomR;
        toArray[6] = z1_topL;
        toArray[7] = z1_topR;

        return upscale3D(toArray, 8, direction);
    }

    static PermValType upscale(const PermValType* z0Layer, const PermValType* z1Layer, FlowDirection direction = X)
    {
        PermValType toArray[8] = {0};
        for (int i = 0; i <4; ++i) {
            toArray[i] = z0Layer[i];
            toArray[4 + i] = z1Layer[i];
        }
        return upscale3D(toArray, 8, direction);
    }

    static PermValType upscale2D(const PermValType* permVals, int sz, FlowDirection direction = X)
    {
        assert(sz == 4);
        if (sz != 4) {
            return -1;
        }
        double renormVal = 0;
        if (direction == X) {
            renormVal = execute(permVals,sz);
        }
        if (direction == Y) {
            PermValType toRenom[4];
            for (int i = 0; i < 4; ++i) {
                if ((i + 1) & 1) {
                    toRenom[i / 2] = permVals[i];
                }
                else {
                    toRenom[(4 + i) >> 1] = permVals[i];
                }
            }
            renormVal = execute(toRenom, sz);
        }
        return renormVal;
    }

    static PermValType upscale3D(const PermValType* permVals, int sz, FlowDirection direction = X)
    {
        assert(sz == 8);
        if (sz != 8) {
            return -1;
        }
        double renormVal = 0;
        if (direction == X) {
            renormVal = execute(permVals,sz);
        }

        return renormVal;
    }

private:

    static PermValType execute (const PermValType* permVals, int sz)
    {
        bool order[3] = {true, false, false};
        PermValType toRenom[8] = {0};
        int orderSz = -1;
        int blockSz = sz;
        int logval = sz;
        while (logval != 0) {logval >>= 1; ++orderSz;}

        double permMin = 0;
        double permMax = 0;

        for (int i = 0; i < blockSz ; ++i) {
            assert(permVals[i] != 0);
            toRenom[i] = permVals[i];
        }

        if (!renormalize(toRenom, blockSz, order, orderSz)) {
            return -1;
        }
        permMin = toRenom[0];

        // Transpose
        for (int i = 0; i < blockSz; ++i) {
            if ((i + 1) & 1) {
                toRenom[i / 2] = permVals[i];
            }
            else {
                toRenom[(blockSz + i) >> 1] = permVals[i];
            }
        }
        for (int i = 0; i < orderSz; ++ i) {
            order[i] = false;
        }
        order[orderSz - 1] = true;
        if (!renormalize(toRenom, blockSz, order, orderSz)) {
            return -1;
        }

        permMax = toRenom[0];
        return geometricMean(permMin, permMax);

    }

    static bool renormalize (PermValType* permVals, int sz, const bool* order, int orderSz)
    {
        if (sz == 1 && !orderSz) {
            return true;
        }
        // Invalid states
        else if ((sz == 1 && orderSz) || (sz != 1 && !orderSz)) {
            return false;
        }
        else {
            if (*order) {
                int j = 0;
                for (int i = 0; i < sz; i += 2, ++j) {
                    permVals[j] = harmonicMean(permVals[i], permVals[i + 1]);
                }
                return renormalize(permVals, sz >> 1, ++order, --orderSz);
            }
            else if (!(*order)) {
                int j = 0;
                for (int i = 0; i < sz; i += 2, ++j) {
                    permVals[j] = arithmeticMean(permVals[i], permVals[i + 1]);
                }
                return renormalize(permVals, sz >> 1, ++order, --orderSz);
            }
        }
        return false;
    }

    static PermValType harmonicMean(const PermValType v1, const PermValType v2)
    {
        return (2 * (v1 * v2)) / (v1 + v2);
    }

    static PermValType arithmeticMean(const PermValType v1, const PermValType v2)
    {
        return (v1 + v2) / 2;
    }

    static PermValType geometricMean(const PermValType v1, const PermValType v2)
    {
        return  std::sqrt(v1 * v2);
    }
};


//Export types
typedef RVERenormalization<double> SimplifiedRenorm;

#endif //__SRUpscale__
