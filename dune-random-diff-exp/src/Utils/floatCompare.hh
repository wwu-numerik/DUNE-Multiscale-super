#ifndef __FLOAT_COMPARE__
#define __FLOAT_COMPARE__

#include<climits>
#include<limits>
template<typename floatT>
bool floatCompare(const floatT& left, const floatT& right, const floatT& maxTolError)
{
    floatT diff = std::abs(left - right);
    if (left <= right && left >= right) {
        return true;
    }
    if (diff < maxTolError) {
        return true;
    }
    return false;

}

#endif // __FLOAT_COMPARE
