#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>

//This header overloads add, subtract and cout for std::vectors
std::vector<int> operator+(const std::vector<int>& a, const std::vector<int>& b)
{
    assert(a.size() == b.size());

    std::vector<int> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::plus<int>());
    return result;
}
std::vector<long double> operator+(const std::vector<long double>& a, const std::vector<long double>& b)
{
    assert(a.size() == b.size());

    std::vector<long double> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::plus<long double>());
    return result;
}
std::vector<long double> operator-(const std::vector<long double>& a, const std::vector<long double>& b)
{
    assert(a.size() == b.size());

    std::vector<long double> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::minus<long double>());
    return result;
}
std::vector<int> operator-(const std::vector<int>& a, const std::vector<int>& b)
{
    assert(a.size() == b.size());

    std::vector<int> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::minus<int>());
    return result;
}
std::vector<int> operator*(const std::vector<int>& a, const std::vector<int>& b)
{
    assert(a.size() == b.size());

    std::vector<int> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::multiplies<int>());
    return result;
}
inline std::ostream& operator << (std::ostream& os, const std::vector<int>& v)
{
    os << "[";
    for (std::vector<int>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii;
    }
    os << " ]";
    return os;
}
inline std::ostream& operator << (std::ostream& os, const std::vector<long double>& v)
{
    os << "[";
    for (std::vector<long double>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii;
    }
    os << " ]";
    return os;
}
inline std::ostream& operator << (std::ostream& os, const std::vector<double>& v)
{
    os << "[";
    for (std::vector<double>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii;
    }
    os << " ]";
    return os;
}

