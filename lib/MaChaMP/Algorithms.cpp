#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>

#include <tuple>
#include  <random>
#include <numeric>
#include <unordered_map>


// This function returns the Perc percentile of an iterator (begin,end)
template <typename It>
auto Percentile(It begin, It end,double Perc)
{
    using T = typename std::iterator_traits<It>::value_type;
    std::vector<T> data(begin, end);

    if (data.size() > 1){
        //This is in accord with the bosst version: Always round down
        std::nth_element(data.begin(),data.begin()+ (int) (std::distance(data.begin(),data.end())*Perc-1.), data.end());
        return data[(int) (std::distance(data.begin(),data.end())*Perc-1.)];
    }
    else{
        return std::numeric_limits<long double>::quiet_NaN();
    }
}

// Returns list of indices of the k smallest elements in the vector values
template <typename T>
std::vector<size_t> topk_index(std::vector<T> const& values, int k) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::nth_element (indices.begin(), indices.begin()+k, indices.end(),[&](size_t a, size_t b) { return values[a] < values[b]; });
    // std::sort( begin(indices), end(indices), [&](size_t a, size_t b) { return values[a] < values[b]; });

    if (indices.size() > k){
        indices.resize(k);
    }
    return indices;
}

// N choose k
unsigned long nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    unsigned long result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

// returns mean of vector values
auto mean(std::vector<long double> values){
    long double sum = 0;
    for(auto const& el: values){
        sum += el;
    }
    if(values.size() > 0){
        return sum/values.size();
    }
    else{
        return sum;
    }
}

// returns variance of vector values
auto variance(std::vector<long double> values){
    long double val_mean = mean(values);
    long double val_var = 0;
    for(auto const& el: values){
        val_var += std::pow(el-val_mean,2);
    }
    if(values.size()>1){
        return val_var/(values.size()); //biased
    }
    else{
        return std::numeric_limits<long double>::quiet_NaN();
    }
}

// Fisher-Yates shuffle on inices of length sequence. This is a randomization, O(n)
std::vector<int> FisherYates(std::vector<double> sequence){
    std::mt19937 gen(sequence.size());
    std::uniform_int_distribution<int> uni(0,sequence.size());
    std::vector<int> index(sequence.size());
    std::iota(index.begin(),index.end(),0);
    for (int k = 0; k <  sequence.size(); k++) {
        int r = k + uni(gen) % (sequence.size() - k); // careful here!
        std::swap(index[k], index[r]);
    }
    return index;
}

