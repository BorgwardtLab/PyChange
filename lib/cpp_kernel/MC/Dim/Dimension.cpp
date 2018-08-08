#include "../Algorithms.cpp"



std::tuple<double, double> boxlen(int k, int n, double u, int big_k) {
    double level = std::pow(2,k-big_k);
    return std::make_tuple(level*n,level*u);
}

double k_to_logepsilon(int k, int big_k, int n){
    return std::log2(n * std::pow(2,k - big_k));
}


std::tuple<int,int> Box_index(int i, double j, int x, double y, double m){
    /*
    INPUT:
    i: timepoint currently under observation
    j: seq[i] value at that point
    x: x-size of box
    y: y-size of box
    m: minimum of time series (lower bound on j)

    OUTPUT:
    x,y-index of Box corresponding to i,j
    */
    double ep = 1e-10;
    return std::make_tuple((int) ((1.*i) / x) , (int) (((j-m)/y)-ep )  );
}

std::tuple<long double,long double> get_stats(std::vector<long double> logseq){
    auto logseq_min = *std::min_element(logseq.begin(),logseq.end());
    auto logseq_max = *std::max_element(logseq.begin(),logseq.end());

    auto m = logseq_min;
    auto u = logseq_max-logseq_min;

    /* Alternative: Use expected mean and variances to calculate spread. In general, this reduces the BCD for long sequences
    auto logseq_mean = mean(logseq);
    auto logseq_std = std::sqrt(variance(logseq));
    m = logseq_mean-1.5*logseq_std;
    if(m > logseq_min){m = logseq_min;}
    u = logseq_max-logseq_min;
    if(u < 3*logseq_std){u=3*logseq_std;}
    */


    return std::make_tuple(u,m);
}



double slope(const std::vector<int>& x, const std::vector<int>& y) {
    const auto n    = x.size();
    const auto s_x  = std::accumulate(x.begin(), x.end(), 0.0);
    const auto s_y  = std::accumulate(y.begin(), y.end(), 0.0);
    const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    return a;
}

int Ne(std::vector<long double> logseq, std::vector<int> index,int k, int n, int big_k, long double u, long double m){
    // double ep = 1e-10;  // Take care of rounding errors in y-dimension
    auto xy = boxlen(k,n,u,big_k);
    double x = std::get<0>(xy);
    double y = std::get<1>(xy);
    int sum_of_boxes = (int) (n/x+1.5);
    if (u != 0.){
        std::unordered_map<int,int> Boxes; //Initialze map of indiactor map of occupied boxes
        for(int i = 0; i < index.size();i++){
            auto ab = Box_index(index[i], logseq[i], (int) x+0.5, y, m); //Extract Box index of element i in logseq with box size x,y
            int a = std::get<0>(ab); //Box index return a tuple, extract first element
            int b = std::get<1>(ab); // a,b are now the indices
            Boxes[n*a+b] = 1; // Label Box at n*a+b as occupied, this is on average O(1) because http://www.cplusplus.com/reference/unordered_map/unordered_map/find/
        }
        sum_of_boxes = Boxes.size(); //count how many boxes are occupied this is O(lgn)
    }
    return sum_of_boxes; //The 2**sum of boxes is missing because it cancels out with log2 at a later point
}

// Pass rand index, sequence and Test, then compute DM values when ist needed
double BCD(std::tuple<std::vector<long double>,std::vector<int>> log_sample, int n){
    std::vector<long double> logseq = std::get<0>(log_sample);
    std::vector<int> index = std::get<1>(log_sample);

    int exclude_exponents = 1;  // Recommended by Liebovic and Toth, exclude two largest and smallest Box sizes
    auto um = get_stats(logseq); //Get spread and minimum SLOW!!!!!
    long double u = std::get<0>(um);
    long double m = std::get<1>(um);
    int big_k = (int) std::log2(n); //Largest box size exponent
    std::vector<int> listk; //List of evaluated box-size exponents
    std::vector<int> logepsilon; //In epsilon notation
    std::vector<int> Nepsilon;  //number of occupied Boxes

    // First iteration
    listk.push_back(big_k-exclude_exponents);
    while ( Nepsilon.size() < big_k - std::pow(2,exclude_exponents) - 1){
        Nepsilon.push_back( Ne(logseq,index,listk.back(),n,big_k,u,m)) ; //We do not have log2 here because it cancels out with 2**x earlier
        logepsilon.push_back(k_to_logepsilon(listk.back(),big_k,n));
        listk.push_back(listk.back()-1);
    }

    return -1.* slope(logepsilon, Nepsilon);
}



