#include "Algorithms.cpp"
#include "./MaChaMP/MaChaMP.cpp"
#include "./CUSUM/CUSUM.cpp"
#include "./EWMA/EWMA.cpp"
#include "./QChart/QChart.cpp"


extern "C" {
    // Changepoint detection on the vector v
    void change(std::vector<double>* v,const char* TEST, double *seqarr, double *timearr, int length){
        // It is possible to change the testing function here, K-S test or Mann-Whitney, for instance
        std::vector<double> seq(seqarr, seqarr + length);
        std::vector<double> time(timearr, timearr + length);

        MCP Result;
        if (strcmp(TEST,"MaChaMP")==0 ){
            Result = MaChaMP(seq,Multiple_TTest);
        }
        else if(strcmp(TEST,"CUSUM")==0 ){
            Result = CUSUM(seq);
        }
        else if(strcmp(TEST,"EWMA")==0 ){
            Result = EWMA(seq,time);
        }
        else if(strcmp(TEST,"QChart")==0 ){
            Result = QChart(seq);
        }
        else {
            std::cout<< "Unrecognized testing framework" << std::endl;
        }

        v->push_back(Result.locations.size()); // v[v.size()] = number of changepoints
        for (auto const& loc : Result.locations){ // v[v.size()+1,v.size()+v[v.size()]] = locations of changepoints
            v->push_back(loc);
        }
        v->push_back(Result.p_value); //v[-3] = p-value
        v->push_back(Result.window); // v[-2] = window size
        v->push_back(Result.threshold); // v[-1] = best candidate
    }

    // Helper fucntions to build the vector containing the time series
    std::vector<double>* new_vector(){
        return new std::vector<double>;
    }
    void delete_vector(std::vector<double>* v){
        // std::cout << "destructor called in C++ for " << v << std::endl;
        delete v;
    }
    int vector_size(std::vector<double>* v){
        return v->size();
    }
    double vector_get(std::vector<double>* v, int i){
        return v->at(i);
    }
    void vector_push_back(std::vector<double>* v, double i){
        v->push_back(i);
    }
}

int main(){

}



/*

#include <gtest/gtest.h>

//nChoosek

TEST(TnChoosek, zerohandle){
    EXPECT_EQ(1,nChoosek(0,0));
}

#include <boost/math/special_functions/factorials.hpp>
TEST(TnChoosek, boost){
    // More than 40 candidates would mean n > 10^10
    for(int i = 1; i < 40; i = i + 3){
        for(int j = 1; j<i; j=j+3){
            double fac = boost::math::factorial<double>(i)/( boost::math::factorial<double>(j) *  boost::math::factorial<double>(i-j));
            EXPECT_NEAR(fac,nChoosek(i,j), 0.0005);
        }
    }
}

//mean and variance

TEST(Stats, zerohandle){
    std::vector<long double> none {};
    EXPECT_EQ(0, mean(none));
    EXPECT_TRUE(std::isnan(variance(none)));
    EXPECT_TRUE(std::isnan(Percentile(none.begin(), none.end(),0.5)));
}

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
TEST(Stats, MeanvarBoost){
    std::mt19937 gen(0);
    std::normal_distribution<double> norm(0,1.);

    for(int i = 10; i < 100000; i = i*2){
        boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> acc;
        std::vector<long double> seq;
        for(int k = 0; k<i;k++){
            double append = norm(gen);
            acc(append);
            seq.push_back(append);
        }
        EXPECT_NEAR(mean(seq), boost::accumulators::mean(acc),0.000005);
        EXPECT_NEAR(variance(seq), boost::accumulators::variance(acc),0.000005);
    }
}

// Test Fisher Yates

TEST(FiYa, simple){
    for(int i = 1; i < 100000; i=i*2){
        std::vector<double> test(i);
        boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> acc;
        std::vector<int> trans=FisherYates(test);
        for(auto const& k:trans){
            acc(k);
        }
        EXPECT_NEAR(i/2., boost::accumulators::mean(acc),0.5);
    }
}


// test Top-k

TEST(Topk,Examples){
    std::vector<double> vec{3.1,5.7,3.2,9.4,10.4,5.3,92.34,3245.234,4567.2345,345.23,345.2345,3456.234,5342.253,23.423};
    std::vector<size_t> topfive{0,2,5,1,3};

    EXPECT_EQ(topfive,topk_index(vec,5));
}

// test percentile: This test is difficult!!!
TEST(Stats,PercNormBoost){
    std::mt19937 gen(0);
    std::normal_distribution<double> norm(0,1.);
    // for sequences longer than 200, this test fails
    for(int i = 10; i < 100000; i = i*2){
        boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::tail_quantile<boost::accumulators::left>>> acc( boost::accumulators::left_tail_cache_size = 1000);
        std::vector<long double> seq;
        for(int k = 0; k<i;k++){
            double append = norm(gen);
            acc(append);
            seq.push_back(append);
        }
        for(double j = 0.01; j < 0.5; j+=0.01){
            if(std::isnan(boost::accumulators::quantile(acc,boost::accumulators::quantile_probability = j)) == false){
                //std::cout << j << " " << seq << std::endl;
                // Precision = 1 std deviation
                double precision =  1./j; //smaller percentiles evaluate to larger error margins
                EXPECT_NEAR(Percentile(seq.begin(),seq.end(),j), boost::accumulators::quantile(acc,boost::accumulators::quantile_probability = j),precision);
            }
        }
    }
}





//The Hausdorff file

TEST(Boxcount, ExampleSlope){
    std::vector<int> x {1,2,3,4,5};
    std::vector<int> yp2 {0,2,4,6,8};
    std::vector<int> yp1 {1,2,3,4,5};
    std::vector<int> yn3 {1,-2,-5,-8,-11};

    EXPECT_EQ(1, slope(x,yp1));
    EXPECT_EQ(2, slope(x,yp2));
    EXPECT_EQ(-3, slope(x,yn3));
}

#include "./polyfit.hpp"

TEST(Boxcount, PolyfitBoost){
    std::mt19937 gen(0);
    std::uniform_int_distribution<> integer(-100,100);

    for(int i = 10; i < 100000; i = i*2){
        std::vector<int> x;
        std::vector<double> xb;
        std::vector<int> y;
        std::vector<double> yb;
        for(int k = 0; k<i;k++){
            int append = integer(gen);
            y.push_back(append);
            yb.push_back((double) append);
            x.push_back(k);
            xb.push_back((double) k);
        }

        //Fit for polynomial of order 1
        EXPECT_NEAR(slope(x,y), polyfit(xb,yb,1)[1],0.05);
    }
}

TEST(Boxcount,ExamplesIndex){
    EXPECT_EQ(std::make_tuple(2,2), Box_index(5,99.923,2,2,95));
    EXPECT_EQ(std::make_tuple(0,0), Box_index(0,0,2,2,0));
}


TEST(Boxcount, Brownian){
    std::mt19937 gen(0);
    std::normal_distribution<long double> norm(0,1.);
    std::vector<long double> BM {0.};
    int length = 100000;
    for(int k = 0; k<length;k++){
        BM.push_back(BM.back()+norm(gen));
    }
    std::vector<int> index(length);
    std::iota(index.begin(),index.end(),0);

    // There is a bias twoards small BCD, therefore 0.3
    EXPECT_NEAR(1.5,BCD(BM,index),0.3);
}


// If variance/drift increases, the BCD odes so, too
TEST(Boxcount, Order){
    std::mt19937 gen(0);
    std::normal_distribution<long double> norm(0,1.);
    std::vector<long double> BMsmall {0.};
    std::vector<long double> BMlarge {0.};
    int length = 1000000;
    for(int k = 0; k<length;k++){
        long double append = norm(gen);
        BMsmall.push_back(BMsmall.back()+append+0.01);
        BMlarge.push_back(BMlarge.back()+5*append+0.01);
    }
    std::vector<int> index(length);
    std::iota(index.begin(),index.end(),0);

    EXPECT_TRUE(BCD(BMsmall,index)< BCD(BMlarge,index));
}



// T-Test and stuff

TEST(TTest, Table){
    for(int seed = 0; seed < 42; seed++){
        std::mt19937 gen(seed);
        std::normal_distribution<double> norm(0,1.);
        std::vector<double> seq;
        int length = 10000;
        for(int k = 0; k<length;k++){
            seq.push_back( norm(gen) );
        }
        std::vector<MeanStd> Tab = Create_Table(seq);

        for(int k = 3; k < length-3; k=k*2){
            std::vector<long double> sub(&seq[0], &seq[k]);
            EXPECT_NEAR(Tab[k].m, mean( sub ),  0.01/k);
            EXPECT_NEAR(Tab[k].s, variance(sub), 0.01/k);
        }
    }

}


TEST(TTest, mean_var){
    for(int seed = 0; seed < 42; seed++){
        std::mt19937 gen(seed);
        std::normal_distribution<double> norm(0,1.);
        std::vector<double> seq;
        int length = 10000;
        for(int k = 0; k<length;k++){
            seq.push_back( norm(gen) );
        }
        std::vector<MeanStd> Tab = Create_Table(seq);

        MeanStd Full = Tab[static_cast<int>(seq.size())];
        std::vector<long double> full(&seq[0], &seq[length]);
        EXPECT_NEAR(mean( full ), Full.m, 0.001);
        EXPECT_NEAR(variance( full ), Full.s, 0.001);

        for (int k = 3; k<length-3; k=k*2){
            std::vector<long double> sub(&seq[k], &seq[length]);
            MeanStd First = Tab[k];
            MeanStd loc = updating(Full.m, First.m, Full.s, First.s ,length, k);
            EXPECT_NEAR(mean( sub ), loc.m, 0.1);
            EXPECT_NEAR(variance( sub ), loc.s, 0.1);
        }
    }
}

*/


