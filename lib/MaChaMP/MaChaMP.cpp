#include "./Dim/Dimension.cpp"
#include "DM.cpp"
#include <math.h>


//Retrun type
struct MCP{
	std::vector<int>  locations = {}; //indices in sequence
	long double p_value =1.; //p-vlaue
            int window = 2;  // Window-size
            long double threshold; // Threshold of candidates
};

// Tree-like data structure for recursive computations of changepoint configuration q-values hallo ich bin leonie
class TopDownCombinations {
    std::vector<long double> Scores; // List of q-values corrseponding to Configurations
    std::vector<std::vector<int>> Configurations; // List of checked configuratioins of candidate changepoints
    int max; // maximal tree depth = number of candidates
    int n; // sequence length
    int runs; // number of times the recursive function call is executed
    long double threshold = 1.; //for returning


public:
    void Init(std::vector<int> Initial,std::vector<double> sequence,long double Score, long double thres);  // Initialization of the tree
    void Fit(std::function<MC(std::vector<int> &,std::vector<MeanStd> &, bool, bool)> Test,std::vector<MeanStd> &Table,std::vector<int> &Configurations_values,long double Master_score, int before); // Recursive computation of all branches
    MCP Transform(); // Search among Configurations for best Score
};

// initial: List of candidates
// sequence: Time series
// Score: Test(sequence,Initial,Table,false).p
void TopDownCombinations::Init(std::vector<int> Initial, std::vector<double> sequence,long double Score, long double thres){


    max = Initial.size();
    n = sequence.size();
    runs = 0;
    threshold = thres;

    Configurations.push_back(Initial);
    Scores.push_back( Score);
}

// This computes the branch
// sequence: Time series
// Test: t-Test
// Table: n precomputed [0,n] means and variances
// Master_score: q-value of mother branch
// before: How many candidates are before the current left-out candidate?
void TopDownCombinations::Fit(std::function<MC(std::vector<int> &,std::vector<MeanStd> &, bool, bool)> Test,std::vector<MeanStd> & Table,std::vector<int> &Configurations_values, long double Master_score, int before){
    for(int left_out = 0; left_out < Configurations_values.size()-before; left_out++){
        runs++;

        // Everything here is O(1)
        // Only compute branch, if [1,0,1,1,.... i ... 1,0,1,1,0] i at location indexing can be changed from 1->0
        //if(Configurations_master.back()[indexing] == 1 && Configurations_values.back().size() > left_out+before ){
            std::vector<int> Current_values = Configurations_values;
            int new_before = before + left_out; // Update how many candidates there have been left untouched
            int c = Current_values[new_before]; //The candidate which will be removed
            Current_values.erase(Current_values.begin()+new_before); //remove candidate value of confirguations
            // Update everything but the q-value
            Configurations.push_back(Current_values);
            double new_Master_score = 0;

            // If there are more than 2 candidates before and after the to be removed candaite, use Dean's update rule
            if ((new_before > 1) && (((static_cast<int>(Current_values.size())-new_before )>1))) {
                // Extract indices of the current candiated and the two changes before and after it
                int _c = Current_values[new_before-1];
                int __c = Current_values[new_before-2];
                int c_ = Current_values[new_before];
                int c__ = Current_values[new_before+1];
                // Construct the configurations needed for recursive computation
                std::vector<int> tp1 {__c,_c,c};
                std::vector<int> tp2 {_c,c,c_};
                std::vector<int> tp3 {c,c_,c__};
                std::vector<int> tn1 {__c,_c,c_};
                std::vector<int> tn2 {_c,c_,c__};
                // remove all tpx and then add tnx to get current configruation q-s ore
                new_Master_score = Master_score + 2*( std::log(Test(tp1,Table,false,true).p) +  std::log(Test(tp2,Table,false,true).p) + std::log(Test(tp3,Table,false,true).p) -  std::log(Test(tn1,Table,false,true).p) - std::log(Test(tn2,Table,false,true).p)); //Deans update rule
                // Update the list of q-scores
                Scores.push_back( new_Master_score );
            } // Otherwise, compute the q-value directly using Test function
            else if (Current_values.size()>0){
                new_Master_score = Test(Current_values,Table,false,false).p;
                Scores.push_back( new_Master_score );
            }

            // Call to remove another changepoint form the current configruation
            if (Current_values.size()>1){
                    Fit(Test,Table,Current_values,new_Master_score,new_before);
            }


    }

}

// Search for best configurtaion after fit.
MCP TopDownCombinations::Transform(){
    // Convert q-scores into p-values
    std::vector<long double> p_values;
    for(int i =0; i < Configurations.size(); i++){
        int m = static_cast<int>(Configurations[i].size()); // degree of freedom = 2*number of candidates
        p_values.push_back( nChoosek(n- (m + 1) * 2 + m, m) * Chisq(Scores[i],  2*m ) ); // Bonferroni correction for each m
    }
    // Index is location in all vectors of configuration with lowest p-value
    int index = std::distance(p_values.begin(), std::min_element( p_values.begin(), p_values.end() ));
    MCP Best;
    Best.locations = Configurations[index];
    Best.p_value = p_values[index];
    Best.threshold = threshold;

    if(Best.p_value < threshold){
        return Best;
    }
    else{
        std::vector<int> none {};
        MCP Null;
        Null.locations = {};
        Null.p_value = 1.;
        Null.threshold = threshold;
        return Null;
    }
}

// Main routine
MCP MaChaMP(std::vector<double> sequence,  std::function<MC(std::vector<int> &,std::vector<MeanStd> &, bool, bool)> Test);

// Subroutine of former Welch-Fisher computing every configuration with Test function.
MCP Combinatorical(std::vector<double> sequence, int m, std::function<MC(std::vector<int> &,std::vector<MeanStd> &, bool, bool)> Test, std::vector<MeanStd> Table, std::vector<Changepoint> Candidates);

// Solves the Welch-Fisher optimization objective, given candidates
MCP Welch_Fisher(std::vector<double> sequence,std::function<MC(std::vector<int> &,std::vector<MeanStd> &, bool, bool)> Test, std::vector<MeanStd> Table, std::vector<Changepoint> Candidates);

// Estimate the window size for extraction of candidates using the Box-counting dimension
int Estimate_Scale(std::vector<double>sequence,std::function<MC(std::vector<int> &,std::vector<MeanStd> &, bool, bool)> Test,std::vector<MeanStd> Table);




MCP MaChaMP(std::vector<double> sequence, std::function<MC(std::vector<int> &,std::vector<MeanStd> &, bool, bool)> Test){
    // Create lookup table
    std::vector<MeanStd> Table = Create_Table(sequence);
    // Initialize return type
    MCP Return_MCP;

    // Estimate window size and p-value threshold using Hausdorff fractal dimension
    int window_size = Estimate_Scale(sequence,Test,Table); //70% of the runtime is spent here
    Return_MCP.window = window_size;
    bool Hack = true; // If overflow occurs at DM?

    // Extract candidates (and threshold)
    std::vector<Changepoint> Candidates = Candidate_extraction( DM(window_size,Test,Table), window_size, Hack);

    // Recombine candidates and use the configuration with minimal p-value
    Return_MCP = Welch_Fisher(sequence, Test,Table,Candidates);

    return Return_MCP;
}



// /*

MCP Welch_Fisher(std::vector<double> sequence,std::function<MC(std::vector<int> &,std::vector<MeanStd> &, bool, bool)> Test, std::vector<MeanStd> Table, std::vector<Changepoint> Candidates){

    // Create tree class
    TopDownCombinations WelchFisher;

    // Compute q-score of configuration taking all candidates
    std::vector<int> Initial;
    long double threshold = 1.;
    for(int i=0;i<Candidates.size();i++){
        Initial.push_back(Candidates[i].i);
        if( Candidates[i].p < threshold){
            threshold = Candidates[i].p;
        }
    }
    long double Score = Test(Initial,Table,false,false).p;

    // Initialize search: Tree class and master node
    WelchFisher.Init(Initial,sequence,Score,threshold);
    std::vector<int> Master(Candidates.size());
    std::fill(Master.begin(),Master.end(),1);

    // Search all configurations
    WelchFisher.Fit(Test,Table,Initial,Score,0);

    // compute multiple hypothesis corrected p-values and take minimum
    return WelchFisher.Transform();
}
// */


// The minimum Box-counting dimension among all dissimilarity manifolds yields the best scale
int Estimate_Scale(std::vector<double> sequence,std::function<MC(std::vector<int>& ,std::vector<MeanStd> &, bool, bool)> Test,std::vector<MeanStd> Table){

    // Initialization of search
    std::vector<int> rand_index = FisherYates(sequence); //The random index permutation, which is accessed be BCD for
    int n = sequence.size();
    // subsampling logn candidates instead of taking all data points
    int window_exp = (int) std::log2(n) -0.5; // the largest possible window size exponent
    std::vector<double> dimension; // Vector containing the result of the BCD
    int delta = 1; //Stepsize of reducing the exponent of the window
    std::vector<int> window; // Vector containing the tested widnow sizes

    // Iterate over all window sizes until size 4
    while(window_exp > 2){
        window.push_back((int) std::pow(2,window_exp)); // Window size=2**window_exp
        //
        // Improve: Only compute DM for log-subsample!!!!
        //
        auto DM_log = log_sample(rand_index,Test,Table,std::pow(2,window_exp));
        dimension.push_back( BCD(  DM_log,n )); // Compute BCD
        window_exp -= delta; // reduce widnow size by 1
    }

    // The optimal window size is the one with lowest Box-counting dimension
    int argmin = std::distance(dimension.begin(),std::min_element(dimension.begin(), dimension.end()));
    int op_window = window[argmin];


    return op_window;
}

