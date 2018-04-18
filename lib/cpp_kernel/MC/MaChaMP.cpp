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
    std::vector<long double> TTest_vals;
    std::vector<std::vector<int>> Configurations; // List of checked configuratioins of candidate changepoints
    std::vector<int> Candidates;
    int max; // maximal tree depth = number of candidates
    int n; // sequence length
    int runs; // number of times the recursive function call is executed
    long double threshold = 1.; //for returning
    MC Test_Return;
    int _c;
    int __c;
    int c_;
    int c__;
    int c;
    int m;
    int msq;
    std::vector<int> loc;
    std::vector<MeanStd> List_of_Stats;


public:
    void Init(std::vector<int> Initial,std::vector<double> sequence,long double Score, long double thres,std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test,std::vector<MeanStd> &Table);  // Initialization of the tree
    void Fit(std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test,std::vector<MeanStd> &Table,std::vector<int> &Configurations_values,long double Master_score, int before); // Recursive computation of all branches
    MCP Transform(); // Search among Configurations for best Score
};

// initial: List of candidates
// sequence: Time series
// Score: Test(sequence,Initial,Table,false).p
void TopDownCombinations::Init(std::vector<int> Initial, std::vector<double> sequence,long double Score, long double thres,std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test,std::vector<MeanStd> &Table){


    max = Initial.size();
    n = sequence.size();
    runs = 0;
    threshold = thres;
    TTest_vals = Test_lookup(Initial, Table,Test); //logn^3
    Candidates = Initial;
    m = Initial.size()+2;
    msq=m*m;

    Configurations.push_back(Initial);
    Scores.push_back( Score);
}

// This computes the branch
// sequence: Time series
// Test: t-Test
// Table: n precomputed [0,n] means and variances
// Master_score: q-value of mother branch
// before: How many candidates are before the current left-out candidate?
void TopDownCombinations::Fit(std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test,std::vector<MeanStd> & Table,std::vector<int> &Configurations_values, long double Master_score, int before){
    for(int left_out = 0; left_out < Configurations_values.size()-before; left_out++){
        runs++;

        // Everything here is O(1)
        // Only compute branch, if [1,0,1,1,.... i ... 1,0,1,1,0] i at location indexing can be changed from 1->0
        //if(Configurations_master.back()[indexing] == 1 && Configurations_values.back().size() > left_out+before ){
            std::vector<int> Current_values = Configurations_values; //Instead have single int with different indices at the correct digit?
            int new_before = before + left_out; // Update how many candidates there have been left untouched
            c = Current_values[new_before]; //The candidate which will be removed
            Current_values.erase(Current_values.begin()+new_before); //remove candidate value of confirguations
            // Update everything but the q-value
            Configurations.push_back(Current_values);
            long double new_Master_score = 0;


            // If there are more than 2 candidates before and after the to be removed candaite, use Dean's update rule
            if ((new_before > 1) && (((static_cast<int>(Current_values.size())-new_before )>1))) {
                //std::cout << "There";
                // Extract indices of the current candiated and the two changes before and after it
                _c = Current_values[new_before-1];
                __c = Current_values[new_before-2];
                c_ = Current_values[new_before];
                c__ = Current_values[new_before+1];
                new_Master_score += Master_score - 2*((TTest_vals[__c + _c*m + c*msq ]) -(TTest_vals[_c + c*m + c_*msq ]) -(TTest_vals[c + c_*m + c__*msq ]) +(TTest_vals[__c + _c*m + c_*msq ])  +(TTest_vals[_c + c_*m + c__*msq ]) );
                Scores.push_back( new_Master_score );
            } // Otherwise, compute the q-value directly using
            else if (Current_values.size()>0){
                if (Current_values.size() > 2){
                    new_Master_score += 2*(TTest_vals[0+Current_values[0]*m+Current_values[1]*msq]) ;//0-1-2
                    for(int i = 1; i < Current_values.size()-1; i++){
                        new_Master_score += 2*(TTest_vals[Current_values[i-1]+Current_values[i]*m+Current_values[i+1]*msq]); //-2-1-n
                    }
                    new_Master_score += 2*(TTest_vals[Current_values[Current_values.size()-2]+Current_values[Current_values.size()-1]*m+(m-1)*msq]); //-2-1-n
                }
                else{
                    loc.resize(Current_values.size());
                    for(int i = 0; i < Current_values.size(); i++){
                        loc[i] = Candidates[Current_values[i]-1];
                    }
                    Test(loc,Table,false,false,Test_Return,List_of_Stats);
                    new_Master_score = Test_Return.p;
                }
                // std::cout << new_Master_score << std::endl;
                Scores.push_back( new_Master_score );
            }

        std::cout << Current_values << " " << Scores[-1] << std::endl;
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
    double best_9 = 1;
    double best_3 = 1.;
    for(int i =0; i < Configurations.size(); i++){
        int number = static_cast<int>(Configurations[i].size()); // degree of freedom = 2*number of candidates
        if(number == 2 && best_3 > Chisq(Scores[i],  2*number ) ){
            best_3 = Chisq(Scores[i],  2*number );}
        if(number == 3 &&best_9 > Chisq(Scores[i],  2*number ) ){
            best_9 = Chisq(Scores[i],  2*number );}
        p_values.push_back( Chisq(Scores[i],  2*number ) ); // Bonferroni correction for each m
    }
    std::cout << best_3 << " " << best_9 << std::endl;
    for(int i =0; i < Configurations.size(); i++){
        int number = static_cast<int>(Configurations[i].size()); // degree of freedom = 2*number of candidates
        p_values[i] =  nChoosek(n- (number + 1) * 2 + number, number) * p_values[i]; // Bonferroni correction for each m
    }
    // Index is location in all vectors of configuration with lowest p-value
    int index = std::distance(p_values.begin(), std::min_element( p_values.begin(), p_values.end() ));
    MCP Best;
    Best.locations.resize(Configurations[index].size());
    for(int i = 0; i < Configurations[index].size(); i++){
        Best.locations[i] = Candidates[Configurations[index][i]-1];
    }
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
MCP MaChaMP(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);

// Subroutine of former Welch-Fisher computing every configuration with Test function.
MCP Combinatorical(std::vector<double> sequence, int m, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test, std::vector<MeanStd> Table, std::vector<Changepoint> Candidates);

// Solves the Welch-Fisher optimization objective, given candidates
MCP Welch_Fisher(std::vector<double> sequence,std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test, std::vector<MeanStd> &Table, std::vector<Changepoint> Candidates);

// Estimate the window size for extraction of candidates using the Box-counting dimension
int Estimate_Scale(std::vector<double>sequence,std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test,std::vector<MeanStd> &Table);




MCP MaChaMP(std::vector<double> sequence, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
    // Create lookup table
    std::vector<MeanStd> Table = Create_Table(sequence);
    // Initialize return type
    MCP Return_MCP;

    // Estimate window size and p-value threshold using Hausdorff fractal dimension
    int window_size = Estimate_Scale(sequence,Test,Table); //70% of the runtime is spent here
    //std::cout << "Scale : " << window_size << std::endl;
    bool Hack = true; // If overflow occurs at DM?

    // Extract candidates (and threshold)
    std::vector<Changepoint> Candidates = Candidate_extraction( DM(window_size,Test,Table), window_size, Hack);
    //std::cout<< "Candidates : " << counter << std::endl;

    // Recombine candidates and use the configuration with minimal p-value
    Return_MCP = Welch_Fisher(sequence, Test,Table,Candidates);
    //std::cout << "Welch-Fisher : " << counter << std::endl;

    return Return_MCP;
}



// /*

MCP Welch_Fisher(std::vector<double> sequence,std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test, std::vector<MeanStd> &Table, std::vector<Changepoint> Candidates){

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
    MC Test_Return;
    std::vector<MeanStd> List_of_Stats(Initial.size());
    Test(Initial,Table,false,false,Test_Return,List_of_Stats);
    long double Score = Test_Return.p;
    std::cout << Initial << std::endl;
    //long double Score = Test(Initial,Table,false,false).p;

    // Initialize search: Tree class and master node
    WelchFisher.Init(Initial,sequence,Score,threshold,Test,Table);
    std::vector<int> Master(Candidates.size());
    std::fill(Master.begin(),Master.end(),1);

    // Search all configurations
    std::vector<int> indexes(Initial.size());
    std::iota(indexes.begin(), indexes.end(),1);
    WelchFisher.Fit(Test,Table,indexes,Score,0);

    // compute multiple hypothesis corrected p-values and take minimum
    return WelchFisher.Transform();
}
// */


// The minimum Box-counting dimension among all dissimilarity manifolds yields the best scale
int Estimate_Scale(std::vector<double> sequence,std::function<void(std::vector<int>& ,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test,std::vector<MeanStd> &Table){

    // Initialization of search
    std::vector<int> rand_index = FisherYates(sequence); //The random index permutation, which is accessed be BCD for
    int n = sequence.size();
    // subsampling logn candidates instead of taking all data points
    int window_exp = (int) (std::log2(n) -0.5) - 1.; // the largest possible window size exponent
    std::vector<double> dimension; // Vector containing the result of the BCD
    int delta = 1; //Stepsize of reducing the exponent of the window
    std::vector<int> window; // Vector containing the tested widnow sizes
    std::tuple<std::vector<long double>,std::vector<int>> DM_log;

    // Iterate over all window sizes until size 4
    while(window_exp > 2){
        window.push_back((int) std::pow(2,window_exp)); // Window size=2**window_exp
        DM_log = log_sample(rand_index,Test,Table,std::pow(2,window_exp));
        //std::cout << window_exp << " " << (int) std::pow(2,window_exp) << " " << BCD(  DM_log,n ) << std::endl;
        dimension.push_back( BCD(  DM_log,n )); // Compute BCD
        window_exp -= delta; // reduce widnow size by 1
    }

    // The optimal window size is the one with lowest Box-counting dimension
    int argmin = std::distance(dimension.begin(),std::min_element(dimension.begin(), dimension.end()));
    int op_window = window[argmin];
    //std::cout << op_window << " ";


    return op_window;
}

