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
                new_Master_score += Master_score + 2*(std::log(TTest_vals[__c + _c*m + c*msq ]) +std::log(TTest_vals[_c + c*m + c_*msq ]) +std::log(TTest_vals[c + c_*m + c__*msq ]) - std::log(TTest_vals[__c + _c*m + c_*msq ])  - std::log(TTest_vals[_c + c_*m + c__*msq ]) );
                Scores.push_back( new_Master_score );
            } // Otherwise, compute the q-value directly using
            else if (Current_values.size()>0){
                if (Current_values.size() > 2){
                    new_Master_score -= 2*std::log(TTest_vals[0+Current_values[0]*m+Current_values[1]*msq]) ;//0-1-2
                    for(int i = 1; i < Current_values.size()-1; i++){
                        new_Master_score -= 2*std::log(TTest_vals[Current_values[i-1]+Current_values[i]*m+Current_values[i+1]*msq]); //-2-1-n
                    }
                    new_Master_score -= 2*std::log(TTest_vals[Current_values[Current_values.size()-2]+Current_values[Current_values.size()-1]*m+(m-1)*msq]); //-2-1-n
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
    int count_tests = 0;
    for(int i =0; i < Configurations.size(); i++){
        int number = static_cast<int>(Configurations[i].size()); // degree of freedom = 2*number of candidates
        p_values.push_back( nChoosek(n- (number + 1) * 2 + number, number) * Chisq(Scores[i],  2*number ) ); // Bonferroni correction for each m
        //p_values.push_back(Chisq(Scores[i],  2*number ) ); // Bonferroni correction for each m
        count_tests += number; //Bonferroni for conducted tests
    }
    // Index is location in all vectors of configuration with lowest p-value
    int index = std::distance(p_values.begin(), std::min_element( p_values.begin(), p_values.end() ));
    MCP Best;
    Best.locations.resize(Configurations[index].size());
    if (index == 0){
        Best.locations = Candidates;
    }
    else{
        for(int i = 0; i < Configurations[index].size(); i++){
            if (Candidates[Configurations[index][i]-1] < 3 || Candidates[Configurations[index][i]-1] > n){
                //std::cout << Configurations[index][i] << " ";
                Best.locations[i] = Configurations[index][i];
            }
            else{
                Best.locations[i] = Candidates[Configurations[index][i]-1];
            }
        }
    }

    Best.p_value = p_values[index];//*count_tests;
    Best.threshold = threshold;

    return Best;
}

// Main routine
MCP MaChaMP(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);

// Subroutine of former Welch-Fisher computing every configuration with Test function.
MCP Combinatorical(std::vector<double> sequence, int m, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test, std::vector<MeanStd> Table, std::vector<Changepoint> Candidates);

// Solves the Welch-Fisher optimization objective, given candidates
MCP Welch_Fisher(std::vector<double> sequence,std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test, std::vector<MeanStd> &Table, std::vector<Changepoint> Candidates);

// Estimate the window size for extraction of candidates using the Box-counting dimension
int Estimate_Scale(std::vector<double>sequence,std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test,std::vector<MeanStd> &Table);

// Stage 0: Random lg n point from sequence as changes
MCP Stage_random(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);
// Stage 1: Random lg n points recombination
MCP Stage_random_recomb(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);
// Stage 2: Fixed window size candidates
MCP Stage_window(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);
// Stage 3: Fixed window size recombination good window
MCP Stage_window_recomb(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);
// Stage 3: Fixed window size recombination bad window
MCP Stage_window_recomb_bad(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);
//Stage 4: Dimensionality analysis window plain candidates
MCP Stage_dimension(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);
//Stage 5: MaChaMP - DImensionality analysis recombination with threshold 0.01
MCP Stage_fixed_alpha(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);

MCP Stage_window_recomb_random(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);




MCP Single_Test(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test);


MCP MaChaMP(std::vector<double> sequence, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
    // Create lookup table
    std::vector<MeanStd> Table = Create_Table(sequence);
    // Initialize return type
    MCP Return_MCP;

    // Estimate window size and p-value threshold using Hausdorff fractal dimension
    int window_size = Estimate_Scale(sequence,Test,Table); //70% of the runtime is spent here
    // std::cout << "Scale : " << window_size << std::endl;
    bool Hack = true; // If overflow occurs at DM?

    // Extract candidates (and threshold)
    std::vector<Changepoint> Candidates = Candidate_extraction( DM(window_size,Test,Table), window_size, Hack);

    // Recombine candidates and use the configuration with minimal p-value
    Return_MCP = Welch_Fisher(sequence, Test,Table,Candidates);
    Return_MCP.window = window_size;

    //Additional Bonferroni on threshold according to window size
    //long double Bonferroni_p = 1.; // * ((std::log2(sequence.size()-4)*std::log2(sequence.size()))*((sequence.size()-2*window_size-6)));
    long double Bonferroni_p = 1. * (((sequence.size()-2*window_size-6)));
    long double Bonferroni_a = 1. * (((sequence.size()-2*window_size-6))); //Number test for window size, number tests for candidates

    Return_MCP.p_value *= Bonferroni_p;
    Return_MCP.threshold *= Bonferroni_a;
    if (Return_MCP.p_value > Return_MCP.threshold){
        Return_MCP.locations = {};
    }
    //if( Return_MCP.p_value > 1.){
    //    Return_MCP.p_value = 1.;
    //}

    return Return_MCP;
}



// /*

MCP Welch_Fisher(std::vector<double> sequence,std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test, std::vector<MeanStd> &Table, std::vector<Changepoint> Candidates){

    // Create tree class
    TopDownCombinations WelchFisher;

    // Compute q-score of configuration taking all candidates
    std::vector<int> Initial;
    long double threshold = 1.;
    //std::cout << "Candidates:" << std::endl;
    for(int i=0;i<Candidates.size();i++){
        Initial.push_back(Candidates[i].i);
        //std::cout << "[" << Candidates[i].i << "," << Candidates[i].p << "]" << std::endl;
        if( Candidates[i].p < threshold){
            threshold = Candidates[i].p;
        }
    }
    MC Test_Return;
    std::vector<MeanStd> List_of_Stats(Initial.size());
    Test(Initial,Table,false,false,Test_Return,List_of_Stats);
    long double Score = Test_Return.p;
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
    int window_exp = (int) std::log2(n) -1.5; // the largest possible window size exponent
    std::vector<double> dimension; // Vector containing the result of the BCD
    int delta = 1; //Stepsize of reducing the exponent of the window
    std::vector<int> window; // Vector containing the tested widnow sizes
    std::tuple<std::vector<long double>,std::vector<int>> DM_log;

    // Iterate over all window sizes until size 4
    while(window_exp > 4){
        window.push_back((int) std::pow(2,window_exp)); // Window size=2**window_exp
        DM_log = log_sample(rand_index,Test,Table,std::pow(2,window_exp));
        dimension.push_back( BCD(  DM_log,n )); // Compute BCD
        window_exp -= delta; // reduce widnow size by 1
    }
    // std::cout << dimension << std::endl;

    // The optimal window size is the one with lowest Box-counting dimension
    int argmin = std::distance(dimension.begin(),std::min_element(dimension.begin(), dimension.end()));
    int op_window = window[argmin];


    return op_window;
}

MCP Stage_fixed_alpha(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
    // Create lookup table
    std::vector<MeanStd> Table = Create_Table(sequence);
    // Initialize return type
    MCP Return_MCP;

    // Estimate window size and p-value threshold using Hausdorff fractal dimension
    int window_size = Estimate_Scale(sequence,Test,Table); //70% of the runtime is spent here
    // std::cout << "Scale : " << window_size << std::endl;
    bool Hack = true; // If overflow occurs at DM?

    // Extract candidates (and threshold)
    std::vector<Changepoint> Candidates = Candidate_extraction( DM(window_size,Test,Table), window_size, Hack);

    // Recombine candidates and use the configuration with minimal p-value
    Return_MCP = Welch_Fisher(sequence, Test,Table,Candidates);
    Return_MCP.window = window_size;

    //Additional Bonferroni on threshold according to window size
    //long double Bonferroni_p = 1.; // * ((std::log2(sequence.size()-4)*std::log2(sequence.size()))*((sequence.size()-2*window_size-6)));
    long double Bonferroni_p = 1. * (((sequence.size()-2*window_size-6)));
    long double Bonferroni_a = 1. * (((sequence.size()-2*window_size-6))); //Number test for window size, number tests for candidates

    Return_MCP.p_value *= Bonferroni_p;
    Return_MCP.threshold *= Bonferroni_a;
    if (Return_MCP.p_value > 0.01){
        Return_MCP.locations = {};
    }
    //if( Return_MCP.p_value > 1.){
    //    Return_MCP.p_value = 1.;
    //}

    return Return_MCP;
}


MCP Stage_random(std::vector<double> sequence, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
    int number_candidates = std::log2(sequence.size());

    std::mt19937 rng(std::random_device{}());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(3,sequence.size()-3);; // guaranteed unbiased

    std::vector<int> Candidates;

    for(int i = 0; i < number_candidates; i++){
        Candidates.push_back(uni(rng));
    }

    std::sort(Candidates.begin(),Candidates.end());
    MCP Return_MCP;
    Return_MCP.window = 0;
    Return_MCP.p_value = 1.; // Compute real p_value
    Return_MCP.threshold = 0.;
    Return_MCP.locations = Candidates;

    return Return_MCP;
}

MCP Stage_window(std::vector<double> sequence, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
    int window_size = (int) sequence.size()/10.; //std::log2(sequence.size());
    // Create lookup table
    std::vector<MeanStd> Table = Create_Table(sequence);
    bool Hack = true; // Log candidates?

    // Extract candidates (and threshold)
    std::vector<Changepoint> Candidates = Candidate_extraction( DM(window_size,Test,Table), window_size, Hack);

    std::vector<int> locations = {};
    double threshold =1.;
    for(auto cp: Candidates){
        locations.push_back(cp.i);
        if(cp.p < threshold){
            threshold = cp.p;
        }
    }

    MCP Return_MCP;
    Return_MCP.window = 0;
    Return_MCP.p_value = 1.; // Compute real p_value
    Return_MCP.threshold = threshold;
    Return_MCP.locations = locations;

    return Return_MCP;
}

MCP Stage_dimension(std::vector<double> sequence, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
    // Create lookup table
    std::vector<MeanStd> Table = Create_Table(sequence);
    // Estimate window size and p-value threshold using Hausdorff fractal dimension
    int window_size = Estimate_Scale(sequence,Test,Table); //70% of the runtime is spent here
    //std::cout << "Scale : " << counter << std::endl;
    bool Hack = true; // Take the log candidates?

    // Extract candidates (and threshold)
    std::vector<Changepoint> Candidates = Candidate_extraction( DM(window_size,Test,Table), window_size, Hack);

    std::vector<int> locations = {};
    double threshold =1.;
    for(auto cp: Candidates){
        locations.push_back(cp.i);
        if(cp.p < threshold){
            threshold = cp.p;
        }
    }

    MCP Return_MCP;
    Return_MCP.window = 0;
    Return_MCP.p_value = 1.; // Compute real p_value
    Return_MCP.threshold = threshold;
    Return_MCP.locations = locations;

    return Return_MCP;
}

MCP Stage_random_recomb(std::vector<double> sequence, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
    int number_candidates = std::log2(sequence.size());

    std::mt19937 rng(std::random_device{}());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(3,sequence.size()-3); // guaranteed unbiased

    std::vector<Changepoint> Candidates = {};
    std::vector<int> loc = {};

    for(int i = 0; i < number_candidates; i++){
        loc.push_back(uni(rng));
    }

    std::sort(loc.begin(),loc.end());

    for(int i = 0; i < number_candidates; i++){
        Changepoint cp;
        if(i>0){
            if(loc[i]-loc[i-1]>3){
                cp.i = loc[i];
                cp.p = 1.;
                Candidates.push_back(cp);
            }
        }
        else{
            cp.i = loc[i];
            cp.p = 1.;
            Candidates.push_back(cp);
        }
    }

    // Create lookup table
    std::vector<MeanStd> Table = Create_Table(sequence);
    MCP Return_MCP = Welch_Fisher(sequence, Test,Table,Candidates);

    return Return_MCP;
}

MCP Stage_window_recomb(std::vector<double> sequence, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
    // Create lookup table
    std::vector<MeanStd> Table = Create_Table(sequence);
    // Initialize return type
    MCP Return_MCP;

    // Estimate window size and p-value threshold using Hausdorff fractal dimension
    int window_size = (int) (sequence.size() / 10.);
    //std::cout << "Scale : " << counter << std::endl;
    bool Hack = true; // If overflow occurs at DM?

    // Extract candidates (and threshold)
    std::vector<Changepoint> Candidates = Candidate_extraction( DM(window_size,Test,Table), window_size, Hack);
    //std::cout<< "Candidates : " << counter << std::endl;

    // Recombine candidates and use the configuration with minimal p-value
    Return_MCP = Welch_Fisher(sequence, Test,Table,Candidates);
    Return_MCP.window = window_size;
    //std::cout << "Welch-Fisher : " << counter << std::endl;

    return Return_MCP;
}

MCP Stage_window_recomb_bad(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
    // Create lookup table
    std::vector<MeanStd> Table = Create_Table(sequence);
    // Initialize return type
    MCP Return_MCP;

    // Estimate window size and p-value threshold using Hausdorff fractal dimension
    int window_size = (int) std::log2(sequence.size());
    //std::cout << "Scale : " << counter << std::endl;
    bool Hack = true; // If overflow occurs at DM?

    // Extract candidates (and threshold)
    std::vector<Changepoint> Candidates = Candidate_extraction( DM(window_size,Test,Table), window_size, Hack);
    //std::cout<< "Candidates : " << counter << std::endl;

    // Recombine candidates and use the configuration with minimal p-value
    Return_MCP = Welch_Fisher(sequence, Test,Table,Candidates);
    Return_MCP.window = window_size;
    //std::cout << "Welch-Fisher : " << counter << std::endl;

    return Return_MCP;
}

MCP Stage_window_recomb_random(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
    // Create lookup table
    std::vector<MeanStd> Table = Create_Table(sequence);
    // Initialize return type
    MCP Return_MCP;

    std::mt19937 gen( (int) (sequence[0]*sequence.size()));
    std::uniform_int_distribution<int> uni(4,(int) (sequence.size()/3));
    int window_size = uni(gen);
    std::cout << window_size << std::endl;
    // Estimate window size and p-value threshold using Hausdorff fractal dimension
    //int window_size = (int) std::log2(sequence.size());
    //std::cout << "Scale : " << counter << std::endl;
    bool Hack = true; // If overflow occurs at DM?

    // Extract candidates (and threshold)
    std::vector<Changepoint> Candidates = Candidate_extraction( DM(window_size,Test,Table), window_size, Hack);
    //std::cout<< "Candidates : " << counter << std::endl;

    // Recombine candidates and use the configuration with minimal p-value
    Return_MCP = Welch_Fisher(sequence, Test,Table,Candidates);
    Return_MCP.window = window_size;
    //std::cout << "Welch-Fisher : " << counter << std::endl;

    return Return_MCP;
}

MCP Single_Test(std::vector<double> sequence,  std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){

    std::vector<MeanStd> Table = Create_Table(sequence);


    int window_size = Estimate_Scale(sequence,Test,Table); //70% of the runtime is spent here
    // std::cout << "Scale : " << window_size << std::endl;
    bool Hack = true; // If overflow occurs at DM?

    // Extract candidates (and threshold)
    std::vector<Changepoint> Candidates = Candidate_extraction( DM(window_size,Test,Table), window_size, Hack);

    MCP Return_MCP;
    // Recombine candidates and use the configuration with minimal p-value
    Return_MCP = Welch_Fisher(sequence, Test,Table,Candidates);
    Return_MCP.window = window_size;


    Return_MCP.locations = {300,700,1000};
    MC Return_MC;
    std::vector<MeanStd> list;
    Test(Return_MCP.locations,Table,true,false,Return_MC,list);
    int n = sequence.size();
    int number = 3;
    Return_MCP.p_value = nChoosek(n- (number + 1) * 2 + number, number)*Return_MC.p;
    //std::cout << Return_MCP.p_value ;

    return Return_MCP;




}
