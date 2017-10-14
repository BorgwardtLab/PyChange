
//Retrun type, loaded in ./MaChaMP/MaChaMP.cpp
/*
struct MCP{
	std::vector<int>  locations = {}; //indices in sequence
	long double p_value =1.; //p-vlaue
            int window = 2;  // Window-size
            long double threshold; // Threshold of candidates
};
*/


MCP EWMA(std::vector<double> sequence);






MCP EWMA(std::vector<double> sequence){
    /*
    def EMWA(seq, parameters=[0.11, 3.74], onset=30): #does not work
    l = parameters[0] #recommended 0.05-0.2
    L = parameters[1] # recommended 3 for l<0.1, 2.6-2.8 for l>0.1
    k = onset #onset
    if len(seq) < k:
        print "Sequence too short EWMA"
        return 0
    hit = 0
    mu0 = np.mean(seq[:k])
    z = mu0
    test = [z]
    sig = np.std(seq[:k])
    diff = L * sig * np.sqrt(l/(2-l)) #statistical quality control 9.27 9.28
    for i in xrange(0,k):
        z = l*seq[i] + (1-l)*z #statistical quality control 9.22
        test.append(z)
    for i in xrange(k,len(seq)):
        z = l*seq[i] + (1-l)*z #statistical quality control 9.22
        test.append(z)
        if np.abs(z-mu0)>diff:
            hit = i
            #check = diff-np.mean(seq[:i])-np.abs(z)
            break
    return hit#, test,mu0,diff
    */

    MCP Return_MCP;
    double l = 0.05;
    double L = 3.7;
    int burn_in_length = 30;

    double mean = 0.;
    double mean0 = 0;
    double M2 = 0.;
    double diff = 0;
    double z = 0.;

    int i = 0;
    int k = i;
    int run_length = 0;
    int skip = burn_in_length;
    while (i < sequence.size()){
        run_length++;
        Welford(mean,M2,run_length,sequence[i-1]);
        if(i>skip){
            z = l*sequence[i] + (1-l)*z;
            if(std::abs(z-mean0)>diff){
                Return_MCP.locations.push_back(i);
                skip = i+burn_in_length;

                mean = 0;
                M2 = 0;
                run_length = 0;
            }
        }
        else if(i==skip){
            mean0 = mean;
            z = mean;
            k = i-burn_in_length;
            while (k<=burn_in_length){
                z = l*sequence[k] + (1-l)*z;
                k++;
            }

            diff = L * std::sqrt(l/(2-l)) *M2/(run_length) ; //stat quality control 9.27&9.28
        }

        i++;
    }

    Return_MCP.threshold = burn_in_length;

    return Return_MCP;
}

