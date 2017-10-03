
//Retrun type, loaded in ./MaChaMP/MaChaMP.cpp
/*
struct MCP{
	std::vector<int>  locations = {}; //indices in sequence
	long double p_value =1.; //p-vlaue
            int window = 2;  // Window-size
            long double threshold; // Threshold of candidates
};
*/


MCP CUSUM(std::vector<double> sequence);






MCP CUSUM(std::vector<double> sequence){

    //Init CUSUM std. control parameters. (CUMULATIVE SUM CONTROL CHARTING: AN UNDERUTILIZED SPC TOOL; DOUGLAS M. HAWKINS). No life or death choice
    double k = 0.25; //  0.75;  // 1.0; // 0.5; //     1.5;  //    0.25; // 0.75; 1.0; 1.25; 1.5;  // \mu \pm k * \sigma
    double h = 8.01;  //  3.34; // 2.52; // 4.77; //    1.61; //   3.34; 2.52; 1.99; 1.61;

    int skip = 3;
    int min_dist = skip;
    int run_length = 0;

    MCP Return_MCP;

    double mean = 0.;
    double M2 = 0;
    double stat_u = 0.;
    double stat_l = 0.;
    double prev_stat_u = 0;
    double prev_stat_l = 0;

    int i = 1;
    while (i<sequence.size()){
        run_length++;
        Welford(mean,M2,run_length,sequence[i-1]);

        if(i>skip){
            // p.418, Introduction to statistical quality control (Montgomery)
            stat_u =  (sequence[i-1] - mean - (k/(M2/run_length))) +prev_stat_u;
            stat_l =  (mean - (k/(M2/run_length)) - sequence[i-1] ) +prev_stat_l;
            if (stat_u < 0){
                stat_u = 0;
            }
            if (stat_l < 0){
                stat_l = 0;
            }

            if (stat_u> h || stat_l > h){
                Return_MCP.locations.push_back(i);
                skip = i+min_dist;

                // Restart: Set values to zero again
                mean = 0.;
                M2 = 0.;
                stat_u = 0.;
                stat_l = 0.;
                run_length = 0;
            }
            prev_stat_u=stat_u;
            prev_stat_l=stat_l;
        }

        i++;
    }

    Return_MCP.threshold = h;

    return Return_MCP;
}

