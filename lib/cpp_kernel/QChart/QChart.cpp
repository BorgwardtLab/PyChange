
//Retrun type, loaded in ./MaChaMP/MaChaMP.cpp
/*
struct MCP{
    std::vector<int>  locations = {}; //indices in sequence
    long double p_value =1.; //p-vlaue
            int window = 2;  // Window-size
            long double threshold; // Threshold of candidates
};
*/


MCP QChart(std::vector<double> sequence);






MCP QChart(std::vector<double> sequence){
    /*
    #QChart, 4of5 test
def QChart(seq,maxlike=5. (my analysis maxlike=6.5)):
    if len(seq) < 7:
        print "Sequence too short QChart"
        return 0

    mean = float(seq[0]+seq[1])/2
    variance = (seq[0]-mean)**2 + (seq[1]-mean)**2
    stat = []
    chart=[0.,0,0] #sum[0] and score[1] of #elements[2]
    cp = 0

    stat_score = [0.]
    for r in xrange(3,len(seq)+1):
        stat.append(stats.norm.ppf(stats.t.cdf(np.sqrt(float(r-1)/r) * (float(seq[r-1]-mean)/np.sqrt(variance)),r-2))) #eq 7, Q Statsistics
        variance = float(r-2)*variance/(r-1) + 1./r*(seq[r-1]-mean)**2 #eq3 in 1991 paper
        mean = float((r-1)*mean+seq[r-1])/r
        chart[0] = chart[0]+ stat[-1]
        chart[1] = chart[1] +int(stat[-1]>0)
        stat_score.append(chart[0])
        if chart[2] < 5:
            chart[2] = chart[2] + 1 #first two: no checks, only increase
        else:
            chart[0] = chart[0]-stat[-5]
            chart[1] = chart[1]-int(stat[-5]>0)
            if (chart[1] < 1 or chart[1] > 4) and np.abs(chart[0]) > maxlike: #resturn most probable location for cp
                cp = r - 6
                maxlike = np.abs(chart[0])
                break #leave this out for making it offline and fpr emwa, enable it for online
    */

    MCP Return_MCP;
    double maxlike = 2.1;
    int burn_in_length = 3;

    double mean = 0.;
    double M2 = 0.;
    double stat = 0;
    double stat_1 = 0;
    double stat_2 = 0;
    double stat_3 = 0;
    double stat_4 = 0;
    double stat_5 = 0;
    double chart_0 = 0.;
    double chart_1 = 0.;

    int i = 0;
    int k = i;
    int run_length = 0;
    int skip = burn_in_length;
    while (i < sequence.size()){
        run_length++;
        if(i>skip){
            stat = std::sqrt((float)(run_length-1)/run_length) * (sequence[i-1]-mean)/std::sqrt(M2/run_length);
            //std::cout << i << " " << std::sqrt(float(run_length-1)/run_length)* (sequence[i-1]-mean)/std::sqrt(M2/run_length) << " "; //here is non zero
            boost::math::students_t_distribution<double> distst(run_length-2);
            boost::math::normal_distribution<double> distnorm;
            if (std::isnan(stat)==false && std::isinf(stat) == false){
                stat = boost::math::cdf(boost::math::complement(distst, fabs(stat)));
                //std::cout << stat << " ";
                stat = boost::math::quantile(boost::math::complement(distnorm,fabs(stat)));
                //std::cout << stat << " "; //here is zero
            }
            else{
                if(std::isnan(stat)==true){
                    //std::cout << i << " " << M2 << " " << run_length << " " << (float)(run_length-1)/run_length << std::endl;
                    stat=stat_1;
                }
                else{
                    stat = stat_1+maxlike;
                }
            }


            // 4 out of 5 test : too many FP
            chart_0 = stat;
            if(stat > maxlike){
                chart_1 = chart_1 +1;
            }
            if(stat_5> maxlike){
                chart_1 = chart_1 -1;
            }
            //if(stat_3 > maxlike){
            //  chart_1 = chart_1 -1;
            //}
            //std::cout << chart_1 << " " << chart_0 << std::endl;
            if ((chart_1 > 3)){//} && (std::abs(chart_0) > maxlike)){
            //if ((chart_1 > 3) && (std::abs(chart_0) > maxlike)){
                Return_MCP.locations.push_back(i);
                skip = i+burn_in_length;

                mean = 0;
                M2 = 0;
                run_length = 1;
                chart_0 = 0;
                chart_1 = 1;
                stat = 0;
                stat_1 = 0;
                stat_2 = 0;
                stat_3 = 0;
            }
        }

        Welford(mean,M2,run_length,sequence[i-1]);
        stat_5 = stat_4;
        stat_4 = stat_3;
        stat_3 = stat_2;
        stat_2 = stat_1;
        stat_1 = stat;
        i++;
    }

    Return_MCP.threshold = maxlike;

    return Return_MCP;
}

