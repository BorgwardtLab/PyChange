#include "MaChaMP.cpp"

int main(){
    std::mt19937 gen(0);
    std::normal_distribution<double> norm(0,1.);
    std::vector<double> seq;
    for(int k = 0; k<1000000;k++){
        double append = norm(gen);
        seq.push_back(append);
    }

    MCP Result = MaChaMP(seq,Multiple_TTest);

    return 0;

}

