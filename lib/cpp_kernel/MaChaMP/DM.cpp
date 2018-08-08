#include "Multiple_TTest.cpp"


std::vector<long double> DM( int window, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test, std::vector<MeanStd> Table);



std::vector<Changepoint> Candidate_extraction(std::vector<long double> dissimilarity_manifold, int window, bool Hack);
std::tuple<std::vector<long double>,std::vector<int>> log_sample(std::vector<int> rand_index, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test, std::vector<MeanStd> Table, int window);


//What we want is a random subsequence of DM
std::tuple<std::vector<long double>,std::vector<int>> log_sample(std::vector<int> rand_index, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test, std::vector<MeanStd> Table, int window){
    std::mt19937 gen(Table.size());
    std::uniform_int_distribution<int> uni(0,Table.size());
    int start = uni(gen) % rand_index.size();

    std::vector<int> index;
    int length = ((int) std::log2(Table.size())+0.5)+1.;
    for (int k = 0; k < length; k++) {
        int current = rand_index[start+k];
        if (current < Table.size() - window && current > window){
            index.push_back(current);
        }
        else{
            length++;
        }
    }
    index.resize( ((int) std::log2(Table.size())+0.5)+1.);
    std::sort(index.begin(),index.end());

    std::vector<long double> DM_sample;
    MC Test_Return;
    std::vector<MeanStd> List_of_Stats(4);
    for(auto const& it : index){
        std::vector<int> loc = {it-window, it, it+window};
        Test(loc,Table,true,true,Test_Return,List_of_Stats);
        DM_sample.push_back(std::log(Test_Return.p));
        //DM_sample.push_back(std::log(Test(loc,Table,true,true).p));
    }
    return std::make_tuple(DM_sample,index);
}


std::vector<long double> DM(int window, std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test, std::vector<MeanStd> Table){
    std::vector<long double> dissimilarity_manifold(Table.size()-2*window-4);
    MC Test_Return;
    std::vector<MeanStd> List_of_Stats(4);
    std::vector<int> loc(3);

    for(int l = window+2; l < static_cast<int>(Table.size())-window-2; l++){
        loc = {l-window, l, l+window};
        Test(loc,Table,true,true,Test_Return,List_of_Stats);
        dissimilarity_manifold[l-window-2] = std::log(Test_Return.p);
        //dissimilarity_manifold[l-window-2] = std::log(Test(loc,Table,true,true).p);
    }

    return dissimilarity_manifold;
}




std::vector<Changepoint> Candidate_extraction(std::vector<long double> dissimilarity_manifold, int window, bool Hack){
    std::vector<Changepoint> Candidates;

    long double Reference;
    if (Hack == false){
        Reference = Percentile(dissimilarity_manifold.begin(), dissimilarity_manifold.end(), 0.01);
    }
    else{
        Reference = Percentile(dissimilarity_manifold.begin(), dissimilarity_manifold.end(), 0.5);
    }

    std::vector<int> cand_ind {};
    std::vector<long double> cand_val {};

    int ilocal_max = 0;
    long double local_max = Reference;

    for(int pos = 0; pos < dissimilarity_manifold.size(); pos++){
        if (dissimilarity_manifold[pos] < local_max){
            local_max = dissimilarity_manifold[pos];
            ilocal_max = pos;
        }

        if ((dissimilarity_manifold[pos] >= Reference) || (pos == dissimilarity_manifold.size()-1)  ){
            if (local_max < Reference){
                cand_ind.push_back(ilocal_max+window+2);
                cand_val.push_back(Chisq(local_max,2));
            }
            ilocal_max = pos;
            local_max = Reference;
        }

    }


    std::vector<size_t> index_cand_val = topk_index(cand_val, (int) std::log2(static_cast<int>(dissimilarity_manifold.size())));

    std::sort(index_cand_val.begin(),index_cand_val.end());

    for(int pos = 0; pos < index_cand_val.size(); pos++){
        Changepoint next;
        next.i = cand_ind[index_cand_val[pos]];
        next.p = cand_val[index_cand_val[pos]];
        Candidates.push_back(next);

    }


    return  Candidates;
}
