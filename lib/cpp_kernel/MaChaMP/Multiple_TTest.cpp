#include "Fisher.cpp"


void Multiple_TTest(std::vector<int> &locations,std::vector<MeanStd> &Table, bool Fish, bool Mid, MC& Returning, std::vector<MeanStd>& List_of_Stats);
bool check_setting(std::vector<int> locations, std::vector<double> seq);
std::vector<MeanStd> Create_Table(std::vector<double>seq);
MeanStd updating(double m, double m1, double v, double v1,double n, double n1);

//Are locations sufficiently far away spaced and are these valid idices?
bool check_setting(std::vector<int> locations, std::vector<double> seq,unsigned offset){

	for (std::vector<int>::size_type i = 1; i != locations.size(); i++){
		if (locations[i]-locations[i-1] < offset){
			return false;
		}
	}


	if (locations[locations.size()-1] >= static_cast<int>(seq.size())-offset){
		return false;
	}
	if (locations[0] <= offset){
		return false;
	}

	return true;
}

//One pass to create Table
std::vector<MeanStd> Create_Table(std::vector<double>seq){
	std::vector<MeanStd> Table;
	double mean = 0;
	double M2 = 0;
	double delta = 0;
	double delta2 = 0;

	MeanStd zero;
	Table.push_back(zero);

	for(std::vector<int>::size_type i = 0; i < seq.size(); i++) {
		delta = seq[i] - mean;
		mean += delta/(i+1);
		delta2 = seq[i] - mean;
		M2 += delta*delta2;
		MeanStd at_i;
		at_i.m = mean;
		at_i.s = M2/(i+1); //biased VARIACNE

		Table.push_back(at_i);
	}
	return Table;
}

//Equation 1.5 from Chan, Golub, LeVeque
MeanStd updating(double m, double m1, double v, double v1,double n, double n1){
	double T = n*m;
	double T1 = n1*m1;
	double S = v*(n); //biased
	double S1 = v1*(n1); //biased
	double T2 = T - T1;
	double S2 = S - S1 - std::pow(T1*n/n1 - T,2) *n1/(n*(n-n1));
	MeanStd MS2;
	MS2.m = T2/(n-n1);
	MS2.s = S2/(n-n1);
	return MS2;
}

std::vector<long double> Test_lookup(std::vector<int> Candidates, std::vector<MeanStd> &Table,std::function<void(std::vector<int> &,std::vector<MeanStd> &, bool, bool,MC&, std::vector<MeanStd>& )> Test){
	Candidates.insert(Candidates.begin(),0);
	Candidates.push_back(Table.size()-1);
	std::vector<long double> TTest_vals(Candidates.size()*Candidates.size()*Candidates.size());
	MC Testresult;
	std::vector<MeanStd> List_of_Stats;
	std::vector<int> loc(3);
	for(int i = 0; i<Candidates.size()-2; i++){
		for (int j = i+1; j< Candidates.size()-1; j++){
			for (int k=j+1; k < Candidates.size(); k++){
				if (i != 0 && k != Candidates.size()-1){
					loc = {Candidates[i],Candidates[j],Candidates[k]};
					Test(loc,Table,false,true,Testresult,List_of_Stats);
					TTest_vals[i+j*Candidates.size()+k*Candidates.size()*Candidates.size()]=Testresult.p;
				}
				else if (i == 0 && k == Candidates.size()-1){
					loc = {Candidates[j]};
					Test(loc,Table,false,false,Testresult,List_of_Stats);
					TTest_vals[i+j*Candidates.size()+k*Candidates.size()*Candidates.size()]=Testresult.Change[0].p;
				}
				else if (i == 0){
					loc = {Candidates[j],Candidates[k]};
					Test(loc,Table,false,false,Testresult,List_of_Stats);
					TTest_vals[i+j*Candidates.size()+k*Candidates.size()*Candidates.size()]=Testresult.Change[0].p;
				}
				else if (k == Candidates.size()-1){
					loc = {Candidates[i],Candidates[j]};
					Test(loc,Table,false,false,Testresult,List_of_Stats);
					TTest_vals[i+j*Candidates.size()+k*Candidates.size()*Candidates.size()]=Testresult.Change[1].p;
				}
				else{
					std::cout << "This should not have happened" << i << " " << j << " " << k << std::endl;
				}
			}
		}
	}

	return TTest_vals;
}

//This routine returns a vector of p-values corresponding to the tested locations
void Multiple_TTest(std::vector<int> &locations, std::vector<MeanStd> &Table, bool Fish, bool Mid, MC& Returning, std::vector<MeanStd>& List_of_Stats){
	List_of_Stats.resize(locations.size()+1);
	MeanStd First;
	MeanStd Full;

	int before,after,division;
	Changepoint next_changepoint;
	Returning.Change.resize(locations.size()); //Linear on the number of elements inserted/erased
	//Iterate through locations
	for(int i = 0; i != locations.size(); i++) {
		if (i == 0){
			List_of_Stats[0] = Table[locations[i]];
			before = 0;
		}
		else {
			before = locations[i-1];
		}
		if (i == locations.size()-1){
			after = static_cast<int>(Table.size())-1;
		}
		else{
			after = locations[i+1];
		}
		division = locations[i];

		Full = Table[after];
		First = Table[division];
		List_of_Stats[i+1] = updating(Full.m, First.m, Full.s, First.s, after, division);

		next_changepoint.i = division;
		if(Mid == false || i==1){
			next_changepoint.p = two_samples_t_test_unequal_sd( List_of_Stats[i].m, std::sqrt(List_of_Stats[i].s), division-before, List_of_Stats[i+1].m, std::sqrt(List_of_Stats[i+1].s), after-division);
		}
		else{
			next_changepoint.p=1.;
		}
		Returning.Change[i]  = next_changepoint;
	}

	//MC Returning;
	if (Mid == false){
		if (Fish==true){
			Returning.p = Fisher(Returning.Change);
		}
		else{
			Returning.p = LogComb(Returning.Change);
		}
	}
	else{
		Returning.p = Returning.Change[1].p;
	}

	//return Returning;
}

