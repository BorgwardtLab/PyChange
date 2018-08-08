#include <iomanip>
#include <iostream>
#include "student.hpp" //dependency difficult for windows

//int counter = 0;

//boost::math::students_t_distribution<long double> distst(1);

long double two_samples_t_test_unequal_sd(double Sm1,double Sd1,unsigned Sn1,double Sm2,double Sd2,unsigned Sn2);
// http://www.boost.org/doc/libs/1_36_0/libs/math/doc/sf_and_dist/html/math_toolkit/dist/stat_tut/weg/st_eg/two_sample_students_t.html
// also known as welchs ttest
long double two_samples_t_test_unequal_sd(
        double Sm1,       	// Sm1 = Sample 1 Mean.
        double Sd1,       	// Sd1 = Sample 1 Standard Deviation.
        unsigned Sn1,     		// Sn1 = Sample 1 Size.
        double Sm2,       	// Sm2 = Sample 2 Mean.
        double Sd2,       	// Sd2 = Sample 2 Standard Deviation.
        unsigned Sn2)     		// Sn2 = Sample 2 Size.
{
	// Degrees of freedom:
    if (Sn1 == 0 || Sn2 == 0){
        return 1.;
    }
    else{
	double v = Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2;
	v *= v;
	double t1 = Sd1 * Sd1 / Sn1;
	t1 *= t1;
	t1 /=  (Sn1 - 1);
	double t2 = Sd2 * Sd2 / Sn2;
	t2 *= t2;
	t2 /= (Sn2 - 1);
	v /= (t1 + t2);
            boost::math::students_t_distribution<long double> distst(v);
            //distst.df_ = v;
	// t-statistic:
	double t_stat = (Sm1 - Sm2) / sqrt(Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2);	// boost::math::students_t_distribution<long double,boost::math::policies::policy<boost::math::policies::digits2<0> >> dist(v);
            //counter++;
	long double q = 2* boost::math::cdf(boost::math::complement(distst, fabs(t_stat)));
	// q is a p-value
	return q;
    }
}

