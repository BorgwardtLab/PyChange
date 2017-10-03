#include <iomanip>
#include <iostream>
#include <boost/math/distributions/students_t.hpp> //dependency difficult for windows


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
	double v = Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2;
	v *= v;
	double t1 = Sd1 * Sd1 / Sn1;
	t1 *= t1;
	t1 /=  (Sn1 - 1);
	double t2 = Sd2 * Sd2 / Sn2;
	t2 *= t2;
	t2 /= (Sn2 - 1);
	v /= (t1 + t2);
	// t-statistic:
	double t_stat = (Sm1 - Sm2) / sqrt(Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2);
	boost::math::students_t dist(v);
	long double q = 2* boost::math::cdf(boost::math::complement(dist, fabs(t_stat))); //this is the bottleneck
	// q is a p-value
	return q;
}

