#include "Statistics.h"
Statistics::Statistics(vector<float> dat)
{
	dat_=dat;		
	float sum=0;
	float min=INT_MAX;
	float max=-INT_MAX;
	
	// statistics
	for(int i=0;i<dat.size();i++)
	{
		sum+=dat[i];
		min=min<dat[i] ? min:dat[i];
		max=max>dat[i] ? max:dat[i];
	}
	sum_=sum;
	min_=min;
	max_=max;
	mean_=sum_/dat.size();
	
	sum=0;
	for(int i=0;i<dat.size();i++)
	{
		sum+=pow(dat[i]-mean_,2);				
	}
	stdevp_=sqrt(sum/dat.size());
	stdev_=sqrt(sum/(dat.size()-1));
}

double GaussErrorFunction(double x)
{	
	double denominator=1+0.278393*x+0.230389*x*x+0.000972*pow(x,3)+0.078108*pow(x,4);
	double rst=1-1.0/pow(denominator,4);
	// double rst=1-1.0/denominator;
	return rst;
}

double LossFunc(double x)
{
	return max((double)0.0f,GaussErrorFunction(x/sqrt(2)));
}

double GaussianKernel(double x)
{
	return 1.0/sqrt(2.0*M_PI)*exp(-x*x/2.0);
}
double Erf(double x)
{
	double rst=0;
	if(x>=0){		
		double p = 0.3275911;
		double a1 = 0.254829592;
		double a2 = -0.284496736;
		double a3 = 1.421413741;
		double a4 = -1.453152027;
		double a5 = 1.061405429;
		double t=1.0/(1+p*x);
		rst=1-(a1*t+a2*t*t+a3*pow(t,3)+a4*pow(t,4)+a5*pow(t,5))*exp(-x*x);
	}
	else{
		rst=-1.0*erf(-x);
	}
	return rst;
}

double Quantile(vector<int>& dat,double ratio)
{
	vector<int> dat2;
	for(int i=0;i<dat.size();i++) dat2.push_back(dat[i]);
	sort(dat2.begin(),dat2.end());

	double Q_idx=(dat.size()+1)*ratio-1;
    int Q_idx_integer=(int)Q_idx;
    double Q_idx_decimal=Q_idx-Q_idx_integer;
    double Q=dat2[Q_idx_integer]+(dat2[Q_idx_integer+1]-dat2[Q_idx_integer])*Q_idx_decimal;    
    return Q;
}