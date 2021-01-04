#include "SignalProcessing.h"
void DaubechiesWavelet(vector<double>& dat,vector<double>& cA, vector<double>& cD)
{
    double h0=(1+sqrt(3))/(4*sqrt(2));
    double h1=(3+sqrt(3))/(4*sqrt(2));
    double h2=(3-sqrt(3))/(4*sqrt(2));
    double h3=(1-sqrt(3))/(4*sqrt(2));
    double g0=h3;
    double g1=-h2;
    double g2=h1;
    double g3=-h0;
    cA.resize(dat.size()/2);
    cD.resize(dat.size()/2);
    for(int i=0;i<cA.size();i++){
        cA[i]=h0*dat[2*i]+h1*dat[2*i+1]+h2*dat[2*i+2]+h3*dat[2*i+3];
        cD[i]=g0*dat[2*i]+g1*dat[2*i+1]+g2*dat[2*i+2]+g3*dat[2*i+3];
    }
}

void DaubechiesWavelet(vector<float>& dat,vector<double>& cA, vector<double>& cD)
{
    double h0=(1+sqrt(3))/(4*sqrt(2));
    double h1=(3+sqrt(3))/(4*sqrt(2));
    double h2=(3-sqrt(3))/(4*sqrt(2));
    double h3=(1-sqrt(3))/(4*sqrt(2));
    // double h0=0.6830127;
    // double h1=1.1830127;
    // double h2=0.3169873;
    // double h3=-0.1830127;
    double g0=h3;
    double g1=-h2;
    double g2=h1;
    double g3=-h0;
    int N=dat.size();
    int N_half=N/2;

    cA.resize(N_half-1);
    cD.resize(N_half);
    int i=0,j=0;
    for(;i<N_half-1;j+=2,i++){
        cA[i]=h0*dat[j]+h1*dat[j+1]+h2*dat[j+2]+h3*dat[j+3];
        cD[i]=g0*dat[j]+g1*dat[j+1]+g2*dat[j+2]+g3*dat[j+3];
    }
    // cA[i]=h0*dat[N-2]+h1*dat[N-1]+h2*dat[0]+h3*dat[1];
    // cD[i]=g0*dat[N-2]+g1*dat[N-1]+g2*dat[0]+g3*dat[1];
}