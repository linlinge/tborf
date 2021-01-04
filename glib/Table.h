#pragma once
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <string>
#include "Statistics.h"
#include "PCLExtend.h"
#include "VectorExtend.h"
#include "Color.h"
using namespace std;
class Cluster
{
    public:
        vector<int> c_;
};
class Rrd1
{
	public:
		int id_;
		double item1_;
		Rrd1(){}
		Rrd1(int id,double item1){
			id_=id;
			item1_=item1;
		}
        double GetItem(int index=1){ return item1_;};
        void SetItem(double val, int index=1){ item1_=val;}       
};
class Rrd2
{
    public:
        int id_;
        double item1_;
        double item2_;
        Rrd2(){id_=0;item1_=0;item2_=0;};
        Rrd2(int id,double item1,double item2){
            id_=id;
            item1_=item1_;
            item2_=item2_;
        }
        double GetItem(int index=1){
            if(1==index)
                return item1_;
            else
                return item2_;
        }
        void SetItem(double val,int item_id=1){
            if(1==item_id)
                item1_=val;
            else if(2==item_id)
                item2_=val;
        }
};

template <class Rrd>
class Table
{
    private:
        vector<Rrd> records_sorted_;
        double mean_,sigma_;		
		double min_=INT_MAX;
		double max_=-INT_MAX;
        double median_;	

	public:
		vector<Rrd> records_;
        vector<int> aidx_; // active indices
		vector<V3> color_;

        /**/
        Table(){}
        Table(int n){ Resize(n);}
        /**/
        void Clear();
		void Resize(int n);
		void push_back(Rrd e);
		void Ascending(int ItemID=1);
		void Descending(int ItemID=1);
		void SortBackup(int ItemID=1);		
		void Print(int ItemID=1);
		double PDF(double t,double n,double h);
		double ReversePDF(double P);
		void Find(vector<int>& out,string str=">",double thresh1=0, double thresh2=0);
		void Standardize_Zscore(int ItemID=1);
		void Normalize_Tanh(double scale=1.0, int ItemID=1);		
		void Normalize_Min_Max(int ItemID=1);
		void Write(string path);
        void Multiply(double val, int ItemID=1);
		void Overturn(int ItemID=1);
		void LocalFilter(string str,pcl::PointCloud<PointType>::Ptr cloud,int K);
		void nPLOF(pcl::PointCloud<PointType>::Ptr cloud,int K);
		void SetActiveIndex(string erea_color, string status);
		// Get 
        void GetNormalDistributionError(int ItemID=1);
        void GetCorrespondingColor(int ItemID=1);
		void GetHistogram(int k=30);
        void GetMeanAndVariance(int ItemID=1);
		void GetMinimumAndMaximum(int ItemID=1);
        double GetMaximum(int ItemID=1);
        double GetMinimum(int ItemID=1);
        double GetMean(int ItemID=1);
        double GetQuantile(double p);
		double GetMedian();
		double GetStd();		
        int GetSize();

		// Fast Remove 
		void FastRemove(int index);
		void FastRemove(vector<int> indices);
		
};
void TableHierarchy();


template <class Rrd>
void Table<Rrd>::Resize(int n)
{
    records_.resize(n);
}
template <class Rrd>
void Table<Rrd>::Clear()
{
    aidx_.clear();
    records_.clear();
    records_sorted_.clear();
    color_.clear();
}
template <class Rrd>
int Table<Rrd>::GetSize()
{
    return records_.size();
}
template <class Rrd>
void Table<Rrd>::GetMeanAndVariance(int ItemID)
{    
    double sum=0;
    // calculate average
    for(int i=0;i<records_.size();i++){
        sum+=records_[i].GetItem(ItemID);
    }
    mean_=sum/records_.size();

    // calculate variance
    sum=0;
    for(int i=0;i<records_.size();i++){
        sum+=pow(records_[i].GetItem(ItemID)-mean_,2);
    }
    sigma_=sqrt(sum/(records_.size()-1));

}

template <class Rrd>
void Table<Rrd>::GetMinimumAndMaximum(int ItemID)
{
    min_=INT_MAX;
    max_=-INT_MAX;    
    for(int i=0;i<records_.size();i++){
        double item_tmp=records_[i].GetItem(ItemID);
        if(min_>item_tmp)
            min_=item_tmp;
        if(max_<item_tmp)
            max_=item_tmp;
    }    
}

template<class Rrd>
double Table<Rrd>::GetMinimum(int ItemID)
{
    min_=INT_MAX;    
    for(int i=0;i<records_.size();i++){
        double item_tmp=records_[i].GetItem(ItemID);
        if(min_>item_tmp)
            min_=item_tmp;
    }
    return min_;
}

template<class Rrd>
double Table<Rrd>::GetMean(int ItemID)
{
    double sum=0;
    // calculate average
    for(int i=0;i<records_.size();i++){
        sum+=records_[i].GetItem(ItemID);
    }
    mean_=sum/records_.size();
    return mean_;
}

template <class Rrd>
void Table<Rrd>::Standardize_Zscore(int ItemID)
{
    GetMeanAndVariance();
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++){
        double new_dat_tmp=(records_[i].GetItem(ItemID)-mean_)/sigma_;
        records_[i].SetItem(new_dat_tmp,ItemID);
    }    
}
template <class Rrd>
void Table<Rrd>::Normalize_Min_Max(int ItemID)
{        
    GetMinimumAndMaximum();
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++){
        double new_dat_tmp=(records_[i].GetItem(ItemID)-min_)/(max_-min_);
        records_[i].SetItem(new_dat_tmp,ItemID);
    }        
}
template <class Rrd>
void Table<Rrd>::Normalize_Tanh(double scale,int ItemID)
{     
    GetMeanAndVariance();
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++){
        double new_dat_tmp=tanh(scale*records_[i].item1_);
        records_[i].SetItem(new_dat_tmp,ItemID);
    } 
}
template <class Rrd>
void Table<Rrd>::Ascending(int ItemID)
{
    if(1==ItemID)
        sort(records_.begin(),records_.end(),[](Rrd& e1,Rrd& e2){ return e1.GetItem()<e2.GetItem();});
    else if(2==ItemID)
        sort(records_.begin(),records_.end(),[](Rrd& e1,Rrd& e2){ return e1.GetItem(2)<e2.GetItem(2);});
}
template <class Rrd>
void Table<Rrd>::Descending(int ItemID)
{
    if(1==ItemID)
        sort(records_.begin(),records_.end(),[](Rrd& e1,Rrd& e2){ return e1.GetItem()>e2.GetItem();});
    else if(2==ItemID)
        sort(records_.begin(),records_.end(),[](Rrd& e1,Rrd& e2){ return e1.GetItem(2)>e2.GetItem(2);});
}
template <class Rrd>
void Table<Rrd>::Write(string path)
{
    ofstream fout(path);
    for(int i=0;i<records_.size();i++){
        fout<<records_[i].item1_<<endl;
    }    
    fout.close();
}

template <class Rrd>
void Table<Rrd>::Multiply(double val, int ItemID)
{
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++){
        double new_dat_tmp=records_[i].GetItem(ItemID)*val;
        records_[i].SetItem(new_dat_tmp);
    }
        
}
template <class Rrd>
void Table<Rrd>::Overturn(int ItemID)
{
    GetMaximum();
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++){  
        double new_dat_tmp=max_-records_[i].GetItem(ItemID);
        records_[i].SetItem(new_dat_tmp);      
    }
}
template <class Rrd>
void Table<Rrd>::GetNormalDistributionError(int ItemID)
{
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++){
        double new_dat_tmp=Erf(records_[i].GetItem(ItemID));
        records_[i].SetItem(new_dat_tmp);     
    }
}
template <class Rrd>
void Table<Rrd>::LocalFilter(string str, pcl::PointCloud<PointType>::Ptr cloud,int K)
{
    vector<double> lf;
    lf.resize(cloud->points.size());
    pcl::search::KdTree<PointType>::Ptr kdtree(new pcl::search::KdTree<PointType>());
    kdtree->setInputCloud(cloud);

    if("average"==str){
         #pragma omp parallel for
        for(int i=0;i<cloud->points.size();i++){
            vector<int> idx(K+1);
            vector<float> dist(K+1);        
            kdtree->nearestKSearch(cloud->points[i], K+1, idx, dist);
            double sum=0;
            for(int j=1;j<K+1;j++){
                sum+=records_[idx[j]].item1_;
            }
            sum/=K;
            lf[i]=records_[i].item1_/sum-1.0f;
        }

        for(int i=0;i<cloud->points.size();i++){
            records_[i].item1_=lf[i];
        }
    }
    else if("std"==str){
        for(int i=0;i<cloud->points.size();i++){
            vector<int> idx(K+1);
            vector<float> dist(K+1);        
            kdtree->nearestKSearch(cloud->points[i], K+1, idx, dist);
            vector<double> tmp;
            for(int j=1;j<K+1;j++){
                tmp.push_back(abs(records_[i].item1_-records_[idx[j]].item1_));
            }
            lf[i]=VectorSum(tmp);
            // lf[i]=VectorStd(tmp);
        }

        for(int i=0;i<cloud->points.size();i++){
            records_[i].item1_=lf[i];
        }
    }
   
}


template <class Rrd>
double Table<Rrd>::GetQuantile(double p)
{   
    SortBackup(1);
    double Q_idx=1+(records_sorted_.size()-1)*p;
    int Q_idx_integer=(int)Q_idx;
    double Q_idx_decimal=Q_idx-Q_idx_integer;
    double Q=records_sorted_[Q_idx_integer-1].item1_+(records_sorted_[Q_idx_integer].item1_-records_sorted_[Q_idx_integer-1].item1_)*Q_idx_decimal;    
    return Q;
}
template <class Rrd>
double Table<Rrd>::GetMedian()
{
    SortBackup(1);
    double median=0;
    if(records_sorted_.size()%2==1){// odd
        median= records_sorted_[(records_sorted_.size()+1.0)/2.0-1].item1_;
    }
    else{ // even        
        median=(records_sorted_[records_sorted_.size()/2-1].item1_+records_sorted_[records_sorted_.size()/2].item1_)/2.0f;
    }
    median_=median;
    return median;
}
template <class Rrd>
void Table<Rrd>::SortBackup(int item_index)
{
    // Step 1: Does it have backup ?
    if(records_sorted_.size()==0){
        records_sorted_.resize(records_.size());
        for(int i=0;i<records_.size();i++){
            records_sorted_[i].id_=records_[i].id_;
            records_sorted_[i].item1_=records_[i].item1_;
        }
    }

    // Step 2: Sort Backup Records
    if(item_index==1){
        sort(records_sorted_.begin(),records_sorted_.end(),[](Rrd& e1, Rrd& e2){ return e1.item1_<e2.item1_;});
    }
    else if(item_index==0){
        sort(records_sorted_.begin(),records_sorted_.end(),[](Rrd& e1, Rrd& e2){ return e1.id_<e2.id_;});
    }
}

template <class Rrd>
void Table<Rrd>::push_back(Rrd e)
{
    records_.push_back(e);
}
template <class Rrd>
void Table<Rrd>::Print(int ItemID)
{
    for(int i=0;i<records_.size();i++){
        cout<<records_[i].GetItem(ItemID)<<" ";
    }
    cout<<endl;
}
template <class Rrd>
void Table<Rrd>::GetCorrespondingColor(int ItemID)
{
    color_.resize(records_.size());
    GetMinimumAndMaximum();
    for(int i=0;i<records_.size();i++){
        V3 ctmp=get_color(min_,max_,records_[i].GetItem(ItemID));
        color_[i].r=ctmp.r;
        color_[i].g=ctmp.g;
        color_[i].b=ctmp.b;
    }
}
template <class Rrd>
void Table<Rrd>::nPLOF(pcl::PointCloud<PointType>::Ptr cloud,int K)
{
    pcl::search::KdTree<PointType>::Ptr kdtree(new pcl::search::KdTree<PointType>());
    kdtree->setInputCloud(cloud);
    
    vector<double> plof;    
    plof.resize(cloud->points.size());
	// #pragma omp parallel for
	for (int i = 0; i < cloud->points.size(); i++){        
        vector<int> idx(K+1);
		vector<float> dist(K+1);
		kdtree->nearestKSearch (cloud->points[i], K+1, idx, dist);
        double sum = 0;
        for (int j = 1; j < K+1; j++)
          sum += records_[idx[j]].item1_;
        sum /= K;
        plof[i] = records_[i].item1_ / sum  - 1.0f;
    }

    for(int i=0;i<cloud->points.size();i++){
        records_[i].item1_=plof[i];
    }
}

template <class Rrd>
double Table<Rrd>::PDF(double t,double n,double h)
{
    double rst=0;
    for(int i=0;i<records_.size();i++){
        double tmp=sqrt(2)/2.0*(t-records_[i].item1_)/h;
        rst+=Erf(tmp);
    }
    rst=rst/(2.0*n)+1/2.0;
    return rst;
}
template <class Rrd>
double Table<Rrd>::GetStd()
{
    double rst=0;
    double vmean=GetMean();
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++) 
        rst+=pow(records_[i].item1_-vmean,2);
    return sqrt(rst/(records_.size()-1));
}
template <class Rrd>
double Table<Rrd>::ReversePDF(double P)
{
    // double IQR=Quantile(0.75)-Quantile(0.25);
    // vector<int> new_indices;
    // new_indices=GetIndex(">=",Quantile(0.75)+3.0*IQR);
    // FastRemove(new_indices);

    // update parameters
    double n=records_.size();
    double IQR=GetQuantile(0.75)-GetQuantile(0.25);
    double h=0.9*pow(n,-0.2)*min(GetStd(),IQR/1.34);

    double t0=GetMinimum();
    double t1=GetMean();
    double t2=GetMaximum();
    double ptmp=PDF(t1,n,h);
    while(abs(ptmp-P)>0.001){
        if(ptmp>P)
            t2=t1;
        else
            t0=t1;
        t1=(t0+t2)/2.0;
        ptmp=PDF(t1,n,h);
    }
    return t1;
}
template <class Rrd>
void Table<Rrd>::Find(vector<int>& out, string str,double thresh1, double thresh2)
{
    if(">"==str){
       double t1=ReversePDF(thresh1);
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_>t1)
                out.push_back(i);
        }
    }
    else if("<"==str){
        double t1=ReversePDF(thresh1);
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_<t1)
                out.push_back(i);
        }
    }
    else if("t1<x<t2"==str){
        double t1=ReversePDF(thresh1);
        double t2=ReversePDF(thresh2);
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_>t1 && records_[i].item1_<t2){
                out.push_back(i);
            }
        }
    }
    else if("x<t1 || x>t2"==str){
        double t1=ReversePDF(thresh1);
        double t2=ReversePDF(thresh2);
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_<t1 || records_[i].item1_>t2){
                out.push_back(i);
            }
        }
    }
}

template <class Rrd>
double Table<Rrd>::GetMaximum(int ItemID)
{
    max_=-INT_MAX;
    for(int i=0;i<records_.size();i++){
        if(max_<records_[i].GetItem(ItemID))
            max_=records_[i].GetItem(ItemID);
    }        
    return max_;
}

template <class Rrd>
void Table<Rrd>::FastRemove(int index)
{
    int tail_idx=records_.size()-1;
    records_[index].id_=records_[tail_idx].id_;
    records_[index].item1_=records_[tail_idx].item1_;
    records_.pop_back();
}
template <class Rrd>
void Table<Rrd>::FastRemove(vector<int> indices)
{
    for(int i=0;i<indices.size();i++){
        FastRemove(indices[i]);
    }
}