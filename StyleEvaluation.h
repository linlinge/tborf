/*
    Name:   Outliers Level Evaluation Module
    Author: linlinge
    Date:   2020.04.02
*/
#pragma once
#include "Table.h"
#include "PCLExtend.h"
#include <stdlib.h>
#include <time.h>
#include <numeric>
class StyleEvaluation
{
    public:
		pcl::PointCloud<PointType>::Ptr cloud_;
		pcl::search::KdTree<PointType>::Ptr kdtree_;
		Table<Rrd1> rst_dmean;
		Table<Rrd1> rst_meval;
		Table<Rrd1> rst_db2;
		Table<Rrd1> rst_singularity;
		double dnst_,dmean_,meval_,db2_;
		double dmin_;
		double db2_Q3_,db2_IQR_;
		double meval_Q3_,meval_IQR_;
		double og_;
		int n_,N_;
		/* Style Evaluation Metrics */
		void OutlierGrade(string path, int K=30);
		void Homogeneity(string path, int K=30);
		void Sigularity(string path, int K1, int K2);

        double ApplyEigenvalue(int K=32);
		double GetMinorEval(int K=32, string str="Common");
		double GetDmean(int K=32,string str="Common");
		double GetDB2(int K=100, string str="Common");
		StyleEvaluation(pcl::PointCloud<PointType>::Ptr cloud,pcl::search::KdTree<PointType>::Ptr kdtree=NULL);
};