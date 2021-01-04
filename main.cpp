/*
	Author: linlinge
	Date: 2020.04.13
*/
#include <iostream>	
#include "PCLExtend.h"
#include "V3.hpp"
#include <Eigen/Dense>
#include <omp.h>
#include "Table.h"
#include <algorithm>
#include "VectorExtend.h"
#include "SignalProcessing.h"
#include "StyleEvaluation.h"
#include <time.h>
#include "stdlib.h"
#include "time.h"   
#include "HybridMethods.h" 
#include "Statistics.h"
#include "StringExtend.h"
#include <random>

double gaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
             
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
         
    phase = 1 - phase;
 
    return X;
}

double F(double x)
{
	double y=0;
	if(x<=80 && x>0){
		y=0.009375;
	}
	else if(x>=80 && x<=240){
		y=0.0015625;
	}
	return y;
}

int main(int argc,char** argv)
{
	string Mtype,input_path,output_path,json_path;
	string filename;
	string fileNoType;

	/* parameters */
	for(int i=1;i<argc;i++){
		string param=argv[i];
		if("-Mtype"==param || "-M"==param){
			Mtype=argv[i+1];
		}
		else if("-i"==param){
			input_path=argv[i+1];
		}
		else if("-o"==param){
			output_path=argv[i+1];
		}
		else if("-json"==param || "-J"==param){
			json_path=argv[i+1];
		}
		else if("-show"==param){
			is_show_progress=atoi(argv[i+1]);
		}
		else if("--help"==param || "-H"== param){
			cout<<"MGS -i path/to/input -o path/to/output -Mtype Methods_Type -json path/of/json/file"<<endl;
			cout<<"Mtype: OGEM UEM SEM BSC1_1 BSC1_2 BSC1_3 MGS1 TBCM2"<<endl;
			return 0;
		}
	}
	/* is show information*/
	if(is_show_progress==1){
		vector<string> str_tmp;
		StrSplit(input_path,"/",str_tmp);
		filename=str_tmp[str_tmp.size()-1];
		str_tmp.clear();
		StrSplit(input_path,".",str_tmp);
		fileNoType=str_tmp[0];
		cout<<"**********************************"<<endl;
		cout<<filename<<endl;
	}		

	/* Load Model */
	pcl::PointCloud<PointType>::Ptr cloud(new pcl::PointCloud<PointType>);	
	if (pcl::io::loadPLYFile<PointType>(input_path, *cloud) == -1){
		PCL_ERROR("Couldn't read file test_pcd.pcd \n");
		return (-1);
	}	
	pcl::search::KdTree<PointType>::Ptr kdtree(new pcl::search::KdTree<PointType>());
	kdtree->setInputCloud(cloud);


	/* Hybrid Methods Generation */	
	StyleEvaluation ole(cloud,kdtree);
	HybridMethods hrd01(cloud,kdtree,&ole);	
	if("OutlierGrade"==Mtype){		
		ole.OutlierGrade(output_path,50);
	}
	else if("meval"==Mtype){
		hrd01.FM_MEval(50,3.0);
		hrd01.DemonstrateResult_Color("Result/color.ply",hrd01.rst_meval_,"nonlinear");
	}
	else if("Homogeneity"==Mtype){		
		ole.Homogeneity(output_path);
	}
	else if("Singularity"==Mtype){
		ole.Sigularity(output_path,60,2000);
	}
	else if("TBCM2"==Mtype){
		// cout<<"TBCM2"<<endl;	
		hrd01.FM_RegionGrowth(9,40);	
		hrd01.FM_NID(100,0.75);		
		hrd01.FM_LoOP("1");		
		hrd01.DemonstrateResult(output_path);
	}
	else if("TBCM3"==Mtype){
		// cout<<"TBCM3"<<endl;
		/* Tanks and Temples */		
		hrd01.FM_Slope(150,5.5);
		hrd01.FM_RegionGrowth(8,200.0,"-1");
		// hrd01.FM_NID(100,0.99,"1,2");			
		hrd01.FM_NID(100,0.9,"1,2");	
		hrd01.DemonstrateResult(output_path);

		/* DTU scan*/
		// hrd01.store_path_=output_path;
		// hrd01.FM_Slope(800,3);
		// hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		// hrd01.FM_NID(100,0.8,"1,2");
		// hrd01.DemonstrateResult(output_path);


		/* DTU furu*/				
		// hrd01.FM_Slope(150,3);
		// hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		// hrd01.FM_NID(100,0.75,"1,2");
		// hrd01.DemonstrateResult(output_path);

		/* THU */		
		// hrd01.store_path_=output_path;
		// hrd01.FM_Slope(150,3);
		// hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		// hrd01.FM_NID(100,0.999,"1,2");
		// hrd01.DemonstrateResult(output_path);
	}
	else if("TBCM4"==Mtype){
		/* THU */		
		hrd01.FM_MEval(400,1.5);
		hrd01.FM_MajorityVote(4000,0.01,"-1");
		hrd01.FM_D4(500,30,"1,2");
		hrd01.FM_RegionGrowth(9,200.0,"-3");
		hrd01.FM_NID(100,0.7,"1,2,3,4");
		hrd01.DemonstrateResult(output_path);

		/* Tanks and Temples*/
		// HybridMethods hrd01(cloud,kdtree,&ole); 
		// hrd01.store_path_=output_path;
		// hrd01.FM_MEval(100,9);
		// hrd01.FM_MajorityVote(4000,"-1");
		// hrd01.FM_D4(500,30,"1,2");
		// hrd01.FM_RegionGrowth(2.0,80.0,200.0,"-3");
		// hrd01.FM_NID(30,0.7,"1,2,3,4");
		// hrd01.DemonstrateResult(fileNoType+"_TBCM2");
	}
	else if("sampling"==Mtype){		
		double range=240;
		double n_of_period=2.0;
		std::random_device rd;
    	std::mt19937 gen(rd());
    	std::uniform_real_distribution<> dis(0, 1);
		vector<double> dist_remove;
		int n_remove=2000000;
		dist_remove.resize(n_remove);
		#pragma omp parallel for
		for(int i=0;i<n_remove;i++){
			int accept=0;
			while(accept==0){
				double u=dis(gen)*range;
				double z=dis(gen)*range;
				
				if(u<=F(z)/(1.0/range)){
					accept=1;
					dist_remove[i]=z;
				}
			}			
		}			

		pcl::PointCloud<PointType>::Ptr pseed(new pcl::PointCloud<PointType>);
		// stl024_cropped
		// pseed->points.push_back(PointType(-79.613190,72.122627,498.238525,0,0,0,0.0,0.0,0.0));
		// pseed->points.push_back(PointType(-121.923134,159.563110,676.015442,0,0,0,0.0,0.0,0.0));
		// pseed->points.push_back(PointType(34.580273,-52.410431,659.988037,0,0,0,0.0,0.0,0.0));
		// pseed->points.push_back(PointType(-47.711945,93.107819,611.769836,0,0,0,0.0,0.0,0.0));

		// stl006_cropped		
		// pseed->points.push_back(PointType(28.291155,-30.487371,652.664856,0,0,0,0.0,0.0,0.0));
		// pseed->points.push_back(PointType(14.794615,-110.490906,637.100525,0,0,0,0.0,0.0,0.0));

		// stl001_cropped		
		pseed->points.push_back(PointType(53.114029,-1.290610,611.314880,0,0,0,0.0,0.0,0.0));		

		pcl::search::KdTree<PointType>::Ptr kdtree (new pcl::search::KdTree<PointType>());
		kdtree->setInputCloud(cloud);
		vector<int> idx_remove;
		idx_remove.resize(n_remove*pseed->points.size());

		
		for(int j=0;j<pseed->points.size();j++){
			cout<<j<<endl;
			int base=j*n_remove;
			#pragma omp parallel for
			for(int i=0;i<dist_remove.size();i++){
				vector<int> idx;
				vector<float> dist;
				kdtree->radiusSearch (pseed->points[j], dist_remove[i], idx, dist);
				if(idx.size()!=0)
					idx_remove[base+i]=idx[idx.size()-1];
			}
		}		
 
		sort(idx_remove.begin(),idx_remove.end());
		idx_remove.erase(unique(idx_remove.begin(), idx_remove.end()), idx_remove.end());

		// vector<int> idx_psv;
		// for(int i=0;i<cloud->points.size();i++)
		// 		idx_psv.push_back(i);
		 
		// vector<int> result;
    	// set_difference(idx_psv.begin(), idx_psv.end(),idx_remove.begin(), idx_remove.end(),back_inserter(result));

		pcl::PointCloud<PointType>::Ptr cloud_out(new pcl::PointCloud<PointType>);
		for(int i=0;i<idx_remove.size();i++)
			cloud_out->points.push_back(cloud->points[idx_remove[i]]);
		pcl::io::savePLYFileBinary("/home/llg/dataset_paper/stl001_cropped_inhm.ply",*cloud_out);

		

	}
	return 0;
}


