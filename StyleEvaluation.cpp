#include "StyleEvaluation.h"
StyleEvaluation::StyleEvaluation(pcl::PointCloud<PointType>::Ptr cloud,
					 pcl::search::KdTree<PointType>::Ptr kdtree)
{
	cloud_=cloud;
	N_=cloud_->points.size();
	if(kdtree==NULL){
		kdtree_=pcl::search::KdTree<PointType>::Ptr(new pcl::search::KdTree<PointType>);
		kdtree_->setInputCloud(cloud);
	}
	else
		kdtree_=kdtree;
	dmean_=ComputeMeanDistance(cloud_);
	dmin_=ComputeMinDistance(cloud_);
	

	// fast mode parameter
	double z=2.58;	// confidence 0.99
	double MOE=0.001;
	double X=z*z*0.5*0.5/(MOE*MOE);
	n_=round((N_*X)/(X+N_-1));
}

// void StyleEvaluation::OutlierGradeMetric(string output_path, int K)
// {
// 	double ratio=0;
// 	vector<double> og;
// 	og.resize(cloud_->points.size());
// 	// #pragma omp parallel for
//     for(int i=0;i<cloud_->points.size();i++){
// 		vector<V3> arrows;
// 		V3 mevec;
// 		double count=0;

// 		/* Get local point cloud */
// 		vector<int> idx(K);
// 		vector<float> dist(K);
// 		pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);
// 		ctmp->points.resize(K);
// 		kdtree_->nearestKSearch(i, K, idx, dist);        
// 		// #pragma omp parallel for
// 		for(int i=0;i<idx.size();i++){		
// 			ctmp->points[i]=cloud_->points[idx[i]];
// 		}
        
// 		/* ratio */
// 		mevec=GetMEvec(ctmp);		
// 		arrows.resize(ctmp->points.size()-1);
// 		for(int j=1;j<ctmp->points.size();j++){
// 			// cout<<ctmp->points[0]<<endl;
// 			// cout<<ctmp->points[j]<<endl;;
// 			double dx=ctmp->points[j].x-ctmp->points[0].x;
// 			double dy=ctmp->points[j].y-ctmp->points[0].y;
// 			double dz=ctmp->points[j].z-ctmp->points[0].z;
// 			// V3 tmp(dx,dy,dz);
// 			// cout<<tmp<<endl<<endl;;
// 			arrows[j-1].x=dx;
// 			arrows[j-1].y=dy;
// 			arrows[j-1].z=dz;
// 		}

// 		double thresh=70.0/180.0*M_PI;
// 		for(int j=0;j<arrows.size();j++){
// 			// double plen=Dot(arrows[j],mevec)/mevec.GetLength()/arrows[j].GetLength();
// 			double plen=Dot(arrows[j],mevec)/mevec.GetLength()/arrows[j].GetLength();
// 			double arc=0;
// 			if(plen<0)
// 				arc=M_PI-acos(plen);
// 			else
// 			{
// 				arc=acos(plen);
// 			}
// 			if(arrows[j].GetLength()>10*dmean_)
// 				arc=M_PI/2.0;
			
// 			// if(plen<0){
// 			// 	cout<<plen<<endl;
// 			// 	cout<<Dot(arrows[j],mevec)<<endl;
// 			// 	cout<<mevec.GetLength()<<endl;
// 			// 	cout<<arrows[j].GetLength()<<endl<<endl;
// 			// 	return ;
// 			// }
				
// 			if(arc<thresh)
// 				count++;
// 		}		
// 		og[i]=count;


// 		// double plen=0;
// 		// for(int j=0;j<arrows.size();j++){
// 		// 	plen+=Dot(arrows[j],mevec)/mevec.GetLength()/arrows[j].GetLength();
// 		// }
// 		// plen=plen/arrows.size();
// 		// // #pragma omp atomic
// 		// // ratio+=count*1.0/arrows.size();		
// 		// og[i]=plen;
//     }

// 	VectorWrite(output_path,og,"cover");

// 	// ratio=ratio/cloud_->points.size();
// 	// ofstream fout(output_path);
//     // fout<<ratio<<endl;
// 	// fout<<ratio<<endl;
// 	// fout.close();
// }


void StyleEvaluation::OutlierGrade(string output_path, int K)
{
	double ratio=0;
	vector<double> og;
	vector<double> rnum;
	og.resize(cloud_->points.size());
	rnum.resize(cloud_->points.size());
	#pragma omp parallel for
    for(int i=0;i<cloud_->points.size();i++){
		vector<V3> arrows;
		V3 mevec;
		double count=0;

		/* Get local point cloud */
		vector<int> idx;
		vector<float> dist;
		pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);
		
		kdtree_->nearestKSearch(i, K, idx, dist);        
		ctmp->points.resize(idx.size());
		// #pragma omp parallel for
		for(int i=0;i<idx.size();i++){		
			ctmp->points[i]=cloud_->points[idx[i]];
		}
        
		/* ratio */
		mevec=GetMEvec(ctmp);		
		arrows.resize(ctmp->points.size()-1);
		for(int j=1;j<ctmp->points.size();j++){
			arrows[j-1].x=ctmp->points[j].x-ctmp->points[0].x;
			arrows[j-1].y=ctmp->points[j].y-ctmp->points[0].y;
			arrows[j-1].z=ctmp->points[j].z-ctmp->points[0].z;
		}

		double thresh=70.0/180.0*M_PI;
		for(int j=0;j<arrows.size();j++){
			// double plen=Dot(arrows[j],mevec)/mevec.GetLength()/arrows[j].GetLength();
			double plen=Dot(arrows[j],mevec)/mevec.GetLength()/arrows[j].GetLength();
			double arc=0;
			if(plen<0)
				arc=M_PI-acos(plen);
			else
				arc=acos(plen);			

			if(arc<thresh)
				count++;
		}		
		og[i]=count;
		rnum[i]=dist[dist.size()-1];
    }

	// VectorWrite(output_path,og,"cover,column");
	
	ofstream fout;
	fout.open(output_path,ios::out);
    for(int i=0;i<og.size();i++)
		fout<<og[i]<<","<<rnum[i]<<endl;
	fout.close();
}

void StyleEvaluation::Homogeneity(string output_path, int K)
{
	vector<double> density;
	density.resize(cloud_->points.size());
	#pragma omp parallel for
	for(int i=0;i<cloud_->points.size();i++){
		/* Get local point cloud */
		vector<int> idx(K);
		vector<float> dist(K);
		pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);
		ctmp->points.resize(K);
		kdtree_->nearestKSearch(i, K, idx, dist);        
		density[i]=K/dist[K-1];
	}
	VectorNormalize(density);
	VectorWrite(output_path,density,"cover,column");
}

void StyleEvaluation::Sigularity(string output_path, int K1, int K2)
{	
	int N=cloud_->points.size();
	rst_meval.Resize(N);
	vector<int> oidx_meval;
	vector<int> st(N);
	

	#pragma omp parallel for
	for(int i=0;i<cloud_->points.size();i++){
		/* Get local point cloud */
		vector<int> idx(K1);
		vector<float> dist(K1);
		pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);
		ctmp->points.resize(K1);
		kdtree_->nearestKSearch(i, K1, idx, dist);        		
		for(int i=0;i<idx.size();i++){		
			ctmp->points[i]=cloud_->points[idx[i]];
		}
		double meval=GetMEval(ctmp);
		// sg[i]=N*pow(meval/dist[K1-1],2);
		rst_meval.records_[i].id_=i;
		rst_meval.records_[i].item1_=meval;
	}
	double IQR=rst_meval.GetQuantile(0.75)-rst_meval.GetQuantile(0.25);
	double thresh1=1.5*IQR+rst_meval.GetQuantile(0.75);
	for(int i=0;i<cloud_->points.size();i++){
		if(rst_meval.records_[i].item1_>thresh1){
			oidx_meval.push_back(i);
			st[i]=1;
		}
		else
			st[i]=0;
	}
	
	Table<Rrd1> rst_mv(oidx_meval.size());
	for(int i=0;i<oidx_meval.size();i++){
		double count=0;
		vector<int> idx(K2);
		vector<float> dist(K2);
		kdtree_->nearestKSearch(oidx_meval[i], K2, idx, dist);
		for(int j=0;j<K2;j++){
			if(st[idx[j]]==0)
				count++;
		}
		rst_mv.records_[i].id_=oidx_meval[i];
		rst_mv.records_[i].item1_=count/K2;
	}
	IQR=rst_meval.GetQuantile(0.75)-rst_meval.GetQuantile(0.25);
	double thresh2=-1.5*IQR+rst_meval.GetQuantile(0.25);
	thresh2=K1*1.0/K2*8.0;
	vector<int> oidx_mv;
	for(int i=0;i<rst_mv.records_.size();i++){
		if(rst_mv.records_[i].item1_<thresh2){
			int itmp=rst_mv.records_[i].id_;
			cloud_->points[itmp].r=255;
			cloud_->points[itmp].g=0;
			cloud_->points[itmp].b=0;
			oidx_mv.push_back(itmp);
		}
		else{
			int itmp=rst_mv.records_[i].id_;
			cloud_->points[itmp].r=0;
			cloud_->points[itmp].g=255;
			cloud_->points[itmp].b=0;
		}
	}	
	double og=log10(oidx_mv.size());
	cout<<oidx_mv.size()<<endl;
	cout<<"Order of magnitude: "<<og <<endl;
	pcl::io::savePLYFileBinary("Result/1.ply",*cloud_);

	ofstream fout(output_path);
	fout<<og<<endl;
	fout<<og<<endl;
	fout.close();
}

double StyleEvaluation::GetDmean(int K,string str)
{
	if("Common"==str){
		rst_dmean.Resize(cloud_->points.size());
		for(int i=0;i<cloud_->points.size();i++){
			double dmean_tmp=ComputeNearestDistance(kdtree_,i);
			rst_dmean.records_[i].id_=i;
			rst_dmean.records_[i].item1_=dmean_tmp;
		}
		dmean_=rst_dmean.GetMean();
		return dmean_;
	}
	else if("Fast"==str){
		srand((int)time(0));
		rst_dmean.Resize(n_);
		for(int i=0;i<n_;i++){
			int itmp=rand()%N_;
			double dmean_tmp=ComputeNearestDistance(kdtree_,itmp);
			rst_dmean.records_[i].id_=itmp;
			rst_dmean.records_[i].item1_=dmean_tmp;
		}
		dmean_=rst_dmean.GetMean();
		return dmean_;
	}
	else{
		cout<<"Dmean Mode Error!"<<endl;
		return -1;
	}
}

double StyleEvaluation::GetMinorEval(int K,string str)
{
	
	if("Common"==str){
		if(rst_meval.GetSize()==0)
			rst_meval.Resize(cloud_->points.size());
	
		#pragma omp parallel for
		for(int i=0;i<cloud_->points.size();i++){
			double meval_tmp=GetMEvalRadius(cloud_,dmean_*10.0,kdtree_,i);
			rst_meval.records_[i].id_=i;
			rst_meval.records_[i].item1_=meval_tmp;
		}
		rst_meval.Write("Result/meval.csv");
		meval_=rst_meval.GetMean();
		double meval_Q1=rst_meval.GetQuantile(0.25);
		meval_Q3_=rst_meval.GetQuantile(0.75);
		meval_IQR_=meval_Q3_-meval_Q1;
		og_=meval_/(dmean_*dmean_);
		return meval_;
	}
	else if("Fast"==str){
		if(rst_meval.GetSize()!=0){
			rst_meval.Clear();
		}
		rst_meval.Resize(n_);
		srand((int)time(0)); 

		#pragma omp parallel for
		for(int i=0;i<n_;i++){
			int itmp=rand()%N_;
			double meval_tmp=GetMEvalRadius(cloud_,dmean_*10.0,kdtree_,itmp);
			double dmean_tmp=ComputeNearestDistance(kdtree_,itmp);
			rst_meval.records_[i].id_=itmp;
			rst_meval.records_[i].item1_=meval_tmp;
		}
		meval_=rst_meval.GetMean();
		return meval_;
	}
	else
		cout<<"Minor Eigenvalue Mode Error!"<<endl;
}

// double StyleEvaluation::GetDB2(int K, string str)
// {	
// 	if("Common"==str){
// 		if(rst_db2.GetSize()==0)
// 			rst_db2.Resize(cloud_->points.size());
		
// 		#pragma omp parallel for
// 		for(int i=0;i<cloud_->points.size();i++){
// 			double db2_tmp=ComputeDB2(kdtree_,100,i);
// 			rst_db2.records_[i].id_=i;
// 			rst_db2.records_[i].item1_=db2_tmp;
// 		}
// 		meval_=rst_db2.GetMean();


// 		double db2_Q1=rst_db2.GetQuantile(0.25);
// 		db2_Q3_=rst_db2.GetQuantile(0.75);
// 		db2_IQR_=db2_Q3_-db2_Q1;

// 		return meval_;
// 	}
// 	else if("Fast"==str){
// 		if(rst_db2.GetSize()!=0)
// 			rst_db2.Clear();
// 		rst_db2.Resize(n_);

// 		#pragma omp parallel for
// 		for(int i=0;i<n_;i++){
// 			int itmp=rand()%N_;
// 			double db2_tmp=ComputeDB2(kdtree_,100,itmp);
// 			rst_db2.records_[i].id_=i;
// 			rst_db2.records_[i].item1_=db2_tmp;
// 		}
// 		meval_=rst_db2.GetMean();
// 		return meval_;	
// 	}
// 	else{
// 		cout<<"DB2 Mode Error!"<<endl;
// 	}	
// }