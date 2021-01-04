#include "HybridMethods.h"
int is_show_progress=0;
/* Private Functions and Variables */
double cthresh=INT_MAX;
bool customRegionGrowing(const PointType& point_a, const PointType& point_b, float squared_distance)
{
  if (squared_distance < cthresh)
    return true;
  else
    return false;
}


/* Public Function */
void HybridMethods::GetScopeIndices(string str,vector<int>& cIdx)
{
    if("full"==str){
        cIdx.resize(cloud_->points.size());
        #pragma omp parallel for
        for(int i=0;i<cIdx.size();i++)
            cIdx[i]=i;
    }
    else{
        vector<int> emnt; // elements
        Str2Vec(str,",",emnt);
        VecFindPos(status_,emnt,cIdx);
    }
}


void HybridMethods::FM_Slope(int K, double kIQR,string domain)
{
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    if(is_show_progress)
    cout<<"Prox start!"<<endl; 

    vector<int> scope;
    GetScopeIndices(domain,scope);

    rst_slope_.Clear();
    rst_slope_.Resize(scope.size());
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    kdtree_tmp->setInputCloud(cloud_tmp);

    Eigen::MatrixXd xtmp(K,1);
    for(int i=0;i<K;i++) xtmp(i,0)=i;
    // calculate y
    #pragma omp parallel for
    for(int i=0;i<scope.size();i++){
        vector<int> idx(K);
        vector<float> dist(K);
        Eigen::MatrixXd ytmp(K,1);
        int itmp=scope[i];
        kdtree_tmp->nearestKSearch(cloud_->points[itmp], K, idx, dist);

        for(int j=0;j<dist.size();j++) ytmp(j,0)=dist[j];
        Eigen::MatrixXd atmp=(xtmp.transpose()*xtmp).inverse()*xtmp.transpose()*ytmp;
        rst_slope_.records_[i].id_=itmp;
        rst_slope_.records_[i].item1_=atmp(0,0);
    }
    // double t=rst_slope_.ReversePDF(p);
    // for(int i=0;i<scope.size();i++){
    //     if(rst_slope_.records_[i].item1_>t){
    //         status_[scope[i]]+=1;
    //     }
    // }

    double IQR=rst_slope_.GetQuantile(0.75)-rst_slope_.GetQuantile(0.25);
    double thresh=rst_slope_.GetQuantile(0.75)+IQR*kIQR;
    if(is_show_progress==1){       
        for(int i=0;i<scope.size();i++){
            if(rst_slope_.records_[i].item1_>thresh){
                status_[scope[i]]=-accumulator_; 
                rst1->points.push_back(cloud_->points[scope[i]]);
            }
            else{
                status_[scope[i]]=accumulator_; 
                rst0->points.push_back(cloud_->points[scope[i]]);
            }
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_regular.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_irregular.ply",*rst1);
        cout<<"Prox end!"<<endl;
    }
    else{
        #pragma omp parallel for
        for(int i=0;i<scope.size();i++){
            if(rst_slope_.records_[i].item1_>thresh)
                status_[scope[i]]=-accumulator_;            
            else
                status_[scope[i]]=accumulator_;
        }
    }
    
    accumulator_++;
}

// void HybridMethods::FM_Gradient(string output_path, int K, double kIQR,string domain)
// {
//     pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
//     pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
//     if(is_show_progress)
//         cout<<"Gradient start!"<<endl; 

//     vector<int> scope;
//     GetScopeIndices(domain,scope);

//     rst_gradient_.Clear();
//     rst_gradient_.Resize(scope.size());
//     pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
//     pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
//     for(int i=0;i<scope.size();i++)
//         cloud_tmp->points.push_back(cloud_->points[scope[i]]);
//     kdtree_tmp->setInputCloud(cloud_tmp);

//     /* Step 01: Calculate Density */
//     vector<double> dty;
//     dty.resize(scope.size());
//     int k2=40;
//     for(int i=0;i<scope.size();i++){
//         vector<int> idx(k2);
//         vector<float> dist(k2);        
//         kdtree_tmp->nearestKSearch(i,k2, idx, dist);
//         // dty[i]=dist[k2-1]/k2;
//         dty[i]=k2/dist[k2-1];
//     }
//     VectorNormalize(dty);

//     /* Step 02: Insensitive to density */
//     Eigen::MatrixXd xtmp(K,1);
//     for(int i=0;i<K;i++) xtmp(i,0)=i;
//     // calculate y
//     // #pragma omp parallel for
//     for(int i=0;i<scope.size();i++){
//         vector<int> idx(K);
//         vector<float> dist(K);
//         Eigen::MatrixXd ytmp(K,1);
//         kdtree_tmp->nearestKSearch(i, K, idx, dist);

//         /* Get BS */ 
//         vector<double> db;
//         DaubechiesWavelet(dist,db);
//         double d4=VectorMaximum(db);
//         // cout<<d4<<endl;
//         // return ;
    
//         rst_gradient_.records_[i].id_=i;
//         rst_gradient_.records_[i].item1_=d4*dty[i];
//         // rst_gradient_.records_[i].item1_=dty[i];
//     }

//     ofstream fout;
// 	fout.open(output_path,ios::out);
//     // rst_gradient_.Standardize_Zscore();
//     // rst_gradient_.Normalize_Tanh();
//     rst_gradient_.GetCorrespondingColor();
//     for(int i=0;i<cloud_->points.size();i++){        
//         fout<<rst_gradient_.records_[i].item1_<<endl;
//         V3 ctmp=rst_gradient_.color_[i];
//         cloud_->points[i].r=(int)ctmp.r;
//         cloud_->points[i].g=(int)ctmp.g;
//         cloud_->points[i].b=(int)ctmp.b;
//     }
// 	pcl::io::savePLYFileBinary("Result/1.ply",*cloud_);	
// 	fout.close();

//     // double IQR=rst_gradient_.GetQuantile(0.75)-rst_gradient_.GetQuantile(0.25);
//     // double thresh=rst_gradient_.GetQuantile(0.75)+IQR*1.5;
//     // for(int i=0;i<scope.size();i++){
//     //     if(rst_gradient_.records_[i].item1_>thresh){
//     //         status_[scope[i]]=-accumulator_; 
//     //         rst1->points.push_back(cloud_->points[scope[i]]);
//     //     }
//     //     else{
//     //         status_[scope[i]]=accumulator_; 
//     //         rst0->points.push_back(cloud_->points[scope[i]]);
//     //     }
//     // }
//     // pcl::io::savePLYFileBinary("Result/0.ply",*rst0);
//     // pcl::io::savePLYFileBinary("Result/1.ply",*rst1);
//     // cout<<"Gradient end!"<<endl;

//     // double IQR=rst_gradient_.GetQuantile(0.75)-rst_gradient_.GetQuantile(0.25);
//     // double thresh=rst_gradient_.GetQuantile(0.75)+IQR*kIQR;
//     // if(is_show_progress==1){       
//     //     for(int i=0;i<scope.size();i++){
//     //         if(rst_gradient_.records_[i].item1_>thresh){
//     //             status_[scope[i]]=-accumulator_; 
//     //             rst1->points.push_back(cloud_->points[scope[i]]);
//     //         }
//     //         else{
//     //             status_[scope[i]]=accumulator_; 
//     //             rst0->points.push_back(cloud_->points[scope[i]]);
//     //         }
//     //     }
//     //     pcl::io::savePLYFileBinary(store_path_+"rst("+to_string(accumulator_)+").ply",*rst0);
//     //     pcl::io::savePLYFileBinary(store_path_+"rst("+to_string(accumulator_)+"-).ply",*rst1);
//     //     cout<<"Gradient end!"<<endl;
//     // }
//     // else{
//     //     #pragma omp parallel for
//     //     for(int i=0;i<scope.size();i++){
//     //         if(rst_gradient_.records_[i].item1_>thresh)
//     //             status_[scope[i]]=-accumulator_;            
//     //         else
//     //             status_[scope[i]]=accumulator_;
//     //     }
//     // }
    
//     accumulator_++;
// }

void HybridMethods::FM_DensityInsensitive(string output_path, int K)
{
    double ratio=0;
	vector<double> og;
	vector<double> rnum;
	og.resize(cloud_->points.size());
	rnum.resize(cloud_->points.size());
    Table<Rrd1> rst_DI;
    rst_DI.Resize(cloud_->points.size());

	// #pragma omp parallel for
    for(int i=0;i<cloud_->points.size();i++){
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
		og[i]=count;
		rnum[i]=dist[dist.size()-1];
    }

	// VectorWrite(output_path,og,"cover");
    double rnum_sum=exp(VectorSum(rnum));
	for(int i=0;i<cloud_->points.size();i++){
        rst_DI.records_[i].id_=i;
        rst_DI.records_[i].item1_=og[i]*exp(rnum[i])/rnum_sum;
    }

	ofstream fout;
	fout.open(output_path,ios::out);
    rst_DI.GetCorrespondingColor();
    for(int i=0;i<og.size();i++){        
        fout<<rst_DI.records_[i].item1_<<endl;
        V3 ctmp=rst_DI.color_[i];
        cloud_->points[i].r=ctmp.r;
        cloud_->points[i].g=ctmp.g;
        cloud_->points[i].b=ctmp.b;
    }
	pcl::io::savePLYFileBinary("Result/1.ply",*cloud_);	
	fout.close();
}

void HybridMethods::FM_MEval(int K, double alpha,string domain)
{
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    if(is_show_progress==1){
        cout<<"MEval start!"<<endl;
    }

    vector<int> scope;
    GetScopeIndices(domain,scope);
    rst_meval_.Clear();
    rst_meval_.Resize(scope.size());
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    kdtree_tmp->setInputCloud(cloud_tmp);

    #pragma omp parallel for
    for(int i=0;i<scope.size();i++){
        int itmp=scope[i];
        // double eval_tmp=GetMEvalRadius(cloud_tmp,K*pOle_->dmean_,kdtree_tmp,i);
        double eval_tmp=GetMEvalKNeighours(cloud_tmp,K,kdtree_tmp,i);
        rst_meval_.records_[i].id_=itmp;
        rst_meval_.records_[i].item1_=eval_tmp;
    }
    // rst_meval_.Standardize_Zscore();
    // rst_meval_.Normalize_Tanh();
    double IQR=rst_meval_.GetQuantile(0.75)-rst_meval_.GetQuantile(0.25);
    double thresh=rst_meval_.GetQuantile(0.75)+IQR*alpha;        
    if(is_show_progress==1){
        for(int i=0;i<scope.size();i++){
            if(rst_meval_.records_[i].item1_>thresh){
                status_[scope[i]]=-accumulator_;
                rst1->points.push_back(cloud_->points[scope[i]]);
            }
            else{
                status_[scope[i]]=accumulator_;
                rst0->points.push_back(cloud_->points[scope[i]]);
            }
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_regular.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_irregular.ply",*rst1);
        cout<<"MEval end!"<<endl;
    }
    else{
        #pragma omp parallel for
        for(int i=0;i<scope.size();i++){
        if(rst_meval_.records_[i].item1_>thresh)
            status_[scope[i]]=-accumulator_;        
        else
            status_[scope[i]]=accumulator_;            
        }
    }
    accumulator_++;
}

void HybridMethods::FM_NID(int K, double p, string domain)
{    
    if(is_show_progress==1){
        cout<<"NID start!"<<endl;
    }

    vector<int> scope;
    GetScopeIndices(domain,scope);
    rst_nid_.Clear();
    rst_nid_.Resize(scope.size());
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    kdtree_tmp->setInputCloud(cloud_tmp);

    Eigen::MatrixXd xtmp(K,1);
    for(int i=0;i<K;i++) xtmp(i,0)=i;
    // calculate y
    #pragma omp parallel for
    for(int i=0;i<scope.size();i++){
        vector<int> idx(K);
        vector<float> dist(K);
        Eigen::MatrixXd ytmp(K,1);
        int itmp=scope[i];
        kdtree_tmp->nearestKSearch(cloud_->points[itmp], K, idx, dist);
        for(int j=0;j<dist.size();j++) ytmp(j,0)=dist[j];
        Eigen::MatrixXd atmp=(xtmp.transpose()*xtmp).inverse()*xtmp.transpose()*ytmp;
        rst_nid_.records_[i].id_=itmp;
        rst_nid_.records_[i].item1_=atmp(0,0);
    }

    rst_nid_.LocalFilter("average",cloud_tmp,80);
    rst_nid_.Standardize_Zscore();

    if(is_show_progress==1){
        pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
        pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
        for(int i=0;i<rst_nid_.GetSize();i++){
            double etmp=GaussErrorFunction(rst_nid_.records_[i].item1_);
            if(etmp>p){
                status_[scope[i]]=-accumulator_;
                rst1->points.push_back(cloud_->points[scope[i]]);
            }
            else{
                status_[scope[i]]=accumulator_;
                rst0->points.push_back(cloud_->points[scope[i]]);
            }            
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+").ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"-).ply",*rst1);
        cout<<"NID end!"<<endl;
    }
    else{
        #pragma omp parallel for
        for(int i=0;i<rst_nid_.GetSize();i++){
            double etmp=GaussErrorFunction(rst_nid_.records_[i].item1_);
            if(etmp>p)
                status_[scope[i]]=-accumulator_;                
            else
                status_[scope[i]]=accumulator_;                         
        } 
    }
    
    accumulator_++;     
}

void HybridMethods::FM_LoOP(string domain, int K, double thresh)
{    
    if(is_show_progress==1){
        cout<<"LoOP start!"<<endl;
    }    
   vector<int> scope;
    if("0"==domain){
        scope.resize(cloud_->points.size());
        #pragma omp parallel for
        for(int i=0;i<scope.size();i++)
            scope[i]=i;
    }
    else
        GetScopeIndices(domain,scope);
    
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    kdtree_tmp->setInputCloud(cloud_tmp);

	vector<double> sigma;
	vector<double> plof;
    double nplof=0;
    // Resize Scores	
	sigma.resize(cloud_tmp->points.size());
	plof.resize(cloud_tmp->points.size());
	rst_LoOP_.Resize(cloud_tmp->points.size());

	// Step 01: Calculate sigma
	#pragma omp parallel for
	for(int i=0;i<cloud_tmp->points.size();i++){
		// find k-nearest neighours
		vector<int> idx(K+1);
		vector<float> dist(K+1);
		kdtree_tmp->nearestKSearch (cloud_tmp->points[i], K+1, idx, dist);
		// cout<<cloud->points[i]<<endl;
		double sum=0;
		for(int j=1;j<K+1;j++){
			sum+=dist[j];
		}
		sum=sum/K;
		sigma[i]=sqrt(sum);
	}
	
	// Step 02: calculate mean
	double mean=0;
	#pragma omp parallel for
	for (int i = 0; i < cloud_tmp->points.size(); i++){        
        vector<int> idx(K+1);
		vector<float> dist(K+1);
		kdtree_tmp->nearestKSearch (cloud_tmp->points[i], K+1, idx, dist);
        double sum = 0;
        for (int j = 1; j < K+1; j++)
          sum += sigma[idx[j]];
        sum /= K;
        plof[i] = sigma[i] / sum  - 1.0f;				
        mean += plof[i] * plof[i];
    }
	nplof=sqrt(mean/cloud_tmp->points.size());	

	// Step 03: caculate score
	#pragma omp parallel for
	for(int i=0;i<cloud_tmp->points.size();i++){
		double value = plof[i] / (nplof * sqrt(2.0f));
		// rst_.records_[i].item1_=value;

        double dem = 1.0 + 0.278393 * value;
        dem += 0.230389 * value * value;
        dem += 0.000972 * value * value * value;
        dem += 0.078108 * value * value * value * value;
        double op = std::max(0.0, 1.0 - 1.0 / dem);
        rst_LoOP_.records_[i].id_=i;
        rst_LoOP_.records_[i].item1_ = op;
	}

	// #pragma omp parallel for
    if(is_show_progress==1){
        pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
        pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
        for(int i=0;i<scope.size();i++){
            if(rst_LoOP_.records_[i].item1_>thresh){
                status_[scope[i]]=-accumulator_;
                rst1->points.push_back(cloud_tmp->points[scope[i]]);
            }
            else{
                status_[scope[i]]=accumulator_;
                rst0->points.push_back(cloud_tmp->points[scope[i]]);
            }        
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_regular.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_irregular.ply",*rst1);
        cout<<"LoOP end!"<<endl;
    }
    else{
        #pragma omp parallel for
       	for(int i=0;i<cloud_tmp->points.size();i++){
            int itmp=GetIndex(kdtree_,cloud_tmp->points[i]);		
            if( rst_LoOP_.records_[i].item1_>thresh)
                status_[itmp]=-accumulator_;		
            else        
                status_[itmp]=accumulator_;                
	    }
    }
    accumulator_++;	
}
void HybridMethods::FM_D4(int K, double P,string domain)
{
    // Initialization
    rst_d4_.Clear();
    if(is_show_progress==1){
        cout<<"D4 start!"<<endl;
    }
    vector<int> scope;
    GetScopeIndices(domain,scope);
    rst_d4_.Resize(scope.size());
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    kdtree_tmp->setInputCloud(cloud_tmp);
    rst_d4_.Resize(cloud_tmp->points.size());

    /* calculate d4 */
    #pragma omp parallel for
    for(int i=0;i<cloud_tmp->points.size();i++){
        vector<int> idx(K);
        vector<float> dist(K);        
        vector<double> ca,cd;        
        kdtree_tmp->nearestKSearch(cloud_tmp->points[i], K, idx, dist);                
        DaubechiesWavelet(dist,ca,cd);        
        rst_d4_.records_[i].id_=i;
        rst_d4_.records_[i].item1_=VectorMaximum(cd);
    }

    /* Manage Status List */
    double IQR=rst_d4_.GetQuantile(0.75)-rst_d4_.GetQuantile(0.25);
    double thresh=rst_d4_.GetQuantile(0.75)+IQR*P;
    if(is_show_progress==1){
        pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
        pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
        for(int i=0;i<scope.size();i++){
            if(rst_d4_.records_[i].item1_>thresh){
                status_[scope[i]]=-accumulator_;
                rst1->points.push_back(cloud_->points[scope[i]]);
            }
            else{
                status_[scope[i]]=accumulator_;
                rst0->points.push_back(cloud_->points[scope[i]]);
            }        
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_regular.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_outlier.ply",*rst1);
        cout<<"D4 end!"<<endl;
    }
    else{
        #pragma omp parallel for
        for(int i=0;i<scope.size();i++){
            if(rst_d4_.records_[i].item1_>thresh)
                status_[scope[i]]=-accumulator_;
            else
                status_[scope[i]]=accumulator_;                     
        }
    }
    
    accumulator_++;
}

void HybridMethods::FM_LNFS(int K, double P,string domain)
{
    // Initialization
    rst_lnfs_.Clear();
    if(is_show_progress==1){
        cout<<"LNFS start!"<<endl;
    }
    vector<int> scope;
    GetScopeIndices(domain,scope);
    rst_lnfs_.Resize(scope.size());
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    kdtree_tmp->setInputCloud(cloud_tmp);

    // find all k-nearest neighbours
    vector<vector<int>> buf_nei_idx(cloud_tmp->points.size(),vector<int>(K+1));
    vector<vector<float>> buf_nei_dist(cloud_tmp->points.size(),vector<float>(K+1));
    #pragma omp parallel for
    for(int i=0;i<cloud_tmp->points.size();i++){
        vector<int> idx;
        vector<float> dist;
        kdtree_tmp->nearestKSearch(cloud_tmp->points[i], K+1, idx, dist);
        for(int j=0;j<K+1;j++){
            buf_nei_idx[i][j]=idx[j];
            buf_nei_dist[i][j]=dist[j];
        }            
    }
    
    // Init ytmp
    #pragma omp parallel for
    for(int i=0;i<scope.size();i++){   
        vector<double> gap(K);
        vector<int> active_idx(K);  
        for(int j=0;j<K;j++)
            active_idx[j]=buf_nei_idx[i][j];

        for(int j=0;j<K;j++){
            // find nearest neighbour which is not belong to local k-nearest neighbours            
            int idx_tmp=VectorFindFirstNot(active_idx,buf_nei_idx[active_idx[j]]);
            double dtmp=pcl::geometry::squaredDistance(cloud_tmp->points[i],cloud_tmp->points[idx_tmp]);
            gap[j]=dtmp;
        }
        rst_lnfs_.records_[i].id_=i;
        rst_lnfs_.records_[i].item1_=VectorMaximum(gap);
    }    
    
    /* Manage Status List */
    double IQR=rst_lnfs_.GetQuantile(0.75)-rst_lnfs_.GetQuantile(0.25);
    double thresh=rst_lnfs_.GetQuantile(0.75)+IQR*P;
    if(is_show_progress==1){
        pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
        pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
        for(int i=0;i<scope.size();i++){
            if(rst_lnfs_.records_[i].item1_>thresh){
                status_[scope[i]]=-accumulator_;
                rst1->points.push_back(cloud_->points[scope[i]]);
            }
            else{
                status_[scope[i]]=accumulator_;
                rst0->points.push_back(cloud_->points[scope[i]]);
            }        
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_regular.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_outlier.ply",*rst1);
        cout<<"LNFS end!"<<endl;
    }
    else{
        #pragma omp parallel for
        for(int i=0;i<scope.size();i++){
            if(rst_lnfs_.records_[i].item1_>thresh)
                status_[scope[i]]=-accumulator_;
            else
                status_[scope[i]]=accumulator_;                     
        }
    }
    
    accumulator_++;
}

void HybridMethods::FM_Density(int K, double alpha,string domain)
{
    if(is_show_progress==1){
        cout<<"Desisty 1 start!"<<endl;
    }

    /* Get domain */
    vector<int> scope;
    if("0"==domain){
        scope.resize(cloud_->points.size());
        #pragma omp parallel for
        for(int i=0;i<scope.size();i++)
            scope[i]=i;
    }
    else
        GetScopeIndices(domain,scope);

    /* kdtree for point cloud in domain */
    rst_density_.Clear();
    rst_density_.Resize(scope.size());
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    kdtree_tmp->setInputCloud(cloud_tmp);

    /* calculate density for each point in domain */
    for(int i=0;i<scope.size();i++){
        vector<int> idx(K);
        vector<float> dist(K);        
        kdtree_tmp->nearestKSearch(i,K, idx, dist);
        rst_density_.records_[i].id_=i;
        rst_density_.records_[i].item1_=K/dist[K-1];
    }

    // double IQR=rst_density_.GetQuantile(0.75)-rst_density_.GetQuantile(0.25);
    // double thresh=rst_density_.GetQuantile(0.25)-IQR*3.0;
    // double max_tmp=rst_density_.GetMaximum();
    // for(int i=0;i<scope.size();i++){
    //     if(rst_density_.records_[i].item1_<thresh){         
    //         rst_density_.records_[i].item1_=max_tmp;
    //     }
    // }
    // if(is_show_progress==1){
    //     pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    //     pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    //     for(int i=0;i<scope.size();i++){
    //         if(rst_density_.records_[i].item1_>thresh){
    //             status_[scope[i]]=-accumulator_;
    //             rst1->points.push_back(cloud_->points[scope[i]]);
    //         }
    //         else{
    //             status_[scope[i]]=accumulator_;
    //             rst0->points.push_back(cloud_->points[scope[i]]);
    //         }
    //     }
    //     pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_regular.ply",*rst0);
    //     pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_irregular.ply",*rst1);
    //     cout<<"Desisty 1 end!"<<endl;
    // }
    // else{
    //     #pragma omp parallel for
    //     for(int i=0;i<N_;i++){
    //         if(rst_density_.records_[i].item1_>thresh)
    //             status_[i]=-accumulator_;
    //         else
    //             status_[i]=accumulator_;
    //     }
    // }
    
    accumulator_++;
}
void HybridMethods::FM_PolyFitting(string domain)
{

}

/*
    This module is based pcl library, however, its performance is not good. 
    So the quantity of input point cloud should be as smaller as possible.
*/
void HybridMethods::FM_RegionGrowth(int level,double thresh_kIQR, string domain)
{
    if(is_show_progress==1){
        cout<<"RegionGrowth start!"<<endl;
    }

    vector<int> scope;
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);    
    pcl::IndicesClustersPtr clusters (new pcl::IndicesClusters);
    pcl::IndicesClustersPtr small_clusters (new pcl::IndicesClusters);
    pcl::IndicesClustersPtr large_clusters (new pcl::IndicesClusters);
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    GetScopeIndices(domain,scope);

    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    double cellsize=GetCellSize(cloud_tmp,level);

    // Set up a Conditional Euclidean Clustering class
    pcl::ConditionalEuclideanClustering<PointType> cec (true);
    cec.setInputCloud (cloud_tmp);
    cec.setConditionFunction (&customRegionGrowing);
    cec.setClusterTolerance (cellsize);
    cec.segment (*clusters);
    vector<int> cluster_size_set;
    for(int i=0;i<clusters->size();i++){
       cluster_size_set.push_back((*clusters)[i].indices.size());
    }


    double IQR=Quantile(cluster_size_set,0.75)-Quantile(cluster_size_set,0.25);
    double thresh=Quantile(cluster_size_set,0.75)+IQR*thresh_kIQR;
    if(is_show_progress==1){
        cout<<"quantity of cluster:"<<cluster_size_set.size()<<endl;
        cout<<"cluster thresh:"<<thresh<<endl;

        for(int i=0;i<clusters->size();i++){
            int current_cluster_size=(*clusters)[i].indices.size();
            if(current_cluster_size<=thresh){            
                for(int j=0;j<(*clusters)[i].indices.size();j++){
                    int itmp=GetIndex(kdtree_,cloud_tmp->points[(*clusters)[i].indices[j]]);
                    status_[itmp]=-accumulator_;
                    rst1->points.push_back(cloud_->points[itmp]);
                }
            }
            else{
                for(int j=0;j<(*clusters)[i].indices.size();j++){
                    int itmp=GetIndex(kdtree_,cloud_tmp->points[(*clusters)[i].indices[j]]);
                    status_[itmp]=accumulator_;
                    rst0->points.push_back(cloud_->points[itmp]);
                }
            }     
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_regular.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"_irregular.ply",*rst1);
        cout<<"RegionGrowth end!"<<endl;
    }
    else{
        #pragma omp parallel for
        for(int i=0;i<clusters->size();i++){
            int current_cluster_size=(*clusters)[i].indices.size();
            if(current_cluster_size<=thresh){            
                for(int j=0;j<(*clusters)[i].indices.size();j++){
                    int itmp=GetIndex(kdtree_,cloud_tmp->points[(*clusters)[i].indices[j]]);
                    status_[itmp]=-accumulator_;                    
                }
            }
            else{
                for(int j=0;j<(*clusters)[i].indices.size();j++){
                    int itmp=GetIndex(kdtree_,cloud_tmp->points[(*clusters)[i].indices[j]]);
                    status_[itmp]=accumulator_;
                }
            }     
        }
    }
    
    accumulator_++;
}

void HybridMethods::FM_MajorityVote(int K,double thresh, string domain)
{
    if(is_show_progress==1){
        cout<<"MajorityVote start!"<<endl;
    }

    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    vector<int> scope;
    GetScopeIndices(domain,scope);
    Table<Rrd1> rst_mv_;
    rst_mv_.Resize(scope.size());

    #pragma omp parallel for
    for(int i=0;i<scope.size();i++){
        int itmp=scope[i];
        vector<int> idx;
        vector<float> dist;
        kdtree_->nearestKSearch(itmp,K,idx,dist);
        // kdtree_->radiusSearch(itmp,pOle_->dmean_*200,idx,dist);

        double conf=0.0;
        double T1=5*pOle_->dmean_;
        for(int j=0;j<idx.size();j++){
            if(status_[idx[j]]>0 && dist[j]<T1)
                conf+=1.0;
        }
        conf=conf/K;
        rst_mv_.records_[i].id_=scope[i];
        rst_mv_.records_[i].item1_=conf;
    }

    // double IQR=rst_mv_.GetQuantile(0.75)-rst_mv_.GetQuantile(0.25);
    // double thresh=rst_mv_.GetQuantile(0.75)+IQR*1.5;
    // double thresh=(rst_mv_.GetMaximum()-rst_mv_.GetMinimum())*0.25+rst_mv_.GetMinimum();
    // double thresh=0.01;   
    if(is_show_progress==1){
        for(int i=0;i<scope.size();i++){
            if(rst_mv_.records_[i].item1_<thresh){
                status_[scope[i]]=-accumulator_;
                rst1->points.push_back(cloud_->points[scope[i]]);
            }
            else{
                status_[scope[i]]=accumulator_;
                rst0->points.push_back(cloud_->points[scope[i]]);
            }      
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+").ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(accumulator_)+"-).ply",*rst1);
        cout<<"MajorityVote end!"<<endl;
    }
    else{
        #pragma omp parallel for
        for(int i=0;i<scope.size();i++){
            if(rst_mv_.records_[i].item1_<thresh)
                status_[scope[i]]=-accumulator_;
            else
                status_[scope[i]]=accumulator_;                   
        }
    }
    
    accumulator_++;
}

void HybridMethods::DemonstrateResult(string path)
{
    endTime_=omp_get_wtime();
    pcl::PointCloud<PointType>::Ptr rst_regular_cloud(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst_irregular_cloud(new pcl::PointCloud<PointType>);
    double max_tmp=*max_element(status_.begin(),status_.end());
    VecUnique(status_);

    if(max_tmp>0){
            for(int i=0;i<status_.size();i++){
            if(status_[i]>0){                    
                rst_regular_cloud->points.push_back(cloud_->points[i]);
            }
            else{
                rst_irregular_cloud->points.push_back(cloud_->points[i]);
            }
        }
        

        // pcl::io::savePLYFileBinary("Result/rst_regular.ply",*rst_regular_cloud);        
        // pcl::io::savePLYFileBinary("Result/rst_irregular.ply",*rst_irregular_cloud);
        // cout<<path<<endl;
        pcl::io::savePLYFileBinary(path,*rst_regular_cloud);
    }
    else
        cout<<"Status Error!"<<endl;
    
    cout<<"Elapse: "<<(endTime_-startTime_)<<endl;
}
void HybridMethods::DemonstrateResult_Color(string path, Table<Rrd1>& tb, string mode)
{
    endTime_=omp_get_wtime();
    if(mode=="linear"){
        tb.GetCorrespondingColor();
        cout<<tb.GetSize()<<endl;
        cout<<cloud_->points.size()<<endl;

        for(int i=0;i<cloud_->points.size();i++){
            V3 ctmp=tb.color_[i];
            cloud_->points[i].r=ctmp.r;
            cloud_->points[i].g=ctmp.g;
            cloud_->points[i].b=ctmp.b;        
       }
       pcl::io::savePLYFileBinary(path,*cloud_);
    }
    else if("nonlinear"==mode){
        tb.Standardize_Zscore();
        tb.Normalize_Tanh();
        tb.GetCorrespondingColor();        
        for(int i=0;i<tb.GetSize();i++){
            V3 ctmp=tb.color_[i];
            cloud_->points[i].r=ctmp.r;
            cloud_->points[i].g=ctmp.g;
            cloud_->points[i].b=ctmp.b;        
       }
       pcl::io::savePLYFileBinary(path,*cloud_);
    }
   
    cout<<"Elapse: "<<(endTime_-startTime_)<<endl;
}