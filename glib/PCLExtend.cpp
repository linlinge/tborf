#include "PCLExtend.h"
#include "SignalProcessing.h"
#include "VectorExtend.h"
#define MAX2(A,B) ((A)>(B) ? (A):(B))
#define MAX3(A,B,C) MAX2(MAX2(A,B),C)
#define MIN2(A,B) ((A)<(B) ? (A):(B))
#define MIN3(A,B,C) MIN2(MIN2(A,B),C)
double ComputeMeanDistance(const pcl::PointCloud<PointType>::ConstPtr cloud)
{
	double res = 0.0;
	int n_points = 0;
	int nres;
	std::vector<int> idx(2);
	std::vector<float> dist(2);
	pcl::search::KdTree<PointType> tree;
	tree.setInputCloud(cloud);

	#pragma omp parallel for
	for (size_t i = 0; i < cloud->size(); ++i)
	{
		nres = tree.nearestKSearch(i, 2, idx, dist);
		if (nres == 2)
		{
			res += sqrt(dist[1]);
			++n_points;
		}
	}
	if (n_points != 0)
	{
		res /= n_points;
	}
	return res;
}




double ComputeMaxDistance(const pcl::PointCloud<PointType>::ConstPtr cloud)
{
	double rst = 0.0;
	int nres;
	std::vector<int> idx(2);
	std::vector<float> dist(2);
	pcl::search::KdTree<PointType> tree;
	tree.setInputCloud(cloud);

	for (size_t i = 0; i < cloud->size(); ++i)
	{
		if (!std::isfinite((*cloud)[i].x))
		{
			continue;
		}
		nres = tree.nearestKSearch(i, 2, idx, dist);
		if (nres == 2)
		{			
			rst = rst>dist[1] ? rst:dist[1];
		}
	}

	return rst;
}

double ComputeMinDistance(const pcl::PointCloud<PointType>::ConstPtr cloud)
{
	double rst = INT_MAX;
	int nres;
	std::vector<int> idx(2);
	std::vector<float> dist(2);
	pcl::search::KdTree<PointType> tree;
	tree.setInputCloud(cloud);

	for (size_t i = 0; i < cloud->size(); ++i)
	{
		if (!std::isfinite((*cloud)[i].x))
		{
			continue;
		}
		nres = tree.nearestKSearch(i, 2, idx, dist);
		if (nres == 2)
		{			
			rst = rst<dist[1] ? rst:dist[1];
		}
	}

	return rst;
}

// compute centroid of point cloud
PointType ComputeCentroid(const pcl::PointCloud<PointType>::ConstPtr cloud)
{
   PointType centroid;
   centroid.x=0;
   centroid.y=0;
   centroid.z=0;
   for(int i=0;i<cloud->points.size();i++){
	   centroid.x+=cloud->points[i].x;
	   centroid.y+=cloud->points[i].y;
	   centroid.z+=cloud->points[i].z;
   }
   centroid.x=centroid.x*1.0/cloud->points.size();
   centroid.y=centroid.y*1.0/cloud->points.size();
   centroid.z=centroid.z*1.0/cloud->points.size();
   return centroid;
}

void TransformPointCloud(pcl::PointCloud<PointType>::Ptr cloud, pcl::PointCloud<PointType>::Ptr cloud_tf,Eigen::Affine3f tf)
{
	cloud_tf->points.resize(cloud->points.size());	
	#pragma omp parallel for
	for(int i=0;i<cloud->points.size();i++)
	{
		Eigen::Vector3f v1(cloud->points[i].x,cloud->points[i].y,cloud->points[i].z);		
		Eigen::Vector3f v2=tf*v1;
		cloud_tf->points[i].x=v2(0,0);
		cloud_tf->points[i].y=v2(1,0); 
		cloud_tf->points[i].z=v2(2,0);
	}
}

void GetEvalAndEvec(pcl::PointCloud<PointType>::Ptr cloud,int k,
					pcl::search::KdTree<PointType>::Ptr kdtree,int index,
                    vector<double>& eval,vector<V3>& evec)
{
	// Init
	eval.resize(3);
    evec.resize(3);

	// local cloud
	vector<int> idx(k);
	vector<float> dist(k);		
	kdtree->nearestKSearch(index, k, idx, dist);
	pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);     
	for(int j=0;j<idx.size();j++)
		ctmp->points.push_back(cloud->points[idx[j]]);

	// Calculate Evec and Eval
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covariance;
	pcl::compute3DCentroid(*ctmp, centroid);
	pcl::computeCovarianceMatrixNormalized(*ctmp, centroid, covariance);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
	Eigen::Matrix3f eig_vec = eigen_solver.eigenvectors();
	Eigen::Vector3f eig_val = eigen_solver.eigenvalues();
	
	// eigenvalue
	eval[0]=eig_val(0);
	eval[1]=eig_val(1);
	eval[2]=eig_val(2);

	// eigenvector
	eig_vec.col(2) = eig_vec.col(0).cross(eig_vec.col(1));
	eig_vec.col(0) = eig_vec.col(1).cross(eig_vec.col(2));
	eig_vec.col(1) = eig_vec.col(2).cross(eig_vec.col(0));
	evec[0].x=eig_vec(0,0);evec[0].y=eig_vec(1,0);evec[0].z=eig_vec(2,0);
	evec[1].x=eig_vec(0,1);evec[1].y=eig_vec(1,1);evec[1].z=eig_vec(2,1);
	evec[2].x=eig_vec(0,2);evec[2].y=eig_vec(1,2);evec[2].z=eig_vec(2,2);
}


V3 GetMEvec(pcl::PointCloud<PointType>::Ptr local_cloud)
{
	V3 mevec;
	// Calculate Evec and Eval
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covariance;
	pcl::compute3DCentroid(*local_cloud, centroid);
	pcl::computeCovarianceMatrixNormalized(*local_cloud, centroid, covariance);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
	Eigen::Matrix3f eig_vec = eigen_solver.eigenvectors();
	Eigen::Vector3f eig_val = eigen_solver.eigenvalues();
	// eigenvector
	eig_vec.col(2) = eig_vec.col(0).cross(eig_vec.col(1));
	eig_vec.col(0) = eig_vec.col(1).cross(eig_vec.col(2));
	eig_vec.col(1) = eig_vec.col(2).cross(eig_vec.col(0));
	mevec.x=eig_vec(0,0); mevec.y=eig_vec(1,0); mevec.z=eig_vec(2,0);
	return mevec;
}

double GetMEval(pcl::PointCloud<PointType>::Ptr local_cloud)
{
	V3 mevec;
	// Calculate Evec and Eval
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covariance;
	pcl::compute3DCentroid(*local_cloud, centroid);
	pcl::computeCovarianceMatrixNormalized(*local_cloud, centroid, covariance);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);	
	Eigen::Vector3f eig_val = eigen_solver.eigenvalues();
	return eig_val(0);
}

void GetLocalCloud(pcl::PointCloud<PointType>::Ptr whole_cloud, 
				   pcl::search::KdTree<PointType>::Ptr kdtree,
				   int index, int K, 
				   pcl::PointCloud<PointType>::Ptr local_cloud)
{
	vector<int> idx(K);
	vector<float> dist(K);
	local_cloud=pcl::PointCloud<PointType>::Ptr(new pcl::PointCloud<PointType>);
	pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);
	ctmp->resize(K);
	kdtree->nearestKSearch (index, K, idx, dist);        
	// #pragma omp parallel for
	for(int i=0;i<idx.size();i++){		
		ctmp->points[i]=whole_cloud->points[idx[i]];
	}        
}



void GetEval(pcl::PointCloud<PointType>::Ptr cloud,int k,
			 pcl::search::KdTree<PointType>::Ptr kdtree,int index,vector<int>& eval)
{
	// Init
	eval.resize(3);

	// local cloud
	vector<int> idx(k);
	vector<float> dist(k);		
	kdtree->nearestKSearch(index, k, idx, dist);
	pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);     
	for(int j=0;j<idx.size();j++)
		ctmp->points.push_back(cloud->points[idx[j]]);
	
	// Eval
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covariance;
	pcl::compute3DCentroid(*ctmp, centroid);
	pcl::computeCovarianceMatrixNormalized(*ctmp, centroid, covariance);
	double a=covariance(0,0); double b=covariance(0,1); double c=covariance(0,2);
	double d=covariance(1,1); double e=covariance(1,2); double f=covariance(2,2);
	double e2=e*e; double b2=b*b;double c2=c*c;
	double B=-(d+f+a); double C=d*f+a*d+a*f-e2-b2-c2;
	double D=-a*d*f+a*e2+b2*f-2*b*e*c+c2*d;
	double p=(3*C-B*B)/3.0; double q=(27.0*D-9*B*C+2*B*B*B)/27.0;
	double r=sqrt(-p*p*p/27.0);
	double theta=1.0/3*acos(-q/(2.0*r));
	double tmp1=2*pow(r,1/3.0); double tmp2=B/3.0;
	eval[0]=tmp1*cos(theta)-tmp2;
	eval[1]=tmp1*cos(theta+2.0/3.0*M_PI)-tmp2;
	eval[2]=tmp1*cos(theta+4.0/3.0*M_PI)-tmp2;
	sort(eval.begin(),eval.end());
}

/*
	Not Successed
*/
double GetMEvalKDistance(pcl::PointCloud<PointType>::Ptr cloud,int k,
				pcl::search::KdTree<PointType>::Ptr kdtree,int index)
{
	vector<int> KIdx,RIdx;
	vector<float> KDist,RDist;

	kdtree->nearestKSearch(index, k, KIdx, KDist);
	cout<<k-1<<endl;
	cout<<KDist[k-1]<<endl;
	kdtree->radiusSearch(cloud->points[index],KDist[k-1],RIdx,RDist);
	pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);     
	for(int j=0;j<RIdx.size();j++)
		ctmp->points.push_back(cloud->points[RIdx[j]]);

	// Eval
	vector<double> eval;
	eval.resize(3);
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covariance;
	pcl::compute3DCentroid(*ctmp, centroid);
	pcl::computeCovarianceMatrixNormalized(*ctmp, centroid, covariance);
	double a=covariance(0,0); double b=covariance(0,1); double c=covariance(0,2);
	double d=covariance(1,1); double e=covariance(1,2); double f=covariance(2,2);
	double e2=e*e; double b2=b*b;double c2=c*c;
	double B=-(d+f+a); double C=d*f+a*d+a*f-e2-b2-c2;
	double D=-a*d*f+a*e2+b2*f-2*b*e*c+c2*d;
	double p=(3*C-B*B)/3.0; double q=(27.0*D-9*B*C+2*B*B*B)/27.0;
	double r=sqrt(-p*p*p/27.0);
	double theta=1.0/3*acos(-q/(2.0*r));
	double tmp1=2*pow(r,1/3.0); double tmp2=B/3.0;
	eval[0]=tmp1*cos(theta)-tmp2;
	eval[1]=tmp1*cos(theta+2.0/3.0*M_PI)-tmp2;
	eval[2]=tmp1*cos(theta+4.0/3.0*M_PI)-tmp2;
	sort(eval.begin(),eval.end());

	return eval[0];
}

double GetMEvalKNeighours(pcl::PointCloud<PointType>::Ptr cloud,int k,
				pcl::search::KdTree<PointType>::Ptr kdtree,int index)
{
	vector<int> idx(k);
	vector<float> dist(k);
	kdtree->nearestKSearch(cloud->points[index], k, idx, dist);
	pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);     
	for(int j=0;j<idx.size();j++)
		ctmp->points.push_back(cloud->points[idx[j]]);

	// Eval
	vector<double> eval;
	eval.resize(3);
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covariance;
	pcl::compute3DCentroid(*ctmp, centroid);
	pcl::computeCovarianceMatrixNormalized(*ctmp, centroid, covariance);
	double a=covariance(0,0); double b=covariance(0,1); double c=covariance(0,2);
	double d=covariance(1,1); double e=covariance(1,2); double f=covariance(2,2);
	double e2=e*e; double b2=b*b;double c2=c*c;
	double B=-(d+f+a); double C=d*f+a*d+a*f-e2-b2-c2;
	double D=-a*d*f+a*e2+b2*f-2*b*e*c+c2*d;
	double p=(3*C-B*B)/3.0; double q=(27.0*D-9*B*C+2*B*B*B)/27.0;
	double r=sqrt(-p*p*p/27.0);
	double theta=1.0/3*acos(-q/(2.0*r));
	double tmp1=2*pow(r,1/3.0); double tmp2=B/3.0;
	eval[0]=tmp1*cos(theta)-tmp2;
	eval[1]=tmp1*cos(theta+2.0/3.0*M_PI)-tmp2;
	eval[2]=tmp1*cos(theta+4.0/3.0*M_PI)-tmp2;
	sort(eval.begin(),eval.end());

	return eval[0];
}

double GetMEvalRadius(pcl::PointCloud<PointType>::Ptr cloud,double radius,
				pcl::search::KdTree<PointType>::Ptr kdtree,int index)
{
	vector<int> idx;
	vector<float> dist;
	kdtree->radiusSearch(index, radius, idx, dist);
	pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);     
	for(int j=0;j<idx.size();j++)
		ctmp->points.push_back(cloud->points[idx[j]]);

	// Eval
	vector<double> eval;
	eval.resize(3);
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covariance;
	pcl::compute3DCentroid(*ctmp, centroid);
	pcl::computeCovarianceMatrixNormalized(*ctmp, centroid, covariance);
	double a=covariance(0,0); double b=covariance(0,1); double c=covariance(0,2);
	double d=covariance(1,1); double e=covariance(1,2); double f=covariance(2,2);
	double e2=e*e; double b2=b*b;double c2=c*c;
	double B=-(d+f+a); double C=d*f+a*d+a*f-e2-b2-c2;
	double D=-a*d*f+a*e2+b2*f-2*b*e*c+c2*d;
	double p=(3*C-B*B)/3.0; double q=(27.0*D-9*B*C+2*B*B*B)/27.0;
	double r=sqrt(-p*p*p/27.0);
	double theta=1.0/3*acos(-q/(2.0*r));
	double tmp1=2*pow(r,1/3.0); double tmp2=B/3.0;
	eval[0]=tmp1*cos(theta)-tmp2;
	eval[1]=tmp1*cos(theta+2.0/3.0*M_PI)-tmp2;
	eval[2]=tmp1*cos(theta+4.0/3.0*M_PI)-tmp2;
	sort(eval.begin(),eval.end());

	return eval[0];
}

double GetCounterAmongRadius(pcl::PointCloud<PointType>::Ptr cloud,double radius,
				pcl::search::KdTree<PointType>::Ptr kdtree,int index)
{
	vector<int> idx;
	vector<float> dist;
	kdtree->radiusSearch(index, radius, idx, dist);
	return idx.size()/pow(radius,3);
}				

double GetMEval2(pcl::PointCloud<PointType>::Ptr cloud,int k,
				 pcl::search::KdTree<PointType>::Ptr kdtree,int index)
{
	vector<int> idx(k);
	vector<float> dist(k);
	kdtree->nearestKSearch(cloud->points[index], k, idx, dist);
	pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);     
	for(int j=0;j<idx.size();j++)
		ctmp->points.push_back(cloud->points[idx[j]]);

	Eigen::Vector4f centroid;
	Eigen::Matrix3f covariance;
	pcl::compute3DCentroid(*ctmp, centroid);
	pcl::computeCovarianceMatrixNormalized(*ctmp, centroid, covariance);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
	// Eigen::Matrix3f eig_vec = eigen_solver.eigenvectors();
	Eigen::Vector3f eig_val = eigen_solver.eigenvalues();
	
	// eigenvalue
	return eig_val(0);
}
void GetEvalAndEvec(pcl::PointCloud<PointType>::Ptr ctmp,vector<double>& eval,vector<V3>& evec)
{
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covariance;
	pcl::compute3DCentroid(*ctmp, centroid);
	pcl::computeCovarianceMatrixNormalized(*ctmp, centroid, covariance);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
	Eigen::Vector3f eig_val = eigen_solver.eigenvalues();
	Eigen::Matrix3f eig_vec = eigen_solver.eigenvectors();
	// eigenvalue
	eval[0]=eig_val(0);
	eval[1]=eig_val(1);
	eval[2]=eig_val(2);

	// eigenvector
	eig_vec.col(2) = eig_vec.col(0).cross(eig_vec.col(1));
	eig_vec.col(0) = eig_vec.col(1).cross(eig_vec.col(2));
	eig_vec.col(1) = eig_vec.col(2).cross(eig_vec.col(0));
	evec[0].x=eig_vec(0,0);evec[0].y=eig_vec(1,0);evec[0].z=eig_vec(2,0);
	evec[1].x=eig_vec(0,1);evec[1].y=eig_vec(1,1);evec[1].z=eig_vec(2,1);
	evec[2].x=eig_vec(0,2);evec[2].y=eig_vec(1,2);evec[2].z=eig_vec(2,2);
}


double ComputeNearestDistance(pcl::search::KdTree<PointType>::Ptr kdtree,int index)
{
	std::vector<int> idx(2);
	std::vector<float> dist(2);
	kdtree->nearestKSearch(index, 2, idx, dist);
	return dist[1];
}


// double ComputeDB2(pcl::search::KdTree<PointType>::Ptr kdtree,int K,int index)
// {
// 	vector<int> idx(K);
// 	vector<float> dist(K);
// 	kdtree->nearestKSearch(index, K, idx, dist);
// 	vector<double> db;
// 	DaubechiesWavelet(dist,db);
// 	db.pop_back();
// 	double db2_max=VectorMaximum(db);
// 	return db2_max;
// }

int GetIndex(pcl::search::KdTree<PointType>::Ptr kdtree,PointType ptmp)
{
	vector<int> idx(1);
	vector<float> dist(1);
	kdtree->nearestKSearch(ptmp,1,idx,dist);
	return idx[0];
}

double GetCellSize(pcl::PointCloud<PointType>::Ptr cloud, int level)
{
	PointType pmin,pmax;
	pcl::getMinMax3D(*cloud,pmin,pmax);
	double cellsize=MAX3(abs(pmax.x-pmin.x),abs(pmax.y-pmin.y),abs(pmax.z-pmax.z));
	return cellsize/pow(2,level);
}
double GetBoxMin(pcl::PointCloud<PointType>::Ptr cloud)
{
	PointType pmin,pmax;
    pcl::getMinMax3D(*cloud,pmin,pmax);
	double vmin=MIN3(abs(pmax.x-pmin.x),abs(pmax.y-pmin.y),abs(pmax.z-pmin.z));
	return vmin;
}

double GetBoxMax(pcl::PointCloud<PointType>::Ptr cloud)
{
	PointType pmin,pmax;
    pcl::getMinMax3D(*cloud,pmin,pmax);
	double vmax=MAX3(abs(pmax.x-pmin.x),abs(pmax.y-pmin.y),abs(pmax.z-pmin.z));
	return vmax;
}