#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <string>
#include <iostream>


// PCL includes
#include <pcl/point_cloud.h>
#include <pcl/octree/octree_search.h>

#include <simarrays.hh>


using namespace std;

class Tools{
public:
  
  Tools(int,string,double);
  void SetVerbosity(int);
  void SetExpSetup(string setup){exp_setup=setup;};
  void GetBoundaries(Boundaries);
  double EvalEspeed(double);
  
  void SetResolution(){resolution = (zmax-zmin)/8.;};
  void SetEOctree();
  void SetWOctrees();
  void SetEField(Electron &, bool);
  void SetWFields(int, Electron &);
  
  bool InDetector(double, double, double);
  int  IsCollected(double, double, double);
  
  void PrintStatus(Electron);
  void PrintEndPoint(Electron);
  int  GetDriftTimeOffset(double, vector<bool>, vector<Electron>);
  vector<int> GetNumberOfStrips(Boundaries,Boundaries);
  int GetStripID(int,Boundaries,Boundaries,int,Electron);
  void ERconvolution(vector<vector<double> >,Strips &);
  void ComputePCBoffset(double,double,double,double,Offset &);
  void PrintHelp();
  
  vector<Field> fields;
  vector<pcl::PointCloud<pcl::PointXYZ>::Ptr > clouds;
  vector<pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> > octrees;

private:
  
  int debug;
  string exp_setup;
  float  resolution;              // octree resolution (see PCL library)
  double xmin, xmax;              // Geometry boundary
  double ymin, ymax;              // -
  double zmin, zmax;              // -
  double hole_r,hole_h,hole_eps;  // PCB hole dimensions for collection stop
  double hole_xs, hole_ys, hole_v;// -
  double zexcmin, zexcmax;        // Inner volume exclusion (like the
  double xc, yc, zc;              // Center point
  double T_lar;                   // Liquid argon temperature
};
#endif

