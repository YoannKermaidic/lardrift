#ifndef READER_H
#define READER_H

// System includes
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// ROOT includes
#include <TFile.h>
#include <TTree.h>

// PCL includes
#include <pcl/point_cloud.h>
#include <pcl/octree/octree_search.h>

// JSON includes
#include <nlohmann/json.hpp>

// LArDrift
#include <simarrays.hh>

using namespace std;

class Reader{
public:
  
  Reader(int);
  void SetVerbosity(int db){debug=db;};
  vector<Electron> GetData(string,double,double,double);
  vector<Electron> GetLarDrift(string);
  vector<Electron> GetGeant4(string,Boundaries,bool*,double &,int);
  vector<vector<double> > GetER(string,double);
  void GetComsol(string, pcl::PointCloud<pcl::PointXYZ>::Ptr &, Field &,bool &,bool*);
  void GetBoundaries(Boundaries &,string,string);
  
  
private:
  int debug;
  double mmTOcm;               // mm to cm conversion
  double mTOcm;                // m  to cm conversion
  double mVToADC;              // mV to ADC TDE conversion
  double degTOrad;             // rad/deg conversion
};
#endif
