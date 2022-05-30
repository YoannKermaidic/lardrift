#ifndef PLOTTER_H
#define PLOTTER_H

#include <vector>
#include <string>
#include <iostream>

#include <TLegend.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TAxis.h>
#include <TAxis3D.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TView.h>
#include <TColor.h>
#include <TObject.h>
#include <TStyle.h>

#include <simarrays.hh>
#include <tools.hh>

using namespace std;

class Plotter{
public:
  
  Plotter(int,Boundaries, double);
  virtual ~Plotter();
  
  void SetBoundaries(Boundaries);
  TCanvas* TwoDshift(vector<bool>,vector<Electron>);
  TCanvas* OneDshift(vector<bool>,vector<Electron>);
  TCanvas* OneDdrift(vector<bool>,vector<Electron>);
  TCanvas* TwoDTracks(vector<bool>,vector<Electron>);
  TCanvas* TwoDProj(vector<bool>,vector<Electron>, vector<Electron>);
  TCanvas* ThreeDProj(vector<bool>,vector<Electron>, vector<Electron>);
  TCanvas* ThreeDTracks(vector<bool>,vector<Electron>);
//TCanvas* Current(string,vector<Electron>, bool *);
//TCanvas* Wpot(vector<Electron>, bool *);
  TCanvas* ER(vector<vector<double> >);
  TCanvas* SpeedModel(TF1*);
  TCanvas* SpeedModel(Tools);
  TCanvas* EventDisplay(string,int, vector<Strips>);
  TCanvas* FullEventDisplay(string, vector<vector<Strips> >);
  vector<TGraph*> Box();
  
private:
  
  int debug;
  // Plane for plotting projection
  vector<string> axaxis;
  vector<string> planes;
  vector<string> pxaxis;
  vector<string> pyaxis;
  double dt_start, xmin, xmax, ymin, ymax, zmin, zmax;
};

#endif
