#ifndef SIMARRAYS_H
#define SIMARRAYS_H

#include <vector>
#include <iostream>
#include <TRandom3.h>

// ROOT includes
#include <TH3.h>
#include <TF3.h>

// PCL includes
#include <pcl/point_cloud.h>
#include <pcl/octree/octree_search.h>

using namespace std;

class Field{
public:
  
  Field();

  vector<double> pot;           // V    - Potential
  vector<double> e;             // V/cm - Norm of Efield
  vector<double> ex;            // V/cm - Efield x-comp.
  vector<double> ey;            // V/cm - Efield y-comp.
  vector<double> ez;            // V/cm - Efield z-comp.
};

class Boundaries{
public:
  
  Boundaries();
  void SetInitialPoint(double, double, double, double);
  void SetInitialPoint(int, int, double, double, double, double, double, double);
  void SetVertexPoint(int, int, double, double, double, double, double);
  void PrintBoundaries();
  
  bool set;                       // Flag to know whether it has been properly set
  double xmin, xmax;              // Geometry boundary
  double ymin, ymax;              // -
  double zmin, zmax;              // -
  double hole_r,hole_h,hole_eps;  // PCB hole dimensions for collection stop
  double hole_xs, hole_ys, hole_v;// -
  double zexcmin, zexcmax;        // Inner volume exclusion (like the
  double xc, yc, zc;              // Center point
  double ti, xi, yi, zi;          // Initial point
  
  vector<double> strip_p;         // cm  - strips pitch
  vector<double> strip_a;         // deg - strips angle

};

class Offset{
public:
  
  Offset();
  bool set;
  double t, x, y, z;
};


class Electron{
public:

  Electron();
  void SetVerbosity(int db){debug = db;}
  void Clear();
  void Resize(bool*);
  void Resize(int,bool*);
  void Copy(int,bool*,Electron);
  void Copy(bool*,Electron);
  void SetSpeed();
  void Propagate(double);
  void SetCurrent(int);
  void ExtrapolateToSurface(double,bool*,Boundaries);
  void ApplyBackOffset(Offset);
  void SetFinalShifts(Boundaries,Offset);
  
  int debug;
  int maxstep;
  int step;                         // #     - Number of calculation points
  int n;                            // #     - Number of electrons in the cloud
  int coll_flag;                    // #     - Flag to check for collection
  
  vector<double> t;                 // us    - Time
  vector<double> v;                 // cm/us - Velocity (v, vx, vy ,vz)
  vector<double> x,  y,  z;         // cm    - Position
  vector<vector<double> > cur;      // au    - Current
  vector<double> i1c,i2c,i3c;       // au    - Convoluted current
  vector<double> w;                 // au    - W-pot from view 1, 2, 3
  vector<double> wf, wfx, wfy, wfz; // au    - W-field from view 1, 2, 3
  
  double ion;                // MeV   - ionisation energy
  double qe;                 // fC    - charge of the electron
  double qc;                 // fC    - charge of the electron cloud
  double ds;                 // cm    - initial cloud extension
  
  double ex, ey, ez, e;      // V/cm  - Efield
  double sx, sy, sz, st;     // cm,us - shift in time and space w.r.t starting position
};

class Noise{
public:
  
  Noise(int);
  void Clear();
  void Resize();
  void SetWhiteNoise(double,double);
  
  int debug;
  int maxstep;
  vector<double> t;
  vector<double> cur;
  vector<double> curc;
};

class Strips{
public:
  
  Strips(int);
  void Clear();
  void Resize();
  void SetPosition(Boundaries,Boundaries,int,int);
  void SetTime(Electron);
  void DumpSignal();
  void AddSignal(int,double,int,Electron);
  void AddNoise(Noise);

  bool set_time, set_pos;     // False to set it once
  int debug;
  int maxstep;
  int id;                     // #   - Strip ID
  double x, y;                // cm  - Strip position
  
  vector<double> t;           // us  - Time
  vector<double> cur;           // fC  - Charge of the signal
  vector<double> curc;          // ADC - Convolved signal
};


#endif
