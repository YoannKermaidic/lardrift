#include "simarrays.hh"

const int MAX_STEPS= 500; // #  - maximum number of calculated points.
double Eion = 23.6e-6;    // MeV - ionisation energy
double Qe   = 1.602e-4;   // fC - electron charge

Field::Field(){}

Electron::Electron(){
  debug     = 0;
  maxstep   = MAX_STEPS;
  ion       = Eion;
  qe        = Qe;
  sx        = 0;
  sy        = 0;
  sz        = 0;
  st        = 0;
  coll_flag = 0;
}

void Electron::Clear(){
  t.clear(); x.clear(); y.clear(); z.clear(); v.clear();
  w.clear(); wf.clear(); wfx.clear(); wfy.clear(); wfz.clear();
  cur.clear(); i1c.clear(); i2c.clear(); i3c.clear();
}

void Electron::Resize(int elSize, bool* isw){
  t.resize(elSize,0); v.resize(4,0);
  x.resize(elSize,0);  y.resize(elSize,0); z.resize(elSize,0);
  if(isw[1] || isw[2] || isw[3]){
    cur.resize(4);
    for(int nw=0;nw<4;nw++) cur[nw].resize(elSize,0);
    w.resize(4,0);        wf.resize(4,0);       wfx.resize(4,0);      wfy.resize(4,0);      wfz.resize(4,0);
    i1c.resize(elSize,0); i2c.resize(elSize,0); i3c.resize(elSize,0);
  }
}

void Electron::Resize(bool* isw){
  t.resize(MAX_STEPS,0);  v.resize(4,0);
  x.resize(MAX_STEPS,0);  y.resize(MAX_STEPS,0); z.resize(MAX_STEPS,0);
  if(isw[1] || isw[2] || isw[3]){
    cur.resize(4);
    for(int nw=0;nw<4;nw++) cur[nw].resize(MAX_STEPS,0);
    w.resize(4,0);           wf.resize(4,0);          wfx.resize(4,0);         wfy.resize(4,0);      wfz.resize(4,0);
    i1c.resize(MAX_STEPS,0); i2c.resize(MAX_STEPS,0); i3c.resize(MAX_STEPS,0);
  }
}

void Electron::Copy(int i,bool* isw,Electron el){
  if(i==0){
    step    = el.step;
    maxstep = el.maxstep;
    n       = el.n;
    ds      = el.ds;
    qe      = el.qe;
    qc      = el.qc;
    ion     = el.ion;
    
    sx      = el.sx;
    sy      = el.sy;
    sz      = el.sz;
    st      = el.st;
  }
  
  t[i] = el.t[i];
  x[i] = el.x[i];
  y[i] = el.y[i];
  z[i] = el.z[i];
  if(isw[1] || isw[2] || isw[3]){
    cur[1][i]  = el.cur[1][i];
    cur[2][i]  = el.cur[2][i];
    cur[3][i]  = el.cur[3][i];
    i1c[i] = el.i1c[i];
    i2c[i] = el.i2c[i];
    i3c[i] = el.i3c[i];
  }
}

void Electron::Copy(bool* isw,Electron el){
  step    = el.step;
  maxstep = el.maxstep;
  n       = el.n;
  ds      = el.ds;
  qe      = el.qe;
  qc      = el.qc;
  ion     = el.ion;
  
  sx      = el.sx;
  sy      = el.sy;
  sz      = el.sz;
  st      = el.st;
  
  for(int i=0;i<el.maxstep;i++){
    t[i] = el.t[i];
    x[i] = el.x[i];
    y[i] = el.y[i];
    z[i] = el.z[i];
    if(isw[1] || isw[2] || isw[3]){
      cur[1][i]  = el.cur[1][i];
      cur[2][i]  = el.cur[2][i];
      cur[3][i]  = el.cur[3][i];
      i1c[i] = el.i1c[i];
      i2c[i] = el.i2c[i];
      i3c[i] = el.i3c[i];
    }
  }
}

void Electron::SetSpeed(){ 
  // Compute speed vector according to E-field and speed norm
  v[1] = -ex/e * v[0];
  v[2] = -ey/e * v[0];
  v[3] = -ez/e * v[0];
}

void Electron::SetCurrent(int view){cur[view][step] = -1. * (wfx[view]*v[1] + wfy[view]*v[2] + wfz[view]*v[3]);}

void Electron::Propagate(double dt){
  t[step] = t[step-1]+dt;
  x[step] = x[step-1]+v[1]*dt;
  y[step] = y[step-1]+v[2]*dt;
  z[step] = z[step-1]+v[3]*dt;
}

void Electron::ExtrapolateToSurface(double coll_flag,bool* isw,Boundaries bounds){
  if(coll_flag==1){
    cout << "LArDrift::Electron has reached the bottom collection plane -> perform extrapolation." << endl;
    x[step] = x[step-1] + (x[step]-x[step-1])/(z[step]-z[step-1])*(bounds.hole_h-z[step-1]);
    y[step] = y[step-1] + (y[step]-y[step-1])/(z[step]-z[step-1])*(bounds.hole_v-z[step-1]);
    z[step] = bounds.hole_h;
    if(isw[3]) w[3] = 1;
  }
  else if(coll_flag==2){
    cout << "LArDrift::Electron has reached the top collection plane -> perform extrapolation." << endl;
    x[step] = x[step-1] + (x[step]-x[step-1])/(z[step]-z[step-1])*(2*bounds.hole_h+bounds.hole_v-z[step-1]);
    y[step] = y[step-1] + (y[step]-y[step-1])/(z[step]-z[step-1])*(2*bounds.hole_h+bounds.hole_v-z[step-1]);
    z[step] = 2*bounds.hole_h+bounds.hole_v;
    if(isw[3]) w[3] = 1;
  }
}

void Electron::ApplyBackOffset(Offset offset){
  for(int i=0;i<maxstep;i++){
    t[i] += offset.t;
    x[i] += offset.x;
    y[i] += offset.y;
    z[i] += offset.z;
  }
}

void Electron::SetFinalShifts(Boundaries bounds,Offset offset){
  sx = x[step] - (bounds.xi - offset.x);
  sy = y[step] - (bounds.yi - offset.y);
  sz = z[step] - (bounds.zi - offset.z);
  st = t[step] - bounds.ti;
}

Boundaries::Boundaries(){
  xmin = 1e6; xmax = -1e6;
  ymin = 1e6; ymax = -1e6;
  zmin = 1e6; zmax = -1e6;

  zexcmin = 1e6; zexcmax = -1e6;

  xi=-1000, yi=-1000, zi=-1000;
  xc = 0; yc = 0; zc = 0;

  hole_r  = 9999; hole_h = 9999; hole_v = 9999; hole_eps = 0.001;
  hole_xs = 9999; hole_ys = 9999;
  
  set = false;
}

void Boundaries::SetInitialPoint(int i, int Nsims, double x1, double x2, double y1, double y2, double z1, double z2){
  xi = x1+(x2-x1)*(i+1)/Nsims;
  yi = y1+(y2-y1)*(i+1)/Nsims;
  zi = z1+(z2-z1)*(i+1)/Nsims;
}

void Boundaries::SetInitialPoint(double t, double x, double y, double z){
  ti = t, xi = x, yi = y, zi = z;
}

void Boundaries::SetVertexPoint(int i, int Nsims, double phi, double theta, double L1, double L2, double z0){
  int split = (L2 > 0) ? Nsims - int(Nsims / (L1/L2+1)) : Nsims;
  
  ti = 0;
  if(i<split){
    xi = L1 * i/double(split) * cos(phi);
    yi = L1 * i/double(split) * sin(phi);
    zi = z0 + L1 * i/double(split) * sin(theta);
  }
  else{
    xi =  L2 * (i-split)/double(Nsims-1-split) * cos(phi);
    yi =  L2 * (i-split)/double(Nsims-1-split) * sin(phi);
    zi = z0 -L2 * (i-split)/double(Nsims-1-split) * sin(theta);
  }
}

void Boundaries::PrintBoundaries(){
  cout << "PrintBoundaries::Center" << endl;
  cout << "  (" << xc << "," << yc << "," << zc << ")" << endl;
  cout << "PrintBoundaries::Initial point" << endl;
  cout << "  (" << xi << "," << yi << "," << zi << ")" << endl;
  cout << "PrintBoundaries::[min,max]" << endl;
  cout << "  x -> [" << xmin << "," << xmax << "]" << endl;
  cout << "  y -> [" << ymin << "," << ymax << "]" << endl;
  cout << "  z -> [" << zmin << "," << zmax << "]" << endl;
  cout << "PrintBoundaries::Exclusion" << endl;
  cout << "  z -> [" << zexcmin << "," << zexcmax << "]" << endl;
  if(hole_r != 9999){
    cout << "PrintBoundaries::PCB hole" << endl;
    cout << "  radius    -> "  << hole_r  << " cm" << endl;
    cout << "  height    -> "  << hole_h  << " cm" << endl;
    cout << "  x spacing -> "  << hole_xs << " cm" << endl;
    cout << "  y spacing -> "  << hole_ys << " cm" << endl;
    cout << "  z spacing -> "  << hole_v  << " cm" << endl;
  }
}

Offset::Offset(){
  set = false;
  x   = 0;
  y   = 0;
  z   = 0;
  t   = 0;
}

Strips::Strips(int db){
  debug     = db;
  maxstep = MAX_STEPS;
  set_pos   = false;
  set_time  = false;
  
  id = 0; x = 0; y = 0;
}

void Strips::Resize(){
  t.resize(maxstep,0);
  cur.resize(maxstep,0);
  curc.resize(maxstep,0);
}

void Strips::Clear(){
  t.clear();
  cur.clear();
  curc.clear();
}

void Strips::SetPosition(Boundaries det, Boundaries anode, int vw, int st){
  set_pos = true;
  id = st;
  
  if(vw == 0){     x = det.xmin + st*anode.strip_p[vw]/sin(anode.strip_a[vw]); y = det.ymin;}
  else if(vw == 1){x = det.xmin;                                               y = det.ymin + st*anode.strip_p[vw];}
  else if(vw == 2){x = det.xmin + st*anode.strip_p[vw];                        y = det.ymin;}
}

void Strips::SetTime(Electron el){
  set_time = true;
  
  for(int i=0;i<maxstep;i++) t[i] += el.t[i];
}

void Strips::DumpSignal(){
  cout << "DumpSignal::t current";
  for(int i=0;i<maxstep;i++) cout << "  " << i << " " << t[i] << " " << cur[i] << endl;
  cout << "DumpSignal::Done";
}

void Strips::AddSignal(int offset_min, double dt, int vw, Electron el){
  if(debug) cout << "AddSignal::Start ... ";
  
  // Move signal in time according to the arrival time from detector volume drift
  int offset = (el.t[0]/dt) - offset_min;
  
  if(debug)cout << "AddSignal::Offset: " << offset << endl;
  
  for(int i=0;i<offset;i++) cur[i] = 0;
  for(int i=maxstep-offset;i<maxstep;i++) cur[i] = 0;
  
  for(int i=0;i<maxstep-offset;i++){
    if(vw==0)      cur[i+offset] += el.cur[1][i];
    else if(vw==1) cur[i+offset] += el.cur[2][i];
    else if(vw==2) cur[i+offset] += el.cur[3][i];
  }
  if(debug) cout << "done" << endl;
}

void Strips::AddNoise(Noise noise){
  if(debug) cout << "AddNoise::Start ... ";
  for(int i=0;i<maxstep;i++) cur[i] += noise.cur[i];
  if(debug) cout << "done" << endl;
}

Noise::Noise(int db){
  debug = db;
  maxstep = MAX_STEPS;
}

void Noise::Clear(){
  t.clear();
  cur.clear();
  curc.clear();
}

void Noise::Resize(){
  t.resize(maxstep,0);
  cur.resize(maxstep,0);
  curc.resize(maxstep,0);
}

void Noise::SetWhiteNoise(double mean, double rms){
  if(debug) cout << "SetWhiteNoise::Gaussian noise with ("<< mean << "," << rms << ") parameters." << endl;
  for(int i=0;i<maxstep;i++) cur[i] = gRandom->Gaus(mean,rms);
}
