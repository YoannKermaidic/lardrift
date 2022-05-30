#include "tools.hh"

Tools::Tools(int db, string setup, double T){
  debug = db;
  exp_setup = setup;
  
  T_lar   = T;

  octrees.resize(4,8);
  clouds.resize(4,pcl::PointCloud<pcl::PointXYZ>::Ptr (new pcl::PointCloud<pcl::PointXYZ>));
  fields.resize(4,Field());
}

void Tools::SetVerbosity(int db){debug=db;};

void Tools::GetBoundaries(Boundaries bdr){
  xmin = bdr.xmin; xmax = bdr.xmax;
  ymin = bdr.ymin; ymax = bdr.ymax;
  zmin = bdr.zmin; zmax = bdr.zmax;
  zexcmin = bdr.zexcmin; zexcmax = bdr.zexcmax;
  xc = bdr.xc; yc = bdr.yc; zc = bdr.zc;
  hole_r = bdr.hole_r; hole_h = bdr.hole_h; hole_eps = bdr.hole_eps;
  hole_xs = bdr.hole_xs; hole_ys = bdr.hole_ys; hole_v = bdr.hole_v;
};

double Tools::EvalEspeed(double x){
  double vd = 0;
  double E = x /1000.; // Conversion to kV/m
  
  // Walkowiak measurement temperature / parameters
  double T_walk = 90.371;
  double dT = T_lar-T_walk;
  
  vector<double> pars = {-0.01481, -0.0075, 0.141, 12.4, 1.627, 0.317,
    -0.03229, 6.231, -10.62, 12.74, -9.112, 2.83};
  
  if(E > 0.5) vd = (pars[0]*dT+1) * (pars[2]*E*TMath::Log(1+pars[3]/E)+pars[4]*pow(E,pars[5])) + pars[1]*dT;
  else{
    for(int i=6;i<12;i++)
      vd += pars[i]*pow(E,i-6);
    
    double Etmp = 0.5;
    double tmp1 = 0;
    
    for(int i=6;i<12;i++)
      tmp1 += pars[i]*pow(Etmp,i-6);
    
    double tmp2 = (pars[0]*dT+1) * (pars[2]*Etmp*TMath::Log(1+pars[3]/Etmp)+pars[4]*pow(Etmp,pars[5])) + pars[1]*dT;
    
    vd = vd * tmp2/tmp1;
  }
  
  return vd*0.1;
}

bool Tools::InDetector(double posx, double posy, double posz){
    
  if(exp_setup=="coldbox"){
    if(posz<zmin || posz>zmax || posx<1.1*xmin || posx>1.1*xmax || posy<1.1*ymin || posy>1.1*ymax){
      if(debug) cout << "InDetector::Point (" << posx << "," << posy << "," << posz << ") not within boundaries." << endl;
      return false;
    }
  }
  else{
    if(posx<xmin || posx>xmax || posy<ymin || posy>ymax || posz<zmin || posz>zmax){
      if(debug) cout << "InDetector::Point (" << posx << "," << posy << "," << posz << ") not within boundaries." << endl;
      return false;
    }
  }
  if(posz>=zexcmin && posz<=zexcmax){
    if(debug) cout << "InDetector::Point (" << posx << "," << posy << "," << posz << ") within excluded region." << endl;
    return false;
  }
  return true;
}

int Tools::IsCollected(double posx, double posy, double posz){
  double radius = sqrt((posx-xc)*(posx-xc) + (posy-yc)*(posy-yc));
  double rx     = sqrt((fabs(posx-xc)-hole_xs)*(fabs(posx-xc)-hole_xs) + (posy-yc)*(posy-yc));
  double rxy    = sqrt((fabs(posx-xc)-hole_xs/2.)*(fabs(posx-xc)-hole_xs/2.) + (fabs(posy-yc)-hole_ys)*(fabs(posy-yc)-hole_ys));
  
  // check collection PCB plane
  if(radius < hole_r || rx < hole_r || rxy < hole_r)                          return 0; // within central/neighbour hole
  else if(posz <= hole_h+hole_eps && posz >= -hole_eps)                       return 1; // within shield     plane
  else if(posz <= 2*hole_h+hole_v+hole_eps && posz >= hole_h+hole_v-hole_eps) return 2; // within collection plane
  return 0;
}

void Tools::SetEOctree(){
  octrees[0] = pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> (resolution);
  octrees[0].setInputCloud(clouds[0]);
  octrees[0].addPointsFromInputCloud();
}

void Tools::SetWOctrees(){
  for(int nw=1;nw<4;nw++){
    octrees[nw] = pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> (resolution);
    octrees[nw].setInputCloud(clouds[nw]);
    octrees[nw].addPointsFromInputCloud();
  }
}

void Tools::SetEField(Electron &el, bool ise){
  
  pcl::PointXYZ searchPoint;
  searchPoint.x = el.x[el.step];
  searchPoint.y = el.y[el.step];
  searchPoint.z = el.z[el.step];
  
  if(debug>1) std::cout << "SetEField::Search neighbors at ("
  << searchPoint.x << " " << searchPoint.y << " " << searchPoint.z << ")" << std::endl;
  
  // K nearest neighbor search
  int K = 20;
  double dx = 1, dy = 1, dz = 1;
  TH3D* hPot;
  std::vector<int> pointIdxNKNSearch;
  std::vector<float> pointNKNSquaredDistance;
  
  if (octrees[0].nearestKSearch (searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
  {
    if(ise){
      double ptmp = 0, ftmp = 0, fxtmp = 0, fytmp = 0, fztmp = 0;
      double sum_dist = 0;
      
      for (std::size_t i = 0; i < pointIdxNKNSearch.size (); ++i){
        if(debug>2){ cout << "      position: ("  << (*clouds[0])[ pointIdxNKNSearch[i] ].x
          << "," << (*clouds[0])[ pointIdxNKNSearch[i] ].y
          << "," << (*clouds[0])[ pointIdxNKNSearch[i] ].z
          << ") -> (distance: " << sqrt(pointNKNSquaredDistance[i]) << " cm) " << endl;
          cout << "      V = "  << fields[0].pot[pointIdxNKNSearch[i]] << " V"    << endl;
          cout << "      E = "  << fields[0].e[pointIdxNKNSearch[i]]   << " V/cm" << endl;
          cout << "      Ex = " << fields[0].ex[pointIdxNKNSearch[i]]  << " V/cm" << endl;
          cout << "      Ey = " << fields[0].ey[pointIdxNKNSearch[i]]  << " V/cm" << endl;
          cout << "      Ez = " << fields[0].ez[pointIdxNKNSearch[i]]  << " V/cm" << endl;
          if(!InDetector((*clouds[0])[ pointIdxNKNSearch[i] ].x,(*clouds[0])[ pointIdxNKNSearch[i] ].y,(*clouds[0])[ pointIdxNKNSearch[i] ].z))
            cout << "      -> Skipped (on boundary) "                         << endl;
          cout << "      "                                                    << endl;
        }
        if(!InDetector((*clouds[0])[ pointIdxNKNSearch[i] ].x,(*clouds[0])[ pointIdxNKNSearch[i] ].y,(*clouds[0])[ pointIdxNKNSearch[i] ].z)) continue;
        if(IsCollected((*clouds[0])[ pointIdxNKNSearch[i] ].x,(*clouds[0])[ pointIdxNKNSearch[i] ].y,(*clouds[0])[ pointIdxNKNSearch[i] ].z)) continue;
        
        ptmp  += fields[0].pot[pointIdxNKNSearch[i]]/ pointNKNSquaredDistance[i];
        ftmp  += fields[0].e[pointIdxNKNSearch[i]]  / pointNKNSquaredDistance[i];
        fxtmp += fields[0].ex[pointIdxNKNSearch[i]] / pointNKNSquaredDistance[i];
        fytmp += fields[0].ey[pointIdxNKNSearch[i]] / pointNKNSquaredDistance[i];
        fztmp += fields[0].ez[pointIdxNKNSearch[i]] / pointNKNSquaredDistance[i];
        sum_dist += 1. / pointNKNSquaredDistance[i];
      }
      el.e  = ftmp  / sum_dist;
      el.ex = fxtmp / sum_dist;
      el.ey = fytmp / sum_dist;
      el.ez = fztmp / sum_dist;
      
    }
    else{
      hPot = new TH3D("hPot","hPot",int((xmax-xmin)/dx),xmin,xmax,int((ymax-ymin)/dy),ymin,ymax,int((zmax-zmin)/dz),zmin,zmax);

      if(debug>1) std::cout << "SetEfield::Neighbors within K nearest search:" << std::endl;
      for (std::size_t i = 0; i < pointIdxNKNSearch.size (); ++i){
        double xnear = (*clouds[0])[ pointIdxNKNSearch[i] ].x;
        double ynear = (*clouds[0])[ pointIdxNKNSearch[i] ].y;
        double znear = (*clouds[0])[ pointIdxNKNSearch[i] ].z;
        double vnear = fields[0].pot[pointIdxNKNSearch[i]];
        
        //if(sqrt(pointNKNSquaredDistance[i]) < boxsize)
        if(debug>1) std::cout << "    "  << xnear
          << " " << ynear
          << " " << znear
          << " (squared distance: " << sqrt(pointNKNSquaredDistance[i]) << ")"
          << " V = " << vnear << " V" << std::endl;
        
        hPot->Fill(xnear,ynear,znear,vnear);
      }
    }
  }
  
  if(!ise){
    // If no Efield found in COMSOL file -> fit electric potential
    TF3* fPot = new TF3("fPot","[0]*(x-[3])+[1]*(y-[4])+[2]*([5]-z)",xmin,xmax,ymin,ymax,zmin,zmax);
    fPot->SetParameter(0,0.1);
    fPot->SetParameter(1,0.1);
    fPot->SetParameter(2,400);
    fPot->FixParameter(3,xc);
    fPot->FixParameter(4,yc);
    fPot->FixParameter(5,zc);
    hPot->Fit("fPot","QR");
    
    el.ex = fPot->GetParameter(0);
    el.ey = fPot->GetParameter(1);
    el.ez = fPot->GetParameter(2); // V/cm
    el.e  = sqrt(pow(el.ex,2) + pow(el.ey,2) + pow(el.ez,2));
    
    if(debug>1) cout << "SetField::Best Efield results: (" << el.e << "," << el.ex << "," << el.ey << "," << el.ez << ") V/cm" << endl;

    delete hPot;
    delete fPot;
  }
}

void Tools::SetWFields(int view, Electron &el){
  
  pcl::PointXYZ searchPoint;
  searchPoint.x = el.x[el.step];
  searchPoint.y = el.y[el.step];
  searchPoint.z = el.z[el.step];
  
  if(debug>1) std::cout << "SetField::Search view " << view << " neighbors at ("
  << searchPoint.x << " " << searchPoint.y << " " << searchPoint.z << ")" << std::endl;
  
  // K nearest neighbor search
  int K = 20;
  double dx = 1, dy = 1, dz = 1;
  std::vector<int> pointIdxNKNSearch;
  std::vector<float> pointNKNSquaredDistance;
  
  if (octrees[view].nearestKSearch (searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0){
    double ptmp = 0, ftmp = 0, fxtmp = 0, fytmp = 0, fztmp = 0;
    double sum_dist = 0;
    
    for (std::size_t i = 0; i < pointIdxNKNSearch.size (); ++i){
      if(!InDetector((*clouds[view])[ pointIdxNKNSearch[i] ].x,(*clouds[view])[ pointIdxNKNSearch[i] ].y,(*clouds[view])[ pointIdxNKNSearch[i] ].z)) continue;
      if(IsCollected((*clouds[view])[ pointIdxNKNSearch[i] ].x,(*clouds[view])[ pointIdxNKNSearch[i] ].y,(*clouds[view])[ pointIdxNKNSearch[i] ].z)) continue;
      
      ptmp  += fields[view].pot[pointIdxNKNSearch[i]]/ pointNKNSquaredDistance[i];
      ftmp  += fields[view].e[pointIdxNKNSearch[i]]  / pointNKNSquaredDistance[i];
      fxtmp += fields[view].ex[pointIdxNKNSearch[i]] / pointNKNSquaredDistance[i];
      fytmp += fields[view].ey[pointIdxNKNSearch[i]] / pointNKNSquaredDistance[i];
      fztmp += fields[view].ez[pointIdxNKNSearch[i]] / pointNKNSquaredDistance[i];
      sum_dist += 1. / pointNKNSquaredDistance[i];
    }
    el.w[view] = ptmp / sum_dist; el.wf[view] = ftmp / sum_dist;
    el.wfx[view] = fxtmp / sum_dist; el.wfy[view] = fytmp / sum_dist; el.wfz[view] = fztmp / sum_dist;
  }
}

int Tools::GetDriftTimeOffset(double dt, vector<bool> isvalid, vector<Electron> els){
  if(debug) cout << "GetDriftTimeOffset:Start" << endl;
  int offset = 1e6;
  
  for(int sim=0;sim<els.size();sim++){
    if(!isvalid[sim]) continue;
    if(offset > els[sim].t[0]/dt)
      offset = int(els[sim].t[0]/dt);
  }
  cout << "GetDriftTimeOffset::Offset: " << offset << endl;
  return offset;
}

void Tools::PrintStatus(Electron el){
  cout << "PrintStatus::Information at time t = " << el.t[el.step] << " us:" << endl;
  cout << "      position: (" << el.x[el.step]  << "," << el.y[el.step]  << "," << el.z[el.step]  << ") cm"                                           << endl;
  cout << "      Efield:   (" << el.ex          << "," << el.ey          << "," << el.ez          << ") V/cm  -> |E| = " << el.e << " V/cm"           << endl;
  cout << "      velocity: (" << el.v[1] << "," << el.v[2] << "," << el.v[3] << ") cm/us -> |v| = " << el.v[0] << " cm/us" << endl;
}

void Tools::PrintEndPoint(Electron el){
  cout << "PrintEndPoint::Tracks end points:" << endl;
  cout << "  Start: " << el.t[0]       << " " << el.x[0]       << " " << el.y[0]       << " " << el.z[0]       << " " << endl;
  cout << "  End:   " << el.t[el.step] << " " << el.x[el.step] << " " << el.y[el.step] << " " << el.z[el.step] << " " << endl;
  cout << "  Delta: " << el.st         << " " << el.sx         << " " << el.sy         << " " << el.sz         << " " << endl;
}

void Tools::ERconvolution(vector<vector<double> > vER, Strips &strip){

  int const nf = vER[0].size();
  int const ng = strip.maxstep;
  int const n  = nf + ng - 1;
  
  if(debug>2){
    cout << "ERconvolution::Size:"   << endl;
    cout << "  ER func:     " << nf << endl;
    cout << "  Current:     " << ng << endl;
    cout << "  Convolution: " << n  << endl;
  }
  
  double norm = 1;
  vector<double> out(n,0);
  
  for(auto i(0); i < n; ++i) {
    int const jmn = (i >= ng - 1)? i - (ng - 1) : 0;
    int const jmx = (i <  nf - 1)? i            : nf - 1;
    for(auto j(jmn); j <= jmx; ++j) out[i] += vER[1][j] * strip.cur[i - j];
  }
  for(std::size_t i = 0; i < ng; ++i) strip.curc[i] = out[i];
}

void Tools::ComputePCBoffset(double t,double x,double y,double z,Offset &offset){
  int nx = round(x/hole_xs);
  int ny = round(y/hole_ys);
  
  offset.set = true;
  offset.x   = hole_xs * nx;
  offset.y   = hole_ys * ny;
  offset.z   = z+2.;
  offset.t   = t;
  
  if(debug) cout << "ComputePCBoffset:(" << offset.x << "," << offset.y << "," << offset.z << ")" << endl;
}

vector<int> Tools::GetNumberOfStrips(Boundaries det, Boundaries anode){
  int nv1 = (det.ymax-det.ymin)/(anode.strip_p[0]/cos(anode.strip_a[0]));
  int nv2 = (det.xmax-det.xmin)/(anode.strip_p[1]);
  int nv3 = (det.ymax-det.ymin)/(anode.strip_p[2]);
  
  cout << "GetNumberOfStrips::(" << nv1 << "," << nv2 << "," << nv3 << ")" << endl;
  
  vector<int> nv={nv1,nv2,nv3};
  return nv;
}

int Tools::GetStripID(int vw, Boundaries det, Boundaries anode, int nStrip, Electron el){
  int id = 0;
    
  if(vw == 0){
    double dy = det.ymax - el.y[0];
    double dx = dy / tan(anode.strip_a[vw]);
    double L  = fabs(el.x[0]) + dx;
        
    double x_pitch = anode.strip_p[vw]/sin(anode.strip_a[vw]);
    
    id = L / x_pitch;
  }
  else if(vw == 1) id = (el.x[0]-det.xmin)/anode.strip_p[vw];
  else if(vw == 2) id = (el.y[0]-det.ymin)/anode.strip_p[vw];

  return id;
}

void Tools::PrintHelp(){
  cout << "                                                     "            << endl;
  cout << "LArDrift commands help:"                                          << endl;
  cout << "  Usage: lardrift <opt> <value>"                                  << endl;
  cout << "  Options:"                                                       << endl;
  cout << "    -ef    input COMSOL E-field map"                              << endl;
  cout << "    -wv1   input COMSOL W-field map of induction view 1"          << endl;
  cout << "    -wv2   input COMSOL W-field map of induction view 2"          << endl;
  cout << "    -wv3   input COMSOL W-field map of collection view"           << endl;
  cout << "    -er    input electronic response function"                    << endl;
  cout << "    -o     output ROOT file"                                      << endl;
  cout << "    -data  data points from Lardon"                               << endl;
  cout << "    -ldf   data points from LArDrift output (pd-vd -> pcb)"       << endl;
  cout << "    -gf    data points from Geant4 output (g4 -> pd-vd -> pcb)"   << endl;
  cout << "    -setup experimental setup (colbox / pd-vd)"                   << endl;
  cout << "    -s     save output file"                                      << endl;
  cout << "    -d     debug mode (0=quiet / 1=minimal / >1=chatty"           << endl;
  cout << "    -n     number of simulated drifting electrons"                << endl;
  cout << "    -sc    scan direction (1=x - 2=y - 3=z - 4/5=spot cart/polar)"<< endl;
  cout << "    -sp    initial spot size (x0,y0,dx,dy) or (x0,y0,r,theta)    "<< endl;
  cout << "    -dm    drift model (1=Lardon / 2=Wolkowiak)"                  << endl;
  cout << "    -p     display plots"                                         << endl;
  cout << "    -dt    calculation timestep in us"                            << endl;
  cout << "    -do    add offset to Lardon data point"                       << endl;
  cout << "    -i     initial position of single simulation"                 << endl;
  cout << "    -l     simulate track along a line with '-n' hits"            << endl;
  cout << "                                                     "            << endl;
}
