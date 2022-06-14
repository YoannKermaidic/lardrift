// PCL includes
#include <pcl/point_cloud.h>
#include <pcl/octree/octree_search.h>

// System includes
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <thread>

// ROOT includes
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

// LArDrift includes
#include <reader.hh>
#include <plotter.hh>
#include <simarrays.hh>
#include <tools.hh>

using namespace std;

bool isw[4] = {false,false,false,false};  // "true" if weighting field recorded in COMSOL file
bool ise    = false;                      // "true" if electric field recorded in COMSOL file
bool iser   = false;                      // "true" if the electronic response file is provided
bool isg4   = false;                      // "true" if the LArSoft Geant4 file is provided
bool isn    = false;                      // "true" if noise properties given
bool isline = false;                      // "true" if simulated points taken from -l <line>

int debug = 0;                       // Enable debugging
int doplot= 0;                       // Enable plotting
int save  = 0;                       // Enable writing output to ROOT file
int Nsims = 1;                       // Number of simulated eletron
int scan  = 0;                       // Geometrical scan scheme of initial position (1=x, 2=y, 3=z, 4=xy)
int drift_model=1;                   // Electron drift speed model (1=LARDON - 2=Walkowiak)
int nViews= 3;                       // Number of view planes (ind - ind - coll)

double v_phi   = 0, v_theta = 10;    // Toy V vertex angles (if scan == 6)
double v_L1 = 20,v_L2 = 12,v_z0 = 10;//              length

double xli=-1000,yli=-1000,zli=-1000;// Starting line point
double xle=-1000,yle=-1000,zle=-1000;// Ending line point
double xdo     = 0, ydo = 0, zdo = 0;// Artificial data offset
double boxsize = 5;                  // cm     - Nearest point box
double dt_start= 40;                 // mus    - Increment simulation time
double q       = -1;                 // Sign of electron charge
double e_thrs  = 0.001;              // V/cm   - Minimum Efield allowed before stopping the drift
double dt_thrs = 0.1;                // us     - Minimum time increment allowed before stopping the drift
double T_lar   = 87;                 // K      - Liquid argon temperature
double lardon_offset[3] = {0,0,0};   //{-168.5,-150,-15};
double dEdS    = 2.1;                // MeV/cm - energy deposit
double dS      = 0.3;                // cm - discretization length
double spot[4] = {0,0,1,0};          // cm - initial spot size for random start: -sc 4 -> cart / -sc 5 -> polar

int nType      = 0;                  // Noise type (0 = white noise)
double nMean   = 0.;                 // Noise baseline offset
double nRMS    = 0.;                 // Noise amplitude

string exp_setup = "";       // Experimental setup ID for boundaries

// Read input 3D track file from LARDON or LArDrift or LArSoft
vector<Electron> inp;

// Store electron drift information
vector<Electron> els;
vector<bool> isvalid;

// Store strip signals if weighting potential is provided
vector<vector<Strips> > strips(nViews);
vector<int> nStrips;

// Multithreading management
std::vector<std::thread> threads;
size_t number_of_threads = 1;

int runsimulation(int tid, Boundaries bounds, Boundaries det, Tools tool, Offset offset){

  int sim_start = tid    * int(Nsims/number_of_threads);
  int sim_end   = (tid+1)* int(Nsims/number_of_threads);

  // Start looping over simulated electrons
  for(int sim=sim_start;sim<sim_end;sim++){

    // Some timing information on the drift calculation
    clock_t start_clock, end_clock;
    double cpu_time_total;
    start_clock = clock();
    
    // Create electrons
    Electron el; Electron el_corr;
    el.SetVerbosity(debug);
    
    // Setup the initial time/position of the electron
    if(isline) bounds.SetInitialPoint(sim,Nsims,xli,xle,yli,yle,zli,zle);
    else if(inp.size() > 0){
      //if(isg4) offset.z = -bounds.zmax*0.5;
      bounds.SetInitialPoint(inp[sim].t[0],inp[sim].x[0],inp[sim].y[0],inp[sim].z[0]);
    }
    if(bounds.xi==-1000)             bounds.SetInitialPoint(0,0,0,-10);
    if(scan==1)                      bounds.SetInitialPoint(0,bounds.xmin + (bounds.xmax-bounds.xmin)*(sim+0.5)/Nsims,bounds.yi,bounds.zi);
    else if(scan==2)                 bounds.SetInitialPoint(0,bounds.xi,bounds.ymin + (bounds.ymax-bounds.ymin)*(sim+0.5)/Nsims,bounds.zi);
    else if(scan==3)                 bounds.SetInitialPoint(0,bounds.xi,bounds.yi,bounds.zmin + (bounds.zmax-bounds.zmin)*(sim+0.5)/Nsims);
    else if(scan==4)                 bounds.SetInitialPoint(0,spot[0] + spot[2]*gRandom->Rndm(),spot[1] + spot[3]*gRandom->Rndm(),bounds.zi);
    else if(scan==5){
      double r     = spot[2]*gRandom->Rndm();
      double theta = 2*M_PI*gRandom->Rndm();
      bounds.SetInitialPoint(0,spot[0] + r*cos(theta),spot[1] + r*sin(theta),bounds.zi);
    }
    else if(scan==6) bounds.SetVertexPoint(sim,Nsims,v_phi*M_PI/180.,v_theta*M_PI/180.,v_L1,v_L2,v_z0);
    
    if(exp_setup=="pcb") tool.ComputePCBoffset(bounds.ti,bounds.xi,bounds.yi,bounds.zi,offset); // If signal sims, operates referential tranformation
    if(zli != zle) offset.z = 0;
    
    // Check if initial position lies within the detector boundaries
    if(!tool.InDetector(bounds.xi-offset.x,bounds.yi-offset.y,bounds.zi-offset.z)){
      cout << "LArDrift::WARNING:" << sim+1 << " / " << Nsims << " -> Initial point (" << bounds.xi-offset.x << " " << bounds.yi-offset.y << " " << bounds.zi-offset.z << ") not within boundaries" << endl;
      bounds.PrintBoundaries();
      continue;
      //return 1;
    }
    
    // Reset step time to standard value
    double dt = dt_start;
    
    // Create a new electron
    el.Clear();
    el.Resize(isw);
    
    // Initialize it with starting parameters
    el.step = 0;
    el.t[0] = bounds.ti-offset.t;
    el.x[0] = bounds.xi-offset.x; el.y[0] = bounds.yi-offset.y; el.z[0] = bounds.zi-offset.z;
    
    if(!isg4){
      el.n    = dEdS*dS/el.ion; if(scan==4 || scan==5) el.n /= Nsims;
      el.ds   = dS;
      el.qc   = el.n*el.qe;
    }
    else{
      el.n    = inp[sim].n;
      el.ds   = inp[sim].ds;
      el.qc   = el.n*el.qe;
    }
    
    // Store electron info with the input offset for plotting
    if(offset.set){el_corr.Clear(); el_corr.Resize(isw); el_corr.Copy(0,isw,el);}
    
    // Start the drift calculation
    cout << " " << endl;
    cout << "LArDrift::Start drifting the electron #" << sim+1 << " / " << Nsims << " from (" << el.x[0] << " " << el.y[0] << " " << el.z[0] << ")" << endl;
    
    int coll_flag   = 0;     // >0 if electron has reached the collection plane
    bool stop_drift = false; // true when electron reaches drift boundaries
    
    // Drift electron until it reaches boundaries
    while(1){
      if(debug>1) cout << " " << endl;
      // Get the local electric field and compute the corresponding velocity
      tool.SetEField(el,ise);
      el.v[0] = tool.EvalEspeed(el.e);
      el.SetSpeed();
      
      if(isw[0]){
        for(int view=1;view<4;view++){
          tool.SetWFields(view,el);
          el.SetCurrent(view);
        }
      }

      if(debug>1) tool.PrintStatus(el);
      
      el.step += 1;

      // Check if next electron position will reach detector boundaries. If yes, refine delta t step and continue the propagation
      if(!isw[0] && !tool.InDetector(el.x[el.step-1]+el.v[1]*dt,el.y[el.step-1]+el.v[2]*dt,el.z[el.step-1]+el.v[3]*dt)) dt /= 10.;
      
      // Move the electron in time by delta t
      el.Propagate(dt);
            
      // Check if electron has reached the collection plane - extrapolate on the surface.
      if((coll_flag=tool.IsCollected(el.x[el.step],el.y[el.step],el.z[el.step]))>0) el.ExtrapolateToSurface(coll_flag,isw,bounds);

      // Stop drift if the electron is collected
      if(coll_flag>0) stop_drift = true;
      // Stop drift if the electron is out of boundaries
      if(!tool.InDetector(el.x[el.step],el.y[el.step],el.z[el.step])) stop_drift = true;
      // Stop drift if electric field is not defined or is below threshold (prevent electron at rest) or other thresholds exceeded
      if(isnan(el.e) || el.e < e_thrs || dt < dt_thrs || el.step == el.maxstep){ el.step -= 1; stop_drift = true;}
      
      if(stop_drift){
        if(debug) cout << "LArDrift::End of the drift calculation after " << el.step << " steps -> z = " << el.z[el.step] << " cm, E = " << el.e << " V/cm, dt = " << dt << " us" << endl;
        // Store shifts to initial position/time before exiting
        el.SetFinalShifts(bounds,offset);
        el.coll_flag = coll_flag;
        break;
      }
    } // end of while loop
    
    // Simulation timing information
    end_clock = clock();
    double cpu_time = ((double) (end_clock - start_clock)) / CLOCKS_PER_SEC;
    cout << "LArDrift::Drift calculation time: " << cpu_time << " sec " << endl;
    
    // Simulation summary
    if(debug) tool.PrintEndPoint(el);
    
    // Bring back calculated positions into initial referential for storage
    if(offset.set){
      cout << "LArDrift::Shift the electron by (" << offset.x << " " << offset.y << " " << offset.z << ")" << endl;
      el_corr.Copy(isw,el);
      el.ApplyBackOffset(offset);
    }

//    if(isw[0]){if(coll_flag==2){isvalid[sim] = true; if(offset.set) els[sim] = el_corr; else els[sim] = el;}}
//    else{ isvalid[sim] = true; if(offset.set) els[sim] = el_corr; else els[sim] = el;}
    if(isw[0]){if(coll_flag==2){isvalid[sim] = true; els[sim] = el;}}
    else{ isvalid[sim] = true; els[sim] = el;}
  }
  
  return 0;
}

// Setup the thread sequence
void spawn(Boundaries bounds, Boundaries det, Tools tool, Offset offset) {
    for(int i=0; i<number_of_threads; ++i) threads.emplace_back(std::thread(runsimulation, i, bounds, det, tool, offset));
    for(int i=0; i<threads.size(); ++i)    threads[i].join();
}

int main (int argc, char** argv){
  // Default input parameters
  string geoFile     = "../config/geometry.json";// geometry file
  string efieldFile  = "./field.txt";   // COMSOL E-field file
  string detector    = "";              // Detector geometry ('pd-vd')
  string w1potFile   = "";              // COMSOL W-pot file
  string w2potFile   = "";              // COMSOL W-pot file
  string w3potFile   = "";              // COMSOL W-pot file
  string erFile      = "";              // Electronic response file
  string outFile     = "./output.root"; // ROOT     output file
  string dataFile    = "";              // LARDON   output file
  string lardriftFile= "";              // LARDRIFT output file
  string g4File      = "";              // GEANT4   output file
  
  cout << " " << endl;
  
  
  // LArDrift utilities
  Tools tool(debug, exp_setup,T_lar);
  Boundaries bounds;    // Detector/anode boundaries
  Boundaries det;       // Detector boundaries
  Reader reader(debug);
  Offset offset;
  
  // Read input arguments
  if(argc==1){tool.PrintHelp(); exit(0);}
  int i=0; while(++i<argc){
    const char *opt=argv[i];
    if(strcmp(opt,"--efield")==0       || strcmp(opt,"-ef")==0)   efieldFile  =argv[++i];
    else if(strcmp(opt,"--w1pot")==0   || strcmp(opt,"-det")==0){ detector    =argv[++i]; reader.GetBoundaries(det,geoFile,detector);}
    else if(strcmp(opt,"--w1pot")==0   || strcmp(opt,"-wv1")==0){ w1potFile   =argv[++i]; isw[1] = true;}
    else if(strcmp(opt,"--w2pot")==0   || strcmp(opt,"-wv2")==0){ w2potFile   =argv[++i]; isw[2] = true;}
    else if(strcmp(opt,"--w3pot")==0   || strcmp(opt,"-wv3")==0){ w3potFile   =argv[++i]; isw[3] = true;}
    else if(strcmp(opt,"--elec")==0    || strcmp(opt,"-er")==0) { erFile      =argv[++i]; iser   = true;}
    else if(strcmp(opt,"--config")==0  || strcmp(opt,"-cf")==0)   geoFile     =argv[++i];
    else if(strcmp(opt,"--output")==0  || strcmp(opt,"-o")==0)    outFile     =argv[++i];
    else if(strcmp(opt,"--datafile")==0|| strcmp(opt,"-data")==0){dataFile    =argv[++i]; inp = reader.GetData(dataFile,xdo, ydo, zdo);       Nsims = inp.size();}
    else if(strcmp(opt,"--simfile")==0 || strcmp(opt,"-ldf")==0) {lardriftFile=argv[++i]; inp = reader.GetLarDrift(lardriftFile);             Nsims = inp.size();}
    else if(strcmp(opt,"--g4file")==0  || strcmp(opt,"-gf")==0)  {g4File      =argv[++i]; inp = reader.GetGeant4(g4File,bounds,isw,dS,Nsims); Nsims = inp.size(); isg4 = true;}
    else if(strcmp(opt,"--expsetup")==0|| strcmp(opt,"-setup")==0){exp_setup=argv[++i];  reader.GetBoundaries(bounds,geoFile,exp_setup); tool.SetExpSetup(exp_setup);}
    else if(strcmp(opt,"--save")==0    || strcmp(opt,"-s")==0)    sscanf(argv[++i],"%i",&save);
    else if(strcmp(opt,"--debug")==0   || strcmp(opt,"-d")==0)   {sscanf(argv[++i],"%i",&debug); reader.SetVerbosity(debug); tool.SetVerbosity(debug);}
    else if(strcmp(opt,"--nsimu")==0   || strcmp(opt,"-n")==0)    sscanf(argv[++i],"%i",&Nsims);
    else if(strcmp(opt,"--scan")==0    || strcmp(opt,"-sc")==0)   sscanf(argv[++i],"%i",&scan);
    else if(strcmp(opt,"--plot")==0    || strcmp(opt,"-p")==0)    sscanf(argv[++i],"%i",&doplot);
    else if(strcmp(opt,"--dmodel")==0  || strcmp(opt,"-dm")==0)   sscanf(argv[++i],"%d",&drift_model);
    else if(strcmp(opt,"--spot")==0    || strcmp(opt,"-sp")==0)   sscanf(argv[++i],"%lf,%lf,%lf,%lf,%lf",&spot[0],&spot[1],&spot[2],&spot[3],&bounds.zi);
    else if(strcmp(opt,"--deltat")==0  || strcmp(opt,"-dt")==0)   sscanf(argv[++i],"%lf",&dt_start);
    else if(strcmp(opt,"--eloss")==0   || strcmp(opt,"-deds")==0) sscanf(argv[++i],"%lf",&dEdS);
    else if(strcmp(opt,"--strip")==0   || strcmp(opt,"-ds")==0)   sscanf(argv[++i],"%lf",&dS);
    else if(strcmp(opt,"--dataoffset")==0|| strcmp(opt,"-do")==0){sscanf(argv[++i],"%lf,%lf,%lf",&offset.x,&offset.y,&offset.z); offset.set = true;}
    else if(strcmp(opt,"--line")==0    || strcmp(opt,"-l")==0)   {sscanf(argv[++i],"%lf,%lf,%lf,%lf,%lf,%lf",&xli,&xle,&yli,&yle,&zli,&zle); isline = true;}
    else if(strcmp(opt,"--vertex")==0  || strcmp(opt,"-vt")==0)   sscanf(argv[++i],"%lf,%lf,%lf,%lf,%lf",&v_phi,&v_theta,&v_L1,&v_L2,&v_z0);
    else if(strcmp(opt,"--noise")==0   || strcmp(opt,"-noise")==0)sscanf(argv[++i],"%d,%lf,%lf",&nType,&nMean,&nRMS);
    else if(strcmp(opt,"--nthread")==0 || strcmp(opt,"-nt")==0)   sscanf(argv[++i],"%zd",&number_of_threads);
    else if(strcmp(opt,"--help")==0    || strcmp(opt,"-h")==0){   tool.PrintHelp(); exit(0);}
    else if(strcmp(opt,"--init")==0    || strcmp(opt,"-i")==0)    {
      if(exp_setup==""){cout<<"LArDrift::Must first provide an experimental setup. Exit."<<endl; exit(1);}
      sscanf(argv[++i],"%lf,%lf,%lf",&bounds.xi,&bounds.yi,&bounds.zi);}
    else {
      cout << "LArDrift::Unknown passed argument: " << argv[i] << endl;
      tool.PrintHelp();
      exit(1);
    }
  }
  
  if(!bounds.set){
    cerr << "LArDrift::Must provide an experimental setup '-setup <setup>'. Exit." << endl;
    exit(1);
  }

  if(number_of_threads>1) cout << "LArDrift::Multithreading enabled with " << number_of_threads << " cores" << endl;
  
  // Read COMSOL fields calculation output
            reader.GetComsol(efieldFile,tool.clouds[0],tool.fields[0],ise,isw);
  if(isw[1]) reader.GetComsol(w1potFile,tool.clouds[1],tool.fields[1],ise,isw);
  if(isw[2]) reader.GetComsol(w2potFile,tool.clouds[2],tool.fields[2],ise,isw);
  if(isw[3]) reader.GetComsol(w3potFile,tool.clouds[3],tool.fields[3],ise,isw);
  
  tool.GetBoundaries(bounds);
  if(debug) bounds.PrintBoundaries();
  
  tool.SetResolution();
  tool.SetEOctree();
  if(isw[0]) tool.SetWOctrees();
    
  if(inp.size() > 0){
    if(dataFile != "")          cout << "LArDrift::Will simulate according to LARDON data point"         << endl;
    else if(lardriftFile != "") cout << "LArDrift::Will simulate according to LArDrift simulation point" << endl;
    else if(isg4)               cout << "LArDrift::Will simulate according to LArSoft simulation point"  << endl;
  }
  
  if(isw[0]){
    if(!det.set){
      cerr << "LArDrift::Must provide a detector setup '-det <detector>' to initialize and compute the number of strips. Exit." << endl;
      exit(1);
    }
    // Get number of strips and initialize them all
    nStrips = tool.GetNumberOfStrips(det,bounds);
    cout << "LArDrift::Initialize strips" << endl;
    for(int vw=0;vw<strips.size();vw++){
      for(int st=0;st<nStrips[vw];st++){
        Strips strip(debug);
        strips[vw].push_back(strip);
        strips[vw][st].Resize();
        // Set strip position
        if(!strips[vw][st].set_pos) strips[vw][st].SetPosition(det,bounds,vw,st);
        if(debug>2) cout << "  -> view: " << vw << " strip: " << st << " (" << strips[vw][st].id << "," << strips[vw][st].x << "," << strips[vw][st].y << ")" << endl;
      }
    }
  }
  
  // Some timing information on the drift calculation
  clock_t start_global_clock, end_global_clock;
  start_global_clock = clock();
  
  cout << "LArDrift::Will perform " << Nsims << " simulations" << endl;
  els.resize(Nsims);
  isvalid.resize(Nsims,false);
   
  // Start multi-threading
  // see https://stackoverflow.com/questions/15716903/pthread-passing-object-as-argument-to-pthread-create
  
  threads.reserve(number_of_threads);
  spawn(bounds,det,tool,offset);

  // Exit if no valid electron recorded
  bool exit_flag = true;
  for(int sim=0;sim<Nsims;sim++) if(isvalid[sim]) {exit_flag = false; break;}
  if(exit_flag) {cout << "LArDrift::No valid electron recorded. Strange termination." << endl; return EXIT_SUCCESS;}
  
  if(isw[0]){
    cout << "LArDrift::Fill in strips with electron induced signal" << endl;
    int drift_offset = tool.GetDriftTimeOffset(dt_start,isvalid,els);

    for(int sim=0;sim<Nsims;sim++){
      if(!isvalid[sim]) continue;
    
      for(int vw=0;vw<nViews;vw++){
        // Get the strip ID based on the electron position
        int st_id = tool.GetStripID(vw,det,bounds,nStrips[vw],els[sim]);
        cout << "  -> view: " << vw << " electron position (" << els[sim].x[0] << "," << els[sim].y[0] << ") strip ID: " << st_id << endl;
        
        if(st_id >= strips[vw].size()){
          cerr << "LArDrift::Error: strip ID " << st_id << " larger than strips size: " << strips[vw].size() << " for view: " << vw+1 << ". Skipped." << endl;
          continue;}

        // Set signal time
        if(!strips[vw][st_id].set_time) strips[vw][st_id].SetTime(els[sim]);
        // Add signal
        strips[vw][st_id].AddSignal(drift_offset,dt_start,vw,els[sim]);
        if(debug>2) strips[vw][st_id].DumpSignal();
        // Add noise
        if(isn){
          Noise noise(debug); noise.Clear(); noise.Resize();
          if(nType==0) noise.SetWhiteNoise(nMean,nRMS);
          strips[vw][st_id].AddNoise(noise);
        }
      }
    }
  }
  
  end_global_clock = clock();
  double cpu_global_time    = ((double) (end_global_clock - start_global_clock)) / double(CLOCKS_PER_SEC) / double(number_of_threads);
  cout << "LArDrift::Global time: " << cpu_global_time << " sec " << endl;
    
  ///
  // Enter plot area
  ///

  TApplication* rootapp;
  if(doplot) rootapp = new TApplication("rootapp",0,0);

  // Plot simulated drift path
  vector<TCanvas*> canvases;
  Plotter plot(debug,bounds,dt_start);
  if(det.set) plot.SetBoundaries(det);
  // Plot electron drift speed model
  canvases.push_back(plot.SpeedModel(tool));
  // (xy) plane position shift plot
  canvases.push_back(plot.OneDshift(isvalid,els));
  // drift time 1D plot
  canvases.push_back(plot.OneDdrift(isvalid,els));
    // Plot displacement maps
  if(scan==4) canvases.push_back(plot.TwoDshift(isvalid,els));
  // Plot electron 2D tracks
  canvases.push_back(plot.TwoDTracks(isvalid,els));
  // Plot electron 3D tracks
  canvases.push_back(plot.ThreeDTracks(isvalid,els));
  // Plot (xy) projections
  if(xli>-1000 || inp.size()>0){
    canvases.push_back(plot.TwoDProj(isvalid,els,inp));
    canvases.push_back(plot.ThreeDProj(isvalid,els,inp));
  }
  
  // Plot event display map
  if(isw[0]){
    vector<vector<double> > vER;
    if(iser){vER = reader.GetER(erFile,dt_start);  canvases.push_back(plot.ER(vER));}
    for(int vw=0;vw<nViews;vw++) for(int st=0;st<strips[vw].size();st++) tool.ERconvolution(vER,strips[vw][st]);
    canvases.push_back(plot.FullEventDisplay("Raw",strips));
    canvases.push_back(plot.FullEventDisplay("Conv",strips));
  }
  
  // Write output to ROOT file
  if(save){
    TFile* rootfile = new TFile(outFile.c_str(),"RECREATE");
    rootfile->cd();
    for(int can=0;can<canvases.size();can++) canvases[can]->Write();
    
    // Output ROOT file settings
    Electron el;
    TTree* roottree;
    roottree = new TTree("tree","electrons tracks information");
    roottree->Branch("n", &el.step);
    roottree->Branch("t", &el.t);
    roottree->Branch("x", &el.x);
    roottree->Branch("y", &el.y);
    roottree->Branch("z", &el.z);

    for(int sim=0;sim<Nsims;sim++){
      if(!isvalid[sim]) continue;
      el.Clear();
      el.Resize(isw);
      el.Copy(isw,els[sim]);
      roottree->Fill();
    }
    roottree->Write();
    rootfile->Close();
    cout << "LArDrift::Output stored in " << outFile << " ROOT file" << endl;
  }
  
  cout << "LArDrift::Normal termination" << endl;
  
  if(doplot)    rootapp->Run();

  return EXIT_SUCCESS;
}

