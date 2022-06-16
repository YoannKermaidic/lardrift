#include "reader.hh"

Reader::Reader(int db){
  debug    = db;
  mmTOcm   = 0.1;                // mm to cm conversion
  mTOcm    = 100;                // m  to cm conversion
  mVToADC  = 2.048;              // mV to ADC TDE conversion
  degTOrad = M_PI/180.;          // rad/deg
}

vector<Electron> Reader::GetData(string filename,double xdo,double ydo,double zdo){
  vector<Electron> els;
  
  string header, x, y, z;
  ifstream File(filename.c_str());

  int numSkip = std::count(std::istreambuf_iterator<char>(File),
                           std::istreambuf_iterator<char>(), '%');
  
  File.clear();
  File.seekg(0);
  
  int numEntries = std::count(std::istreambuf_iterator<char>(File),
                              std::istreambuf_iterator<char>(), '\n');
  
  File.clear();
  File.seekg(0);
  
  if(debug) cout << "GetData::Header:  " << numSkip    << " lines to be skipped" << endl;
  if(debug) cout << "GetData::Entries: " << numEntries << " lines to be readout" << endl;

  if(File){
    cout << "GetData::" << filename << " opened ..." << endl;
    for(int i=0;i<numSkip;i++) getline(File,header);
    
    els.resize(numEntries);
    
    cout << "GetData::Track with " << numEntries << " hits found" << endl;
    if(debug) cout << "GetData::Hits location:" << endl;
    for(std::size_t i = 0; i < numEntries; ++i){
      File >> x >> y >> z;

      Electron data;
      data.step = 1;
      data.x.resize(1);
      data.y.resize(1);
      data.z.resize(1);

      data.x[0] = atof(x.c_str()) + xdo;
      data.y[0] = atof(y.c_str()) + ydo;
      data.z[0] = atof(z.c_str()) + zdo;
      
      if(debug) cout << "  " << i << " " << data.x[0] << " " << data.y[0] << " " << data.z[0] << endl;
      els[i] = data;
    }
  }
  else{
    cerr << "GetData::" << filename << " not found. Exit." << std::endl;
    exit(1);
  }
  File.close();  // Close file
  cout << "GetData::closed. " << endl;

  return els;
}

vector<vector<double> > Reader::GetER(string filename,double dt){
  if(debug) std::cout << "GetER::Start -> " << filename << std::endl;

  vector<vector<double> > data(2);
  string buf, t, x;
  ifstream File(filename.c_str());
  
  if(File){
    cout << "GetER::" << filename << " opened"<< endl;
    
    int numEntries = std::count(std::istreambuf_iterator<char>(File),
                                std::istreambuf_iterator<char>(), '\n');
    
    File.clear();
    File.seekg(0);
  
    getline(File,buf);
    File >> t >> x;
    
    double response = 0;
    double ERdt = atof(t.c_str());
    int iter = 0;
    int ds = dt / ERdt; // Down sampling rate
    
    cout << "GetER::Down-sampling rate: " << ds << endl;
    
    File.clear();
    File.seekg(0);

    data[0].push_back(0);
    data[1].push_back(0);

    for(std::size_t i = 0; i < numEntries; ++i){
      File >> t >> x;
      if(iter<ds){response += atof(x.c_str()); iter++;}
      else{
        data[0].push_back(atof(t.c_str()));
        data[1].push_back(mVToADC * response/double(ds));
        response = atof(x.c_str());
        iter = 1;
      }
    }
  }
  else{
    cerr << "GetER::" << filename << " not found. Exit." << std::endl;
    exit(1);
  }
  
  File.close();  // Close file
  return data;
}

void Reader::GetComsol(string filename, pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, Field &field,bool &ise,bool* isw){
  if(debug) std::cout << "GetComsol::Start -> " << filename << std::endl;
  
  cloud = pcl::PointCloud<pcl::PointXYZ>::Ptr (new pcl::PointCloud<pcl::PointXYZ>);

  bool isEfield = true;
  if(filename.find("wfield") != string::npos) isEfield = false;
  
  double lenght_unit = 1., e_unit = 1;
  bool set_unit = false;
  string flag_unit = "Length unit:";
  char unit[2];
  string header, x, y, z ,v, e, ex, ey, ez;
  ifstream File(filename.c_str());
  
  int numSkip = std::count(std::istreambuf_iterator<char>(File),
                           std::istreambuf_iterator<char>(), '%');
  
  File.clear();
  File.seekg(0);
  
  int numEntries = std::count(std::istreambuf_iterator<char>(File),
                              std::istreambuf_iterator<char>(), '\n');
  
  File.clear();
  File.seekg(0);
  
  if(debug) cout << "GetComsol::Header:  " << numSkip    << " lines to be skipped" << endl;
  if(debug) cout << "GetComsol::Entries: " << numEntries << " lines to be readout" << endl;
  
  // Generate pointcloud data
  cloud->width = numEntries;
  cloud->height = 1;
  cloud->points.resize (cloud->width * cloud->height);
  field.pot.resize(numEntries);
  
  if(File){
    cout << "GetComsol::" << filename << " opened"<< endl;
    for(int i=0;i<numSkip;i++){
      getline(File,header);
      // get length unit
      if(!set_unit && header.find(flag_unit) != string::npos){
        set_unit = true;
        sscanf(header.c_str(),"%% Length unit:        %s",unit);
        if(!strcmp(unit,"m"))       lenght_unit = mTOcm;
        else if(!strcmp(unit,"cm")) lenght_unit = 1.;
        else if(!strcmp(unit,"mm")) lenght_unit = mmTOcm;
        if(debug) cout << "GetComsol::Found length unit: " << unit << " -> conversion factor to cm: " << lenght_unit << endl;
      }
    }
    if(isEfield  && header.find("es.normE") != string::npos) { if(debug) cout << "ReadCOMSOL::Header (electric) -> " << header << endl; ise    = true;}
    if(!isEfield && header.find("es.normE") != string::npos) { if(debug) cout << "ReadCOMSOL::Header (weighting)-> " << header << endl; isw[0] = true;}
    if((isEfield && ise) || (!isEfield && isw[0])){
      cout << "GetComsol::Field found." << endl;
      field.e.resize(numEntries);
      field.ex.resize(numEntries);
      field.ey.resize(numEntries);
      field.ez.resize(numEntries);
      // get electric field unit
      if(header.find("(V/m)") != std::string::npos)       e_unit = mTOcm;
      else if(header.find("(V/cm)") != std::string::npos) e_unit = 1.;
      else if(header.find("(V/mm)") != std::string::npos) e_unit = mmTOcm;
      if(debug) cout << "GetComsol::Found Field unit: -> conversion factor to V/cm: " << 1./e_unit << endl;
    }
    
    for(std::size_t i = 0; i < cloud->size (); ++i){
      File >> x >> y >> z >> v;
      (*cloud)[i].x = lenght_unit*atof(x.c_str());
      (*cloud)[i].y = lenght_unit*atof(y.c_str());
      (*cloud)[i].z = lenght_unit*atof(z.c_str());
      field.pot[i]  = atof(v.c_str());
      if((isEfield && ise) || (!isEfield && isw[0])){
        File >> e >> ex >> ey >> ez;
        field.e[i]  = atof(e.c_str())/e_unit;
        field.ex[i] = atof(ex.c_str())/e_unit;
        field.ey[i] = atof(ey.c_str())/e_unit;
        field.ez[i] = atof(ez.c_str())/e_unit;
      }
    }
  }
  else{
    cerr << "GetComsol::" << filename << " not found. Exit." << std::endl;
    exit(1);
  }
  File.close();  // Close file
  cout << "GetComsol::" << filename << " closed"<< endl;
  if(debug) cout << "GetComsol::Done" << endl;
}

void Reader::GetBoundaries(Boundaries &bdr, string filename, string exp_setup){
  nlohmann::json j;
  
  // read channel map JSON file
  ifstream File(filename.c_str());
  if(File) cout << "GetBoundaries::Geometry settings opened ... ";
  else { cout << "GetBoundaries::Geometry settings not found. Exit. " << endl; exit(1);}
    
  File >> j;
  File.close();
  cout << " closed. " << endl;
  
  if(debug>1) cout << j.dump(2) << endl;
  
  if(!j.contains(exp_setup.c_str())){
    cout << "GetBoundaries::" << exp_setup << " not found." << endl; j.dump(2);
    cout << "GetBoundaries::Exit." << endl; exit(1);
  }
  if(!j[exp_setup.c_str()].contains("boundaries") || !j[exp_setup.c_str()].contains("center")){
    cout << "GetBoundaries:: 'boundaries' or 'center' information not found." << endl; j.dump(2);
    cout << "GetBoundaries::Exit." << endl; exit(1);
  }
  
  bdr.xmin = j[exp_setup.c_str()]["boundaries"]["x"][0];
  bdr.xmax = j[exp_setup.c_str()]["boundaries"]["x"][1];
  bdr.ymin = j[exp_setup.c_str()]["boundaries"]["y"][0];
  bdr.ymax = j[exp_setup.c_str()]["boundaries"]["y"][1];
  bdr.zmin = j[exp_setup.c_str()]["boundaries"]["z"][0];
  bdr.zmax = j[exp_setup.c_str()]["boundaries"]["z"][1];
  
  if(j[exp_setup.c_str()]["boundaries"].contains("hole")){
    bdr.hole_r = j[exp_setup.c_str()]["boundaries"]["hole"]["radius"];
    bdr.hole_h = j[exp_setup.c_str()]["boundaries"]["hole"]["height"];
    bdr.hole_xs= j[exp_setup.c_str()]["boundaries"]["hole"]["x_spacing"];
    bdr.hole_ys= j[exp_setup.c_str()]["boundaries"]["hole"]["y_spacing"];
    bdr.hole_v = j[exp_setup.c_str()]["boundaries"]["hole"]["v_spacing"];
  }
  if(j[exp_setup.c_str()]["boundaries"].contains("strips")){
    bdr.strip_p.resize(3); bdr.strip_a.resize(3);
    for(int i=0;i<3;i++){
      bdr.strip_p[i] = j[exp_setup.c_str()]["boundaries"]["strips"]["pitches"][i];
      bdr.strip_a[i] = j[exp_setup.c_str()]["boundaries"]["strips"]["angles"][i];
      bdr.strip_a[i] *= degTOrad;
    }
  }
  
  bdr.xc   = j[exp_setup.c_str()]["center"][0];
  bdr.yc   = j[exp_setup.c_str()]["center"][1];
  bdr.zc   = j[exp_setup.c_str()]["center"][2];
  
  if(j[exp_setup.c_str()].contains("exclude")){
    if(j[exp_setup.c_str()]["exclude"].contains("z")){
      bdr.zexcmin = j[exp_setup.c_str()]["exclude"]["z"][0];
      bdr.zexcmax = j[exp_setup.c_str()]["exclude"]["z"][1];
    }
  }
  
  bdr.set = true;
}

vector<Electron> Reader::GetLarDrift(string filename){
  vector<Electron> els;
  
  TFile* rootfile = new TFile(filename.c_str(),"READ");
  if(rootfile){
    cout << "GetLarDrift::LArDrift ROOT file opened ... ";
    
    int v_n             = 0;
    vector<double> *v_t = 0;
    vector<double> *v_x = 0;
    vector<double> *v_y = 0;
    vector<double> *v_z = 0;
    
    TTree* roottree = (TTree*)rootfile->Get("tree");
    roottree->SetBranchAddress("n", &v_n);
    roottree->SetBranchAddress("t", &v_t);
    roottree->SetBranchAddress("x", &v_x);
    roottree->SetBranchAddress("y", &v_y);
    roottree->SetBranchAddress("z", &v_z);
    
    cout << "entries: " << roottree->GetEntries();

    els.resize(roottree->GetEntries());
    
    for(int i=0;i<roottree->GetEntries();i++){
      roottree->GetEntry(i);
      
      Electron el;
      el.step = 1;
      el.t.push_back(v_t->at(v_n-1));
      el.x.push_back(v_x->at(v_n-1));
      el.y.push_back(v_y->at(v_n-1));
      el.z.push_back(v_z->at(v_n-1));
      
      els[i] = el;
    }
    
    rootfile->Close();
    cout << " closed. " << endl;
  }
  else { cout << "GetLarDrift::LArDrift ROOT file not found. Exit. " << endl; exit(1);}
  
  return els;
}

vector<Electron> Reader::GetGeant4(string filename, Boundaries bdr, bool *isw, double &ds, int nEvt){
  vector<Electron> els;
  
  if(debug) std::cout << "GetGeant4::Start -> " << filename << std::endl;
  
  int run, subrun, event, numHits, hID, tID, nPhotons, nSlow, eln;
  double hEnergy, elt, elx, ely, elz, elstep;
  string line;
  char buffer[256];
  ifstream File(filename.c_str());
  
  if(File){
    cout << "GetGeant4::" << filename << " opened"<< endl;

//    int numEvents = std::count(std::istreambuf_iterator<char>(File),
//                               std::istreambuf_iterator<char>(), 'Event');

    int numEvents = 0;
    
    while ( std::getline(File, line) )
      if(sscanf(line.c_str(),"Event run: %d %*s",&run) == 1) numEvents++;
    
    File.clear();
    File.seekg(0);
  
    cout << "GetGeant4::Number of events:  " << numEvents    << " to process" << endl;
    
    if(nEvt>numEvents) {cerr << "GetGeant4::Number of events:  " << numEvents    << " smaller than requested event: " << nEvt << ". Exit." << endl; exit(1);}
    
    for(int evt=nEvt-1;evt<nEvt;evt++){
      getline(File,line);
      sscanf(line.c_str(),"Event run: %d subRun: %d event: %d contains %d %*s",&run,&subrun,&event,&numHits);
    
      printf("   Event run: %d subRun: %d event: %d contains %d hits \n",run,subrun,event,numHits);
      els.resize(numHits);
      
      int numClusters = 0;
      int nsteps = 0;
      double steps = 0;
      
      Electron el;

      for(int hit=0;hit<numHits;hit++){

        if(nsteps==0){el.Clear(); el.Resize(1,isw); el.n = 0; el.t[0] = 0; el.x[0] = 0; el.y[0] = 0; el.z[0] = 0;}
        
        getline(File,line);
        if(debug) cout << "   " << line << endl;
        sscanf(line.c_str(),"[#%d]  TrkID=%d (mu-): %lf MeV on %lf ns at (%lf,%lf,%lf) (step: %lf cm); electrons: %d; photons: %d (fast), %d (slow)",
               &hID,&tID,&hEnergy,&elt,&elz,&elx,&ely,&elstep,&eln,&nPhotons,&nSlow);

        nsteps++;
        if(ds==0) ds = elstep;
        steps   += elstep;
        
        el.step = 1;
        el.n    += eln;
        el.t[0] += elt;
        el.x[0] += elx;
        el.y[0] += ely + bdr.ymin;
        el.z[0] += elz;
        
        if(steps >= ds){
          el.t[0] /= double(nsteps);
          el.x[0] /= double(nsteps);
          el.y[0] /= double(nsteps);
          el.z[0] /= double(nsteps);
          
          els[numClusters] = el;
          numClusters++;

          steps  = 0;
          nsteps = 0;
        }
        
        //if(tID>1) break; // simulate only muons
        
      }
      els.resize(numClusters);
      cout << "   Recorded " << numClusters << " clusters" << endl;

    }
    cout << " closed. " << endl;
  }
  else { cout << "GetGeant4::LArSoft LOG file not found. Exit. " << endl; exit(1);}
   
  return els;
}
