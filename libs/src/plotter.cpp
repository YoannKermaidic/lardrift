#include "plotter.hh"

Plotter::~Plotter(){}

Plotter::Plotter(int db, Boundaries bdr, double dt){
  debug = db;
  
  // Plane for plotting projection
  axaxis = {"x","y","z"};
  planes = {"xy","zy","xz"};
  pxaxis = {"x (cm)","z (cm)","x (cm)"};
  pyaxis = {"y (cm)","y (cm)","z (cm)"};
  
  xmin = bdr.xmin; xmax = bdr.xmax;
  ymin = bdr.ymin; ymax = bdr.ymax;
  zmin = bdr.zmin; zmax = bdr.zmax;
  dt_start = dt;
}

void Plotter::SetBoundaries(Boundaries bdr){
  xmin = bdr.xmin; xmax = bdr.xmax;
  ymin = bdr.ymin; ymax = bdr.ymax;
  zmin = bdr.zmin; zmax = bdr.zmax;
}

TCanvas* Plotter::TwoDshift(vector<bool> isvalid, vector<Electron> els){
  TCanvas* c2DS = new TCanvas("c2DS","Electron drift (xy) displacement",700,500);
  c2DS->Divide(3,1);
  
  vector<string> zaxis = {"x displacement (cm)","y displacement (cm)","r displacement (cm)"};
  vector<TGraph2D*> gS2D(3);
  for(int i=0;i<gS2D.size();i++) gS2D[i] = new TGraph2D();
  
  for(int sim=0;sim<els.size();sim++){
    if(!isvalid[sim]) continue;
    gS2D[0]->SetPoint(gS2D[0]->GetN(),els[sim].x[0],els[sim].y[0],els[sim].sx);
    gS2D[1]->SetPoint(gS2D[1]->GetN(),els[sim].x[0],els[sim].y[0],els[sim].sy);
    gS2D[2]->SetPoint(gS2D[2]->GetN(),els[sim].x[0],els[sim].y[0],TMath::Sqrt(pow(els[sim].sx,2)+pow(els[sim].sy,2)));
  }
  
  for(int i=0;i<gS2D.size();i++){
    c2DS->cd(i+1);
    string title = zaxis[i] + ";x (cm);y (cm)";
    gS2D[i]->SetTitle(title.c_str());
    gS2D[i]->GetZaxis()->SetTitle(zaxis[i].c_str());
    gS2D[i]->Draw("A COLZ");
  }
  
  return c2DS;
}

TCanvas* Plotter::OneDshift(vector<bool> isvalid, vector<Electron> els){
  
  vector<string> sdir  = {"x-shift","y-shift","z-shift"};
  vector<vector<TGraph*> > gS1D(axaxis.size());
  vector<TH1D*> hAxis(axaxis.size());
  vector<TLegend*> leg(axaxis.size());
  Color_t colors[3] = {kBlue+1,kRed+1,kGreen+1};
  
  TCanvas* cShift = new TCanvas("cShift","Electron drift displacement",700,500);
  cShift->Divide(axaxis.size(),1);
  
  for(int ax=0;ax<axaxis.size();ax++){
    cShift->cd(ax+1);
    double axmin = (ax == 0) ? xmin : ((ax == 1) ? ymin : zmin);
    double axmax = (ax == 0) ? xmax : ((ax == 1) ? ymax : zmax);
    double aymin = 1e6, aymax = -1e6;
    
    hAxis[ax] = new TH1D(Form("hAxis_%s",axaxis[ax].c_str()),"",100,axmin,axmax);
    gS1D[ax].resize(3);
    for(int i=0;i<gS1D[ax].size();i++){
      gS1D[ax][i] = new TGraph();
      gS1D[ax][i]->SetName(Form("gS1D_%s_%s",axaxis[ax].c_str(),sdir[i].c_str()));
      gS1D[ax][i]->SetMarkerSize(0.5);        gS1D[ax][i]->SetMarkerStyle(20);
      gS1D[ax][i]->SetMarkerColor(colors[i]); gS1D[ax][i]->SetLineColor(colors[i]);
    }
    
    for(int sim=0;sim<els.size();sim++){
      if(!isvalid[sim]) continue;
      if(aymin > els[sim].sx) aymin = els[sim].sx;
      if(aymax < els[sim].sx) aymax = els[sim].sx;
      if(aymin > els[sim].sy) aymin = els[sim].sy;
      if(aymax < els[sim].sy) aymax = els[sim].sy;
      if(aymin > els[sim].sz) aymin = els[sim].sz;
      if(aymax < els[sim].sz) aymax = els[sim].sz;
      
      if(ax==0){
        gS1D[ax][0]->SetPoint(gS1D[ax][0]->GetN(),els[sim].x[0],els[sim].sx);
        gS1D[ax][1]->SetPoint(gS1D[ax][1]->GetN(),els[sim].x[0],els[sim].sy);
        gS1D[ax][2]->SetPoint(gS1D[ax][2]->GetN(),els[sim].x[0],els[sim].sz);
      }
      else if(ax==1){
        gS1D[ax][0]->SetPoint(gS1D[ax][0]->GetN(),els[sim].y[0],els[sim].sx);
        gS1D[ax][1]->SetPoint(gS1D[ax][1]->GetN(),els[sim].y[0],els[sim].sy);
        gS1D[ax][2]->SetPoint(gS1D[ax][2]->GetN(),els[sim].y[0],els[sim].sz);
      }
      else if(ax==2){
        gS1D[ax][0]->SetPoint(gS1D[ax][0]->GetN(),els[sim].z[0],els[sim].sx);
        gS1D[ax][1]->SetPoint(gS1D[ax][1]->GetN(),els[sim].z[0],els[sim].sy);
        gS1D[ax][2]->SetPoint(gS1D[ax][2]->GetN(),els[sim].z[0],els[sim].sz);
      }
    }
    
    hAxis[ax]->GetXaxis()->SetTitle(Form("%s (cm)",axaxis[ax].c_str()));
    hAxis[ax]->GetYaxis()->SetTitle("Displacement (cm)");
    hAxis[ax]->GetXaxis()->SetNdivisions(409); hAxis[ax]->GetYaxis()->SetNdivisions(409);
    hAxis[ax]->SetStats(0);
    hAxis[ax]->SetMinimum(0.9*aymin); hAxis[ax]->SetMaximum(1.1*aymax);
    hAxis[ax]->Draw("axis");

    leg[ax] = new TLegend(0.2,0.6,0.4,0.89);
    leg[ax]->SetLineWidth(0); leg[ax]->SetHeader("Along axis:");
    leg[ax]->SetTextSize(0.05);
    
    for(int i=0;i<gS1D.size();i++){
      leg[ax]->AddEntry(gS1D[ax][i],sdir[i].c_str(),"LP");
      gS1D[ax][i]->Draw("P same");
    }
    leg[ax]->Draw();
  }
  
  return cShift;
}

vector<TGraph*> Plotter::Box(){
  vector<TGraph*> gBox(planes.size());

  for(int it = 0; it<planes.size();it++){
    gBox[it] = new TGraph();
    gBox[it]->SetName(Form("gBox_%s",planes[it].c_str()));
    gBox[it]->GetXaxis()->SetTitle(pxaxis[it].c_str());
    gBox[it]->GetYaxis()->SetTitle(pyaxis[it].c_str());
    gBox[it]->SetFillColor(kCyan);  gBox[it]->SetFillStyle(3002);
    gBox[it]->SetLineColor(kBlack); gBox[it]->SetLineStyle(7);
    if(planes[it] == "xz"){
      gBox[it]->SetPoint(gBox[it]->GetN(),xmin,zmin);
      gBox[it]->SetPoint(gBox[it]->GetN(),xmin,zmax);
      gBox[it]->SetPoint(gBox[it]->GetN(),xmax,zmax);
      gBox[it]->SetPoint(gBox[it]->GetN(),xmax,zmin);
      gBox[it]->SetPoint(gBox[it]->GetN(),xmin,zmin);
    }
    else if(planes[it] == "zy"){
      gBox[it]->SetPoint(gBox[it]->GetN(),zmin,ymin);
      gBox[it]->SetPoint(gBox[it]->GetN(),zmin,ymax);
      gBox[it]->SetPoint(gBox[it]->GetN(),zmax,ymax);
      gBox[it]->SetPoint(gBox[it]->GetN(),zmax,ymin);
      gBox[it]->SetPoint(gBox[it]->GetN(),zmin,ymin);
    }
    else{
      gBox[it]->SetPoint(gBox[it]->GetN(),xmin,ymin);
      gBox[it]->SetPoint(gBox[it]->GetN(),xmin,ymax);
      gBox[it]->SetPoint(gBox[it]->GetN(),xmax,ymax);
      gBox[it]->SetPoint(gBox[it]->GetN(),xmax,ymin);
      gBox[it]->SetPoint(gBox[it]->GetN(),xmin,ymin);
    }
  }
  return gBox;
}

TCanvas* Plotter::OneDdrift(vector<bool> isvalid, vector<Electron> els){
  
  vector<TGraph*> gDT(axaxis.size());
  vector<TH1D*> hDTAxis(axaxis.size());
  
  TCanvas* cDT = new TCanvas("cDT","Electron drift time",700,500);
  cDT->Divide(2,2);
  
  for(int ax=0;ax<axaxis.size();ax++){
    cDT->cd(ax+1);
    double axmin = (ax == 0) ? xmin : ((ax == 1) ? ymin : zmin);
    double axmax = (ax == 0) ? xmax : ((ax == 1) ? ymax : zmax);
    double aymin = 1e6, aymax = -1e6;
    
    hDTAxis[ax] = new TH1D(Form("hDTAxis_%s",axaxis[ax].c_str()),"",100,axmin,axmax);
    gDT[ax] = new TGraph(); gDT[ax]->SetName(Form("gDT_%s",axaxis[ax].c_str()));
    
    for(int sim=0;sim<els.size();sim++){
      if(!isvalid[sim]) continue;
      if(aymin > els[sim].st) aymin = els[sim].st;
      if(aymax < els[sim].st) aymax = els[sim].st;
      if(ax == 0) gDT[ax]->SetPoint(gDT[ax]->GetN(),els[sim].x[0],els[sim].st);
      if(ax == 1) gDT[ax]->SetPoint(gDT[ax]->GetN(),els[sim].y[0],els[sim].st);
      if(ax == 2) gDT[ax]->SetPoint(gDT[ax]->GetN(),els[sim].z[0],els[sim].st);
    }
    
    hDTAxis[ax]->GetXaxis()->SetRangeUser(axmin,axmax);
    hDTAxis[ax]->SetMinimum(0.9*aymin); hDTAxis[ax]->SetMaximum(1.1*aymax);

    hDTAxis[ax]->GetXaxis()->SetTitle(Form("%s (cm)",axaxis[ax].c_str()));
    hDTAxis[ax]->GetYaxis()->SetTitle("Drift time (us)");
    hDTAxis[ax]->GetXaxis()->SetNdivisions(409); hDTAxis[ax]->GetYaxis()->SetNdivisions(409);
    hDTAxis[ax]->SetStats(0);
    hDTAxis[ax]->Draw("axis");
    
    gDT[ax]->SetMarkerSize(0.5);     gDT[ax]->SetMarkerStyle(20);
    gDT[ax]->SetMarkerColor(kBlack); gDT[ax]->SetLineColor(kBlack);
    gDT[ax]->Draw("P same");
  }

  return cDT;
}

TCanvas* Plotter::ER(vector<vector<double> > vER){
  TCanvas* cER = new TCanvas("cER","Electronic response",700,500);
  cER->cd();
  
  TGraph* gER = new TGraph();
  for(int i=0;i<vER[0].size();i++) gER->SetPoint(gER->GetN(),vER[0][i],vER[1][i]);
  gER->SetName("gER");
  gER->GetXaxis()->SetTitle("Time (#ms)");
  gER->GetYaxis()->SetTitle("Electronic response (ADC/fC)");
  gER->Draw("ALP");
  
  return cER;
}

TCanvas* Plotter::TwoDTracks(vector<bool> isvalid, vector<Electron> els){
  
  TCanvas* cDrift = new TCanvas("cDrift","Electron drift line",700,500);
  cDrift->Divide(2,2);
  
  gStyle->SetPalette(57); // kBird ROOT palette
  
  vector<TGraph*> gBox = Box();             // Draw box defined by boundaries
  vector<vector<TGraph*> > gDrift(planes.size());
  
  for(int pl=0;pl<planes.size();pl++){
    cDrift->cd(pl+1);
    gBox[pl]->Draw("AFL");
    gDrift[pl].resize(els.size());
    for(int sim=0;sim<els.size();sim++){
      if(!isvalid[sim]) continue;
      gDrift[pl][sim] = new TGraph();
      gDrift[pl][sim]->SetName(Form("gTrack_%s_%i",planes[pl].c_str(),sim));
      gDrift[pl][sim]->SetMarkerStyle(20);     gDrift[pl][sim]->SetMarkerSize(0.3);
      
      int icol = sim;
      if(icol>254) icol -= 255;
      gDrift[pl][sim]->SetMarkerColor(gStyle->GetColorPalette(icol)); gDrift[pl][sim]->SetLineColor(gStyle->GetColorPalette(icol));

      for(std::size_t i = 0; i < els[sim].step; ++i){
        if(planes[pl] == "xz")      gDrift[pl][sim]->SetPoint(gDrift[pl][sim]->GetN(),els[sim].x[i],els[sim].z[i]);
        else if(planes[pl] == "zy") gDrift[pl][sim]->SetPoint(gDrift[pl][sim]->GetN(),els[sim].z[i],els[sim].y[i]);
        else                        gDrift[pl][sim]->SetPoint(gDrift[pl][sim]->GetN(),els[sim].x[i],els[sim].y[i]);
      }
      gDrift[pl][sim]->Draw("LP same");
    }
  }
  return cDrift;
}

TCanvas* Plotter::TwoDProj(vector<bool> isvalid, vector<Electron> els, vector<Electron> inp){
  
  TCanvas* cProj = new TCanvas("cProj","Projected electron drift line",700,500);
  cProj->Divide(2,2);
  
  vector<TGraph*> gBox = Box();             // Draw box defined by boundaries
  TLegend* leg = new TLegend(0.2,0.2,0.8,0.8);
  leg->SetLineWidth(0); leg->SetTextSize(0.06); leg->SetFillStyle(0);
  
  for(int pl=0;pl<planes.size();pl++){
    cProj->cd(pl+1);
    gBox[pl]->Draw("AL");

    TGraph* gInit = new TGraph(); gInit->SetName("gInit");
    TGraph* gEnd  = new TGraph(); gEnd->SetName("gEnd");
  
    gEnd->SetMarkerSize(0.5);     gEnd->SetMarkerStyle(20);
    gEnd->SetLineColor(kRed);     gEnd->SetMarkerColor(kRed);
  
    gInit->SetMarkerSize(0.5);    gInit->SetMarkerStyle(20);
    gInit->SetLineColor(kBlue);   gInit->SetMarkerColor(kBlue);

    for(int sim=0;sim<els.size();sim++){
      if(!isvalid[sim]) continue;
      if(els[sim].x[els[sim].step] == 9999) continue;
      if(planes[pl] == "xz"){
        gInit->SetPoint(gInit->GetN(),els[sim].x[0],els[sim].z[0]);
        gEnd->SetPoint(gEnd->GetN(),els[sim].x[els[sim].step],els[sim].z[0]);
      }
      else if(planes[pl] == "zy"){
        gInit->SetPoint(gInit->GetN(),els[sim].z[0],els[sim].y[0]);
        gEnd->SetPoint(gEnd->GetN(),els[sim].z[0],els[sim].y[els[sim].step]);
      }
      else{
        gInit->SetPoint(gInit->GetN(),els[sim].x[0],els[sim].y[0]);
        gEnd->SetPoint(gEnd->GetN(),els[sim].x[els[sim].step],els[sim].y[els[sim].step]);
      }
    }
  
    gBox[pl]->Draw("ALF"); gInit->Draw("PL same"); gEnd->Draw("PL same");
    
    TGraph* gData;
    if(inp.size() > 0){
      gData = new TGraph();         gData->SetName("gData");
      gData->SetMarkerSize(1.3);    gData->SetMarkerStyle(24);
      gData->SetLineColor(kGray+2); gData->SetMarkerColor(kGray+2);
    
      for(int pts=0;pts<inp.size();pts++){
        if(planes[pl] == "xz")      gData->SetPoint(gData->GetN(),inp[pts].x[0],inp[pts].z[0]);
        else if(planes[pl] == "zy") gData->SetPoint(gData->GetN(),inp[pts].z[0],inp[pts].y[0]);
        else                        gData->SetPoint(gData->GetN(),inp[pts].x[0],inp[pts].y[0]);
      }
      gData->Draw("PL same");
    }
    if(pl==0){
      leg->AddEntry(gBox[0],"Drift volume","LF");
      leg->AddEntry(gInit,"Initial position","PL");
      leg->AddEntry(gEnd, "End position","PL");
      if(inp.size() > 0) leg->AddEntry(gData, "Data hits","PL");
    }
  }
  cProj->cd(4);
  leg->Draw();

  return cProj;
}

TCanvas* Plotter::ThreeDProj(vector<bool> isvalid, vector<Electron> els, vector<Electron> inp){
  
  TCanvas* c3DProj = new TCanvas("c3DProj","3D electron drift line",700,500);
  c3DProj->cd();
  
  TView *view1 = TView::CreateView(1);
  view1->SetRange(xmin,ymin,zmin,xmax,ymax,zmax);
  
  const int Nbox = 16;
  double xarr[Nbox], yarr[Nbox], zarr[Nbox]; int iarr = 0;
  xarr[iarr] = xmin; yarr[iarr] = ymin; zarr[iarr] = zmin; iarr++;
  xarr[iarr] = xmin; yarr[iarr] = ymin; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymin; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymin; zarr[iarr] = zmin; iarr++;
  xarr[iarr] = xmin; yarr[iarr] = ymin; zarr[iarr] = zmin; iarr++;
  
  xarr[iarr] = xmin; yarr[iarr] = ymax; zarr[iarr] = zmin; iarr++;
  xarr[iarr] = xmin; yarr[iarr] = ymax; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymax; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymax; zarr[iarr] = zmin; iarr++;
  xarr[iarr] = xmin; yarr[iarr] = ymax; zarr[iarr] = zmin; iarr++;
  
  xarr[iarr] = xmin; yarr[iarr] = ymax; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmin; yarr[iarr] = ymin; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymin; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymax; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymax; zarr[iarr] = zmin; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymin; zarr[iarr] = zmin; iarr++;
  
  TPolyLine3D* gBox    = new TPolyLine3D(Nbox, xarr, yarr, zarr); // Draw box defined by boundaries
  TPolyLine3D* lInit   = new TPolyLine3D(els.size());   TPolyLine3D* lEnd    = new TPolyLine3D(els.size());
  TPolyMarker3D* pInit = new TPolyMarker3D(els.size()); TPolyMarker3D* pEnd  = new TPolyMarker3D(els.size());
  
  for(int sim=0;sim<els.size();sim++){
    if(!isvalid[sim]) continue;
    if(els[sim].x[els[sim].step] == 9999) continue;
    lInit->SetPoint(sim,els[sim].x[0],els[sim].y[0],els[sim].z[0]);
    lEnd->SetPoint(sim,els[sim].x[els[sim].step],els[sim].y[els[sim].step],els[sim].z[0]);
    
    pInit->SetPoint(sim,els[sim].x[0],els[sim].y[0],els[sim].z[0]);
    pEnd->SetPoint(sim,els[sim].x[els[sim].step],els[sim].y[els[sim].step],els[sim].z[0]);
  }
  
  pEnd->SetMarkerSize(0.5);   pEnd->SetMarkerStyle(20);
  pEnd->SetMarkerColor(kRed);
  lEnd->SetLineColor(kRed);
  
  pInit->SetMarkerSize(0.5);  pInit->SetMarkerStyle(20);
  pInit->SetMarkerColor(kBlue);
  lInit->SetLineColor(kBlue);
  
  gBox->Draw();
  lInit->Draw(); lEnd->Draw();
  pInit->Draw(); pEnd->Draw();
  
  TPolyMarker3D* pData; TPolyLine3D* lData;
  if(inp.size() > 0){
    pData = new TPolyMarker3D(inp.size());
    lData = new TPolyLine3D(inp.size());
    
    for(int pts=0;pts<inp.size();pts++){
      pData->SetPoint(pts,inp[pts].x[0],inp[pts].y[0],inp[pts].z[0]);
      lData->SetPoint(pts,inp[pts].x[0],inp[pts].y[0],inp[pts].z[0]);
    }
    pData->SetMarkerSize(1.3);   pData->SetMarkerStyle(24);
    pData->SetMarkerColor(kGray+2);
    lData->SetLineColor(kGray+2);
    
    lData->Draw(); pData->Draw();
  }
  
  TLegend* leg = new TLegend(0.05,0.85,0.25,0.99);
  leg->SetLineWidth(0); leg->SetFillStyle(0);
  leg->AddEntry(lInit,"Initial position","L");
  leg->AddEntry(lEnd, "End position","L");
  if(inp.size() > 0) leg->AddEntry(lData, "Data hits","L");
  leg->Draw();

  return c3DProj;
}

TCanvas* Plotter::ThreeDTracks(vector<bool> isvalid, vector<Electron> els){
  
  TCanvas* c3DTracks = new TCanvas("c3DTracks","3D electron tracks",700,500);
  c3DTracks->cd();
  
  gStyle->SetPalette(57); // kBird ROOT palette

  TView *view1 = TView::CreateView(1);
  view1->SetRange(xmin,ymin,zmin,xmax,ymax,zmax);
  
  const int Nbox = 16;
  double xarr[Nbox], yarr[Nbox], zarr[Nbox]; int iarr = 0;
  xarr[iarr] = xmin; yarr[iarr] = ymin; zarr[iarr] = zmin; iarr++;
  xarr[iarr] = xmin; yarr[iarr] = ymin; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymin; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymin; zarr[iarr] = zmin; iarr++;
  xarr[iarr] = xmin; yarr[iarr] = ymin; zarr[iarr] = zmin; iarr++;
  
  xarr[iarr] = xmin; yarr[iarr] = ymax; zarr[iarr] = zmin; iarr++;
  xarr[iarr] = xmin; yarr[iarr] = ymax; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymax; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymax; zarr[iarr] = zmin; iarr++;
  xarr[iarr] = xmin; yarr[iarr] = ymax; zarr[iarr] = zmin; iarr++;
  
  xarr[iarr] = xmin; yarr[iarr] = ymax; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmin; yarr[iarr] = ymin; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymin; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymax; zarr[iarr] = zmax; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymax; zarr[iarr] = zmin; iarr++;
  xarr[iarr] = xmax; yarr[iarr] = ymin; zarr[iarr] = zmin; iarr++;
  
  TPolyLine3D* gBox    = new TPolyLine3D(Nbox, xarr, yarr, zarr); // Draw box defined by boundaries
  vector<TPolyLine3D*> lTracks(els.size());

  gBox->Draw();

  for(int sim=0;sim<els.size();sim++){
    if(!isvalid[sim]) continue;
    lTracks[sim] = new TPolyLine3D(els[sim].step);
    for(std::size_t i = 0; i < els[sim].step; ++i)
      lTracks[sim]->SetPoint(i,els[sim].x[i],els[sim].y[i],els[sim].z[i]);

    int icol = sim;
    if(icol>254) icol -= 255;
    lTracks[sim]->SetLineColor(gStyle->GetColorPalette(icol));

    lTracks[sim]->Draw();
  }
  


  return c3DTracks;
}

TCanvas* Plotter::SpeedModel(TF1* model){
  TCanvas* cSpeed = new TCanvas("cSpeed","LARDON electron speed model",700,500);
  cSpeed->cd();
  
  int dE = 100;
  TGraph* gSpeed = new TGraph();
  for(int i=0;i<1000/dE;i++)
    gSpeed->SetPoint(gSpeed->GetN(),i*dE,model->Eval(i*dE));
  
  gSpeed->GetXaxis()->SetTitle("Electric field (V/cm)");
  gSpeed->GetYaxis()->SetTitle("Electron speed (cm/us)");
  
  gSpeed->SetMarkerStyle(20); gSpeed->SetMarkerSize(0.5);
  gSpeed->SetMarkerColor(kBlack);
  model->SetLineColor(kRed); model->SetLineWidth(2);
  
  gSpeed->Draw("AP");
  model->Draw("L same");
  
  return cSpeed;
}

TCanvas* Plotter::SpeedModel(Tools tool){
  TCanvas* cSpeed = new TCanvas("cSpeed","LARDON electron speed model",700,500);
  cSpeed->cd();
    
  int dE = 10;
  TGraph* gSpeed = new TGraph();
  if(debug>1 || debug == -1) cout << "GetSpeedModel:: E (V/cm) \t V (cm/us)" << endl;
  for(int i=0;i<1000/dE;i++){
    gSpeed->SetPoint(gSpeed->GetN(),i*dE,tool.EvalEspeed(i*dE));
    if(debug>1 || debug == -1) cout << i*dE << "\t" << tool.EvalEspeed(i*dE) << endl;
  }
  
  gSpeed->GetXaxis()->SetTitle("Electric field (V/cm)");
  gSpeed->GetYaxis()->SetTitle("Electron speed (cm/us)");
  
  gSpeed->SetMarkerStyle(20); gSpeed->SetMarkerSize(0.5);
  gSpeed->SetMarkerColor(kBlack);
  gSpeed->SetLineColor(kRed); gSpeed->SetLineWidth(2);
  gSpeed->Draw("ALP");
  
  return cSpeed;
}

TCanvas* Plotter::EventDisplay(string type, int vw, vector<Strips> strips){
  TH2D* hED = new TH2D(Form("hED_%s_%d",type.c_str(),vw),Form("View %d event display - %s current;Strip ID;Time (0.1 #ms tick)",vw,type.c_str()),strips.size(),0,strips.size(),strips[0].maxstep,0,strips[0].maxstep);
  hED->SetStats(0);
  
  for(int st=0;st<strips.size();st++)
    for(int t=0;t<strips[st].maxstep;t++)
      hED->Fill(st,t,strips[st].cur[t]);
  
  TCanvas* cED = new TCanvas(Form("cED_%s_%d",type.c_str(),vw),Form("View %d event display - %s current",vw,type.c_str()),700,500);
  cED->cd();
  hED->Draw("colz");
  
  return cED;
}

TCanvas* Plotter::FullEventDisplay(string type, vector<vector<Strips> > strips){
  
  TCanvas* cFED = new TCanvas(Form("cED_%s",type.c_str()),Form("Full event display - %s current",type.c_str()),700,500);
  cFED->Divide(3,1);
  
  for(int vw=0;vw<strips.size();vw++){
    cFED->cd(vw+1);

    TH2D* hED = new TH2D(Form("hED_%s_%d",type.c_str(),vw),Form("View %d event display - %s current;Strip ID;Time (0.1 #ms tick)",vw,type.c_str()),strips[vw].size(),0,strips[vw].size(),strips[vw][0].maxstep,0,strips[vw][0].maxstep);
    hED->SetStats(0);
    
    for(int st=0;st<strips[vw].size();st++)
      for(int t=0;t<strips[vw][st].maxstep;t++){
        if(type=="Raw")       hED->Fill(st,t,strips[vw][st].cur[t]);
        else if(type=="Conv") hED->Fill(st,t,strips[vw][st].curc[t]);
      }
  
    hED->Draw("colz");
  }
  
  return cFED;
}
