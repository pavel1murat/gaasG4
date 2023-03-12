/////////////////////T//////////////////////////////////////////////////////////
// different sides - different inefficiencies
// fStop =  1: hit photodiode
// fStop =  2: absorption
// fStop =  3: exited through X side - no internal reflection
// fStop =  4: exited through Y side - no internal reflection
// fStop =  5: exited through the Z top side    - no internal reflection
// fStop =  6: exited through the Z bottom size - no internal reflection
// fStop =  9: too many reflections
// fStop = 10+iplane: exited through iplane - because the angle too small and no internal reflection
// fStop = 20+iplane: exited through iplane - diffuse reflection
//////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"

#include "gaasG4/inc/GaasSimAlg.hh"

// ClassImp(GaasSimAlg)
//-----------------------------------------------------------------------------
// a separate script is required to initialize geometry
//-----------------------------------------------------------------------------
GaasSimAlg::GaasSimAlg(): TNamed() {

					// N(e-h pairs) 
				// by default, fixed X,Y
  fPosMode      = 0;
//-----------------------------------------------------------------------------
// geometry is not initialized in the constructor,
// a separate call to InitGeometry is required
//-----------------------------------------------------------------------------
  fNPhMean         = 1000; // 13000  LYSO:Ce 26000 photons/MeV

  fPhotoEff        = 0.25;		// not used so far
  fMaxNReflections = 190;

  fRefrIndGaAs     = 3.4;
  fRefrIndAir      = 1.0;
  fRefrIndEpoxy    = 1.5;

  fAbsLength       = 2.2;    // 4 mm is too generous
  fReflProb        = 0.974;		// 0.95;

  for (int i=0; i<6; i++) fMirror[i] = 0;

  fRn = new TRandom3();
//-----------------------------------------------------------------------------
// folders
//-----------------------------------------------------------------------------
  TObject* o = gROOT->GetRootFolder()->FindObject("Ana");
  if (o != 0) {
    gROOT->GetRootFolder()->Remove(o);
    o->Delete();
  }

  fAnaFolder  = gROOT->GetRootFolder()->AddFolder("Ana","TStnAna Folder");
  fFolder     = fAnaFolder->AddFolder("GaasSimAlg","GaasSimAlg");
  fHistFolder = fFolder   ->AddFolder("Hist"          ,"ListOfHistograms");

  TH1::AddDirectory(0);

  for (int i=0; i<kNEventHistSets; i++) {
    fHist.fEvent[i] = 0;
  }

  for (int i=0; i<kNPhotonHistSets; i++) {
    fHist.fPhoton[i] = 0;
  }

  fDebugLevel = 0;
  fFirstCall  = 1;
}


//-----------------------------------------------------------------------------
GaasSimAlg::~GaasSimAlg() {
}

//-----------------------------------------------------------------------------
// although this part looks generic, only one detector is currently allowed
//-----------------------------------------------------------------------------
int GaasSimAlg::InitGeometry(double* Sensor, Detector_t* Det, int NDet) {
  
  fGeoManager = new TGeoManager("gm", "gm");
  TGeoElementTable *table = fGeoManager->GetElementTable();
//-----------------------------------------------------------------------------
// just for fun, follow ROOT rules
//-----------------------------------------------------------------------------
  TGeoMixture *air = new TGeoMixture("air",4, 0.00120479);
  air->AddElement(table->FindElement("C" ),0.000124);
  air->AddElement(table->FindElement("N" ),0.755268);
  air->AddElement(table->FindElement("O" ),0.231781);
  air->AddElement(table->FindElement("AR"),0.012827);
  
  TGeoMedium*  fMedAir = new TGeoMedium  ("AIR",1,air);
//-----------------------------------------------------------------------------
// describe the GaAs medium
//-----------------------------------------------------------------------------
  TGeoMixture *gaas = new TGeoMixture("gaas",2, 5.32);
  gaas->AddElement(table->FindElement ("GA"),1);
  gaas->AddElement(table->FindElement ("AS"),1);

  fMedGaAs = new TGeoMedium  ("GAAS",2,gaas);
//-----------------------------------------------------------------------------
// add photodiode medium
//-----------------------------------------------------------------------------
  TGeoMixture *ingaas = new TGeoMixture("ingaas",3, 5.32);
  ingaas->AddElement(table->FindElement ("In"),0.05);
  ingaas->AddElement(table->FindElement ("GA"),0.45);
  ingaas->AddElement(table->FindElement ("AS"),0.50);

  fMedInGaAs = new TGeoMedium  ("INGAAS",3,ingaas);
//-----------------------------------------------------------------------------
// now create volumes
//-----------------------------------------------------------------------------
  TGeoVolume   *top = fGeoManager->MakeBox("TOP",fMedAir,1,1,1);
  fGeoManager->SetTopVolume(top);
//-----------------------------------------------------------------------------
// this is the part which goes into the geometry initialization
//-----------------------------------------------------------------------------
  fCrystal   = fGeoManager->MakeBox("gaas_sensor",fMedGaAs,Sensor[0],Sensor[1],Sensor[2]);
  fCrystal->SetTransparency(60);
  top->AddNode(fCrystal,1);

  fNDetectors = NDet;
  
  char name[100];
  for (int i=0; i<NDet; i++) {
    fDetector[i] = new Detector_t(*(Det+i));
    sprintf(name,"Diode_%02i",i);
    TGeoVolume* diode = fGeoManager->MakeBox(name,fMedInGaAs,Det->fDx,Det->fDy,Det->fDz);
    diode->SetLineColor(kRed+2);
    diode->SetTransparency(20);

    if      (Det->fSide == 0) top->AddNode(diode,1,new TGeoTranslation(-Det->fDx-Sensor[0], Det->fOffsetY     , Det->fOffsetZ     ));
    else if (Det->fSide == 1) top->AddNode(diode,1,new TGeoTranslation( Det->fDx+Sensor[0], Det->fOffsetY     , Det->fOffsetZ     ));
    else if (Det->fSide == 2) top->AddNode(diode,1,new TGeoTranslation( Det->fOffsetX     ,-Det->fDy-Sensor[1], Det->fOffsetZ     ));
    else if (Det->fSide == 3) top->AddNode(diode,1,new TGeoTranslation( Det->fOffsetX     , Det->fDy+Sensor[1], Det->fOffsetZ     ));
    else if (Det->fSide == 4) top->AddNode(diode,1,new TGeoTranslation( Det->fOffsetX     , Det->fOffsetY     ,-Det->fDz-Sensor[2]));
    else if (Det->fSide == 5) top->AddNode(diode,1,new TGeoTranslation( Det->fOffsetX     , Det->fOffsetY     , Det->fDz+Sensor[2]));
  }
  
  fGeoManager->CloseGeometry();
//-----------------------------------------------------------------------------
//  retrieve the detector position from the geometry manager
//-----------------------------------------------------------------------------
  TGeoNode* tnode = fGeoManager->GetTopNode();

  fDiodePos = tnode->GetDaughter(1)->GetMatrix()->GetTranslation();

  // top->Draw();
  //  TView *view = gPad->GetView();
  // view->ShowAxis();
  return 0;
}

//_____________________________________________________________________________
void     GaasSimAlg::AddHistogram(TObject* hist, const char* FolderName) {
  TFolder* fol = (TFolder*) fFolder->FindObject(FolderName);
  fol->Add(hist); 
}

//_____________________________________________________________________________
void GaasSimAlg::HBook1F(TH1F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1F(Name,Title,Nx,XMin,XMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void GaasSimAlg::HBook2F(TH2F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 Int_t Ny, Double_t YMin, Double_t YMax,
			 const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH2F(Name,Title,Nx,XMin,XMax,Ny,YMin,YMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void GaasSimAlg::HProf(TProfile*& Hist, const char* Name, const char* Title,
		       Int_t Nx, Double_t XMin, Double_t XMax,
		       Double_t YMin, Double_t YMax,
		       const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TProfile(Name,Title,Nx,XMin,XMax,YMin,YMax);
  AddHistogram(Hist,FolderName);
}


//_____________________________________________________________________________
void GaasSimAlg::ResetHistograms(TFolder* Folder, const char* Opt) {
  // internal method...

  TObject  *o1;

  if (Folder == 0) Folder = fFolder;

  if (strcmp(Opt,"all") == 0) {

    TIter    it1(Folder->GetListOfFolders());

    while ((o1 = it1.Next())) {
      if (o1->InheritsFrom("TFolder")) {
	ResetHistograms((TFolder*) o1,Opt);
      }
      else if (o1->InheritsFrom("TH1")) {
	((TH1*) o1)->Reset();
      }
    }
  }
}


//-----------------------------------------------------------------------------
int GaasSimAlg::BookEventHistograms(EventHist_t* Hist, const char* Folder) {

  HBook1F(Hist->fNPhotons   ,"nphot" ,Form("%s: nphotons",Folder),1000,0,1.e5,Folder);
  HBook1F(Hist->fNDetPhotons,"ndphot",Form("%s: ndet photons",Folder),1000,0,1.e5,Folder);
  HBook1F(Hist->fEff[0]     ,"eff_0" ,Form("%s: eff[0]",Folder),10001,0,1.0001,Folder);
  HBook1F(Hist->fEff[1]     ,"eff_1" ,Form("%s: eff[1]",Folder),10001,0,0.10001,Folder);

  return 0;
}


//-----------------------------------------------------------------------------
int GaasSimAlg::BookPhotonHistograms(PhotonHist_t* Hist, const char* Folder) {

  int    nbx, nby;
  double dx, dy;

  TGeoBBox* box = (TGeoBBox*) fCrystal->GetShape();
  dx = box->GetDX();
  dy = box->GetDY();

  nbx = 100;
  nby = 100;

  HBook1F(Hist->fStop        ,"stop"  ,Form("%s: stop"        ,Folder),  50,0, 50,Folder);
  HBook1F(Hist->fNReflections,"nref"  ,Form("%s: nreflections",Folder), 500,0,500,Folder);
  HBook1F(Hist->fPath        ,"path"  ,Form("%s: path"        ,Folder),5000,0,  5,Folder);
  HBook2F(Hist->fYVsX        ,"y_vs_x",Form("%s: y vs x"      ,Folder),nbx,-dx,dx,nby,-dy,dy,Folder);

  return 0;
}


//-----------------------------------------------------------------------------
int GaasSimAlg::BookHistograms() {

  char folder_name[200];

  TFolder* hist_folder = (TFolder*) fFolder->FindObject("Hist");

  TH1::SetDefaultSumw2(kTRUE);
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  int  book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[0] = 1; // all 

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      TFolder* fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      if (fHist.fEvent[i] == 0) fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// photon histograms
//-----------------------------------------------------------------------------
  int  book_photon_histset[kNPhotonHistSets];
  for (int i=0; i<kNPhotonHistSets; i++) book_photon_histset[i] = 0;

  book_photon_histset[0] = 1;                      // all 
  book_photon_histset[1] = 1;                      // photons reached the photodiode
  book_photon_histset[2] = 1;                      // photons not reached the photodiode

  for (int i=0; i<kNPhotonHistSets; i++) {
    if (book_photon_histset[i] != 0) {
      sprintf(folder_name,"pho_%i",i);
      TFolder* fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      if (fHist.fPhoton[i] == 0) fHist.fPhoton[i] = new PhotonHist_t;
      BookPhotonHistograms(fHist.fPhoton[i],Form("Hist/%s",folder_name));
    }
  }
    
  return 0;
}

//-----------------------------------------------------------------------------
int GaasSimAlg::FillEventHistograms(EventHist_t* Hist) {
  Hist->fNPhotons->Fill(fNPhotons);
  Hist->fNDetPhotons->Fill(fNDetPhotons);

  double eff = fNDetPhotons/fNPhotons;

  Hist->fEff[0]->Fill(eff);
  Hist->fEff[1]->Fill(eff);

  return 0;
}

//-----------------------------------------------------------------------------
int GaasSimAlg::FillPhotonHistograms(PhotonHist_t* Hist) {

  double x    = fLastPoint.X();
  double y    = fLastPoint.Y();
  double path = fLastPoint.S();

  Hist->fStop->Fill(fStop);
  Hist->fNReflections->Fill(fNReflections);

  Hist->fPath->Fill(path);
  Hist->fYVsX->Fill(x,y);

  return 0;
}

//-----------------------------------------------------------------------------
int GaasSimAlg::FillHistograms(const char* Mode) {

  if (strcmp(Mode,"event") == 0) {
//-----------------------------------------------------------------------------
// fill event histograms
//-----------------------------------------------------------------------------
    FillEventHistograms(fHist.fEvent[0]);
  }
  else if (strcmp(Mode,"photon") == 0) {
//-----------------------------------------------------------------------------
// fill photon histograms
//-----------------------------------------------------------------------------
    FillPhotonHistograms(fHist.fPhoton[0]);

    if (fStop == 1) FillPhotonHistograms(fHist.fPhoton[1]);
    if (fStop != 1) FillPhotonHistograms(fHist.fPhoton[2]);
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
// planes 0,1 : low and high X; 2,3: low and high Y; 4,5: bottom and top (Z)
//-----------------------------------------------------------------------------
int GaasSimAlg::SimulateReflection(TTrajectoryPoint* Point, int IPlane) {
  int stop(0);

  int ipl = IPlane / 2;

  TVector3* p_dir = Point->GetDirection();

  double n0 = (*p_dir)[ipl];

  if (fMirror[IPlane] != 0) {
//-----------------------------------------------------------------------------
// the side is a mirror - re
//-----------------------------------------------------------------------------
    (*p_dir)[ipl] = -n0;
  }
  else {
    double sinalp = sqrt(1-n0*n0);
  
    if (sinalp < fRefrIndAir/fRefrIndGaAs) {
					// photon exits  - angle too close to 90 deg
      stop = 10+IPlane;
    }
    else {
//-----------------------------------------------------------------------------
// simulate reflection	
// a) check it a photon undergoes the geometric reflection
//-----------------------------------------------------------------------------
      double rn    = fRn->Rndm();
      if (rn < fReflProb) {
	(*p_dir)[ipl] = -n0;
      }
      else {
//-----------------------------------------------------------------------------
// diffuse scattering - in half of all cases the photon exits the crystal
//-----------------------------------------------------------------------------
	double dir[3];
	fRn->Sphere(dir[0],dir[1],dir[2],1);
	if (dir[ipl]*n0 > 0) {
					// photon exits after the diffuse scattering
	  stop = 20+IPlane;
	}
	else {
					// photon continues with new direction
	  p_dir->SetXYZ(dir[0],dir[1],dir[2]);
	}
      }
    }
  }
  return stop;
}


//-----------------------------------------------------------------------------
// stop = -1: photon got within acceptance of the fiber, but got reflected
// and didn't get out
int GaasSimAlg::SimulateDetector(TTrajectoryPoint* Point, Detector_t* Det) {
  int stop(0);
  
  double dx, dy, dz;
  
  dx = Point->X()-Det->fOffsetX;
  dy = Point->Y()-Det->fOffsetY;
  dz = Point->Z()-Det->fOffsetZ;

  TVector3* p_dir = Point->GetDirection();

  int    ipl    = Det->fSide / 2;
  double n0     = (*p_dir)[ipl];
  double sinalp = sqrt(1-n0*n0);
  
  
  if (fDebugLevel > 0) {
    printf("GaasSimAlg::SimulateDetector: Det->fType, Det->fSide, dx, dy, dz = %3i %3i %12.5f %12.5f %12.5f\n",
	   Det->fType,Det->fSide,dx,dy,dz);
  }

  if ((Det->fSide == 0) || (Det->fSide == 1)) {

    if (fDebugLevel > 0) {
      printf("GaasSimAlg::SimulateDetector: Det->fDy, Det->fDz = %12.5f %12.5f\n",Det->fDy,Det->fDz);
    }

    if ((fabs(dy) < Det->fDy) && (fabs(dz) < Det->fDz)) {
      if (Det->fType == 0) {
//-----------------------------------------------------------------------------
// photodiode - any photon reaching it stops
//-----------------------------------------------------------------------------
	stop = 1;
      }
      else {
//-----------------------------------------------------------------------------
// fiber - to get detected, a photon needs to actually get out
//-----------------------------------------------------------------------------
	if (sinalp < fRefrIndAir/fRefrIndGaAs) fStop = 1;
	else {
//-----------------------------------------------------------------------------
// diffuse scattering - in half of all cases the photon exits the crystal
//-----------------------------------------------------------------------------
	  double dir[3];
	  fRn->Sphere(dir[0],dir[1],dir[2],1);
	  if (dir[ipl]*n0 > 0) {
					// photon exits after diffuse scattering
	    stop = 1;
	  }
	  else {
					// photon continues with new direction
	    p_dir->SetXYZ(dir[0],dir[1],dir[2]);
	    stop = -1;
	  }
	}
      }
    }
  }
  else if ((Det->fSide == 2) || (Det->fSide == 3)) {
    if ((fabs(dx) < Det->fDx) && (fabs(dz) < Det->fDz)) {
      if (Det->fType == 0) stop = 1;
      else {
//-----------------------------------------------------------------------------
// fiber
//-----------------------------------------------------------------------------
	if (sinalp < fRefrIndAir/fRefrIndGaAs) fStop = 1;
	else {
//-----------------------------------------------------------------------------
// diffuse scattering - in half of all cases the photon exits the crystal
//-----------------------------------------------------------------------------
	  double dir[3];
	  fRn->Sphere(dir[0],dir[1],dir[2],1);
	  if (dir[ipl]*n0 > 0) {
	    stop = 1;
	  }
	  else {
					// photon continues with new direction
	    p_dir->SetXYZ(dir[0],dir[1],dir[2]);
	    stop = -1;
	  }
	}
      }
    }
  }
  else if ((Det->fSide == 4) || (Det->fSide == 5)) {
    if ((fabs(dx) < Det->fDx) && (fabs(dy) < Det->fDy)) {
      if (Det->fType == 0) stop = 1;
      else {
//-----------------------------------------------------------------------------
// fiber
//-----------------------------------------------------------------------------
	if (sinalp < fRefrIndAir/fRefrIndGaAs) fStop = 1;
	else {
//-----------------------------------------------------------------------------
// diffuse scattering - in half of all cases the photon exits the crystal
//-----------------------------------------------------------------------------
	  double dir[3];
	  fRn->Sphere(dir[0],dir[1],dir[2],1);
	  if (dir[ipl]*n0 > 0) {
	    stop = 1;
	  }
	  else {
					// photon continues with new direction
	    p_dir->SetXYZ(dir[0],dir[1],dir[2]);
	    stop = -1;
	  }
	}
      }
    }
  }

  if (fDebugLevel > 0) {
    printf("GaasSimAlg::SimulateDetector: END stop = %3i\n",stop);
  }

  return stop;
}

//-----------------------------------------------------------------------------
// return code:
// ------------
// fStop = 1: photon hits the photodiode
// fStop = 2: photon lost due to reflection inefficiency
//-----------------------------------------------------------------------------
int GaasSimAlg::TracePhoton(TTrajectoryPoint* Start, TGeoBBox* Vol) {

  double            x, y, z, nx, ny, nz, s, sx, sy, sz, smin, p;
  int               reflection_side;
  TTrajectoryPoint* tp;

  tp = &fLastPoint;

  tp->SetPoint(Start->X (),Start->Y (),Start->Z (),
	       Start->Nx(),Start->Ny(),Start->Nz(),
	       Start->S(),Start->GetPTotal());

  double xhsize = Vol->GetDX();
  double yhsize = Vol->GetDY();
  double zhsize = Vol->GetDZ();

  fStop = 0;

  fNReflections = 0;

  if (fDebugLevel > 0) {
    printf("GaasSimAlg::TracePhoton : Start :  fStop = %3i\n",fStop);
    tp->Print();
  }

  while (fStop == 0) {
					// find next relection
					// calculate path till the X-wall
    x  = tp->X ();
    y  = tp->Y ();
    z  = tp->Z ();
    nx = tp->Nx();
    ny = tp->Ny();
    nz = tp->Nz();

    smin = 0;

    double deltax, deltay, deltaz;
    
    sx = 1.e12;
    if (nx > 0) deltax = ( xhsize-x);
    else        deltax = (-xhsize-x);
    if (nx != 0) sx = deltax/nx;                     // always positive
					// calculate path till the Y-wall
    sy = 1.e12;
    if (ny > 0) deltay = ( yhsize-y);
    else        deltay = (-yhsize-y);
    if (ny != 0) sy = deltay/ny;
					// calculate path till the Z-wall
    sz = 1.e12;
    if (nz > 0) deltaz = ( zhsize-z);
    else        deltaz = (-zhsize-z);
    if (nz != 0) sz = deltaz/nz;
    
    if (fDebugLevel > 0) printf(" sx,sy,sz = %12.5e %12.5e %12.5e\n",sx,sy,sz);
//-----------------------------------------------------------------------------
// flag the side reached first
//-----------------------------------------------------------------------------
    if (sx < sy) {
      if (sx < sz) {
	smin = sx;
	if (nx < 0) reflection_side = 0;
	else        reflection_side = 1;
      }
      else         {
	smin = sz;
	if (nz < 0) reflection_side = 4;
	else        reflection_side = 5;
      }
    }
    else {
      if (sy < sz) {
	smin = sy;
	if (ny < 0) reflection_side = 2;
	else        reflection_side = 3;
      }
      else         {
	smin = sz;
	if (nz < 0) reflection_side = 4;
	else        reflection_side = 5;
      }
    }

    if (fDebugLevel > 0) printf(" reflection_side, smin = %3i %12.5e\n",reflection_side,smin);
//-----------------------------------------------------------------------------
// model absorption
//-----------------------------------------------------------------------------
    double abs_len(1.e12);
    
    if (fAbsLength > 0) {
      abs_len = fRn->Exp(fAbsLength);
      if (abs_len < smin) {
	fStop = 2;
	smin  = abs_len;
      }
    }

    if (fDebugLevel > 0) printf(" fStop, abs_len, smin = %3i %12.5e %12.5e\n",fStop, abs_len, smin);
//-----------------------------------------------------------------------------
// update coordinates
//-----------------------------------------------------------------------------
    x  = x+smin*nx;
    y  = y+smin*ny;
    z  = z+smin*nz;
    s  = tp->S()+smin;
    p  = tp->GetPTotal();
    
    tp->SetPoint(x,y,z,nx,ny,nz,s,p);

    if (fDebugLevel > 0) {
      printf("GaasSimAlg::TracePhoton : updated coordinates, fStop = %3i\n",fStop);
      tp->Print();
    }

    if (fStop != 0)                                         goto END_OF_PHOTON_SIM;
//-----------------------------------------------------------------------------
// photon didn't get absorbed, see if it ended up in one of detectors
// sides 0,1: X, 2,3: Y, 4,5: Z
//-----------------------------------------------------------------------------
    for (int i=0; i<fNDetectors; i++) {
      Detector_t* det = fDetector[i];
//-----------------------------------------------------------------------------
// check intersection with the detector, to simplify things,
// the side is stored in the detector 
//-----------------------------------------------------------------------------
      if (det->fSide == reflection_side) {
	fStop = SimulateDetector(tp,det);
      }
      if (fStop > 0)                                       goto END_OF_PHOTON_SIM;
    }
//-----------------------------------------------------------------------------
// the photon didn't get into any of photodetectors, simulate reflection
//-----------------------------------------------------------------------------
    if (fStop == 0) {
      fStop = SimulateReflection(tp,reflection_side);
    }
    else if (fStop < 0) {
					// photon reflected from the fiber
      fStop = 0;
    }
      
  END_OF_PHOTON_SIM:
    if (fStop == 0) fNReflections++;
//-----------------------------------------------------------------------------
// to avoid infinite loops, limite max number of reflections
//-----------------------------------------------------------------------------
    if (fNReflections >= fMaxNReflections) fStop = 9;
  }

  if (fDebugLevel > 0) {
    printf("GaasSimAlg::TracePhoton END: fStop,fNReflections: %3i %3i\n",fStop,fNReflections);
  }

  return 0;
}


//-----------------------------------------------------------------------------
int GaasSimAlg::BeginJob() {
  if (fFirstCall) {
    BookHistograms();
    fFirstCall = 0;
  }
  else {
    ResetHistograms(0,"all");
  }
  return 0;
}


//-----------------------------------------------------------------------------
int GaasSimAlg::simulate(int NPhotons) {
  
  float rn[3];
  double nx, ny, nz; // , x, y;


  for (int iph=0; iph<NPhotons; iph++) {
    if (fDebugLevel > 0) {
      printf(" >>>> GaasQuickSim::Run simulate photon number:  %10i\n",iph);
    }

    fRn->RndmArray(2,rn);

    double phi   = TMath::TwoPi()*rn[0];

    // double theta = TMath::ASin(2*(rn[1]-0.5));
    
    // nx  = TMath::Cos(theta)*TMath::Cos(phi);
    // ny  = TMath::Cos(theta)*TMath::Sin(phi);
    // nz  = TMath::Sin(theta);     
					// debugging
    nx  = TMath::Cos(phi);
    ny  = 0.;
    nz  = TMath::Sin(phi);
//-----------------------------------------------------------------------------
// trace photon till it reaches the photodiode of gets absorbed, reflecting
// it from the other faces
//-----------------------------------------------------------------------------
    fPoint.SetPoint(fX0,fY0,fZ0,nx,ny,nz,0,1);

    TracePhoton(&fPoint,(TGeoBBox*) fCrystal->GetShape());
    
    if (fStop == 1) {
//-----------------------------------------------------------------------------
// photon reached the photodiode
//-----------------------------------------------------------------------------
      fNDetPhotons += 1;
    }

    if (fDebugLevel > 0) {
      printf(">>>> GaasQuickSim::Run done simulating photon number:  %10i\n",iph);
    }
//-----------------------------------------------------------------------------
// end of photon tracing - fill "per photon" histograms
//-----------------------------------------------------------------------------
    FillHistograms("photon");
  }

  return 0;
}


//-----------------------------------------------------------------------------
int GaasSimAlg::Run(double Dist, int NEvents) {

  float  rn[3];  // to generate phi and theta

  BeginJob();
//-----------------------------------------------------------------------------
// simulate events, trace photons, and count the number of detected ones
//-----------------------------------------------------------------------------
  for (int ievent=0; ievent<NEvents; ievent++) {
    fEventNumber      = ievent;
    if (fDebugLevel > 0) {
      printf(" >>>> Start event number %10i\n",ievent);
    }
//-----------------------------------------------------------------------------
// reset profile histograms in the beginning
//-----------------------------------------------------------------------------
    ResetHistograms(0,"event");

    fRn->RndmArray(3,rn);

    TGeoBBox* cPos = (TGeoBBox*) fCrystal->GetShape();
    
    if (fPosMode == 0) {
//-----------------------------------------------------------------------------
// fixed distance in X from the photodiode center , Y=0, Z=0
//-----------------------------------------------------------------------------
      fX0 = fDiodePos[0]+Dist;
      fY0 = 0.; // 2*(rn[1]-0.5)*cPos->GetDY();
      fZ0 = 0.; // 2*(rn[2]-0.5)*cPos->GetDZ();
    }
    else if (fPosMode == 1) {
//-----------------------------------------------------------------------------
// randomly within the 50 um beam - it is an approximation
//-----------------------------------------------------------------------------      
      fX0 = fDiodePos[0]+Dist;
      fY0 = 2*(rn[1]-0.5)*0.0050;
      fZ0 = 2*(rn[2]-0.5)*cPos->GetDZ();
    }
      					// number of photons in this event

    if (fNPhMean > 0) fNPhotons    = fRn->Poisson(fNPhMean);
    else              fNPhotons    = -fNPhMean;
   
    fNDetPhotons = 0;

    if (fDebugLevel > 0) {
      printf(" >>>> GaasQuickSim::Run nphotons:  %12.0f\n",fNPhotons);
    }
//-----------------------------------------------------------------------------
// perform simulation
//-----------------------------------------------------------------------------
    simulate(fNPhotons);
//-----------------------------------------------------------------------------
// end of event processing - fill "per event" histograms
//-----------------------------------------------------------------------------
    FillHistograms("event");
  }

//-----------------------------------------------------------------------------
// on exit, print efficiency
//-----------------------------------------------------------------------------
  double eff0 = fHist.fEvent[0]->fEff[0]->GetMean();
  double eff1 = fHist.fEvent[0]->fEff[1]->GetMean();

  printf(" >> eff[0], eff[1] = %12.5e %12.5e\n",eff0, eff1);
  return 0;
}
  
//_____________________________________________________________________________
int  GaasSimAlg::SaveFolder(TFolder* Folder, TDirectory* Dir) {
  // save Folder into a subdirectory
  // do not write TStnModule's - for each TStnModule save contents of its
  // fFolder

  //  TFolder*     fol;
  TDirectory*  dir;
  TObject*     o;
//-----------------------------------------------------------------------------
// create new subdirectory in Dir to save Folder
//-----------------------------------------------------------------------------
  Dir->cd();
  //  dir = new TDirectory(Folder->GetName(),Folder->GetName(),"");
  dir = Dir->mkdir(Folder->GetName(),Folder->GetName());
  dir->cd();

//   printf(" ------------------- Dir: %s, new dir: %s\n",
// 	 Dir->GetName(),dir->GetName());


  TIter  it(Folder->GetListOfFolders());
  while ((o = it.Next())) {
//     printf(" o->GetName, o->ClassName : %-20s %-20s\n",
// 	   o->GetName(),
// 	   o->ClassName());

    if (strcmp(o->ClassName(),"TFolder") == 0) {
      SaveFolder((TFolder*) o, dir);
      //      dir->cd();
    }
    else if (! o->InheritsFrom("TStnModule")) {
      //      printf("gDirectory->GetPath = %s\n",gDirectory->GetPath());
      o->Write();
      //      gDirectory->GetListOfKeys()->Print();
    }
  }

  Dir->cd();
  return 0;
}

//_____________________________________________________________________________
void GaasSimAlg::SaveHist(const char* Filename) {
  // save histograms booked by all the modules into a file with the given name
  // Mode = 2: save directories

  TFile* f = new TFile(Filename,"recreate");

  SaveFolder(fAnaFolder,f);

  f->Close();
  delete f;
}


