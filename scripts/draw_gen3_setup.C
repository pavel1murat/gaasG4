///////////////////////////////////////////////////////////////////////////////
// draw different parts of Mu2e
//
// first example: 
//
// dg = new DrawGen3Setup("gen3_setup.gdml")
// dg->gm->GetVolume("BoxInTheWorld")->Draw("ogl")
//
// comment: TGeoManager::Import chokes on filenames like "~/mu2e.gdml") 
///////////////////////////////////////////////////////////////////////////////

#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TString.h"

class DrawGen3Setup {
public:
  TGeoManager* gm;

  TGeoNode*    fTop; 

  TGeoNode*    fSourceHolder; 
  TGeoNode*    fSensor; 
  TGeoNode*    fPhotodiode[7]; 
  
  int          fTransp;

  DrawGen3Setup(const char* Fn = "gen3_setup.gdml", int OriginalColors = 0);
  ~DrawGen3Setup();
  
  void SetRecursiveVisibility(TGeoVolume* Vol, int OnOff);
  void SetRecursiveVisibility(TGeoNode*   Vol, int OnOff);

  void SetDefaultColorTransp            ();
  void SetRecursiveColorTransp          (TGeoVolume* Vol , Int_t Color, Int_t Transp);
  void SetRecursiveColorTranspByName    (TGeoNode*   Vol , const char* Name   , Int_t Color, Int_t Transp);
  void SetRecursiveColorTranspByMaterial(TGeoNode*   Node, const char* MatName, Int_t Color, Int_t Transp);

  void SetRecursiveVisibilityColorTranspByNameAndMaterial(TGeoNode*   Top         ,
							  const char* Name        ,
							  const char* MatName     ,
							  int         Visibility  ,
							  int         Color       ,
							  int         Transparency);
    
  void SetRecursiveVisibilityByName    (TGeoNode* Node, const char* NamePattern, int OnOff);
  void SetRecursiveVisibilityByMaterial(TGeoNode* Node, const char* Material   , int OnOff);

				// Mu2e-specific - Node name starts with 'Pattern'
				// assume it is unique
  
  TGeoNode* FindNodeByName      (TGeoNode*   Top, const char* Name      );
  TGeoNode* FindNodeByVolumeName(TGeoNode*   Top, const char* VolumeName);
  TGeoNode* FindNodeByVolumeName(TGeoVolume* Top, const char* VolumeName);
    
  void Draw();
};


//-----------------------------------------------------------------------------
DrawGen3Setup::DrawGen3Setup(const char* Fn, int KeepOriginalColors) {
  gm = new TGeoManager();
  gm->Import(Fn);

  fTop       = gGeoManager->GetTopNode();

  fSourceHolder = FindNodeByVolumeName(fTop,"SourceHolder1");
  fSensor       = FindNodeByVolumeName(fTop,"GaasSensor");

  for (int i=0; i<7; i++) {
    fPhotodiode[i]   = FindNodeByVolumeName(fTop,Form("InGaAsPD%i",i));
  }

  fTransp = 40;

  SetDefaultColorTransp();
}

//-----------------------------------------------------------------------------
DrawGen3Setup::~DrawGen3Setup() {
  if (gm) delete gm;
}

//-----------------------------------------------------------------------------
TGeoNode* DrawGen3Setup::FindNodeByName(TGeoNode* Top, const char* Name) {
  TGeoNode  *top, *found(0);

  if (Top) top = fTop;
  else     top = gm->GetTopNode();
  
  TObjArray* o =  top->GetNodes();

  int n = o->GetEntriesFast();
  
  for (int i=0; i<n; i++) {
    TGeoNode* node = (TGeoNode*) o->UncheckedAt(i);
    const char* name = node->GetName();
    if (strcmp(name,Name) == 0) {
      found = node;
      break;
    }
    else if (node->GetNodes() != NULL) {
      found = FindNodeByName(node,Name);
      if (found) break;
    }
  }
  return found;
}

//-----------------------------------------------------------------------------
TGeoNode* DrawGen3Setup::FindNodeByVolumeName(TGeoNode* Top, const char* VolumeName) {
  TGeoNode  *top, *found(0);

  if (Top) top = Top;
  else     top = fTop;
  
  TObjArray* o =  top->GetNodes();

  int n = o->GetEntriesFast();
  
  for (int i=0; i<n; i++) {
    TGeoNode* node = (TGeoNode*) o->UncheckedAt(i);
    const char* name = node->GetVolume()->GetName();
    if (strcmp(name,VolumeName) == 0) {
      found = node;
      break;
    }
    else if (node->GetNodes() != NULL) {
      found = FindNodeByVolumeName(node,VolumeName);
      if (found) break;
    }
  }
  return found;
}

//-----------------------------------------------------------------------------
// assume that we're looking for one of the daughters
//-----------------------------------------------------------------------------
TGeoNode* DrawGen3Setup::FindNodeByVolumeName(TGeoVolume* Top, const char* VolumeName) {
  TGeoVolume  *top;
  TGeoNode*    found(0);

  if (Top) top = Top;
  else     top = gm->GetTopVolume();

  TObjArray* o =  top->GetNodes();

  int n = o->GetEntriesFast();
  
  for (int i=0; i<n; i++) {
    TGeoNode* node = (TGeoNode*) o->UncheckedAt(i);
    const char* name = node->GetVolume()->GetName();
    if (strcmp(name,VolumeName) == 0) {
      found = node;
      break;
    }
    else if (node->GetNodes() != NULL) {
      found = FindNodeByVolumeName(node,VolumeName);
      if (found) break;
    }
  }
  return found;
}

//-----------------------------------------------------------------------------
void DrawGen3Setup::SetRecursiveVisibilityByName(TGeoNode* Node, const char* Pattern, int OnOff) {

  TString name(Node->GetName());
  
  if (name.Index(Pattern) >= 0) {
    Node->SetVisibility(OnOff);
    //std::cout <<"hiding "<< name << std::endl;
  }
				        // Descend recursively into each daughter TGeoNode.
  int nd = Node->GetNdaughters();
  for (int i=0; i<nd; ++i) {
    TGeoNode * d = Node->GetDaughter(i);
    SetRecursiveVisibilityByName(d,Pattern,OnOff);
  }
}

//-----------------------------------------------------------------------------
void DrawGen3Setup::SetRecursiveVisibility(TGeoNode* Node, int OnOff) {

  TString name(Node->GetName());
  
  Node->SetVisibility(OnOff);
    //std::cout <<"hiding "<< name << std::endl;
				        // Descend recursively into each daughter TGeoNode.
  int nd = Node->GetNdaughters();
  for (int i=0; i<nd; ++i) {
    TGeoNode * d = Node->GetDaughter(i);
    SetRecursiveVisibility(d,OnOff);
  }
}

//-----------------------------------------------------------------------------
void DrawGen3Setup::SetRecursiveVisibilityByMaterial(TGeoNode* Node, const char* Material, int OnOff) {

  TString mat(Node->GetVolume()->GetMaterial()->GetName());
  
  if (mat.Index(Material) >= 0) Node->SetVisibility(OnOff);

				        // Descend recursively into each daughter TGeoNode.
  int ndau = Node->GetNdaughters();
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * d = Node->GetDaughter(i);
    SetRecursiveVisibilityByMaterial(d,Material,OnOff);
  }
}

//-----------------------------------------------------------------------------
void DrawGen3Setup::SetRecursiveColorTranspByName(TGeoNode* Node, const char* Name, Int_t Color, Int_t Transp) {

  
  TString node_name = Node->GetName();
  TGeoVolume*  vol = Node->GetVolume();

  if (node_name.Index(Name) >= 0) {
    vol->SetLineColor   (Color );
    vol->SetTransparency(Transp);
  }
     
  int nd = vol->GetNdaughters();
  for (int i=0; i<nd; i++) {
    TGeoNode* dn = vol->GetNode(i);
    SetRecursiveColorTranspByName(dn, Name, Color, Transp);
  }
}

//-----------------------------------------------------------------------------
void DrawGen3Setup::SetRecursiveColorTranspByMaterial(TGeoNode* Node, const char* MatName, Int_t Color, Int_t Transp) {
  
  TGeoVolume*  vol = Node->GetVolume();
  TString mat_name = vol->GetMaterial()->GetName();

  if (mat_name.Index(MatName) >= 0) {
    vol->SetLineColor   (Color );
    vol->SetTransparency(Transp);
  }
     
  int nd = vol->GetNdaughters();
  for (int i=0; i<nd; i++) {
    TGeoNode* dn = vol->GetNode(i);
    SetRecursiveColorTranspByMaterial(dn, MatName, Color, Transp);
  }
}

//-----------------------------------------------------------------------------
void DrawGen3Setup::SetRecursiveVisibilityColorTranspByNameAndMaterial(TGeoNode*   Top    ,
								       const char* Name   ,
								       const char* MatName,
								       int         Visibility,
								       int         Color  ,
								       int         Transp) {
  TGeoVolume*  vol = Top->GetVolume();
  TString name     = vol->GetName();
  TString mat_name = vol->GetMaterial()->GetName();

  if ((name.Index(Name) >= 0) && (mat_name.Index(MatName) >= 0)) {
    Top->SetVisibility  (Visibility);
    vol->SetLineColor   (Color );
    vol->SetTransparency(Transp);
  }
     
  int nd = vol->GetNdaughters();
  for (int i=0; i<nd; i++) {
    TGeoNode* node = vol->GetNode(i);
    SetRecursiveVisibilityColorTranspByNameAndMaterial(node,Name,MatName,Visibility,Color,Transp);
  }
}

//-----------------------------------------------------------------------------
void DrawGen3Setup::SetRecursiveColorTransp(TGeoVolume *Vol, Int_t Color, Int_t Transp) {

  TString name = Vol->GetName();
  
  int col    = Color;
  int transp = Transp;

  if      (name.Index("TargetFoil") >= 0) { col = kBlue+4;  transp = 10; }
  else if (name.Index("CaloPipe"  ) >= 0) { col = kOrange+7; }
    
  if (col    >=0 ) Vol->SetLineColor   (col   );
  if (Transp >=0 ) Vol->SetTransparency(transp);
     
  int nd = Vol->GetNdaughters();
  for (int i=0; i<nd; i++) {
    TGeoVolume* vd = Vol->GetNode(i)->GetVolume();
    SetRecursiveColorTransp(vd, Color, transp);
  }
}

//-----------------------------------------------------------------------------
// everything is kCyan by default
// production tsrget is kRed+3
//-----------------------------------------------------------------------------
void DrawGen3Setup::SetDefaultColorTransp() {
  int default_color = kCyan-10;
  SetRecursiveColorTransp(fTop->GetVolume(),default_color,fTransp);

//-----------------------------------------------------------------------------
// aluminum rod
//-----------------------------------------------------------------------------
  TGeoVolume* sh = gm->GetVolume("ARod");
  int col        = 844;
  SetRecursiveColorTransp(sh,col,fTransp);
//-----------------------------------------------------------------------------
// wall
//-----------------------------------------------------------------------------
  TGeoVolume* src1 = gm->GetVolume("Wall");
  col              = 844;
  SetRecursiveColorTransp(src1,col,20);
//-----------------------------------------------------------------------------
// color coating
//-----------------------------------------------------------------------------
  TGeoVolume* coat = gm->GetVolume("Coat");
  col             = kRed+2;
  SetRecursiveColorTransp(coat,col,fTransp);
//-----------------------------------------------------------------------------
// color Collimator
//-----------------------------------------------------------------------------
  TGeoVolume* src2 = gm->GetVolume("Coll");
  col = 755;
  SetRecursiveColorTransp(src2,col,20);
//-----------------------------------------------------------------------------
// color virtual detector
//-----------------------------------------------------------------------------
  TGeoVolume* vd1 = gm->GetVolume("Vd1");

  col             = kGray;
  SetRecursiveColorTransp(vd1,col,fTransp);
//-----------------------------------------------------------------------------
// color proton absorbers, start from TS1
//-----------------------------------------------------------------------------
  TGeoVolume* sensor = gm->GetVolume("GaasSensor");
  col             = kRed+2;
  SetRecursiveColorTransp(sensor,col,10);
//-----------------------------------------------------------------------------
// proceed with the photodiodes
//-----------------------------------------------------------------------------
  col             = kBlue-2;

  TGeoVolume* pd0 = gm->GetVolume("InGaAsPD0");
  SetRecursiveColorTransp(pd0,col,40);

  TGeoVolume* pd1 = gm->GetVolume("InGaAsPD1");
  if (pd1) SetRecursiveColorTransp(pd1,col,40);

  TGeoVolume* pd2 = gm->GetVolume("InGaAsPD2");
  if (pd2) SetRecursiveColorTransp(pd2,col,40);

  TGeoVolume* pd3 = gm->GetVolume("InGaAsPD3");
  if (pd3) SetRecursiveColorTransp(pd3,col,40);

  TGeoVolume* pd4 = gm->GetVolume("InGaAsPD4");
  if (pd4) SetRecursiveColorTransp(pd4,col,40);

  TGeoVolume* pd5 = gm->GetVolume("InGaAsPD5");
  if (pd5) SetRecursiveColorTransp(pd5,col,40);

  TGeoVolume* pd6 = gm->GetVolume("InGaAsPD6");
  if (pd6) SetRecursiveColorTransp(pd6,col,40);
}

//-----------------------------------------------------------------------------
// the names are this is for v4_0_6
// clipping: xc = -387.93
//           zlen = 1429
//-----------------------------------------------------------------------------
void DrawGen3Setup::Draw() {

  TGeoVolume* top = fTop->GetVolume();

  TObjArray* list_of_nodes = top->GetNodes();

  int n_nodes = list_of_nodes->GetEntries();

  top->Draw("ogl");
}



//-----------------------------------------------------------------------------

