
#include "TRotation.h"
#include "gaasG4/inc/TTrajectoryPoint.hh"

// ClassImp(TTrajectoryPoint)

//_____________________________________________________________________________
TTrajectoryPoint::TTrajectoryPoint() {
}

//_____________________________________________________________________________
TTrajectoryPoint::TTrajectoryPoint(Double_t  x, Double_t  y, Double_t  z,
				   Double_t nx, Double_t ny, Double_t nz,
				   Double_t  s, Double_t  p) 
{
  SetPoint(x,y,z,nx,ny,nz,s,p);
}

//_____________________________________________________________________________
TTrajectoryPoint::~TTrajectoryPoint() {
}

//_____________________________________________________________________________
void TTrajectoryPoint::SetPoint(Double_t* pin) {
  // pin = x,y,z,nx,ny,nz,ptotal
  fPosition.SetXYZ(pin[0],pin[1],pin[2]);
  fDirection.SetXYZ(pin[3],pin[4],pin[5]);
  fS      = pin[6];
  fPTotal = pin[7];
}

//_____________________________________________________________________________
void TTrajectoryPoint::SetPoint(Double_t  x, Double_t  y, Double_t  z,
				Double_t nx, Double_t ny, Double_t nz,
				Double_t  s, Double_t  p) 
{
  // pin = x,y,z,nx,ny,nz,ptotal

  fPosition.SetXYZ(x,y,z);
  fDirection.SetXYZ(nx,ny,nz);
  fS      = s;
  fPTotal = p;
}

//_____________________________________________________________________________
void TTrajectoryPoint::GetPoint(Double_t* point) {
  // retrieve point = x,y,z,nx,ny,nz,ptotal,s

  fPosition.GetXYZ(point);
  fDirection.GetXYZ(&point[3]);
  point[6] = fPTotal;
  point[7] = fS;
}

//_____________________________________________________________________________
void TTrajectoryPoint::GlobalToLocal(const TVector3&   p0,
				     const TRotation&  rot,
				     TTrajectoryPoint& pout) 
{
  // transform a trajectory point from global to local coordinate system
  // first translate, then rotate

  TVector3* pos = pout.GetPosition();
  TVector3* dir = pout.GetDirection();

  pos->SetXYZ(X()-p0.X(),Y()-p0.Y(),Z()-p0.Z());
  pos->Transform(rot);

  dir->SetXYZ(Nx(),Ny(),Nz());
  dir->Transform(rot);

  pout.fPTotal = fPTotal;
  pout.fS      = fS;
}

//_____________________________________________________________________________
void TTrajectoryPoint::LocalToGlobal(const TVector3&   P0,
				     const TRotation&  Rot,
				     TTrajectoryPoint& Pout) 
{
  // transform a trajectory point from local to global coordinate system
  // first rotate, then translate (p0 and rot are assumed to be defined
  // in the global coordinate system). 
  // `pout' may be the same with `*this'

  TVector3* pos = Pout.GetPosition();
  TVector3* dir = Pout.GetDirection();
					// do it quick and dirty - just need 
					// the code to work...
  TRotation inv_rot = Rot.Inverse();
  pos->Transform(inv_rot);
  pos->SetXYZ(X()+P0.X(),Y()+P0.Y(),Z()+P0.Z());
  
  dir->SetXYZ(Nx(),Ny(),Nz());
  dir->Transform(inv_rot);

  Pout.fPTotal = fPTotal;
  Pout.fS      = fS;
}

//_____________________________________________________________________________
void TTrajectoryPoint::Clear(Option_t* opt) {
}

//_____________________________________________________________________________
void TTrajectoryPoint::Print(Option_t* opt) const {

  if ((strstr(opt,"banner") != 0) || (opt[0] == 0)) {
    printf("         X          Y           Z           Nx         Ny ");
    printf("         Nz       P(total)      S  \n");
  }

  if ((strstr(opt,"data") != 0) || (opt[0] == 0)) {
    printf("%11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n",
	   fPosition.X(),
	   fPosition.Y(),
	   fPosition.Z(),
	   fDirection.X(),
	   fDirection.Y(),
	   fDirection.Z(),
	   fPTotal,
	   fS);
  }
}

