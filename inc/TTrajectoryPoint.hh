#ifndef TTrajectoryPoint_hh
#define TTrajectoryPoint_hh
//-----------------------------------------------------------------------------
//  Dec 23 2000 P.Murat: point on particle's trajectory
//-----------------------------------------------------------------------------

#include "TVector3.h"


class TTrajectoryPoint: public TObject {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TVector3    fPosition;		// position
  TVector3    fDirection;		// direction cosines
  Double_t    fS;			// path length over the trajectory
  Double_t    fPTotal;			// total momentum (GeV) 
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TTrajectoryPoint();
  TTrajectoryPoint(Double_t  x, Double_t  y, Double_t  z,
		   Double_t nx, Double_t ny, Double_t nz,
		   Double_t  s, Double_t  p);

  virtual ~TTrajectoryPoint();
					// ****** accessors

  TVector3* GetPosition () { return &fPosition;  }
  TVector3* GetDirection() { return &fDirection; }
  Double_t  GetPTotal   () { return fPTotal;     }
  Double_t  S           () { return fS;          }
  Double_t  X           () { return fPosition.X(); }
  Double_t  Y           () { return fPosition.Y(); }
  Double_t  Z           () { return fPosition.Z(); }
  Double_t  Nx          () { return fDirection.X(); }
  Double_t  Ny          () { return fDirection.Y(); }
  Double_t  Nz          () { return fDirection.Z(); }

  Double_t  GetS        () { 
    printf("obsolete, use TTrajectoryPoint::S()\n"); 
    return fS;       
  }
					// point is an array of dimension 8...
					// (x,y,z,nx,ny,nz,ptotal,s)
  void     GetPoint (Double_t* point);

					// ****** modifiers

  void     SetPoint    (Double_t* point);
  void     SetPoint    (Double_t  x, Double_t  y, Double_t  z,
			Double_t nx, Double_t ny, Double_t nz,
			Double_t  s, Double_t  p);

  void     SetPosition (const TVector3* mom) { fPosition = *mom; }
  void     SetDirection(const TVector3* dir) { fPosition = *dir; }
  void     SetPTotal   (Double_t       ptot) { fPTotal   = ptot; }
  void     SetS        (Double_t       s   ) { fS        = s;    }

					// ****** other methods

					// transform point from global to local 
					// coordinate system - first translate,
					// then rotate
					// p0: center of the lrs wrt grs
					// rot: rotation from global to local

  void     GlobalToLocal (const TVector3&   p0, 
			  const TRotation&  rot, 
			  TTrajectoryPoint& pout);

					// transform point from local to global
					// coordinate system - first rotate,
					// then translate. 
					// p0: center of the lrs wrt grs
					// rot: rotation from global to local

  void     LocalToGlobal (const TVector3&   p0, 
			  const TRotation&  rot, 
			  TTrajectoryPoint& pout);

					// ****** overloaded methods of 
					// TObject
  void     Clear(Option_t* opt = "");
  void     Print(Option_t* opt = "") const ;

  //   ClassDef(TTrajectoryPoint,1)
};

#endif
