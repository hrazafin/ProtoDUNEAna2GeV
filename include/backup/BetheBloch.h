#ifndef BETHEBLOCH_H
#define BETHEBLOCH_H

#include <map>

using namespace std;

class TSpline3;

class BetheBloch {

 public:

  BetheBloch();
  BetheBloch(int pdg);

  void SetPdgCode(int pdg);

  int GetPdgCode(){ return pdgcode;};

  double meandEdx(double KE);

  double MPVdEdx(double KE, double pitch);

  double IntegratedEdx(double KE0, double KE1, int n = 10000);

  double RangeFromKE(double KE);

  double RangeFromKESpline(double KE);

  double KEFromRangeSpline(double range);

  double KEAtLength(double KE0, double tracklength);

  void CreateSplineAtKE(int iKE);

  
 private:

  int pdgcode;
  double mass;
  int charge;

  TSpline3 *sp_KE_range;
  TSpline3 *sp_range_KE;

  map<int, TSpline3*> spmap;
  
  double densityEffect(double beta, double gamma);

  double betaGamma(double KE);

  void CreateSplines(int np = 1000, double minke = .01, double maxke = 2e5);

};

#endif
