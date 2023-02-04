#ifndef UNFOLD_H
#define UNFOLD_H

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TH2D.h"

class Unfold {
  
 public:

  Unfold(int nb, double xlo, double xhi);

  RooUnfoldResponse response_SliceID_Int;  //Interaction
  RooUnfoldResponse response_SliceID_BeamInt;  //Interaction

  RooUnfoldResponse response_SliceID_Inc;  //Incident
  RooUnfoldResponse response_SliceID_Ini;  //Initial

  RooUnfoldResponse response_SliceID_Pi0KE;
  RooUnfoldResponse response_SliceID_Pi0CosTheta;
  RooUnfoldResponse response_SliceID_Pi0Theta;


  TH1D *eff_num_Int; //Interaction efficiency numerator
  TH1D *eff_den_Int; //Interaction efficiency denominator

  TH1D *eff_num_Inc; //Incident efficiency numerator
  TH1D *eff_den_Inc; //Incident efficiency denominator

  TH1D *eff_num_Ini; //Initial efficiency numerator
  TH1D *eff_den_Ini; //Initial efficiency denominator

  TH1D *eff_num_BeamInt; //Interaction efficiency numerator
  TH1D *eff_den_BeamInt; //Interaction efficiency denominator

  TH1D *eff_num_Pi0KE; //Pi0KE efficiency numerator
  TH1D *eff_den_Pi0KE; //Pi0KE efficiency denominator

  TH1D *eff_num_Pi0CosTheta; //Pi0CosTheta efficiency numerator
  TH1D *eff_den_Pi0CosTheta; //Pi0CosTheta efficiency denominator

  TH1D *eff_num_Pi0Theta; //Pi0Theta efficiency numerator
  TH1D *eff_den_Pi0Theta; //Pi0Theta efficiency denominator



  TH1D *pur_num_Int; //Interaction purity numerator
  TH1D *pur_num_Inc; //Incident purity numerator
  TH1D *pur_den_Int; //Interaction purity denominator
  TH1D *pur_den_Inc; //Incident purity denominator


  TH1D *eff_Int;
  TH1D *eff_Inc;
  TH1D *eff_Ini;
  TH1D *eff_BeamInt;
  TH1D* eff_Pi0KE;
  TH1D* eff_Pi0CosTheta;
  TH1D* eff_Pi0Theta;


  TH1D *pur_Int;
  TH1D *pur_Inc;

  void SaveHistograms();

};

#endif
