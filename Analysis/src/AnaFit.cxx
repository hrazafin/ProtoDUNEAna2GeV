#ifndef _ANAFIT_H_
#define _ANAFIT_H_

namespace AnaFit
{

  void Chi2FCN(int &npars, double *grad, double &value, double *par, int flag)
  {
    value = 0;
    double pi0Mass   = par[0];
    double optimalE1 = par[1];
    double optimalE2 = par[2];
    double openAngle = par[3];
    double sigmaM    = par[4];
    double E1        = par[5];
    double E2        = par[6];
    double sigmaE1   = par[7];
    double sigmaE2   = par[8]; 
    double optimalopenAngle = par[9];
    double sigmaopenAngle = par[10];
    double demo = pow((pi0Mass*pi0Mass - 2*optimalE1*optimalE2*(1-cos(optimalopenAngle))),2);
    double tmp = demo/(sigmaM*sigmaM);
    double penalty = pow((optimalE1 - E1),2)/(sigmaE1*sigmaE1) + pow((optimalE2 - E2),2)/(sigmaE2*sigmaE2) + pow((optimalopenAngle - openAngle),2)/(sigmaopenAngle*sigmaopenAngle);
    value = tmp+penalty;
  }

  void KinematicFitting(double pi0Mass, double sigmaM, double openAngle, double E1, double E2, double sigmaE1, double sigmaE2, double &optE1, double &optE2)
  {
    TMinuit *min = new TMinuit(9);
    min->SetFCN(Chi2FCN);
    
    min->DefineParameter(0, "pi0Mass", pi0Mass, 130, 200, 500);
    min->DefineParameter(1, "optimalE1", E1, 0.1, -100, 300);
    min->DefineParameter(2, "optimalE2", E2, 0.1, -100, 300);
    min->DefineParameter(3, "openAngle", openAngle, 1, 0, 80);
    min->DefineParameter(4, "sigmaM", sigmaM, 0.1, 0, 5);
    min->DefineParameter(5, "E1", E1, 0.5, 0, 3);
    min->DefineParameter(6, "E2", E2, 0.5, 0, 3);
    min->DefineParameter(7, "sigmaE1", sigmaE1, 0.1, 0, 5);
    min->DefineParameter(8, "sigmaE2", sigmaE2, 0.1, 0, 5);
    min->DefineParameter(9, "optimalopenAngle", openAngle, 1, 0, 80);
    min->DefineParameter(10, "sigmaopenAngle", 0.3593, 0.1, 0, 5);

    min->FixParameter(10);
    min->FixParameter(0);
    min->FixParameter(3);
    min->FixParameter(4);
    min->FixParameter(5);
    min->FixParameter(6);
    min->FixParameter(7);
    min->FixParameter(8);
    
    min->Command("MIGRAD");
    cout << "openAngle1: " << openAngle << endl;
    double par[11], parerr[11];
    TString name[11] = {"pi0Mass", "optimalE1", "optimalE2", "openAngle", "sigmaM", "E1", "E2", "sigmaE1", "sigmaE2", "optimalopenAngle", "sigmaopenAngle" }; 
    for (int i = 0; i < 11; i++) {
      min->GetParameter(i, par[i], parerr[i]);
      cout << name[i] << ": " << par[i] << " " << parerr[i] << endl;
    }
    optE1 = par[1]; optE2 = par[2];
    cout << "ini Mass: " << sqrt(2*E1*E2*(1-cos(openAngle))) << endl;
    cout << "Mass: " << sqrt(2*par[1]*par[2]*(1-cos(par[9]))) << endl;
  }
}

#endif
