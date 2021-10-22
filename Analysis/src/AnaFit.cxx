#ifndef _ANAFIT_H_
#define _ANAFIT_H_

#include "../include/AnaIO.h"

namespace AnaFit
{
  // Get CovarianceMatrix in a vector
  vector<double> CovarianceMatrixElements(vector<double> X, vector<double> Y, double X_bar, double Y_bar){
    vector<double> CovarianceMatrix;
    double sum_sigma11 = 0, sum_sigma22 = 0, sum_sigma12 = 0;
    for(unsigned int i=0; i < X.size(); i++){
      sum_sigma11 += pow((X[i]-X_bar),2);
      sum_sigma22 += pow((Y[i]-Y_bar),2);
      sum_sigma12 += (X[i]-X_bar)*(Y[i]-Y_bar);
    }
    
    CovarianceMatrix.push_back(sum_sigma11/X.size());
    CovarianceMatrix.push_back(sum_sigma12/X.size());
    CovarianceMatrix.push_back(sum_sigma12/X.size());
    CovarianceMatrix.push_back(sum_sigma22/X.size());

    return CovarianceMatrix;

  }


  // Fill the CVM bin map
  void FillBinMaps(vector<double> LdShowerEnergyTruth, vector<double> SlShowerEnergyTruth, vector<double> OpenAngleTruth, 
                   vector<double> LdShowerEnergyRaw, vector<double> SlShowerEnergyRaw, vector<double> OpenAngle, 
                   vector<double> BinE1, vector<double> BinE2, std::map<std::pair <int,int>,vector<double>> &BinMapE1, 
                   std::map<std::pair <int,int>,vector<double>> &BinMapE2, std::map<std::pair <int,int>, vector<double>> &BinMapTheta, 
                   std::map<std::pair <int,int>,vector<double>> &BinMapE1Truth, std::map<std::pair <int,int>,vector<double>> &BinMapE2Truth, 
                   std::map<std::pair <int,int>,vector<double>> &BinMapThetaTruth)
  {
    // Create an empty vector
    vector<double> Vect;
    // Set empty vector to maps
    for(unsigned int xx = 1; xx <= BinE1.size() - 1; xx++){
      for(unsigned int yy = 1; yy <= BinE2.size() - 1; yy++){
        std::pair <int,int> bin = std::make_pair (xx,yy);
        BinMapE1[bin] = Vect;
        BinMapE2[bin] = Vect;
        BinMapTheta[bin] = Vect;
        BinMapE1Truth[bin] = Vect;
        BinMapE2Truth[bin] = Vect;
        BinMapThetaTruth[bin] = Vect;
      }
    }
    // Loop over shower energy
    for(unsigned int i = 0; i < LdShowerEnergyTruth.size(); i++){
      for(unsigned int xx = 1; xx <= BinE1.size() - 1; xx++){
        for(unsigned int yy = 1; yy <= BinE2.size() - 1; yy++){
          // Need to use truth info 
          if(LdShowerEnergyTruth[i] > BinE1[xx - 1]  && LdShowerEnergyTruth[i] < BinE1[xx] 
            && SlShowerEnergyTruth[i] > BinE2[yy - 1] && SlShowerEnergyTruth[i] < BinE2[yy])
          {
            // Create the pair for map key
            std::pair <int,int> bin = std::make_pair (xx,yy);
            BinMapE1[bin].push_back(LdShowerEnergyRaw[i]);
            BinMapE2[bin].push_back(SlShowerEnergyRaw[i]);
            BinMapTheta[bin].push_back(OpenAngle[i]);
            BinMapE1Truth[bin].push_back(LdShowerEnergyTruth[i]);
            BinMapE2Truth[bin].push_back(SlShowerEnergyTruth[i]);
            BinMapThetaTruth[bin].push_back(OpenAngleTruth[i]);
          }       
        }
      }
    }
  }

  // Get CVM elements and fill histograms
  vector<double> FillVij(vector<double> BinE1, vector<double> BinE2, std::map<std::pair <int,int>,vector<double>> BinMapE1,
                         std::map<std::pair <int,int>,vector<double>> BinMapE2, std::map<std::pair <int,int>,vector<double>> BinMapE1Truth,
                         std::map<std::pair <int,int>,vector<double>> BinMapE2Truth, TH2D *V11,TH2D *V12,TH2D *V22, TH2D *TruthBin_Entries)
  {
    vector<double> CVM; // Covariance Matrix vector
    // Calculate the true mean value of each bin
    for(unsigned int xx = 1; xx <= BinE1.size() - 1; xx++){
      for(unsigned int yy = 1; yy <= BinE2.size() - 1; yy++){
        std::pair <int,int> bin = std::make_pair (xx,yy);
        double E1TruthMean = 0;
        double E2TruthMean = 0;
        for(unsigned int jj = 0; jj < BinMapE1Truth[bin].size(); jj++){
          E1TruthMean += BinMapE1Truth[bin][jj]/BinMapE1Truth[bin].size();
          E2TruthMean += BinMapE2Truth[bin][jj]/BinMapE2Truth[bin].size();
        }
        // Get the Covariance Matrix Elements in a vector
        CVM = CovarianceMatrixElements(BinMapE1[bin],BinMapE2[bin],E1TruthMean,E2TruthMean);
        
        // Set each truth bin entries
        TruthBin_Entries->SetBinContent(xx,yy,BinMapE1Truth[bin].size());
        V11->SetBinContent(xx,yy,CVM[0]*1000);
        V12->SetBinContent(xx,yy,CVM[1]*1000);
        V22->SetBinContent(xx,yy,CVM[3]*1000);
      }
    }
  return CVM;
}

void DoCVM(vector<double> LdShowerEnergyTruth, vector<double> SlShowerEnergyTruth, vector<double> OpenAngleTruth, 
           vector<double> LdShowerEnergyRaw, vector<double> SlShowerEnergyRaw, vector<double> OpenAngle)
{
  // Now creat the bin sizes for Covariance Matrix calculation
  vector<double> BinE1_1x1={0,1.2};
  vector<double> BinE2_1x1={0,0.8};

  vector<double> BinE1_2x2={0,0.3,1.2};
  vector<double> BinE2_2x2={0,0.15,0.8};

  vector<double> BinE1_3x3={0,0.25,0.35,1.2};
  vector<double> BinE2_3x3={0,0.12,0.2,0.8};

  vector<double> BinE1_4x4={0,0.2,0.3,0.4,1.2};
  vector<double> BinE2_4x4={0,0.1,0.16,0.26,0.8};

  // Reco bin maps
  std::map<std::pair <int,int>,vector<double>> BinMapE1_1x1, BinMapE1_2x2, BinMapE1_3x3, BinMapE1_4x4;
  std::map<std::pair <int,int>,vector<double>> BinMapE2_1x1, BinMapE2_2x2, BinMapE2_3x3, BinMapE2_4x4;
  std::map<std::pair <int,int>,vector<double>> BinMapTheta_1x1, BinMapTheta_2x2, BinMapTheta_3x3, BinMapTheta_4x4;
  // Truth bin maps
  std::map<std::pair <int,int>,vector<double>> BinMapE1Truth_1x1, BinMapE1Truth_2x2, BinMapE1Truth_3x3, BinMapE1Truth_4x4;
  std::map<std::pair <int,int>,vector<double>> BinMapE2Truth_1x1, BinMapE2Truth_2x2, BinMapE2Truth_3x3, BinMapE2Truth_4x4;
  std::map<std::pair <int,int>,vector<double>> BinMapThetaTruth_1x1, BinMapThetaTruth_2x2, BinMapThetaTruth_3x3, BinMapThetaTruth_4x4;

  // Fill each bin maps
  FillBinMaps(LdShowerEnergyTruth,SlShowerEnergyTruth,OpenAngleTruth,LdShowerEnergyRaw,SlShowerEnergyRaw,OpenAngle,BinE1_1x1,BinE2_1x1,BinMapE1_1x1,BinMapE2_1x1,BinMapTheta_1x1,BinMapE1Truth_1x1,BinMapE2Truth_1x1,BinMapThetaTruth_1x1);
  FillBinMaps(LdShowerEnergyTruth,SlShowerEnergyTruth,OpenAngleTruth,LdShowerEnergyRaw,SlShowerEnergyRaw,OpenAngle,BinE1_2x2,BinE2_2x2,BinMapE1_2x2,BinMapE2_2x2,BinMapTheta_2x2,BinMapE1Truth_2x2,BinMapE2Truth_2x2,BinMapThetaTruth_2x2);
  FillBinMaps(LdShowerEnergyTruth,SlShowerEnergyTruth,OpenAngleTruth,LdShowerEnergyRaw,SlShowerEnergyRaw,OpenAngle,BinE1_3x3,BinE2_3x3,BinMapE1_3x3,BinMapE2_3x3,BinMapTheta_3x3,BinMapE1Truth_3x3,BinMapE2Truth_3x3,BinMapThetaTruth_3x3);
  FillBinMaps(LdShowerEnergyTruth,SlShowerEnergyTruth,OpenAngleTruth,LdShowerEnergyRaw,SlShowerEnergyRaw,OpenAngle,BinE1_4x4,BinE2_4x4,BinMapE1_4x4,BinMapE2_4x4,BinMapTheta_4x4,BinMapE1Truth_4x4,BinMapE2Truth_4x4,BinMapThetaTruth_4x4);

  vector<double> CVM_1x1;
  CVM_1x1 = FillVij(BinE1_1x1,BinE2_1x1,BinMapE1_1x1,BinMapE2_1x1,BinMapE1Truth_1x1,BinMapE2Truth_1x1,AnaIO::hV11_1x1,AnaIO::hV12_1x1,AnaIO::hV22_1x1,AnaIO::hBinSize_1x1);
  vector<double> CVM_2x2;
  CVM_2x2 = FillVij(BinE1_2x2,BinE2_2x2,BinMapE1_2x2,BinMapE2_2x2,BinMapE1Truth_2x2,BinMapE2Truth_2x2,AnaIO::hV11_2x2,AnaIO::hV12_2x2,AnaIO::hV22_2x2,AnaIO::hBinSize_2x2);
  vector<double> CVM_3x3;
  CVM_3x3 = FillVij(BinE1_3x3,BinE2_3x3,BinMapE1_3x3,BinMapE2_3x3,BinMapE1Truth_3x3,BinMapE2Truth_3x3,AnaIO::hV11_3x3,AnaIO::hV12_3x3,AnaIO::hV22_3x3,AnaIO::hBinSize_3x3);
  vector<double> CVM_4x4;
  CVM_4x4 = FillVij(BinE1_4x4,BinE2_4x4,BinMapE1_4x4,BinMapE2_4x4,BinMapE1Truth_4x4,BinMapE2Truth_4x4,AnaIO::hV11_4x4,AnaIO::hV12_4x4,AnaIO::hV22_4x4,AnaIO::hBinSize_4x4);

}


}

#endif
