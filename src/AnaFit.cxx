#ifndef _ANAFIT_H_
#define _ANAFIT_H_

#include "../include/AnaIO.h"

namespace AnaFit
{
  // Get Covariance Matrix in a vector
  vector<double> CovarianceMatrixElements(vector<double> X, vector<double> Y, double X_bar, double Y_bar){
    vector<double> CovarianceMatrix;
    double sum_sigma11 = 0, sum_sigma22 = 0, sum_sigma12 = 0;
    for(unsigned int i=0; i < X.size(); i++){
      sum_sigma11 += pow((X[i]-X_bar),2);
      sum_sigma22 += pow((Y[i]-Y_bar),2);
      sum_sigma12 += (X[i]-X_bar)*(Y[i]-Y_bar);
    }
    /*
    CovarianceMatrix.push_back(sum_sigma11/X.size());
    CovarianceMatrix.push_back(sum_sigma12/X.size());
    CovarianceMatrix.push_back(sum_sigma12/X.size());
    CovarianceMatrix.push_back(sum_sigma22/X.size());
    */
/*    
    CovarianceMatrix.push_back(sum_sigma11/X.size());
    CovarianceMatrix.push_back(sum_sigma12/X.size());
    CovarianceMatrix.push_back(0.00000001);
    CovarianceMatrix.push_back(sum_sigma12/X.size());
    CovarianceMatrix.push_back(sum_sigma22/X.size());
    CovarianceMatrix.push_back(0.00000001);
    CovarianceMatrix.push_back(0.00000001);
    CovarianceMatrix.push_back(0.00000001);
    CovarianceMatrix.push_back(1);
*/   

    CovarianceMatrix.push_back(1);
    CovarianceMatrix.push_back(0.00000001);
    CovarianceMatrix.push_back(0.00000001);
    CovarianceMatrix.push_back(0.00000001);
    CovarianceMatrix.push_back(sum_sigma22/X.size());
    CovarianceMatrix.push_back(sum_sigma12/X.size());
    CovarianceMatrix.push_back(0.00000001);
    CovarianceMatrix.push_back(sum_sigma12/X.size());
    CovarianceMatrix.push_back(sum_sigma11/X.size());

    return CovarianceMatrix;

  }

  // Get Correlation Matrix in a vector
  vector<double> CorrelationMatrixElements(vector<double> X, vector<double> Y, double X_bar, double Y_bar){
    vector<double> CorrelationMatrix;
    double sum_sigma11 = 0, sum_sigma22 = 0, sum_sigma12 = 0;
    double sigma1 = 0, sigma2 = 0;
    for(unsigned int i=0; i < X.size(); i++){
      sum_sigma11 += pow((X[i]-X_bar),2);
      sum_sigma22 += pow((Y[i]-Y_bar),2);
      sum_sigma12 += (X[i]-X_bar)*(Y[i]-Y_bar);
    }

    sigma1 = sqrt(sum_sigma11/X.size());
    sigma2 = sqrt(sum_sigma22/Y.size());
    
    double cor11 = (sum_sigma11/X.size())/(sigma1*sigma1);
    double cor12 = (sum_sigma12/X.size())/(sigma1*sigma2);
    double cor21 = (sum_sigma12/X.size())/(sigma2*sigma1);
    double cor22 = (sum_sigma22/Y.size())/(sigma2*sigma2);
    CorrelationMatrix.push_back(cor11);
    CorrelationMatrix.push_back(cor12);
    CorrelationMatrix.push_back(cor21);
    CorrelationMatrix.push_back(cor22);

    return CorrelationMatrix;

  }


  // Fill the CVM bin map
  void FillBinMaps(vector<double> LdShowerEnergyTruth, vector<double> SlShowerEnergyTruth, vector<double> OpenAngleTruth, 
                   vector<double> LdShowerEnergyRaw, vector<double> SlShowerEnergyRaw, vector<double> OpenAngle, 
                   vector<double> Bin1, vector<double> Bin2, std::map<std::pair <int,int>,vector<double>> &BinMapE1, 
                   std::map<std::pair <int,int>,vector<double>> &BinMapE2, std::map<std::pair <int,int>, vector<double>> &BinMapTheta, 
                   std::map<std::pair <int,int>,vector<double>> &BinMapE1Truth, std::map<std::pair <int,int>,vector<double>> &BinMapE2Truth, 
                   std::map<std::pair <int,int>,vector<double>> &BinMapThetaTruth)
  {
    // Create an empty vector
    vector<double> Vect;
    // Set empty vector to maps
    for(unsigned int xx = 1; xx <= Bin1.size() - 1; xx++){
      for(unsigned int yy = 1; yy <= Bin2.size() - 1; yy++){
        std::pair <int,int> bin = std::make_pair (xx,yy);
        BinMapE1[bin] = Vect;
        BinMapE2[bin] = Vect;
        BinMapTheta[bin] = Vect;
        BinMapE1Truth[bin] = Vect;
        BinMapE2Truth[bin] = Vect;
        BinMapThetaTruth[bin] = Vect;
      }
    }
    // Loop over shower sample
    for(unsigned int i = 0; i < LdShowerEnergyTruth.size(); i++){
      for(unsigned int xx = 1; xx <= Bin1.size() - 1; xx++){
        for(unsigned int yy = 1; yy <= Bin2.size() - 1; yy++){
          // Create the pair for map key
          std::pair <int,int> bin = std::make_pair (xx,yy);
          // Need to use truth info 
          if(LdShowerEnergyTruth[i] > Bin1[xx - 1]  && LdShowerEnergyTruth[i] < Bin1[xx] 
            && SlShowerEnergyTruth[i] > Bin2[yy - 1] && SlShowerEnergyTruth[i] < Bin2[yy])
          {
            // Fill in the info for each bin
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

  vector<double> GetVij(std::map<std::pair <int,int>,vector<double>> Bin1Map,
                        std::map<std::pair <int,int>,vector<double>> Bin2Map, std::map<std::pair <int,int>,
                        vector<double>> Bin1MapTruth,std::map<std::pair <int,int>,vector<double>> Bin2MapTruth, std::pair <int,int> bin)
  {
    vector<double> CVM;

    double TruthMean1 = 0;
    double TruthMean2 = 0;
    for(unsigned int jj = 0; jj < Bin1MapTruth[bin].size(); jj++){
      TruthMean1 += Bin1MapTruth[bin][jj]/Bin1MapTruth[bin].size();
      TruthMean2 += Bin2MapTruth[bin][jj]/Bin2MapTruth[bin].size();
    }
        
    CVM = CovarianceMatrixElements(Bin1Map[bin],Bin2Map[bin],TruthMean1,TruthMean2);

    return CVM;
  }

vector<double> GetBinCVM(const double &LdShowerEnergyTruth, const double &SlShowerEnergyTruth, const double &OpenAngleTruth, 
           const double &LdShowerEnergyRaw, const double &SlShowerEnergyRaw, const double &OpenAngle,
           vector<double> LdShowerEnergyTruthVet, vector<double> SlShowerEnergyTruthVet, vector<double> OpenAngleTruthVet, 
           vector<double> LdShowerEnergyRawVet, vector<double> SlShowerEnergyRawVet, vector<double> OpenAngleVet)
{
  // Now creat the bin sizes for Covariance Matrix calculation
  vector<double> BinE1_1x1={0,1.2};
  vector<double> BinE2_1x1={0,0.8};
  vector<double> BinOA_1x1={0,180*TMath::DegToRad()};

  vector<double> BinE1_2x2={0,0.3,1.2};
  vector<double> BinE2_2x2={0,0.15,0.8};
  vector<double> BinOA_2x2={0,40*TMath::DegToRad(),180*TMath::DegToRad()};

  vector<double> BinE1_3x3={0,0.25,0.35,1.2};
  vector<double> BinE2_3x3={0,0.12,0.2,0.8};
  vector<double> BinOA_3x3={0,25*TMath::DegToRad(),60*TMath::DegToRad(),180*TMath::DegToRad()};

  vector<double> BinE1_4x4={0,0.2,0.3,0.4,1.2};
  vector<double> BinE2_4x4={0,0.1,0.16,0.26,0.8};
  vector<double> BinOA_4x4={0,20*TMath::DegToRad(),40*TMath::DegToRad(),80*TMath::DegToRad(),180*TMath::DegToRad()};

  // Reco bin maps
  std::map<std::pair <int,int>,vector<double>> BinMapE1_1x1, BinMapE1_2x2, BinMapE1_3x3, BinMapE1_4x4;
  std::map<std::pair <int,int>,vector<double>> BinMapE2_1x1, BinMapE2_2x2, BinMapE2_3x3, BinMapE2_4x4;
  std::map<std::pair <int,int>,vector<double>> BinMapTheta_1x1, BinMapTheta_2x2, BinMapTheta_3x3, BinMapTheta_4x4;
  // Truth bin maps
  std::map<std::pair <int,int>,vector<double>> BinMapE1Truth_1x1, BinMapE1Truth_2x2, BinMapE1Truth_3x3, BinMapE1Truth_4x4;
  std::map<std::pair <int,int>,vector<double>> BinMapE2Truth_1x1, BinMapE2Truth_2x2, BinMapE2Truth_3x3, BinMapE2Truth_4x4;
  std::map<std::pair <int,int>,vector<double>> BinMapThetaTruth_1x1, BinMapThetaTruth_2x2, BinMapThetaTruth_3x3, BinMapThetaTruth_4x4;

  // Fill each bin maps
  FillBinMaps(LdShowerEnergyTruthVet,SlShowerEnergyTruthVet,OpenAngleTruthVet,LdShowerEnergyRawVet,SlShowerEnergyRawVet,OpenAngleVet,BinE1_1x1,BinE2_1x1,BinMapE1_1x1,BinMapE2_1x1,BinMapTheta_1x1,BinMapE1Truth_1x1,BinMapE2Truth_1x1,BinMapThetaTruth_1x1);
  FillBinMaps(LdShowerEnergyTruthVet,SlShowerEnergyTruthVet,OpenAngleTruthVet,LdShowerEnergyRawVet,SlShowerEnergyRawVet,OpenAngleVet,BinE1_2x2,BinE2_2x2,BinMapE1_2x2,BinMapE2_2x2,BinMapTheta_2x2,BinMapE1Truth_2x2,BinMapE2Truth_2x2,BinMapThetaTruth_2x2);
  FillBinMaps(LdShowerEnergyTruthVet,SlShowerEnergyTruthVet,OpenAngleTruthVet,LdShowerEnergyRawVet,SlShowerEnergyRawVet,OpenAngleVet,BinE1_3x3,BinE2_3x3,BinMapE1_3x3,BinMapE2_3x3,BinMapTheta_3x3,BinMapE1Truth_3x3,BinMapE2Truth_3x3,BinMapThetaTruth_3x3);
  FillBinMaps(LdShowerEnergyTruthVet,SlShowerEnergyTruthVet,OpenAngleTruthVet,LdShowerEnergyRawVet,SlShowerEnergyRawVet,OpenAngleVet,BinE1_4x4,BinE2_4x4,BinMapE1_4x4,BinMapE2_4x4,BinMapTheta_4x4,BinMapE1Truth_4x4,BinMapE2Truth_4x4,BinMapThetaTruth_4x4);

  vector<double> tmp_CVM;
  for(unsigned int xx = 1; xx <= BinOA_3x3.size() - 1; xx++){
    for(unsigned int yy = 1; yy <= BinE2_3x3.size() - 1; yy++){
      if(OpenAngleTruth > BinOA_3x3[xx - 1]  && OpenAngleTruth < BinOA_3x3[xx] 
        && SlShowerEnergyTruth > BinE2_3x3[yy - 1] && SlShowerEnergyTruth < BinE2_3x3[yy]){
        std::pair <int,int> bin = std::make_pair (xx,yy);
        tmp_CVM = GetVij(BinMapTheta_3x3,BinMapE2_3x3,BinMapThetaTruth_3x3,BinMapE2Truth_3x3,bin);
        // E1 & E2 dependent
        //tmp_CVM = GetVij(BinMapE1_3x3,BinMapE2_3x3,BinMapE1Truth_3x3,BinMapE2Truth_3x3,bin);
      }
    }
  }
  
  return tmp_CVM;
}


/*
vector<double> GetCVM(vector<double> LdShowerEnergyTruth, vector<double> SlShowerEnergyTruth, vector<double> OpenAngleTruth, 
                      vector<double> LdShowerEnergyRaw, vector<double> SlShowerEnergyRaw, vector<double> OpenAngle,
                      const double &LdShowerEnergyTruth_val, const double &SlShowerEnergyTruth_val, const double &OpenAngleTruth_val, 
                      const double &LdShowerEnergyRaw_val, const double &SlShowerEnergyRaw_val, const double &OpenAngle_val)
{
  vector<double> CVM;
  double V_11_bin1 = 0, V_12_bin1 = 0, V_22_bin1 = 0, V_13_bin1 = 0, V_23_bin1 = 0, V_33_bin1 = 0;
  int sampleSize = LdShowerEnergyRaw.size();
  cout <<  "sampleSize: " << sampleSize << endl;
  for(int ii = 0; ii < sampleSize; ii++){
    V_11_bin1 += pow((LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii]),2);
    V_22_bin1 += pow((SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii]),2);
    V_33_bin1 += pow((OpenAngle[ii]-OpenAngleTruth[ii]),2);
    V_12_bin1 += (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])*(SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii]);
    V_13_bin1 += (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])*(OpenAngle[ii]-OpenAngleTruth[ii]);
    V_23_bin1 += (SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii])*(OpenAngle[ii]-OpenAngleTruth[ii]);
    
    AnaIO::hsigma11->Fill(LdShowerEnergyRaw[ii],(LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])/LdShowerEnergyTruth[ii]);
    AnaIO::hsigma22->Fill(SlShowerEnergyRaw[ii],(SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii])/SlShowerEnergyTruth[ii]);
    //AnaIO::hsigma11->Fill(LdShowerEnergyRaw[ii],(LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])/LdShowerEnergyTruth[ii]);
    //cout << "LD E: " << LdShowerEnergyRaw[ii] <<" sigma11: " << (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])/LdShowerEnergyTruth[ii] << endl;
  }
  CVM.push_back(V_11_bin1/sampleSize); 
  CVM.push_back(V_12_bin1/sampleSize); 
  CVM.push_back(V_13_bin1/sampleSize); 
  CVM.push_back(V_12_bin1/sampleSize); 
  CVM.push_back(V_22_bin1/sampleSize); 
  CVM.push_back(V_23_bin1/sampleSize); 
  CVM.push_back(V_13_bin1/sampleSize); 
  CVM.push_back(V_23_bin1/sampleSize); 
  CVM.push_back(V_33_bin1/sampleSize);

  return CVM;
}
*/

vector<double> GetCVM(vector<double> LdShowerEnergyTruth, vector<double> SlShowerEnergyTruth, vector<double> OpenAngleTruth, 
                      vector<double> LdShowerEnergyRaw, vector<double> SlShowerEnergyRaw, vector<double> OpenAngle,
                      const double &LdShowerEnergyTruth_val, const double &SlShowerEnergyTruth_val, const double &OpenAngleTruth_val, 
                      const double &LdShowerEnergyRaw_val, const double &SlShowerEnergyRaw_val, const double &OpenAngle_val)
{

  // Energy dependent CVM
  double V_11_bin1 = 0, V_12_bin1 = 0, V_22_bin1 = 0, V_13_bin1 = 0, V_23_bin1 = 0, V_33_bin1 = 0;
  double V_11_bin2 = 0, V_12_bin2 = 0, V_22_bin2 = 0, V_13_bin2 = 0, V_23_bin2 = 0, V_33_bin2 = 0;
  double V_11_bin3 = 0, V_12_bin3 = 0, V_22_bin3 = 0, V_13_bin3 = 0, V_23_bin3 = 0, V_33_bin3 = 0;
  double V_11_bin4 = 0, V_12_bin4 = 0, V_22_bin4 = 0, V_13_bin4 = 0, V_23_bin4 = 0, V_33_bin4 = 0;
  double sample_bin1 = 0, sample_bin2 = 0, sample_bin3 = 0, sample_bin4 = 0;

  int sampleSize = LdShowerEnergyRaw.size();
  for(int ii = 0; ii < sampleSize; ii++){
    //if(SlShowerEnergyRaw[ii] < 0.15 && OpenAngleRaw[ii] < 50*TMath::RadToDeg()){
    if(LdShowerEnergyRaw[ii] < 0.3 && SlShowerEnergyRaw[ii] < 0.15){
      V_11_bin1 += pow((LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii]),2);
      V_22_bin1 += pow((SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii]),2);
      V_33_bin1 += pow((OpenAngle[ii]-OpenAngleTruth[ii]),2);
      V_12_bin1 += (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])*(SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii]);
      V_13_bin1 += (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])*(OpenAngle[ii]-OpenAngleTruth[ii]);
      V_23_bin1 += (SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii])*(OpenAngle[ii]-OpenAngleTruth[ii]);
      sample_bin1++;
    }
    //else if(SlShowerEnergyTruth[ii] > 0.15 && OpenAngleTruth[ii] < 50*TMath::RadToDeg()){
    else if(LdShowerEnergyRaw[ii] < 0.3 && SlShowerEnergyRaw[ii] > 0.15){

      V_11_bin2 += pow((LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii]),2);
      V_22_bin2 += pow((SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii]),2);
      V_33_bin2 += pow((OpenAngle[ii]-OpenAngleTruth[ii]),2);
      V_12_bin2 += (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])*(SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii]);
      V_13_bin2 += (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])*(OpenAngle[ii]-OpenAngleTruth[ii]);
      V_23_bin2 += (SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii])*(OpenAngle[ii]-OpenAngleTruth[ii]);
      sample_bin2++;
    }
    //else ifSlShowerEnergyTruth[ii] < 0.15 && OpenAngleTruth[ii] > 50*TMath::RadToDeg()){
    else if(LdShowerEnergyRaw[ii] > 0.3 && SlShowerEnergyRaw[ii] < 0.15){
      
      V_11_bin3 += pow((LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii]),2);
      V_22_bin3 += pow((SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii]),2);
      V_33_bin3 += pow((OpenAngle[ii]-OpenAngleTruth[ii]),2);
      V_12_bin3 += (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])*(SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii]);
      V_13_bin3 += (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])*(OpenAngle[ii]-OpenAngleTruth[ii]);
      V_23_bin3 += (SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii])*(OpenAngle[ii]-OpenAngleTruth[ii]);
      sample_bin3++;
    }
    else{
      V_11_bin4 += pow((LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii]),2);
      V_22_bin4 += pow((SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii]),2);
      V_33_bin4 += pow((OpenAngle[ii]-OpenAngleTruth[ii]),2);
      V_12_bin4 += (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])*(SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii]);
      V_13_bin4 += (LdShowerEnergyRaw[ii]-LdShowerEnergyTruth[ii])*(OpenAngle[ii]-OpenAngleTruth[ii]);
      V_23_bin4 += (SlShowerEnergyRaw[ii]-SlShowerEnergyTruth[ii])*(OpenAngle[ii]-OpenAngleTruth[ii]);
      sample_bin4++;
    }
  }
  vector<double> CVM;
 
  //if(SlShowerEnergyTruth_val < 0.15 && OpenAngle_val < 50*TMath::RadToDeg()){
  if(LdShowerEnergyRaw_val < 0.3 && SlShowerEnergyRaw_val < 0.15){

    CVM.push_back(V_11_bin1/sample_bin1); 
    CVM.push_back(V_12_bin1/sample_bin1); 
    CVM.push_back(V_13_bin1/sample_bin1); 
    CVM.push_back(V_12_bin1/sample_bin1); 
    CVM.push_back(V_22_bin1/sample_bin1); 
    CVM.push_back(V_23_bin1/sample_bin1); 
    CVM.push_back(V_13_bin1/sample_bin1); 
    CVM.push_back(V_23_bin1/sample_bin1); 
    CVM.push_back(V_33_bin1/sample_bin1);
    cout << "CVM 1: " << endl;
  }
  //else if(SlShowerEnergyTruth_val > 0.15 && OpenAngle_val < 50*TMath::RadToDeg()){
  else if(LdShowerEnergyRaw_val < 0.3 && SlShowerEnergyRaw_val > 0.15){

    CVM.push_back(V_11_bin2/sample_bin2); 
    CVM.push_back(V_12_bin2/sample_bin2); 
    CVM.push_back(V_13_bin2/sample_bin2); 
    CVM.push_back(V_12_bin2/sample_bin2); 
    CVM.push_back(V_22_bin2/sample_bin2); 
    CVM.push_back(V_23_bin2/sample_bin2); 
    CVM.push_back(V_13_bin2/sample_bin2); 
    CVM.push_back(V_23_bin2/sample_bin2); 
    CVM.push_back(V_33_bin2/sample_bin2);
    cout << "CVM 2: " << endl;
  }
  //else if(SlShowerEnergyTruth_val < 0.15 && OpenAngle_val > 50*TMath::RadToDeg()){
  else if(LdShowerEnergyRaw_val > 0.3 && SlShowerEnergyRaw_val < 0.15){

    CVM.push_back(V_11_bin3/sample_bin3); 
    CVM.push_back(V_12_bin3/sample_bin3); 
    CVM.push_back(V_13_bin3/sample_bin3); 
    CVM.push_back(V_12_bin3/sample_bin3); 
    CVM.push_back(V_22_bin3/sample_bin3); 
    CVM.push_back(V_23_bin3/sample_bin3); 
    CVM.push_back(V_13_bin3/sample_bin3); 
    CVM.push_back(V_23_bin3/sample_bin3); 
    CVM.push_back(V_33_bin3/sample_bin3);
    cout << "CVM 3: " << endl;
  }
  else{
    CVM.push_back(V_11_bin4/sample_bin4); 
    CVM.push_back(V_12_bin4/sample_bin4); 
    CVM.push_back(V_13_bin4/sample_bin4); 
    CVM.push_back(V_12_bin4/sample_bin4); 
    CVM.push_back(V_22_bin4/sample_bin4); 
    CVM.push_back(V_23_bin4/sample_bin4); 
    CVM.push_back(V_13_bin4/sample_bin4); 
    CVM.push_back(V_23_bin4/sample_bin4); 
    CVM.push_back(V_33_bin4/sample_bin4);
    cout << "CVM 4: " << endl;
  }

  for(auto i : CVM){
      cout << "ele: " << i << endl;
  }
  
  return CVM;

}


}

#endif
