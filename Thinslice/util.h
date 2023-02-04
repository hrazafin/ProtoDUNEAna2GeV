#ifndef UTIL_H
#define UTIL_H

#include "TH1D.h"
#include "TH2D.h"
#include "EventType.h"

void FillHistVec1D(TH1D *hist[pi::nIntTypes+1], const double &value, const int &partype, bool fill_underflow=true, bool fill_overflow=true);

void FillHistVec2D(TH2D *hist[pi::nIntTypes+1], const double &value1, const double &value2, const int &partype);

void FillHist1D(TH1D *hist, const double &value, const double &wei);

double GetPionKE(double length);

double GetTheta(double x, double y, double z);

double GetPhi(double x, double y, double z);

double GetThetaxz(double x, double y, double z);

double GetThetayz(double x, double y, double z);

#endif
