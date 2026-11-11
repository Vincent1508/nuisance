#ifndef MINERVA_CCINC_XSEC_2DEAVQ3_NU_2022_H
#define MINERVA_CCINC_XSEC_2DEAVQ3_NU_2022_H

#include "Measurement1D.h"
#include "MINERvA_SignalDef.h"

class MINERvA_CCinc_XSec_2DEavq3_nu_2022 : public Measurement1D {

 public:

  MINERvA_CCinc_XSec_2DEavq3_nu_2022(nuiskey samplekey);

  virtual ~MINERvA_CCinc_XSec_2DEavq3_nu_2022() {}


  void FillEventVariables(FitEvent *event);

  bool isSignal(FitEvent *event);

  void ScaleEvents();


 private:

  int GetPublishedBin(
		      double q3,
		      double eavail) const;

  double GetBinArea(
		    int i) const;


  bool splitMEC_PN_NN;
};

#endif
