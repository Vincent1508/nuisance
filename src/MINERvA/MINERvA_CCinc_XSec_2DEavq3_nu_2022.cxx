#include "MINERvA_CCinc_XSec_2DEavq3_nu_2022.h"

#include <cmath>
#include <cstdlib>

namespace {

struct EavQ3Bin {
  double q3_min;
  double q3_max;
  double eav_min;
  double eav_max;
};

// Exact ordering of Table II of the 2022 MINERvA release.
static const int kNBins = 44;

static const EavQ3Bin kBins[kNBins] = {
    {0.0, 0.2, 0.00, 0.04}, // 0
    {0.2, 0.3, 0.00, 0.04}, // 1
    {0.3, 0.4, 0.00, 0.04}, // 2
    {0.4, 0.6, 0.00, 0.04}, // 3
    {0.6, 0.9, 0.00, 0.04}, // 4
    {0.9, 1.2, 0.00, 0.04}, // 5

    {0.0, 0.2, 0.04, 0.08}, // 6
    {0.2, 0.3, 0.04, 0.08}, // 7
    {0.3, 0.4, 0.04, 0.08}, // 8
    {0.4, 0.6, 0.04, 0.08}, // 9
    {0.6, 0.9, 0.04, 0.08}, // 10
    {0.9, 1.2, 0.04, 0.08}, // 11

    {0.0, 0.2, 0.08, 0.12}, // 12
    {0.2, 0.3, 0.08, 0.12}, // 13
    {0.3, 0.4, 0.08, 0.12}, // 14
    {0.4, 0.6, 0.08, 0.12}, // 15
    {0.6, 0.9, 0.08, 0.12}, // 16
    {0.9, 1.2, 0.08, 0.12}, // 17

    {0.0, 0.2, 0.12, 0.24}, // 18

    {0.2, 0.3, 0.12, 0.16}, // 19
    {0.3, 0.4, 0.12, 0.16}, // 20
    {0.4, 0.6, 0.12, 0.16}, // 21
    {0.6, 0.9, 0.12, 0.16}, // 22
    {0.9, 1.2, 0.12, 0.16}, // 23

    {0.2, 0.3, 0.16, 0.24}, // 24
    {0.3, 0.4, 0.16, 0.24}, // 25
    {0.4, 0.6, 0.16, 0.24}, // 26
    {0.6, 0.9, 0.16, 0.24}, // 27
    {0.9, 1.2, 0.16, 0.24}, // 28

    {0.3, 0.4, 0.24, 0.34}, // 29
    {0.4, 0.6, 0.24, 0.32}, // 30
    {0.6, 0.9, 0.24, 0.32}, // 31
    {0.9, 1.2, 0.24, 0.32}, // 32

    {0.3, 0.4, 0.34, 0.40}, // 33
    {0.4, 0.6, 0.32, 0.40}, // 34
    {0.6, 0.9, 0.32, 0.40}, // 35
    {0.9, 1.2, 0.32, 0.40}, // 36

    {0.4, 0.6, 0.40, 0.60}, // 37
    {0.6, 0.9, 0.40, 0.60}, // 38
    {0.9, 1.2, 0.40, 0.60}, // 39

    {0.6, 0.9, 0.60, 1.00}, // 40

    {0.9, 1.2, 0.60, 0.80}, // 41
    {0.9, 1.2, 0.80, 1.00}, // 42
    {0.9, 1.2, 1.00, 1.20}  // 43
};


// -------------------------------------------------------------------------
// Truth-level Eavail.
//
// Units returned: MeV.
//
// Proton and charged-pion contributions are kinetic energies.
// Neutrons and the outgoing muon are excluded.
// Other visible final-state particles contribute their total energy.
//
// Kept local to this measurement so that we do not modify the historical
// MINERvA low-recoil implementation used by the older sample.
// -------------------------------------------------------------------------
double GetEavail2022(FitEvent *event) {

  double eavail = 0.0; // MeV

  // Nucleon mass used in the MINERvA strange-baryon prescription.
  // FitEvent four-vectors are in MeV.
  static const double mN = 939.0; // MeV

  for (UInt_t i = 0; i < event->Npart(); ++i) {

    FitParticle *part = event->PartInfo(i);

    // We only want final-state particles.
    if (!part->fIsAlive)
      continue;

    const int pdg = part->fPID;
    const int apdg = std::abs(pdg);

    // ---------------------------------------------------------------
    // Do not include the outgoing charged lepton.
    // ---------------------------------------------------------------
    if (apdg == 13)
      continue;

    // ---------------------------------------------------------------
    // Neutrons contribute no available energy.
    // ---------------------------------------------------------------
    if (pdg == 2112)
      continue;

    // ---------------------------------------------------------------
    // Protons:
    //
    //     Eavail contribution = kinetic energy
    // ---------------------------------------------------------------
    if (pdg == 2212) {
      eavail += part->KE();
      continue;
    }

    // ---------------------------------------------------------------
    // Charged pions:
    //
    //     Eavail contribution = kinetic energy
    // ---------------------------------------------------------------
    if (apdg == 211) {
      eavail += part->KE();
      continue;
    }

    // ---------------------------------------------------------------
    // Strange baryons:
    //
    // MINERvA removes one nucleon mass from their total energy.
    //
    // Lambda:       3122
    // Sigma-:       3112
    // Sigma0:       3212
    // Sigma+:       3222
    // Xi-:          3312
    // Xi0:          3322
    // Omega-:       3334
    // ---------------------------------------------------------------
    // Strange baryons
    if (pdg == 3122 ||
	pdg == 3112 ||
	pdg == 3212 ||
	pdg == 3222 ||
	pdg == 3312 ||
	pdg == 3322 ||
	pdg == 3334) {

      eavail += part->fP.E() - mN;
      continue;
    }

    // Antibaryons
    if (pdg == -2212 ||
	pdg == -2112 ||
	pdg == -3122 ||
	pdg == -3112 ||
	pdg == -3212 ||
	pdg == -3222 ||
	pdg == -3312 ||
	pdg == -3322 ||
	pdg == -3334) {

      eavail += part->fP.E() + mN;
      continue;
    }

    // ---------------------------------------------------------------
    // All other final-state particles:
    //
    //     Eavail contribution = total energy
    //
    // This automatically includes pi0, photons, electrons, etc.
    // ---------------------------------------------------------------
    eavail += part->fP.E();
  }

  return eavail;
}

} // namespace


// *************************************************************************
MINERvA_CCinc_XSec_2DEavq3_nu_2022::
MINERvA_CCinc_XSec_2DEavq3_nu_2022(nuiskey samplekey) {
// *************************************************************************

  std::string descrip =
      "MINERvA 2022 inclusive CC numu low-q3 measurement on CH. "
      "Measurement represented internally as 44 flattened "
      "(q3,Eavail) bins.";

  fSettings = LoadSampleSettings(samplekey);

  fSettings.SetDescription(descrip);

  // The statistical histogram is flattened. Bin i corresponds to
  // published Table-II bin i.
  fSettings.SetXTitle("Published (q_{3},E_{avail}) bin index");

  fSettings.SetYTitle(
      "d^{2}#sigma/dq_{3}dE_{avail} "
      "(10^{-38} cm^{2}/GeV^{2}/nucleon)");

  // Artificial flattened bins must NOT introduce a bin-width correction.
  this->SetFitOptions("NOWIDTH");

  fSettings.SetAllowedTypes(
      "FIX,FREE,SHAPE/FULL,DIAG/MASK",
      "FIX/FULL");

  // 2022 cross section is flux averaged over the full 0--100 GeV flux.
  fSettings.SetEnuRange(0.0, 100.0);

  fSettings.DefineAllowedTargets("C,H");
  fSettings.DefineAllowedSpecies("numu");

  fSettings.SetTitle(
      "MINERvA_CCinc_XSec_2DEavq3_nu_2022");

  fSettings.SetDataInput(
      FitPar::GetDataBase() +
      "/MINERvA/CCEavq3_2022/data_1D.txt");

  fSettings.SetCovarInput(
      FitPar::GetDataBase() +
      "/MINERvA/CCEavq3_2022/covar_1D.txt");

  splitMEC_PN_NN =
      FitPar::Config().GetParB("Modes_split_PN_NN");

  FinaliseSampleSettings();


  // Standard flux-averaged cross-section normalization.
  //
  // With NOWIDTH, Measurement1D::ScaleEvents() will apply fScaleFactor
  // without dividing by the artificial flattened-bin width.
  fScaleFactor =
      (GetEventHistogram()->Integral("width") *
       1E-38 / (fNEvents + 0.0)) /
      TotalIntegratedFlux();


  // SetDataFromTextFile creates fDataHist itself.
  //
  // data_1D.txt therefore contains 45 rows defining 44 unit-width
  // flattened bins [0,1], [1,2], ..., [43,44].
  SetDataFromTextFile(
      fSettings.GetDataInput());

  ScaleData(1E-38);

  if (fDataHist->GetNbinsX() != kNBins) {
    NUIS_ABORT(
        "MINERvA_CCinc_XSec_2DEavq3_nu_2022: "
        "expected 44 data bins, got "
        << fDataHist->GetNbinsX());
  }


  // This is the covariance matrix itself, NOT a Cholesky decomposition.
  SetCovarFromTextFile(
      fSettings.GetCovarInput(),
      kNBins);

  FinaliseMeasurement();
}


// *************************************************************************
void MINERvA_CCinc_XSec_2DEavq3_nu_2022::
FillEventVariables(FitEvent *event) {
// *************************************************************************

  // Preserve the existing MINERvA MEC PN/NN mode splitting convention.
  if (splitMEC_PN_NN) {

    int npr = 0;
    int nne = 0;

    for (UInt_t j = 0; j < event->Npart(); ++j) {

      // Same convention as the historical implementation:
      // count initial/non-alive nucleons for the MEC mode split.
      if ((event->PartInfo(j))->fIsAlive)
        continue;

      if ((event->PartInfo(j))->fPID == 2212)
        ++npr;
      else if ((event->PartInfo(j))->fPID == 2112)
        ++nne;
    }

    if (event->Mode == 2 && npr == 1 && nne == 1) {
      event->Mode = 2;
      Mode = 2;
    }
    else if (event->Mode == 2 && npr == 0 && nne == 2) {
      event->Mode = 3;
      Mode = 3;
    }
  }


  FitParticle *muon =
      event->GetHMFSParticle(13);

  FitParticle *neutrino =
      event->GetNeutrinoIn();

  if (!muon || !neutrino) {
    fXVar = -999.0;
    return;
  }


  // Truth three-momentum transfer.
  //
  // GENIE/NUISANCE four-vectors are in MeV, hence /1000 -> GeV.
  TLorentzVector q =
      neutrino->fP - muon->fP;

  const double q3 =
      q.Vect().Mag() / 1000.0;


  // Truth Eavail in GeV.
  const double eavail =
      GetEavail2022(event) / 1000.0;


  const int ibin =
      GetPublishedBin(q3, eavail);

  if (ibin < 0) {
    fXVar = -999.0;
    return;
  }


  // Unit-width flattened bin:
  //
  // published bin 0 -> x=0.5
  // published bin 1 -> x=1.5
  // ...
  // published bin 43 -> x=43.5
  fXVar =
      static_cast<double>(ibin) + 0.5;
}


// *************************************************************************
bool MINERvA_CCinc_XSec_2DEavq3_nu_2022::
isSignal(FitEvent *event) {
// *************************************************************************

  // Inclusive CC numu interaction within the configured neutrino-energy
  // range.
  if (!SignalDef::isCCINC(
          event,
          14,
          EnuMin,
          EnuMax))
    return false;


  FitParticle *muon =
      event->GetHMFSParticle(13);

  FitParticle *neutrino =
      event->GetNeutrinoIn();

  if (!muon || !neutrino)
    return false;


  // Muon momentum in GeV/c.
  const double pmu =
      muon->fP.P() / 1000.0;

  if (pmu <= 1.5 || pmu >= 20.0)
    return false;


  // Muon angle relative to the incoming neutrino direction.
  const double cosTheta =
      muon->fP.Vect().Unit().Dot(
          neutrino->fP.Vect().Unit());

  // theta_mu < 20 degrees
  static const double cos20 =
      0.9396926207859084;

  if (cosTheta <= cos20)
    return false;


  // Restrict signal to the actual published phase space.
  TLorentzVector q =
      neutrino->fP - muon->fP;

  const double q3 =
      q.Vect().Mag() / 1000.0;

  const double eavail =
      GetEavail2022(event) / 1000.0;

  if (GetPublishedBin(q3, eavail) < 0)
    return false;


  return true;
}


// *************************************************************************
int MINERvA_CCinc_XSec_2DEavq3_nu_2022::
GetPublishedBin(double q3, double eavail) const {
// *************************************************************************

  for (int i = 0; i < kNBins; ++i) {

    if (q3 >= kBins[i].q3_min &&
        q3 <  kBins[i].q3_max &&
        eavail >= kBins[i].eav_min &&
        eavail <  kBins[i].eav_max) {

      return i;
    }
  }

  return -1;
}


// *************************************************************************
double MINERvA_CCinc_XSec_2DEavq3_nu_2022::
GetBinArea(int i) const {
// *************************************************************************

  if (i < 0 || i >= kNBins)
    return 0.0;

  return
      (kBins[i].q3_max - kBins[i].q3_min) *
      (kBins[i].eav_max - kBins[i].eav_min);
}


// *************************************************************************
void MINERvA_CCinc_XSec_2DEavq3_nu_2022::
ScaleEvents() {
// *************************************************************************

  // Standard Measurement1D scaling.
  //
  // Since NOWIDTH is set, this does NOT divide by the artificial
  // flattened-bin width.
  Measurement1D::ScaleEvents();


  // Convert the integrated cross section accumulated in each physical
  // (q3,Eavail) region into the published double-differential quantity:
  //
  //   d2sigma / dq3 dEavail
  //
  // by dividing by the true physical bin area.
  for (int i = 0; i < kNBins; ++i) {

    const double area =
        GetBinArea(i);

    if (area <= 0.0)
      continue;

    const int rootBin =
        i + 1;

    fMCHist->SetBinContent(
        rootBin,
        fMCHist->GetBinContent(rootBin) / area);

    fMCHist->SetBinError(
        rootBin,
        fMCHist->GetBinError(rootBin) / area);
  }
}
