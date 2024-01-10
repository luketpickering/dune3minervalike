#include "MINERvA_SignalDef.h"
#include "Measurement1D.h"

#include "TH3D.h"

class DUNE_3D_MINERvALike : public Measurement1D {
public:
  std::unique_ptr<TH3D> f3DHist;

  std::unique_ptr<TH3D> f3DHist_CCQE;
  std::unique_ptr<TH3D> f3DHist_CC1pi;
  std::unique_ptr<TH3D> f3DHist_CCDIS;

  //********************************************************************
  DUNE_3D_MINERvALike(nuiskey samplekey) {
    //********************************************************************

    // START boilerplate
    //
    // Setup common settings
    fSettings = LoadSampleSettings(samplekey);
    fSettings.SetAllowedTypes("FIX,FREE,SHAPE/FULL,DIAG/MASK", "FIX/FULL");
    fSettings.SetEnuRange(0.0, 100.0);
    fSettings.DefineAllowedTargets("Ar");
    fSettings.DefineAllowedSpecies("numu");

    fSettings.SetTitle("");
    fSettings.SetXTitle("p_{T} [GeV/#it{c}]");
    fSettings.SetYTitle("p_{z} [GeV/#it{c}]");
    fSettings.SetZTitle("#Sigma T_{prot} [GeV]");

    // END boilerplate

    // 3D binning
    std::vector<double> ptbins = {0, 0.15, 0.25, 0.33, 0.7, 1, 2.5}; // GeV
    std::vector<double> pzbins = {1.5, 3.5, 4.5, 7, 8, 20};          // GeV
    std::vector<double> sumTpbins = {
        0,   40,  80,  120, 160, 200,  240,  280,
        320, 360, 400, 600, 800, 1000, 10000}; // MeV

    // This histogram is just used to help with the binning, we could manually
    // write the bin-mapping function ourselves
    f3DHist =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike", "", ptbins.size() - 1,
                               ptbins.data(), pzbins.size() - 1, pzbins.data(),
                               sumTpbins.size() - 1, sumTpbins.data());

    f3DHist_CCQE = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_CCQE", "", ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), sumTpbins.size() - 1,
        sumTpbins.data());
    f3DHist_CC1pi = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_CC1pi", "", ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), sumTpbins.size() - 1,
        sumTpbins.data());
    f3DHist_CCDIS = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_CCDIS", "", ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), sumTpbins.size() - 1,
        sumTpbins.data());

    // this is just to appease NUISANCE, it is copied to make the tracked 'MC'
    // histogram
    fDataHist = new TH1D("DUNE_3D_MINERvALike_data", "", f3DHist->GetNcells(),
                         0, f3DHist->GetNcells());

    // more boilerplate
    FinaliseSampleSettings();
    fScaleFactor =
        (GetEventHistogram()->Integral("width") * 1E-38 / (fNEvents + 0.)) /
        this->TotalIntegratedFlux();
    FinaliseMeasurement();
  };

  //********************************************************************
  void FillEventVariables(FitEvent *event) {
    //********************************************************************
    // Checking to see if there is a Muon
    if (event->NumFSParticle(13) == 0)
      return;

    // Get the muon kinematics
    TLorentzVector Pmu = event->GetHMFSParticle(13)->fP;
    TLorentzVector Pnu = event->GetNeutrinoIn()->fP;

    Double_t px = Pmu.X() / 1000;
    Double_t py = Pmu.Y() / 1000;
    Double_t pt = sqrt(px * px + py * py);

    // Don't want to assume the event generators all have neutrino coming along
    // z pz is muon momentum projected onto the neutrino direction
    Double_t pz = Pmu.Vect().Dot(Pnu.Vect() * (1.0 / Pnu.Vect().Mag())) / 1000.;
    // Set Hist Variables

    // Sum up kinetic energy of protons
    double sum_TProt = 0.0;
    for (auto prot : event->GetAllFSProton()) {
      sum_TProt += prot->KE();
    }

    // find the bin number along each axis
    int binx = f3DHist->GetXaxis()->FindFixBin(pt);
    int biny = f3DHist->GetYaxis()->FindFixBin(pz);
    int binz = f3DHist->GetZaxis()->FindFixBin(sum_TProt);

    // set this as the global bin number, could also use
    // f3DHist->FindFixBin(pt,pz,sum_TProt)
    fXVar = f3DHist->GetBin(binx, biny, binz);

    f3DHist->Fill(pt, pz, sum_TProt);

    // Example of EventIsCCQE as a lambda
    // auto EventIsCCQE = [](FitEvent *event){ return (event->Mode == 1); };

    // could implement
    // if(EventIsCCQE(event)){
    //   f3DHist_CCQE->Fill(pt,pz,sum_TProt);
    // } else if(EventIsCC1pi(event)){
    //   f3DHist_CC1pi->Fill(pt,pz,sum_TProt);
    // } else if(EventIsCCDIS(event)){
    //   f3DHist_CCDIS->Fill(pt,pz,sum_TProt);
    // }
  }

  //********************************************************************
  bool isSignal(FitEvent *event) {
    //********************************************************************
    return SignalDef::isCC0pi_MINERvAPTPZ(event, 14, EnuMin, EnuMax);
  }

  void Write(std::string drawOpt) {

    f3DHist->Write();
    f3DHist->SetDirectory(nullptr);

    for (int x = 0; x < f3DHist->GetXaxis()->GetNbins(); ++x) {
      for (int y = 0; y < f3DHist->GetYaxis()->GetNbins(); ++y) {
        TH1 *DUNE_3D_MINERvALike_proj =
            f3DHist->ProjectionZ((std::string("DUNE_3D_MINERvALike_SumTp_x") +
                                  std::to_string(x) + "_y" + std::to_string(y))
                                     .c_str(),
                                 x + 1, x + 1, y + 1, y + 1);

        std::stringstream ss;
        ss << f3DHist->GetXaxis()->GetBinLowEdge(x + 1) << " < p_t < "
           << f3DHist->GetXaxis()->GetBinUpEdge(x + 1) << ", "
           << f3DHist->GetYaxis()->GetBinLowEdge(y + 1) << " < p_t < "
           << f3DHist->GetYaxis()->GetBinUpEdge(y + 1);

        DUNE_3D_MINERvALike_proj->SetTitle(ss.str().c_str());

        DUNE_3D_MINERvALike_proj->Write();
        DUNE_3D_MINERvALike_proj->SetDirectory(nullptr);
      }
    }

    // we have to tidy this up in this SO if we don't want horrible crashes on
    // program tear down
    fDataHist->Write();
    fDataHist->SetDirectory(nullptr);
    delete fDataHist;
    fDataHist = nullptr;
    fMCHist->Write();
    fMCHist->SetDirectory(nullptr);
    delete fMCHist;
    fMCHist = nullptr;
  }
};
