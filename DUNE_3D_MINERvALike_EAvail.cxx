#include "InteractionModes.h"
#include "/root/software/nuisance_version2/nuisance/src/MINERvA/MINERvA_SignalDef.h"
#include "Measurement1D.h"

#include "TH3D.h"

class DUNE_3D_MINERvALike_EAvail: public Measurement1D {
public:
  std::unique_ptr<TH3D> f3DHist;

  std::unique_ptr<TH3D> f3DHist_CCQE;
  std::unique_ptr<TH3D> f3DHist_CC2p2h;
  std::unique_ptr<TH3D> f3DHist_CC1pi;
  std::unique_ptr<TH3D> f3DHist_CCDIS;
  std::unique_ptr<TH3D> f3DHist_Other;

  std::unique_ptr<TH3D> f3DHist_EAvail;  //pointers to the new histograms for the new EAvailable binning
  std::unique_ptr<TH3D> f3DHist_EAvail_CCQE;
  std::unique_ptr<TH3D> f3DHist_EAvail_CC2p2h;
  std::unique_ptr<TH3D> f3DHist_EAvail_CC1p1pi;
  std::unique_ptr<TH3D> f3DHist_EAvail_CCDIS;
  std::unique_ptr<TH3D> f3DHist_EAvail_Other;

  //********************************************************************
  DUNE_3D_MINERvALike_EAvail(nuiskey samplekey){
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
    std::vector<double> ptbins = {0,     0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55,  0.7,  0.85, 1.0,   2.5}; // GeV
    std::vector<double> pzbins = {1.5, 3.5, 4.5, 7.0, 8.0, 10.0, 20.0};  // GeV
    std::vector<double> sumTpbins = {0,    0.02, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.4,  0.6,  0.8}; // GeV
    std::vector<double> EAvail_bins = {0.04, 0.08}; // GeV


    // This histogram is just used to help with the binning, we could manually
    // write the bin-mapping function ourselves
    f3DHist =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike", "", ptbins.size() - 1,
                               ptbins.data(), pzbins.size() - 1, pzbins.data(),
                               sumTpbins.size() - 1, sumTpbins.data());


    f3DHist_EAvail =
        std::make_unique<TH3D>("f3DHist_f3DHist_EAvail", "", ptbins.size() - 1,
                               ptbins.data(), pzbins.size() - 1, pzbins.data(),
                               EAvail_bins.size() - 1, EAvail_bins.data());

    f3DHist_CCQE = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_CCQE", "", ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), sumTpbins.size() - 1,
        sumTpbins.data());
    f3DHist_CC2p2h = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_CC2p2h", "", ptbins.size() - 1, ptbins.data(),
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
    f3DHist_Other = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_CCOther", "", ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), sumTpbins.size() - 1,
        sumTpbins.data());


    ///////////////////////
    
    f3DHist_EAvail_CCQE =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CCQE", "",ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data()) ;

    f3DHist_EAvail_CC2p2h =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC2p2h", "", ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());

    f3DHist_EAvail_CC1p1pi =
        std::make_unique<TH3D>("f3DHist_EAvail_CC1p1pi", "", ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());
    
    f3DHist_EAvail_CCDIS = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_EAvail_CCDIS", "", ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());
    f3DHist_EAvail_Other = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_EAvail_Other", "", ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());


    // this is just to appease NUISANCE, it is copied to make the tracked 'MC'
    // histogram
    fDataHist = new TH1D("DUNE_3D_MINERvALike_data", "", f3DHist->GetNcells(),
                         0, f3DHist->GetNcells());

    fDataHist = new TH1D("DUNE_3D_MINERvALike_data_EAvail", "", f3DHist_EAvail->GetNcells(),
                         0, f3DHist_EAvail->GetNcells());

    // more boilerplate
    FinaliseSampleSettings();
    fScaleFactor =
        (GetEventHistogram()->Integral("width") * 1E-38 / (fNEvents + 0.)) /
        this->TotalIntegratedFlux();
    FinaliseMeasurement();
  };

  //********************************************************************
  void FillEventVariables(FitEvent *event){
    //********************************************************************
    // Checking to see if there is a Muon
    if (event->NumFSParticle(13) == 0)
      return;

    auto d3dml_bx = dynamic_cast<D3DMLBox *>(GetBox());
    d3dml_bx->mode = event->Mode;    
    // Get the muon kinematics
    TLorentzVector Pmu = event->GetHMFSParticle(13)->fP;
    TVector3 nudir = event->GetNeutrinoIn()->fP.Vect().Unit();

      // Set Defaults
    double Eav = -999.9;
    double q3 = -999.9;

    // If muon found get kinematics
    FitParticle* muon      = event->GetHMFSParticle(13);
    FitParticle* neutrino  = event->GetNeutrinoIn();
    if (muon && neutrino) {

    // Set Q from Muon
    TLorentzVector q = neutrino->fP - muon->fP;
    double q0 = (q.E()) / 1.E3;
    //double q3_true = (q.Vect().Mag())/1.E3;
    double thmu = muon->fP.Vect().Angle(neutrino->fP.Vect());
    double pmu  = muon->fP.Vect().Mag() / 1.E3;
    double emu  = muon->fP.E() / 1.E3;
    double mmu  = muon->fP.Mag() / 1.E3;

    //// Get Enu Rec
    double enu_rec = emu + q0;

    // Set Q2 QE
    double q2qe = 2 * enu_rec * (emu - pmu * cos(thmu)) - mmu * mmu;


    // Calc Q3 from Q2QE and EnuTree
    q3 = sqrt(q2qe + q0 * q0);

    // Get Eav too
    Eav = FitUtils::GetErecoil_MINERvA_LowRecoil(event) / 1.E3;
    //Put the available energy into a box function
    d3dml_bx->E_Avail = Eav;
    }

    static const double toGeV = 1E-3;

    d3dml_bx->p_para = Pmu.Vect().Dot(nudir) * toGeV;
    d3dml_bx->p_perp = Pmu.Vect().Cross(nudir).Mag() * toGeV;

    // Sum up kinetic energy of protons
    d3dml_bx->sum_TProt = 0.0;
    for (auto prot : event->GetAllFSProton()) {
      d3dml_bx->sum_TProt += prot->KE() * toGeV;
    }
    

    // find the bin number along each axis
    int binx = f3DHist->GetXaxis()->FindFixBin(d3dml_bx->p_perp);
    int biny = f3DHist->GetYaxis()->FindFixBin(d3dml_bx->p_para);
    int binz = f3DHist->GetZaxis()->FindFixBin(d3dml_bx->sum_TProt);
    int binz_EAvail = f3DHist_EAvail->GetZaxis()->FindFixBin(d3dml_bx->E_Avail);

    // set this as the global bin number, could also use
    // f3DHist->FindFixBin(pt,pz,sum_TProt)
    f3DHist_EAvail->GetBin(binx, biny, binz_EAvail);
    fXVar = f3DHist->GetBin(binx, biny, binz);
  }

  struct D3DMLBox : public MeasurementVariableBox1D {
    D3DMLBox()
        : MeasurementVariableBox1D(), mode{0},p_para{0}, p_perp{0}, sum_TProt{0}, E_Avail{0} {}

    MeasurementVariableBox *CloneSignalBox() {
      auto cl = new D3DMLBox();

      cl->fX = fX;
      cl->mode = mode;
      cl->p_para = p_para;
      cl->p_perp = p_perp;
      cl->sum_TProt = sum_TProt;
      cl->E_Avail = E_Avail;
      return cl;
    }

    void Reset() {
      MeasurementVariableBox1D::Reset();
      mode = 0;
      p_para = 0;
      p_perp = 0;
      sum_TProt = 0;
      E_Avail = 0;
    }

    int mode;
    double p_para;
    double p_perp;
    double sum_TProt;
    double E_Avail;
  };
  
  MeasurementVariableBox *CreateBox() { return new D3DMLBox(); };

  void FillExtraHistograms(MeasurementVariableBox *vars, double weight) {
    Measurement1D::FillExtraHistograms(vars, weight);

    auto d3dml_bx = dynamic_cast<D3DMLBox *>(vars);

    f3DHist->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                         d3dml_bx->sum_TProt, weight);
    
    f3DHist_EAvail->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                         d3dml_bx->E_Avail, weight);

    int amode = std::abs(d3dml_bx->mode);

    if (amode == InputHandler::kCCQE) {
      f3DHist_CCQE->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                         d3dml_bx->sum_TProt, weight);
      f3DHist_EAvail_CCQE->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
      d3dml_bx->E_Avail, weight);
    } else if (amode == InputHandler::kCC2p2h) {
      f3DHist_CC2p2h->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                           d3dml_bx->sum_TProt, weight);
      f3DHist_EAvail_CC2p2h->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
      d3dml_bx->E_Avail, weight);
    } else if ((amode == InputHandler::kCC1piponp) ||
               (amode == InputHandler::kCC1pi0onn) ||
               (amode == InputHandler::kCC1piponn)) {
      f3DHist_CC1pi->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                          d3dml_bx->sum_TProt, weight);
      f3DHist_EAvail_CC1p1pi->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
      d3dml_bx->E_Avail, weight);
    } else if ((amode == InputHandler::kCCmultipi) ||
               (amode == InputHandler::kCCDIS)) {
      f3DHist_CCDIS->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                          d3dml_bx->sum_TProt, weight);
      f3DHist_EAvail_CCDIS->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
      d3dml_bx->E_Avail, weight);
    } else {
      f3DHist_Other->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                          d3dml_bx->sum_TProt, weight);
      f3DHist_EAvail_Other->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
      d3dml_bx->E_Avail, weight);
    }
  }
 

  //********************************************************************
  bool isSignal(FitEvent *event) {
    //********************************************************************
    //return SignalDef::isCC0pi_MINERvAPTPZ(event, 14, EnuMin, EnuMax);
    return SignalDef::isCCincLowRecoil_MINERvA(event, EnuMin, EnuMax);
  }

  void SplitWrite3D(std::unique_ptr<TH3D> const &h) {
    h->Write();
    h->SetDirectory(nullptr);

    for (int x = 0; x < h->GetXaxis()->GetNbins(); ++x) {
      for (int y = 0; y < h->GetYaxis()->GetNbins(); ++y) {
        TH1 *proj =
            h->ProjectionZ((std::string(h->GetName()) + "_x" +
                            std::to_string(x) + "_y" + std::to_string(y))
                               .c_str(),
                           x + 1, x + 1, y + 1, y + 1, "e");

        std::stringstream ss;
        ss << h->GetXaxis()->GetBinLowEdge(x + 1) << " < p_t < "
           << h->GetXaxis()->GetBinUpEdge(x + 1) << " [GeV/c], "
           << h->GetYaxis()->GetBinLowEdge(y + 1) << " < p_z < "
           << h->GetYaxis()->GetBinUpEdge(y + 1) << " [GeV/c]";

        double proj_cell_area = (h->GetXaxis()->GetBinUpEdge(x + 1) -
                                 h->GetXaxis()->GetBinLowEdge(x + 1)) *
                                (h->GetYaxis()->GetBinUpEdge(y + 1) -
                                 h->GetYaxis()->GetBinLowEdge(y + 1));

        proj->Scale(fScaleFactor / proj_cell_area, "WIDTH");

        proj->SetTitle(ss.str().c_str());
        proj->GetXaxis()->SetTitle("#Sigma T_{p} [GeV] or Available Energy ");

        proj->Write();
        proj->SetDirectory(nullptr);
        delete proj;
      }
    }
  }

  void Write(std::string drawOpt) {

    SplitWrite3D(f3DHist);
    SplitWrite3D(f3DHist_CCQE);
    SplitWrite3D(f3DHist_CC2p2h);
    SplitWrite3D(f3DHist_CC1pi);
    SplitWrite3D(f3DHist_CCDIS);
    SplitWrite3D(f3DHist_Other);

    SplitWrite3D(f3DHist_EAvail);
    SplitWrite3D(f3DHist_EAvail_CCQE);
    SplitWrite3D(f3DHist_EAvail_CC2p2h);
    SplitWrite3D(f3DHist_EAvail_CC1p1pi);
    SplitWrite3D(f3DHist_EAvail_CCDIS);
    SplitWrite3D(f3DHist_EAvail_Other);

    // we have to tidy this up in this SO if we don't want horrible crashes on
    // program tear down
    fDataHist->Write();
    fDataHist->SetDirectory(nullptr);
    fMCHist->Write();
    fMCHist->SetDirectory(nullptr);
  }

  void ResetAll() {
    Measurement1D::ResetAll();
    f3DHist->Reset();
    f3DHist_CCQE->Reset();
    f3DHist_CC2p2h->Reset();
    f3DHist_CC1pi->Reset();
    f3DHist_CCDIS->Reset();
    f3DHist_Other->Reset();

    f3DHist_EAvail->Reset();
    f3DHist_EAvail_CCQE->Reset();
    f3DHist_EAvail_CC2p2h->Reset();
    f3DHist_EAvail_CC1p1pi->Reset();
    f3DHist_EAvail_CCDIS->Reset();
    f3DHist_EAvail_Other->Reset();
  }
};