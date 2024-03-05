#include "InteractionModes.h"
#include "/root/software/nuisance_version2/nuisance/src/MINERvA/MINERvA_SignalDef.h"
#include "Measurement1D.h"
#include "/root/software/nuisance_version2/nuisance/src/MCStudies/GenericFlux_Vectors.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH2.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;


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

  std::unique_ptr<TH3D> f3DHist_EAvail_lowW;  //pointers to the new histograms for the new EAvailable binning
  std::unique_ptr<TH3D> f3DHist_EAvail_midW;
  std::unique_ptr<TH3D> f3DHist_EAvail_highW;


  /* a function to generate numpy linspace */
//template <typename T>
  template <typename T>
    std::vector<T> linspace(T a, T b, size_t N) {
        T h = (b - a) / static_cast<T>(N-1);
        std::vector<T> xs(N);
        typename std::vector<T>::iterator x;
        T val;
        for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
            *x = val;
        return xs;
    }

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

    std::vector<double> pzbins = linspace(0.5,20.0,399);

    std::vector<double> ptbins =linspace(0.0, 2.5,125) ;


    /*
    std::vector<double> ptbins = {0,     0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55,  0.7,  0.85, 1.0,   2.5}; // GeV
    std::vector<double> pzbins = {1.5, 3.5, 4.5, 7.0, 8.0, 10.0, 20.0};  // GeV
    std::vector<double> sumTpbins = {0,    0.02, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.4,  0.6,  0.8}; // GeV
    std::vector<double> EAvail_bins = {0.04, 0.08}; // GeV
    */
   // std::vector<double> ptbins = {0.4, 0.475}; // GeV

   /*
    std::vector<double> ptbins = {0.00, 0.02,0.04, 0.06,0.08, 
                                  0.1, 0.12,0.14 ,0.16,0.18, 
                                  0.2, 0.22,0.24,0.26, 0.28, 
                                  0.3,0.32,0.34 ,0.36,0.38, 
                                  0.4,0.42, 0.44, 0.46,0.48, 
                                  0.5,0.52, 0.54, 0.56,0.58, 
                                  0.6, 0.62, 0.64, 0.66,0.68, 
                                  0.7,0.72, 0.74, 0.76,0.78, 
                                  0.8,0.82, 0.84, 0.86,0.88, 
                                  0.9,0.92, 0.94, 0.96,0.98, 
                                  1.0, 1.1,1.2,1.3,1.4,1.5, 1.6,1.7,1.8,1.9,
                                  2.0, 2.1,2.2,2.3,2.4,2.5, 2.6,2.7,2.8,2.9,
                                }; // GeV */
    //std::vector<double> pzbins = {1.5, 3.5};  // GeV
    //std::vector<double> pzbins = {0.5,1.0,1.5,2,2.5,3, 3.5, 4, 4.5,5.0,6.0, 7.0, 8.0,9.0,10.0, 12.5, 15.0,17.5, 20.0};  // GeV
    /*
   std::vector<double> pzbins={0.5,  0.55 ,0.6 , 0.65 ,0.7 , 0.75, 0.8 , 0.85, 0.9 , 0.95, 1. ,  1.05 ,1.1,  1.15,
                                1.2,  1.25 ,1.3 , 1.35 ,1.4 , 1.45, 1.5 , 1.55, 1.6 , 1.65, 1.7,  1.75 ,1.8,  1.85,
                                1.9 , 1.95, 2.0,   2.05, 2.1 , 2.15, 2.2 , 2.25, 2.3 , 2.35, 2.4,  2.45, 2.5,  2.55,
                                2.6  ,2.65, 2.7 , 2.75, 2.8 , 2.85, 2.9,  2.95, 3. ,  3.05, 3.1,  3.15, 3.2,  3.25,
                                3.3 , 3.35, 3.4 , 3.45, 3.5 , 3.55, 3.6 , 3.65, 3.7 , 3.75, 3.8,  3.85, 3.9,  3.95,
                                4.   ,4.05, 4.1 , 4.15, 4.2 , 4.25, 4.3,  4.35, 4.4 ,4.45, 4.5,  4.55, 4.6,  4.65,
                                4.7 , 4.75, 4.8 , 4.85, 4.9 , 4.95, 5.0,   5.05, 5.1 , 5.15, 5.2,  5.25, 5.3,  5.35,
                                5.4 , 5.45, 5.5 , 5.55, 5.6 , 5.65, 5.7,  5.75, 5.8 , 5.85, 5.9,  5.95, 6.0,   6.05,
                                6.1 , 6.15, 6.2 , 6.25, 6.3 , 6.35, 6.4 , 6.45, 6.5,  6.55, 6.6,  6.65, 6.7, 6.75,
                                6.8 , 6.85, 6.9 , 6.95, 7.  , 7.05, 7.1  ,7.15, 7.2,  7.25, 7.3,  7.35, 7.4,  7.45,
                                7.5 , 7.55, 7.6 , 7.65, 7.7 , 7.75, 7.8 , 7.85, 7.9 , 7.95, 8.0,  8.05, 8.1,  8.15,
                                8.2 , 8.25, 8.3 , 8.35, 8.4 , 8.45, 8.5 , 8.55, 8.6 , 8.65, 8.7,  8.75, 8.8,  8.85,
                                8.9 , 8.95, 9.0,   9.05, 9.1,  9.15, 9.2,  9.25, 9.3 , 9.35, 9.4,  9.45, 9.5,  9.55,
                                9.6 , 9.65, 9.7 , 9.75, 9.8 , 9.85 ,9.9  ,9.95,10, 10.2,
                                10.4,10.6,10.8,12,12.2,12.4,12.6,12.8,13,14,15,16,17,18,19,20}; */

   // std::vector<double> ptbins = {0,  0.05,   0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,0.55,0.6,0.65 , 0.7,0.75,0.8 ,0.9, 1.0, 1.5,2.0,  2.5}; // GeV
  //  std::vector<double> ptbins = {0, 0.025, 0.05,0.075,   0.1,0.125,  0.15, 0.175, 0.2,0.225, 0.25,0.275, 0.3, 0.325, 0.35,0.375, 0.4,0.425, 0.45,0.475,  0.5, 0.525, 0.55,0.575, 0.6,0.625, 0.65, 0.675, 0.7,0.725,0.750,0.775,0.8,0.825 ,0.85,0.875,0.9, 0.925, 0.95, 0.975, 1.0, 1.05, 1.1, 1.15,  1.2 , 1.25 ,1.3, 1.35,1.4, 1.45,  1.5,1.55, 1.6, 1.65, 1.7, 1.75,1.8, 1.85, 1.9, 1.95,2.0, 2.1,2.2,2.3,2.4, 2.5}; 
    //std::vector<double> pzbins = {0.5,1.0,1.5,2,2.5,3, 3.5, 4, 4.5,5.0,6.0, 7.0, 8.0,9.0,10.0, 12.5, 15.0,17.5, 20.0};  // GeV

   /* std::vector<double> pzbins = {0.5,0.75, 1.0, 1.25, 1.5, 1.75, 2,2.25,2.5,2.75, 3, 3.25, 3.5, 3.75,  4, 4.25, 4.5,4.75, 5.0, 5.25, 
                                  5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75,
                                  10.0, 10.25, 10.5, 10.75, 11, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13, 13.25, 13.5, 13.75, 14.0, 
                                  15.0, 16.0, 17, 18.0, 19.0, 20.0};  // GeV*/


    //std::vector<double> pzbins(100) ; // vector with 100 ints.
    //std::iota (std::begin(v), std::end(v), 0); // Fill with 0, 1, ..., 99.


    std::vector<double> sumTpbins = {0,  0.02, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.4,  0.6}; // GeV
    std::vector<double> EAvail_bins = {0, 0.01,   0.02, 0.04,0.06,  0.08,0.1, 0.12,0.14, 0.16, 0.2, 0.24, 0.28, 0.32, 0.4, 0.5, 0.6,0.8}; // GeV

    // This histogram is just used to help with the binning, we could manually
    // write the bin-mapping function ourselves
    f3DHist =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike", "", ptbins.size() - 1,
                               ptbins.data(), pzbins.size() - 1, pzbins.data(),
                               sumTpbins.size() - 1, sumTpbins.data());


    f3DHist_EAvail =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail", "", ptbins.size() - 1,
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
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC1p1pi", "", ptbins.size() - 1, ptbins.data(),
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

    // Now the invariant mass histograms for available energy
    f3DHist_EAvail_lowW =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_lowW", "",ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data()) ;

    f3DHist_EAvail_midW =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_midW", "",ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data()) ;

    f3DHist_EAvail_highW =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_highW", "",ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data()) ;

    // this is just to appease NUISANCE, it is copied to make the tracked 'MC'
    // histogram
    //fDataHist = new TH1D("DUNE_3D_MINERvALike_data", "", f3DHist->GetNcells(),0, f3DHist->GetNcells());
     fDataHist = new TH1D("DUNE_3D_MINERvALike_data", "", 1, 0, 1);


    // more boilerplate
    FinaliseSampleSettings();
     //fScaleFactor =((GetEventHistogram()->Integral("width") * 1E-38) / (fNEvents + 0.)) /this->TotalIntegratedFlux();
    fScaleFactor =1.0;
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

    // Basic interaction kinematics
    double Q2 = -1 * (neutrino->fP - muon->fP).Mag2() / 1E6;
    //double  q0 = (neutrino->fP - muon->fP).E() / 1E3;
    // Get W_true with assumption of initial state nucleon at rest
    float m_n = (float)PhysConst::mass_proton;
    // Q2 assuming nucleon at rest
    double W_nuc_rest = sqrt(-Q2 + 2 * m_n * q0 + m_n * m_n);
    double  W = W_nuc_rest; // For want of a better thing to do
    //std::cout << " W = "  << W << "   Q2 = " << Q2 << "   q0 = " << q0 << "   m_n =   "  << m_n <<  "  E available  =  " << Eav << std::endl;
    d3dml_bx->W = W_nuc_rest;
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
        : MeasurementVariableBox1D(), mode{0},p_para{0}, p_perp{0}, sum_TProt{0}, E_Avail{0}, W{0} {}

    MeasurementVariableBox *CloneSignalBox() {
      auto cl = new D3DMLBox();

      cl->fX = fX;
      cl->mode = mode;
      cl->p_para = p_para;
      cl->p_perp = p_perp;
      cl->sum_TProt = sum_TProt;
      cl->E_Avail = E_Avail;
      cl->W = W;
      return cl;
    }

    void Reset() {
      MeasurementVariableBox1D::Reset();
      mode = 0;
      p_para = 0;
      p_perp = 0;
      sum_TProt = 0;
      E_Avail = 0;
      W = 0;
    }

    int mode;
    double p_para;
    double p_perp;
    double sum_TProt;
    double E_Avail;
    double W;
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
    double W = std::abs(d3dml_bx->W);

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

     //if statments to fill W hists
    if (W < 1.4) {
      f3DHist_EAvail_lowW->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
      d3dml_bx->E_Avail, weight);
      //std::cout << " W should be less than 1.4 , it is" << W <<std::endl;
    } else if (1.4 < W && W < 2.0) {
      f3DHist_EAvail_midW->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
      d3dml_bx->E_Avail, weight);
     // std::cout << " W should be 1.4 < W < 2.0 , it is" << W <<std::endl;
    } 
    else if (W > 2.0){
      f3DHist_EAvail_highW->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
      d3dml_bx->E_Avail, weight);
      //std::cout << " W should be > 2.0 , it is" << W <<std::endl;
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

    
    auto ptpzproj = std::unique_ptr<TH2>(static_cast<TH2*>(h->Project3D("xy")));
    

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
        ptpzproj->SetBinContent(x+1,y+1, ptpzproj->GetBinContent(x+1,y+1)/proj_cell_area); //divide by the x, y bin widths 
        proj->SetTitle(ss.str().c_str());
        proj->GetXaxis()->SetTitle("#Sigma T_{p} [GeV] or Available Energy ");

        proj->Write();
        proj->SetDirectory(nullptr);
        delete proj;
      }
    }

    ptpzproj->Write();
    //make sure to tell ROOT that it doesn't own this histogram so that we can delete it
    ptpzproj->SetDirectory(nullptr);

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

    SplitWrite3D(f3DHist_EAvail_lowW);
    SplitWrite3D(f3DHist_EAvail_midW);
    SplitWrite3D(f3DHist_EAvail_highW);
    


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

    f3DHist_EAvail_lowW->Reset();
    f3DHist_EAvail_midW->Reset();
    f3DHist_EAvail_highW->Reset();

  }
};