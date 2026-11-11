// #include "InteractionModes.h"
// #include "/home/abipeake/scratch/Formfactorstuff/nuisance/src/MINERvA/MINERvA_SignalDef.h"
// #include "Measurement1D.h"
// #include "/home/abipeake/scratch/Formfactorstuff/nuisance/src/MCStudies/GenericFlux_Vectors.h"
// #include "TH1D.h"
// #include "TH3D.h"
// #include "TH2D.h"
// #include "TH2.h"
// #include <fstream>
// #include <functional>
// #include <iostream>
// #include <numeric>
// #include <stdlib.h>
// #include <string>
// #include <vector>


// using namespace std;


// class DUNE_3D_MINERvALike_EAvail: public Measurement1D {
// public:
//   std::unique_ptr<TH3D> f3DHist;
//   std::unique_ptr<TH1D> f1DHist_EAvail_true_nuenergy;

//   std::unique_ptr<TH3D> f3DHist_CCQE;
//   std::unique_ptr<TH3D> f3DHist_CC2p2h;
//   std::unique_ptr<TH3D> f3DHist_CC1pi;
//   std::unique_ptr<TH3D> f3DHist_CCDIS;
//   std::unique_ptr<TH3D> f3DHist_Other;

//   std::unique_ptr<TH3D> f3DHist_EAvail;  //pointers to the new histograms for the new EAvailable binning
//   std::unique_ptr<TH3D> f3DHist_EAvail_CCQE;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CC2p2h;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CC1p1pi;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CCDIS;
//   std::unique_ptr<TH3D> f3DHist_EAvail_Other;


//   std::unique_ptr<TH2D> f2DHist_EAvail;  //pointers to the new histograms for the new EAvailable binning
//   std::unique_ptr<TH2D> f2DHist_EAvail_CCQE;
//   std::unique_ptr<TH2D> f2DHist_EAvail_CC2p2h;
//   std::unique_ptr<TH2D> f2DHist_EAvail_CC1p1pi;
//   std::unique_ptr<TH2D> f2DHist_EAvail_CCDIS;
//   std::unique_ptr<TH2D> f2DHist_EAvail_Other;

  
//   std::unique_ptr<TH3D> f3DHist_EAvail_pzbins_double;  //pointers to the new histograms for the new EAvailable binning
//   std::unique_ptr<TH3D> f3DHist_EAvail_CCQE_pzbins_double;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CC2p2h_pzbins_double;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CC1p1pi_pzbins_double;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CCDIS_pzbins_double;
//   std::unique_ptr<TH3D> f3DHist_EAvail_Other_pzbins_double;

//   /*
//   std::unique_ptr<TH3D> f3DHist_EAvail_pzbins_uniform; 
//   std::unique_ptr<TH3D> f3DHist_EAvail_pzbins_uniform_double;  //pointers to the new histograms for the new EAvailable binning
//   std::unique_ptr<TH3D> f3DHist_EAvail_CCQE_pzbins_uniform;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CC2p2h_pzbins_uniform;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CC1p1pi_pzbins_uniform;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CCDIS_pzbins_uniform;
//   std::unique_ptr<TH3D> f3DHist_EAvail_Other_pzbins_uniform;

//    std::unique_ptr<TH3D> f3DHist_EAvail_CCQE_pzbins_uniform_double;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CC2p2h_pzbins_uniform_double;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CC1p1pi_pzbins_uniform_double;
//   std::unique_ptr<TH3D> f3DHist_EAvail_CCDIS_pzbins_uniform_double;
//   std::unique_ptr<TH3D> f3DHist_EAvail_Other_pzbins_uniform_double;
//   */
//   std::unique_ptr<TH3D> f3DHist_EAvail_lowW;  //pointers to the new histograms for the new EAvailable binning
//   std::unique_ptr<TH3D> f3DHist_EAvail_midW;
//   std::unique_ptr<TH3D> f3DHist_EAvail_highW;


//   /* a function to generate numpy linspace */
// //template <typename T>
//   template <typename T>
//     std::vector<T> linspace(T a, T b, size_t N) {
//         T h = (b - a) / static_cast<T>(N);
//         std::vector<T> xs(N);
//         typename std::vector<T>::iterator x;
//         T val;
//         for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
//             *x = val;
//         return xs;
//     }

//   void ExtendLinspace(std::vector<double> &in, double low, double high,size_t nbins){
//     for(auto x : linspace(low, high,nbins)){
//         in.push_back(x);
        
//       }
//       auto x =std::unique(in.begin(),in.end());  //put duplicate bins at end
//       in.erase(x, in.end()); //erase duplicates
//   }

//   //********************************************************************
//   DUNE_3D_MINERvALike_EAvail(nuiskey samplekey){
//     //********************************************************************

//     // START boilerplate
//     //
//     // Setup common settings
//     fSettings = LoadSampleSettings(samplekey);
//     fSettings.SetAllowedTypes("FIX,FREE,SHAPE/FULL,DIAG/MASK", "FIX/FULL");
//     fSettings.SetEnuRange(0.0, 100.0);
//     fSettings.DefineAllowedTargets("Ar");
//     fSettings.DefineAllowedSpecies("numu");

//     fSettings.SetTitle("");
//     fSettings.SetXTitle("p_{T} [GeV/#it{c}]");
//     fSettings.SetYTitle("p_{z} [GeV/#it{c}]");
//     fSettings.SetZTitle("#Sigma T_{prot} [GeV]");
//     // END boilerplate

//     // 3D binning
//     //std::vector<double> ptbins =linspace(0.0, 2.5,20) ;

//     std::vector<double> ptbins;
//     ExtendLinspace(ptbins,1e-8,0.2,4);
//     ExtendLinspace(ptbins,0.2,0.4,5);
//     ExtendLinspace(ptbins,0.4,0.7,4);
//     ExtendLinspace(ptbins,0.7,2.5,8);
//   ptbins.push_back(2.5);
  
//     std::cout<<"pt vector is: "<<std::endl;
//     for(int i=0 ; i< 21 ; i++){
//         std::cout<< ptbins[i]<<std::endl;

//     }

//     std::vector<double> pzbins_uniform = linspace(1.0e-8,6.6,22); //linearly spaced bins
//     std::vector<double> pzbins_uniform_doublebins = linspace(1.0e-8,6.6,44); //linearly spaced bins
    
//     std::vector<double> pzbins;
//     ExtendLinspace(pzbins,1e-8,1.5,5);
//     ExtendLinspace(pzbins,1.5,1.8,2);
//     ExtendLinspace(pzbins,1.8,2.3,7);
//     ExtendLinspace(pzbins,2.3,6.6,10);
//   pzbins.push_back(6.6);
    

//     std::vector<double> pzbins_double;
//     ExtendLinspace(pzbins_double,1e-8,1.5,10);
//     ExtendLinspace(pzbins_double,1.5,1.8,4);
//     ExtendLinspace(pzbins_double,1.8,2.3,14);
//     ExtendLinspace(pzbins_double,2.3,6.6,20);
    
//     std::cout<<"pz vector is: "<<std::endl;
//     for(int i=0 ; i< pzbins.size(); i++){
//         std::cout<< pzbins[i]<<std::endl;
//     }

//     std::cout<<"pzbins= double is: "<<std::endl;
//     for(int i=0 ; i< pzbins_double.size() ; i++){
//         std::cout<< pzbins_double[i]<<std::endl;

//     }

   
//     std::cout<<"pzbins_uniform is: "<<std::endl;
//     for(int i=0 ; i<pzbins_uniform.size() ; i++){
//         std::cout<< pzbins_uniform[i]<<std::endl;
//     }
//      /*
//     std::cout<<"pzbins_uniform double is: "<<std::endl;
//     for(int i=0 ; i< pzbins_uniform_doublebins.size() ; i++){
//         std::cout<< pzbins_uniform_doublebins[i]<<std::endl;
//     }
//     */
    
//     std::vector<double> sumTpbins = {1e-8,  0.02, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.4,  0.6}; // GeV
//     std::vector<double> EAvail_bins = {1e-8, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.2, 0.24, 0.28, 0.32, 0.4, 0.5, 0.6, 0.8, 1.2, 1.6, 2.0, 4.0, 8.0}; // GeV
    
//     // This histogram is just used to help with the binning, we could manually
//     // write the bin-mapping function ourselves
//     f3DHist =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike", "", ptbins.size() - 1,
//                                ptbins.data(), pzbins.size() - 1, pzbins.data(),
//                                sumTpbins.size() - 1, sumTpbins.data());


//     f3DHist_EAvail =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail", "", ptbins.size() - 1,
//                                ptbins.data(), pzbins.size() - 1, pzbins.data(),
//                                EAvail_bins.size() - 1, EAvail_bins.data());

    

//     f3DHist_EAvail_pzbins_double =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_pzbins_double", "", ptbins.size() - 1,
//                                ptbins.data(), pzbins_double.size() - 1, pzbins_double.data(),
//                                EAvail_bins.size() - 1, EAvail_bins.data());

//     /*
//     f3DHist_EAvail_pzbins_uniform =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_pzbins_uniform", "", ptbins.size() - 1,
//                                ptbins.data(), pzbins_uniform.size() - 1, pzbins_uniform.data(),
//                                EAvail_bins.size() - 1, EAvail_bins.data());

//     f3DHist_EAvail_pzbins_uniform_double =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_pzbins_uniform_double", "", ptbins.size() - 1,
//                                ptbins.data(), pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(),
//                                EAvail_bins.size() - 1, EAvail_bins.data());
    
//     */

//     f3DHist_CCQE = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_CCQE", "", ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), sumTpbins.size() - 1,
//         sumTpbins.data());
//     f3DHist_CC2p2h = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_CC2p2h", "", ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), sumTpbins.size() - 1,
//         sumTpbins.data());
//     f3DHist_CC1pi = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_CC1pi", "", ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), sumTpbins.size() - 1,
//         sumTpbins.data());
//     f3DHist_CCDIS = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_CCDIS", "", ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), sumTpbins.size() - 1,
//         sumTpbins.data());
//     f3DHist_Other = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_CCOther", "", ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), sumTpbins.size() - 1,
//         sumTpbins.data());


//     ///////////////////////
    
//     f3DHist_EAvail_CCQE =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CCQE", "",ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data()) ;

//     f3DHist_EAvail_CC2p2h =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC2p2h", "", ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());

//     f3DHist_EAvail_CC1p1pi =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC1p1pi", "", ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());
    
//     f3DHist_EAvail_CCDIS = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_EAvail_CCDIS", "", ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());
//     f3DHist_EAvail_Other = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_EAvail_Other", "", ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());


//     /*
//     //uniform binning
//     f3DHist_EAvail_CCQE_pzbins_uniform =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CCQE_pzbins_uniform", "",ptbins.size() - 1, ptbins.data(),
//         pzbins_uniform.size() - 1, pzbins_uniform.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data()) ;

//     f3DHist_EAvail_CC2p2h_pzbins_uniform =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC2p2h_pzbins_uniform", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_uniform.size() - 1, pzbins_uniform.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());

//     f3DHist_EAvail_CC1p1pi_pzbins_uniform =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC1p1pi_pzbins_uniform", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_uniform.size() - 1, pzbins_uniform.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());
    
//     f3DHist_EAvail_CCDIS_pzbins_uniform = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_EAvail_CCDI_pzbins_uniformS", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_uniform.size() - 1, pzbins_uniform.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());
//     f3DHist_EAvail_Other_pzbins_uniform = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_EAvail_Other_pzbins_uniform", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_uniform.size() - 1, pzbins_uniform.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());

//     //double uniform binning
//     //uniform_double binning
//     f3DHist_EAvail_CCQE_pzbins_uniform_double =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CCQE_pzbins_uniform_double", "",ptbins.size() - 1, ptbins.data(),
//         pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data()) ;

//     f3DHist_EAvail_CC2p2h_pzbins_uniform_double =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC2p2h_pzbins_uniform_double", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());

//     f3DHist_EAvail_CC1p1pi_pzbins_uniform_double =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC1p1pi_pzbins_uniform_double", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());
    
//     f3DHist_EAvail_CCDIS_pzbins_uniform_double = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_EAvail_CCDI_pzbins_uniform_doubleS", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());
//     f3DHist_EAvail_Other_pzbins_uniform_double = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_EAvail_Other_pzbins_uniform_double", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());
//         */

//     //1D hist of Enu_true
//     f1DHist_EAvail_true_nuenergy =
//         std::make_unique<TH1D>("DUNE_3D_MINERvALike__EAvail_true_nuenergy", "", EAvail_bins.size() - 1,
//         EAvail_bins.data()) ;


//     //double binning  over peak
//     f3DHist_EAvail_CCQE_pzbins_double =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CCQE_pzbins_double", "",ptbins.size() - 1, ptbins.data(),
//         pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data()) ;


//     f3DHist_EAvail_CC2p2h_pzbins_double =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC2p2h_pzbins_double", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());

//     f3DHist_EAvail_CC1p1pi_pzbins_double =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC1p1pi_pzbins_double", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());
    
//     f3DHist_EAvail_CCDIS_pzbins_double = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_EAvail_CCDI_pzbins_doubleS", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());
//     f3DHist_EAvail_Other_pzbins_double = std::make_unique<TH3D>(
//         "DUNE_3D_MINERvALike_EAvail_Other_pzbins_double", "", ptbins.size() - 1, ptbins.data(),
//         pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data());
      

//     // Now the invariant mass histograms for available energy
//     f3DHist_EAvail_lowW =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_lowW", "",ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data()) ;

//     f3DHist_EAvail_midW =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_midW", "",ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data()) ;

//     f3DHist_EAvail_highW =
//         std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_highW", "",ptbins.size() - 1, ptbins.data(),
//         pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
//         EAvail_bins.data()) ;

//     // this is just to appease NUISANCE, it is copied to make the tracked 'MC'
//     // histogram
//     //fDataHist = new TH1D("DUNE_3D_MINERvALike_data", "", f3DHist->GetNcells(),0, f3DHist->GetNcells());
//      fDataHist = new TH1D("DUNE_3D_MINERvALike_data", "", 1, 0, 1);


//     // more boilerplate
//     FinaliseSampleSettings();
//      //fScaleFactor =((GetEventHistogram()->Integral("width") * 1E-38) / (fNEvents + 0.)) /this->TotalIntegratedFlux();
//     fScaleFactor =1.0;
//     FinaliseMeasurement();
//   };

//   //********************************************************************
//   void FillEventVariables(FitEvent *event){
//     //********************************************************************
//     // Checking to see if there is a Muon
//     if (event->NumFSParticle(13) == 0)
//       return;

//     auto d3dml_bx = dynamic_cast<D3DMLBox *>(GetBox());
//     d3dml_bx->mode = event->Mode;    
//     // Get the muon kinematics
//     TLorentzVector Pmu = event->GetHMFSParticle(13)->fP;
//     TVector3 nudir = event->GetNeutrinoIn()->fP.Vect().Unit();

//       // Set Defaults
//     double Eav = -999.9;
//     double q3 = -999.9;

//     // If muon found get kinematics
//     FitParticle* muon      = event->GetHMFSParticle(13);
//     FitParticle* neutrino  = event->GetNeutrinoIn();
//     if (muon && neutrino) {

//     // Set Q from Muon
//     TLorentzVector q = neutrino->fP - muon->fP;
//     double q0 = (q.E()) / 1.E3;
//     //double q3_true = (q.Vect().Mag())/1.E3;
//     double thmu = muon->fP.Vect().Angle(neutrino->fP.Vect());
//     double pmu  = muon->fP.Vect().Mag() / 1.E3;
//     double emu  = muon->fP.E() / 1.E3;
//     double mmu  = muon->fP.Mag() / 1.E3;


//     double neutrino_energy = (neutrino->fP.E())/1.E3;
//     double E_nu =neutrino_energy;
//     d3dml_bx->E_nu =neutrino_energy;
//     f1DHist_EAvail_true_nuenergy->Fill(neutrino_energy);

//     //// Get Enu Rec
//     double enu_rec = emu + q0;

//     // Set Q2 QE
//     double q2qe = 2 * enu_rec * (emu - pmu * cos(thmu)) - mmu * mmu;

//     // Calc Q3 from Q2QE and EnuTree
//     q3 = sqrt(q2qe + q0 * q0);

//     // Get Eav too
//     Eav = FitUtils::GetErecoil_MINERvA_LowRecoil(event) / 1.E3;
//     //Put the available energy into a box function
//     d3dml_bx->E_Avail = Eav;

//     // Basic interaction kinematics
//     double Q2 = -1 * (neutrino->fP - muon->fP).Mag2() / 1E6;
//     //double  q0 = (neutrino->fP - muon->fP).E() / 1E3;
//     // Get W_true with assumption of initial state nucleon at rest
//     float m_n = (float)PhysConst::mass_proton;
//     // Q2 assuming nucleon at rest
//     double W_nuc_rest = sqrt(-Q2 + 2 * m_n * q0 + m_n * m_n);
//     double  W = W_nuc_rest; // For want of a better thing to do
//     //std::cout << " W = "  << W << "   Q2 = " << Q2 << "   q0 = " << q0 << "   m_n =   "  << m_n <<  "  E available  =  " << Eav << std::endl;
//     d3dml_bx->W = W_nuc_rest;
//     }

//     static const double toGeV = 1E-3;

    
//     d3dml_bx->p_para = Pmu.Vect().Dot(nudir) * toGeV;
//     d3dml_bx->p_perp = Pmu.Vect().Cross(nudir).Mag() * toGeV;

//     // Sum up kinetic energy of protons
//     d3dml_bx->sum_TProt = 0.0;
//     for (auto prot : event->GetAllFSProton()) {
//       d3dml_bx->sum_TProt += prot->KE() * toGeV;
//     }
    
//     // find the bin number along each axis
//     int binx = f3DHist->GetXaxis()->FindFixBin(d3dml_bx->p_perp);
//     int biny = f3DHist->GetYaxis()->FindFixBin(d3dml_bx->p_para);
//     int binz = f3DHist->GetZaxis()->FindFixBin(d3dml_bx->sum_TProt);
//     int binz_EAvail = f3DHist_EAvail->GetZaxis()->FindFixBin(d3dml_bx->E_Avail);

//     // set this as the global bin number, could also use
//     // f3DHist->FindFixBin(pt,pz,sum_TProt)
//     f3DHist_EAvail->GetBin(binx, biny, binz_EAvail);
//     fXVar = f3DHist->GetBin(binx, biny, binz);
//   }

//   struct D3DMLBox : public MeasurementVariableBox1D {
//     D3DMLBox()
//         : MeasurementVariableBox1D(), mode{0},p_para{0}, p_perp{0}, sum_TProt{0}, E_Avail{0}, W{0}, E_nu{0} {}

//     MeasurementVariableBox *CloneSignalBox() {
//       auto cl = new D3DMLBox();

//       cl->fX = fX;
//       cl->mode = mode;
//       cl->p_para = p_para;
//       cl->p_perp = p_perp;
//       cl->sum_TProt = sum_TProt;
//       cl->E_Avail = E_Avail;
//       cl->W = W;
//       cl->E_nu = E_nu;
//       return cl;
//     }

//     void Reset() {
//       MeasurementVariableBox1D::Reset();
//       mode = 0;
//       p_para = 0;
//       p_perp = 0;
//       sum_TProt = 0;
//       E_Avail = 0;
//       W = 0;
//       E_nu = 0;
//     }

//     int mode;
//     double p_para;
//     double p_perp;
//     double sum_TProt;
//     double E_Avail;
//     double W;
//     double E_nu;
//   };
  
//   MeasurementVariableBox *CreateBox() { return new D3DMLBox(); };
//   int lowenergy_nu = 0;
//   int highenergy_nu = 0;
//   void FillExtraHistograms(MeasurementVariableBox *vars, double weight) {
    
//     Measurement1D::FillExtraHistograms(vars, weight);
    
//     auto d3dml_bx = dynamic_cast<D3DMLBox *>(vars);
//     if(d3dml_bx->E_nu < 8.0){
//       lowenergy_nu=lowenergy_nu+1;

//       //std::cout<< "Reco neutrino energy is ok:  "<< d3dml_bx->E_nu <<std::endl;
//       f3DHist->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,d3dml_bx->sum_TProt, weight);
      
//       f3DHist_EAvail->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//                           d3dml_bx->E_Avail, weight);
      
//       f3DHist_EAvail_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//                                          d3dml_bx->E_Avail, weight);

//       /*
//       f3DHist_EAvail_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//                           d3dml_bx->E_Avail, weight);

//       f3DHist_EAvail_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//                           d3dml_bx->E_Avail, weight);
//       */
//       int amode = std::abs(d3dml_bx->mode);
//       double W = std::abs(d3dml_bx->W);

//       if (amode == InputHandler::kCCQE) {
//         f3DHist_CCQE->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//                           d3dml_bx->sum_TProt, weight);
//         f3DHist_EAvail_CCQE->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         f3DHist_EAvail_CCQE_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         /*
//         f3DHist_EAvail_CCQE_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);

//         f3DHist_EAvail_CCQE_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         */
        
//       } else if (amode == InputHandler::kCC2p2h) {
//         f3DHist_CC2p2h->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//                             d3dml_bx->sum_TProt, weight);

//         f3DHist_EAvail_CC2p2h->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         f3DHist_EAvail_CC2p2h_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         /*
//         f3DHist_EAvail_CC2p2h_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         f3DHist_EAvail_CC2p2h_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         */

        
//       } else if ((amode == InputHandler::kCC1piponp) ||
//                 (amode == InputHandler::kCC1pi0onn) ||
//                 (amode == InputHandler::kCC1piponn)) {
//         f3DHist_CC1pi->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//                             d3dml_bx->sum_TProt, weight);
//         f3DHist_EAvail_CC1p1pi->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         f3DHist_EAvail_CC1p1pi_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);

//         /*
//         f3DHist_EAvail_CC1p1pi_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);

//         f3DHist_EAvail_CC1p1pi_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         */

//       } else if ((amode == InputHandler::kCCmultipi) ||
//                 (amode == InputHandler::kCCDIS)) {
//         f3DHist_CCDIS->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//                             d3dml_bx->sum_TProt, weight);
//         f3DHist_EAvail_CCDIS->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         f3DHist_EAvail_CCDIS_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         /*
//         f3DHist_EAvail_CCDIS_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         f3DHist_EAvail_CCDIS_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);

        
//         */


//       } else {
//         f3DHist_Other->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//                             d3dml_bx->sum_TProt, weight);
//         f3DHist_EAvail_Other->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         f3DHist_EAvail_Other_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         /*
//         f3DHist_EAvail_Other_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         f3DHist_EAvail_Other_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
        
//         */
//       }

//       //if statments to fill W hists
//       if (W < 1.4) {
//         f3DHist_EAvail_lowW->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         //std::cout << " W should be less than 1.4 , it is" << W <<std::endl;
//       } else if (1.4 < W && W < 2.0) {
//         f3DHist_EAvail_midW->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//       // std::cout << " W should be 1.4 < W < 2.0 , it is" << W <<std::endl;
//       } 
//       else if (W > 2.0){
//         f3DHist_EAvail_highW->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
//         d3dml_bx->E_Avail, weight);
//         //std::cout << " W should be > 2.0 , it is" << W <<std::endl;
//       } 
//     }
//     else {highenergy_nu=highenergy_nu+1;
//             }
//   }




//   //********************************************************************
//   bool isSignal(FitEvent *event) {
//     //********************************************************************
//     //return SignalDef::isCC0pi_MINERvAPTPZ(event, 14, EnuMin, EnuMax);
//     //return SignalDef::isCCincLowRecoil_MINERvA(event, EnuMin, EnuMax);
//     return SignalDef::isCCINC(event, 14, EnuMin, EnuMax);
//   }

  
  
//   void Create2DHistogram(std::unique_ptr<TH3D> const &h){ //
//         h->Write();
//         h->SetDirectory(nullptr);
//         TH1::SetDefaultSumw2(true);
//         std::cout<<"making 2D hist"<<std::endl;
//         int number_of_pt_bins = h->GetXaxis()->GetNbins();
//         int number_of_pz_bins = h->GetYaxis()->GetNbins();
//         int number_of_EAvail_bins = h->GetZaxis()->GetNbins();
//         TH3* ptslice = (TH3*)h->Clone();
//         //TH2D *ptslice_2D = new TH2D("ptslice_2D", "ptslice_2D", , number_of_pz_bins,number_of_EAvail_bins );
//         for (int x = 0; x < number_of_pt_bins - 1 ; ++x) {
//             std::cout<<"making slice in pt"<< x+1 << std::endl;
//             ptslice->GetXaxis()->SetRange(x+1,x+1); //take slice of width of a bin in pt
//             auto ptslice_2D= std::unique_ptr<TH2>(static_cast<TH2*>(ptslice->Project3D("yz"))); 

//             std::stringstream ss;
//             ss <<  "twodprojection_ptbin"<< x+1 <<  h->GetName()  ;

//             std::stringstream title;
//             title <<  "p_t"<< h->GetXaxis()->GetBinLowEdge(x + 1) << " < p_t < "
//             << h->GetXaxis()->GetBinUpEdge(x + 1) << " [GeV/c], " << "p_z"<< 
//             h->GetYaxis()->GetBinLowEdge(x + 1) << " < p_z < "<< h->GetYaxis()->GetBinUpEdge(x + 1) << " [GeV/c], "
//             <<h->GetZaxis()->GetBinLowEdge(x + 1) << " < E_avail < "<< h->GetYaxis()->GetBinUpEdge(x + 1) << " [GeV/c], " ;
            
//             std::string s = ss.str();
//             std::cout<< "ss for 2D" << s <<std::endl;

//             ptslice_2D->SetName(ss.str().c_str());
//             ptslice_2D->SetTitle(title.str().c_str());
//             ptslice_2D->GetYaxis()->SetTitle("pz (GeV) ");
//             ptslice_2D->GetZaxis()->SetTitle("Available Energy (GeV) ");
            
//             ptslice_2D->Write();
//             ptslice_2D->SetDirectory(nullptr);
          
//       }
      
      
//   }
  

//   void SplitWrite3D(std::unique_ptr<TH3D> const &h) {

//   std::cout<<"number of Enu < 8GeV"<< lowenergy_nu << std::endl;
//   std::cout<<"number of Enu > 8GeV"<< highenergy_nu << std::endl;
//     h->Write();
//     h->SetDirectory(nullptr);
//     TH1::SetDefaultSumw2(true);
//     std::ofstream outputFile1DALL("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/binerrors/all1Dbinerrors.txt"); // create a new output file or overwrite an existing one
       
//     auto ptpzproj = std::unique_ptr<TH2>(static_cast<TH2*>(h->Project3D("yx")));
//     std::string hist_2Dname = ptpzproj->GetName();
//     std::ofstream outputFile(("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/binerrors/2Dbinerrors" +hist_2Dname+ ".txt").c_str()); // create a new output file or overwrite an existing one
//     outputFile << " Bin errors for 2D projection histogram " << h->GetName() << "      \n" ;

//       int number_of_pt_bins = h->GetXaxis()->GetNbins();
//       int number_of_pz_bins = h->GetYaxis()->GetNbins();
//       int number_of_EAvail_bins = h->GetZaxis()->GetNbins();

//       std::cout << "The 3D histogram " << h->GetName() << "has" << number_of_pt_bins << "bins in pt (x axis) and " << number_of_pz_bins << "Pz bins" <<"and" << h->GetZaxis()->GetNbins()<< "energy bins"<<std::endl;
//       for (int x = 0; x < h->GetXaxis()->GetNbins(); ++x) {
  
//         for (int y = 0; y < h->GetYaxis()->GetNbins(); ++y) {
//           TH1 *proj =
//               h->ProjectionZ((std::string(h->GetName()) + "_x" +
//                               std::to_string(x) + "_y" + std::to_string(y))
//                                 .c_str(),
//                                 x+1, x+1 , y+1 , y+1 , "e");
    
        
//         std::stringstream ss;
//         ss << h->GetXaxis()->GetBinLowEdge(x + 1) << " < p_t < "
//            << h->GetXaxis()->GetBinUpEdge(x + 1) << " [GeV/c], "
//            << h->GetYaxis()->GetBinLowEdge(y + 1) << " < p_z < "
//            << h->GetYaxis()->GetBinUpEdge(y + 1) << " [GeV/c]";

//         double proj_cell_area = (h->GetXaxis()->GetBinUpEdge(x + 1) - h->GetXaxis()->GetBinLowEdge(x + 1)) *
//                                 (h->GetYaxis()->GetBinUpEdge(y + 1) - h->GetYaxis()->GetBinLowEdge(y + 1));

        
//         double x_bin_width_pt = h->GetXaxis()->GetBinUpEdge(x + 1) - h->GetXaxis()->GetBinLowEdge(x + 1);
//         double y_bin_width_pz = h->GetYaxis()->GetBinUpEdge(y + 1) - h->GetYaxis()->GetBinLowEdge(y + 1);
      
//         //proj->Scale(1.0); 
//         proj->Scale(1.0/proj_cell_area,"WIDTH"); 
//         //proj->Scale(fScaleFactor/proj_cell_area,"WIDTH"); 
//        // proj->Scale(proj_cell_area/fScaleFactor , "WIDTH"); //, "WIDTH"    fScaleFactor /
//         //std::cout<< "proj_cell_area  =  " << proj_cell_area  <<std::endl;
//         proj->SetTitle(ss.str().c_str());
//         proj->GetXaxis()->SetTitle("#Sigma T_{p} [GeV] or Available Energy ");

//         std::string hist_1Dname = proj->GetName();
//         std::ofstream outputFile1D(("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/binerrors/1Dbinerrors" + hist_1Dname + ".txt").c_str()); // create a new output file or overwrite an existing one
        
//         outputFile1D << " Bin errors for 1D projection histogram " << proj->GetName() << "      \n" ;
//         int binsin1Dhist=proj->GetNbinsX();

//         if (outputFile1D.is_open()) { // check if the file was opened successfully
//           for(int n =1 ; n< binsin1Dhist+1 ; n++){

//           double sqrt_bincontent1D = std::sqrt(proj->GetBinContent(n));
//           double standard_error_bincontent1D = 1.0/sqrt_bincontent1D;
//           outputFile1D << "Hist_name " << proj->GetName()  << "  "  << " Histogram_noofbins " << binsin1Dhist << "  \n";
//           outputFile1D << "Hist_binwidth_x" << proj->GetXaxis()->GetBinWidth(n)<< "  ";
//           //outputFile1D << "Histogram_binwidth_y" << proj->GetYaxis()->GetBinWidth(n)<< "  ";
//           outputFile1D << "Bin_number " << proj->GetBin(n) << "  ";
//           outputFile1D << "proj_cell_area " << proj_cell_area << "  ";
//           outputFile1D << "x_bin_width_pt " << x_bin_width_pt << "  ";
//           outputFile1D << "y_bin_width_pz " << y_bin_width_pz << "  ";

//           outputFile1D << "Bin_Content " << proj->GetBinContent(n) << "  " ;
//           outputFile1D << "Sqrt_Bin_Content " << sqrt_bincontent1D << "  " ;
//           outputFile1D << "Bin_Error_ROOT " << proj->GetBinError(n) << "  \n ";
          
//           outputFile1DALL << "Histogram_name " << proj->GetName() << "  "  << " Histogram_noofbins " << binsin1Dhist << "  \n";
//           outputFile1DALL << "Histogram_binwidth " <<proj->GetXaxis()->GetBinWidth(n)<< "  ";
//           outputFile1DALL << "Bin_number " << proj->GetBin(n) << "  ";
//           outputFile1DALL << "proj_cell_area " << proj_cell_area << "  ";
//           outputFile1DALL << "x_bin_edge " << x_bin_width_pt << "  ";
//           outputFile1DALL << "y_bin_edge " << y_bin_width_pz << "  ";

//           outputFile1DALL << "Bin_Content " << proj->GetBinContent(n) << "  " ;
//           outputFile1DALL << "Sqrt_Bin_Content " << sqrt_bincontent1D << "  " ;
//           outputFile1DALL << "Bin_Error_ROOT " << proj->GetBinError(n) << "  \n";
         

//           }
//         }
//         double sqrt_bincontent = std::sqrt(ptpzproj->GetBinContent(x + 1, y + 1));
//         double standard_error_bincontent = 1.0/sqrt_bincontent;
//         if (outputFile.is_open()) { // check if the file was opened successfully
//           outputFile << " Bin number " << ptpzproj->GetBin(x + 1, y + 1) << "  ";
//           outputFile << " Bin Content = " << ptpzproj->GetBinContent(x + 1, y + 1) << "  " ;
//           outputFile << " Sqrt Bin Content = " << sqrt_bincontent << "  " ;
//           outputFile << " Bin Error = " << ptpzproj->GetBinError(x + 1, y + 1) << "  ";
//           outputFile << " Standard error (1/sqrt(n)) = " << standard_error_bincontent << "      \n";
//         }
      
//         outputFile1D.close();
//         proj->Write();
//         proj->SetDirectory(nullptr);
//         delete proj;
//       }

       
//     }
//     outputFile.close(); // close the file when done
//     outputFile1DALL.close();
//     ptpzproj->Scale(fScaleFactor,"WIDTH");  //get  root to divide out the bin width
//     ptpzproj->Write();
//     //make sure to tell ROOT that it doesn't own this histogram so that we can delete it
//     ptpzproj->SetDirectory(nullptr);

//   }

//   void Write(std::string drawOpt) {    
//     SplitWrite3D(f3DHist);
//     SplitWrite3D(f3DHist_CCQE);
//     SplitWrite3D(f3DHist_CC2p2h);
//     SplitWrite3D(f3DHist_CC1pi);
//     SplitWrite3D(f3DHist_CCDIS);
//     SplitWrite3D(f3DHist_Other);

//     SplitWrite3D(f3DHist_EAvail);
//     SplitWrite3D(f3DHist_EAvail_CCQE);
//     SplitWrite3D(f3DHist_EAvail_CC2p2h);
//     SplitWrite3D(f3DHist_EAvail_CC1p1pi);
//     SplitWrite3D(f3DHist_EAvail_CCDIS);
//     SplitWrite3D(f3DHist_EAvail_Other);

    
    

//     //SplitWrite3D(f3DHist_EAvail_pzbins_double);
//     /*
//     SplitWrite3D(f3DHist_EAvail_pzbins_uniform);
//     SplitWrite3D(f3DHist_EAvail_pzbins_uniform_double);
//     SplitWrite3D(f3DHist_EAvail_CCQE_pzbins_uniform);
//     SplitWrite3D(f3DHist_EAvail_CC2p2h_pzbins_uniform);
//     SplitWrite3D(f3DHist_EAvail_CC1p1pi_pzbins_uniform);
//     SplitWrite3D(f3DHist_EAvail_CCDIS_pzbins_uniform);
//     SplitWrite3D(f3DHist_EAvail_Other_pzbins_uniform);


//     SplitWrite3D(f3DHist_EAvail_CCQE_pzbins_uniform_double);
//     SplitWrite3D(f3DHist_EAvail_CC2p2h_pzbins_uniform_double);
//     SplitWrite3D(f3DHist_EAvail_CC1p1pi_pzbins_uniform_double);
//     SplitWrite3D(f3DHist_EAvail_CCDIS_pzbins_uniform_double);
//     SplitWrite3D(f3DHist_EAvail_Other_pzbins_uniform_double);
//     */
//     SplitWrite3D(f3DHist_EAvail_pzbins_double);
//     SplitWrite3D(f3DHist_EAvail_CCQE_pzbins_double);
//     SplitWrite3D(f3DHist_EAvail_CC2p2h_pzbins_double);
//     SplitWrite3D(f3DHist_EAvail_CC1p1pi_pzbins_double);
//     SplitWrite3D(f3DHist_EAvail_CCDIS_pzbins_double);
//     SplitWrite3D(f3DHist_EAvail_Other_pzbins_double);
    
   

//     SplitWrite3D(f3DHist_EAvail_lowW);
//     SplitWrite3D(f3DHist_EAvail_midW);
//     SplitWrite3D(f3DHist_EAvail_highW);



//     Create2DHistogram(f3DHist_EAvail);
//     Create2DHistogram(f3DHist_EAvail_CCQE);
//     Create2DHistogram(f3DHist_EAvail_CC2p2h);
//     Create2DHistogram(f3DHist_EAvail_CC1p1pi);
//     Create2DHistogram(f3DHist_EAvail_CCDIS);
//     Create2DHistogram(f3DHist_EAvail_Other);
    
//     f3DHist_EAvail->Write();
//     f3DHist_EAvail_CCQE->Write();
//     f3DHist_EAvail_CC2p2h->Write();
//     f3DHist_EAvail_CC1p1pi->Write();
//     f3DHist_EAvail_CCDIS->Write();
//     f3DHist_EAvail_Other->Write();
//     f1DHist_EAvail_true_nuenergy->Write();

//     f3DHist_EAvail_pzbins_double->Write();
//     f3DHist_EAvail_CCQE_pzbins_double->Write();
//     f3DHist_EAvail_CC2p2h_pzbins_double->Write();
//     f3DHist_EAvail_CC1p1pi_pzbins_double->Write();
//     f3DHist_EAvail_CCDIS_pzbins_double->Write();
//     f3DHist_EAvail_Other_pzbins_double->Write();

//     // we have to tidy this up in this SO if we don't want horrible crashes on
//     // program tear down
//     fDataHist->Write();
//     fDataHist->SetDirectory(nullptr);
//     fMCHist->Write();
//     fMCHist->SetDirectory(nullptr);
//   }

//   void ResetAll() {
//     Measurement1D::ResetAll();
//     f3DHist->Reset();
//     f3DHist_CCQE->Reset();
//     f3DHist_CC2p2h->Reset();
//     f3DHist_CC1pi->Reset();
//     f3DHist_CCDIS->Reset();
//     f3DHist_Other->Reset();

//     f3DHist_EAvail->Reset();
//     f3DHist_EAvail_CCQE->Reset();
//     f3DHist_EAvail_CC2p2h->Reset();
//     f3DHist_EAvail_CC1p1pi->Reset();
//     f3DHist_EAvail_CCDIS->Reset();
//     f3DHist_EAvail_Other->Reset();

//     f1DHist_EAvail_true_nuenergy->Reset();

//     /*
//     f3DHist_EAvail_pzbins_uniform->Reset();
//     f3DHist_EAvail_CCQE_pzbins_uniform->Reset();
//     f3DHist_EAvail_CC2p2h_pzbins_uniform->Reset();
//     f3DHist_EAvail_CC1p1pi_pzbins_uniform->Reset();
//     f3DHist_EAvail_CCDIS_pzbins_uniform->Reset();
//     f3DHist_EAvail_Other_pzbins_uniform->Reset();

//     f3DHist_EAvail_pzbins_uniform_double->Reset();
//     f3DHist_EAvail_CCQE_pzbins_uniform_double->Reset();
//     f3DHist_EAvail_CC2p2h_pzbins_uniform_double->Reset();
//     f3DHist_EAvail_CC1p1pi_pzbins_uniform_double->Reset();
//     f3DHist_EAvail_CCDIS_pzbins_uniform_double->Reset();
//     f3DHist_EAvail_Other_pzbins_uniform_double->Reset();

//     f3DHist_EAvail_pzbins_uniform_double->Reset();
//     */

//     f3DHist_EAvail_pzbins_double->Reset();
//     f3DHist_EAvail_CCQE_pzbins_double->Reset();
//     f3DHist_EAvail_CC2p2h_pzbins_double->Reset();
//     f3DHist_EAvail_CC1p1pi_pzbins_double->Reset();
//     f3DHist_EAvail_CCDIS_pzbins_double->Reset();
//     f3DHist_EAvail_Other_pzbins_double->Reset();
    
//     f3DHist_EAvail_lowW->Reset();
//     f3DHist_EAvail_midW->Reset();
//     f3DHist_EAvail_highW->Reset();

//   }
  
// };


#include "InteractionModes.h"
#include "/home/abipeake/scratch/Formfactorstuff/nuisance/src/MINERvA/MINERvA_SignalDef.h"
#include "Measurement1D.h"
#include "/home/abipeake/scratch/Formfactorstuff/nuisance/src/MCStudies/GenericFlux_Vectors.h"
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
  std::unique_ptr<TH1D> f1DHist_EAvail_true_nuenergy;

  std::unique_ptr<TH3D> f3DHist_CCQE;
  std::unique_ptr<TH3D> f3DHist_CC2p2h;
  std::unique_ptr<TH3D> f3DHist_CC1pi;
  std::unique_ptr<TH3D> f3DHist_CCDIS;
  std::unique_ptr<TH3D> f3DHist_Other;
  std::unique_ptr<TH3D> f3DHist_CCCoh;

  std::unique_ptr<TH3D> f3DHist_EAvail;  //pointers to the new histograms for the new EAvailable binning
  std::unique_ptr<TH3D> f3DHist_EAvail_CCQE;
  std::unique_ptr<TH3D> f3DHist_EAvail_CC2p2h;
  std::unique_ptr<TH3D> f3DHist_EAvail_CC1p1pi;
  std::unique_ptr<TH3D> f3DHist_EAvail_CCDIS;
  std::unique_ptr<TH3D> f3DHist_EAvail_Other;
  std::unique_ptr<TH3D> f3DHist_EAvail_CCCoh;


  std::unique_ptr<TH2D> f2DHist_EAvail;  //pointers to the new histograms for the new EAvailable binning
  std::unique_ptr<TH2D> f2DHist_EAvail_CCQE;
  std::unique_ptr<TH2D> f2DHist_EAvail_CC2p2h;
  std::unique_ptr<TH2D> f2DHist_EAvail_CC1p1pi;
  std::unique_ptr<TH2D> f2DHist_EAvail_CCDIS;
  std::unique_ptr<TH2D> f2DHist_EAvail_Other;

  
  std::unique_ptr<TH3D> f3DHist_EAvail_pzbins_double;  //pointers to the new histograms for the new EAvailable binning
  std::unique_ptr<TH3D> f3DHist_EAvail_CCQE_pzbins_double;
  std::unique_ptr<TH3D> f3DHist_EAvail_CC2p2h_pzbins_double;
  std::unique_ptr<TH3D> f3DHist_EAvail_CC1p1pi_pzbins_double;
  std::unique_ptr<TH3D> f3DHist_EAvail_CCDIS_pzbins_double;
  std::unique_ptr<TH3D> f3DHist_EAvail_Other_pzbins_double;
  std::unique_ptr<TH3D> f3DHist_EAvail_CCCoh_pzbins_double;

  /*
  std::unique_ptr<TH3D> f3DHist_EAvail_pzbins_uniform; 
  std::unique_ptr<TH3D> f3DHist_EAvail_pzbins_uniform_double;  //pointers to the new histograms for the new EAvailable binning
  std::unique_ptr<TH3D> f3DHist_EAvail_CCQE_pzbins_uniform;
  std::unique_ptr<TH3D> f3DHist_EAvail_CC2p2h_pzbins_uniform;
  std::unique_ptr<TH3D> f3DHist_EAvail_CC1p1pi_pzbins_uniform;
  std::unique_ptr<TH3D> f3DHist_EAvail_CCDIS_pzbins_uniform;
  std::unique_ptr<TH3D> f3DHist_EAvail_Other_pzbins_uniform;

   std::unique_ptr<TH3D> f3DHist_EAvail_CCQE_pzbins_uniform_double;
  std::unique_ptr<TH3D> f3DHist_EAvail_CC2p2h_pzbins_uniform_double;
  std::unique_ptr<TH3D> f3DHist_EAvail_CC1p1pi_pzbins_uniform_double;
  std::unique_ptr<TH3D> f3DHist_EAvail_CCDIS_pzbins_uniform_double;
  std::unique_ptr<TH3D> f3DHist_EAvail_Other_pzbins_uniform_double;
  */
  std::unique_ptr<TH3D> f3DHist_EAvail_lowW;  //pointers to the new histograms for the new EAvailable binning
  std::unique_ptr<TH3D> f3DHist_EAvail_midW;
  std::unique_ptr<TH3D> f3DHist_EAvail_highW;

  // ---------------------------------------------------------------------
  // Oscillation-window subsample: events whose TRUE neutrino energy lies
  // between the first oscillation minimum (phase pi, ~1.3 GeV) and the first
  // oscillation maximum (phase pi/2, ~2.6 GeV) for DUNE, L = 1300 km,
  // dm^2_32 ~ 2.5e-3 eV^2. Change the limits here if you use a different dm^2.
  //
  // fOscWin index: 0 = all, 1 = CCQE, 2 = CC2p2h, 3 = CC1p1pi, 4 = CCDIS,
  //               5 = Other, 6 = CCCoh
  // Histogram names: DUNE_3D_MINERvALike_EAvail_OscWin[_CCQE|_CC2p2h|...|_CCCoh]
  // ---------------------------------------------------------------------
  const double kEnuOscLo = 1.3; // GeV
  const double kEnuOscHi = 2.6; // GeV
  std::vector<std::unique_ptr<TH3D>> fOscWin;


  /* a function to generate numpy linspace */
//template <typename T>
  template <typename T>
    std::vector<T> linspace(T a, T b, size_t N) {
        T h = (b - a) / static_cast<T>(N);
        std::vector<T> xs(N);
        typename std::vector<T>::iterator x;
        T val;
        for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
            *x = val;
        return xs;
    }

  void ExtendLinspace(std::vector<double> &in, double low, double high,size_t nbins){
    for(auto x : linspace(low, high,nbins)){
        in.push_back(x);
        
      }
      auto x =std::unique(in.begin(),in.end());  //put duplicate bins at end
      in.erase(x, in.end()); //erase duplicates
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
    //std::vector<double> ptbins =linspace(0.0, 2.5,20) ;

    std::vector<double> ptbins;
    ExtendLinspace(ptbins,1e-8,0.2,4);
    ExtendLinspace(ptbins,0.2,0.4,5);
    ExtendLinspace(ptbins,0.4,0.7,4);
    ExtendLinspace(ptbins,0.7,2.5,8);
  ptbins.push_back(2.5);
  
    std::cout<<"pt vector is: "<<std::endl;
    for(int i=0 ; i< 21 ; i++){
        std::cout<< ptbins[i]<<std::endl;

    }

    std::vector<double> pzbins_uniform = linspace(1.0e-8,6.6,22); //linearly spaced bins
    std::vector<double> pzbins_uniform_doublebins = linspace(1.0e-8,6.6,44); //linearly spaced bins
    
    std::vector<double> pzbins;
    ExtendLinspace(pzbins,1e-8,1.5,5);
    ExtendLinspace(pzbins,1.5,1.8,2);
    ExtendLinspace(pzbins,1.8,2.3,7);
    ExtendLinspace(pzbins,2.3,6.6,10);
  pzbins.push_back(6.6);
    

    std::vector<double> pzbins_double;
    ExtendLinspace(pzbins_double,1e-8,1.5,10);
    ExtendLinspace(pzbins_double,1.5,1.8,4);
    ExtendLinspace(pzbins_double,1.8,2.3,14);
    ExtendLinspace(pzbins_double,2.3,6.6,20);
    
    std::cout<<"pz vector is: "<<std::endl;
    for(int i=0 ; i< pzbins.size(); i++){
        std::cout<< pzbins[i]<<std::endl;
    }

    std::cout<<"pzbins= double is: "<<std::endl;
    for(int i=0 ; i< pzbins_double.size() ; i++){
        std::cout<< pzbins_double[i]<<std::endl;

    }

   
    std::cout<<"pzbins_uniform is: "<<std::endl;
    for(int i=0 ; i<pzbins_uniform.size() ; i++){
        std::cout<< pzbins_uniform[i]<<std::endl;
    }
     /*
    std::cout<<"pzbins_uniform double is: "<<std::endl;
    for(int i=0 ; i< pzbins_uniform_doublebins.size() ; i++){
        std::cout<< pzbins_uniform_doublebins[i]<<std::endl;
    }
    */
    
    std::vector<double> sumTpbins = {1e-8,  0.02, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.4,  0.6}; // GeV
    std::vector<double> EAvail_bins = {1e-8, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.2, 0.24, 0.28, 0.32, 0.4, 0.5, 0.6, 0.8, 1.2, 1.6, 2.0, 4.0, 8.0}; // GeV
    
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

    

    f3DHist_EAvail_pzbins_double =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_pzbins_double", "", ptbins.size() - 1,
                               ptbins.data(), pzbins_double.size() - 1, pzbins_double.data(),
                               EAvail_bins.size() - 1, EAvail_bins.data());

    /*
    f3DHist_EAvail_pzbins_uniform =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_pzbins_uniform", "", ptbins.size() - 1,
                               ptbins.data(), pzbins_uniform.size() - 1, pzbins_uniform.data(),
                               EAvail_bins.size() - 1, EAvail_bins.data());

    f3DHist_EAvail_pzbins_uniform_double =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_pzbins_uniform_double", "", ptbins.size() - 1,
                               ptbins.data(), pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(),
                               EAvail_bins.size() - 1, EAvail_bins.data());
    
    */

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
    f3DHist_CCCoh = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_CCCoh", "", ptbins.size() - 1, ptbins.data(),
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
    f3DHist_EAvail_CCCoh = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_EAvail_CCCoh", "", ptbins.size() - 1, ptbins.data(),
        pzbins.size() - 1, pzbins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());


    /*
    //uniform binning
    f3DHist_EAvail_CCQE_pzbins_uniform =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CCQE_pzbins_uniform", "",ptbins.size() - 1, ptbins.data(),
        pzbins_uniform.size() - 1, pzbins_uniform.data(), EAvail_bins.size() - 1,
        EAvail_bins.data()) ;

    f3DHist_EAvail_CC2p2h_pzbins_uniform =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC2p2h_pzbins_uniform", "", ptbins.size() - 1, ptbins.data(),
        pzbins_uniform.size() - 1, pzbins_uniform.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());

    f3DHist_EAvail_CC1p1pi_pzbins_uniform =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC1p1pi_pzbins_uniform", "", ptbins.size() - 1, ptbins.data(),
        pzbins_uniform.size() - 1, pzbins_uniform.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());
    
    f3DHist_EAvail_CCDIS_pzbins_uniform = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_EAvail_CCDI_pzbins_uniformS", "", ptbins.size() - 1, ptbins.data(),
        pzbins_uniform.size() - 1, pzbins_uniform.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());
    f3DHist_EAvail_Other_pzbins_uniform = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_EAvail_Other_pzbins_uniform", "", ptbins.size() - 1, ptbins.data(),
        pzbins_uniform.size() - 1, pzbins_uniform.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());

    //double uniform binning
    //uniform_double binning
    f3DHist_EAvail_CCQE_pzbins_uniform_double =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CCQE_pzbins_uniform_double", "",ptbins.size() - 1, ptbins.data(),
        pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data()) ;

    f3DHist_EAvail_CC2p2h_pzbins_uniform_double =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC2p2h_pzbins_uniform_double", "", ptbins.size() - 1, ptbins.data(),
        pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());

    f3DHist_EAvail_CC1p1pi_pzbins_uniform_double =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC1p1pi_pzbins_uniform_double", "", ptbins.size() - 1, ptbins.data(),
        pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());
    
    f3DHist_EAvail_CCDIS_pzbins_uniform_double = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_EAvail_CCDI_pzbins_uniform_doubleS", "", ptbins.size() - 1, ptbins.data(),
        pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());
    f3DHist_EAvail_Other_pzbins_uniform_double = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_EAvail_Other_pzbins_uniform_double", "", ptbins.size() - 1, ptbins.data(),
        pzbins_uniform_doublebins.size() - 1, pzbins_uniform_doublebins.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());
        */

    //1D hist of Enu_true
    f1DHist_EAvail_true_nuenergy =
        std::make_unique<TH1D>("DUNE_3D_MINERvALike__EAvail_true_nuenergy", "", EAvail_bins.size() - 1,
        EAvail_bins.data()) ;


    //double binning  over peak
    f3DHist_EAvail_CCQE_pzbins_double =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CCQE_pzbins_double", "",ptbins.size() - 1, ptbins.data(),
        pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
        EAvail_bins.data()) ;


    f3DHist_EAvail_CC2p2h_pzbins_double =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC2p2h_pzbins_double", "", ptbins.size() - 1, ptbins.data(),
        pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());

    f3DHist_EAvail_CC1p1pi_pzbins_double =
        std::make_unique<TH3D>("DUNE_3D_MINERvALike_EAvail_CC1p1pi_pzbins_double", "", ptbins.size() - 1, ptbins.data(),
        pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());
    
    f3DHist_EAvail_CCDIS_pzbins_double = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_EAvail_CCDI_pzbins_doubleS", "", ptbins.size() - 1, ptbins.data(),
        pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());
    f3DHist_EAvail_Other_pzbins_double = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_EAvail_Other_pzbins_double", "", ptbins.size() - 1, ptbins.data(),
        pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
        EAvail_bins.data());
    f3DHist_EAvail_CCCoh_pzbins_double = std::make_unique<TH3D>(
        "DUNE_3D_MINERvALike_EAvail_CCCoh_pzbins_double", "", ptbins.size() - 1, ptbins.data(),
        pzbins_double.size() - 1, pzbins_double.data(), EAvail_bins.size() - 1,
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

    // Oscillation-window histograms (same pT, pz, E_avail binning as f3DHist_EAvail)
    {
      std::stringstream wintitle;
      wintitle << kEnuOscLo << " < E_{#nu}^{true} < " << kEnuOscHi << " GeV";
      for (const char *suffix : {"", "_CCQE", "_CC2p2h", "_CC1p1pi", "_CCDIS", "_Other", "_CCCoh"}) {
        fOscWin.emplace_back(std::make_unique<TH3D>(
            (std::string("DUNE_3D_MINERvALike_EAvail_OscWin") + suffix).c_str(),
            wintitle.str().c_str(),
            ptbins.size() - 1, ptbins.data(),
            pzbins.size() - 1, pzbins.data(),
            EAvail_bins.size() - 1, EAvail_bins.data()));
        fOscWin.back()->SetDirectory(nullptr);
      }
    }

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


    double neutrino_energy = (neutrino->fP.E())/1.E3;
    double E_nu =neutrino_energy;
    d3dml_bx->E_nu =neutrino_energy;
    f1DHist_EAvail_true_nuenergy->Fill(neutrino_energy);

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
        : MeasurementVariableBox1D(), mode{0},p_para{0}, p_perp{0}, sum_TProt{0}, E_Avail{0}, W{0}, E_nu{0} {}

    MeasurementVariableBox *CloneSignalBox() {
      auto cl = new D3DMLBox();

      cl->fX = fX;
      cl->mode = mode;
      cl->p_para = p_para;
      cl->p_perp = p_perp;
      cl->sum_TProt = sum_TProt;
      cl->E_Avail = E_Avail;
      cl->W = W;
      cl->E_nu = E_nu;
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
      E_nu = 0;
    }

    int mode;
    double p_para;
    double p_perp;
    double sum_TProt;
    double E_Avail;
    double W;
    double E_nu;
  };
  
  MeasurementVariableBox *CreateBox() { return new D3DMLBox(); };
  int lowenergy_nu = 0;
  int highenergy_nu = 0;
  void FillExtraHistograms(MeasurementVariableBox *vars, double weight) {
    
    Measurement1D::FillExtraHistograms(vars, weight);
    
    auto d3dml_bx = dynamic_cast<D3DMLBox *>(vars);
    if(d3dml_bx->E_nu < 8.0){
      lowenergy_nu=lowenergy_nu+1;

      //std::cout<< "Reco neutrino energy is ok:  "<< d3dml_bx->E_nu <<std::endl;
      f3DHist->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,d3dml_bx->sum_TProt, weight);
      
      f3DHist_EAvail->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                          d3dml_bx->E_Avail, weight);
      
      f3DHist_EAvail_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                                         d3dml_bx->E_Avail, weight);

      /*
      f3DHist_EAvail_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                          d3dml_bx->E_Avail, weight);

      f3DHist_EAvail_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                          d3dml_bx->E_Avail, weight);
      */
      int amode = std::abs(d3dml_bx->mode);
      double W = std::abs(d3dml_bx->W);

      if (amode == InputHandler::kCCQE) {
        f3DHist_CCQE->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                          d3dml_bx->sum_TProt, weight);
        f3DHist_EAvail_CCQE->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        f3DHist_EAvail_CCQE_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        /*
        f3DHist_EAvail_CCQE_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);

        f3DHist_EAvail_CCQE_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        */
        
      } else if (amode == InputHandler::kCC2p2h) {
        f3DHist_CC2p2h->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                            d3dml_bx->sum_TProt, weight);

        f3DHist_EAvail_CC2p2h->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        f3DHist_EAvail_CC2p2h_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        /*
        f3DHist_EAvail_CC2p2h_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        f3DHist_EAvail_CC2p2h_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        */

        
      } else if ((amode == InputHandler::kCC1piponp) ||
                (amode == InputHandler::kCC1pi0onn) ||
                (amode == InputHandler::kCC1piponn)) {
        f3DHist_CC1pi->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                            d3dml_bx->sum_TProt, weight);
        f3DHist_EAvail_CC1p1pi->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        f3DHist_EAvail_CC1p1pi_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);

        /*
        f3DHist_EAvail_CC1p1pi_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);

        f3DHist_EAvail_CC1p1pi_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        */

      } else if ((amode == InputHandler::kCCmultipi) ||
                (amode == InputHandler::kCCDIS)) {
        f3DHist_CCDIS->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                            d3dml_bx->sum_TProt, weight);
        f3DHist_EAvail_CCDIS->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        f3DHist_EAvail_CCDIS_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        /*
        f3DHist_EAvail_CCDIS_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        f3DHist_EAvail_CCDIS_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);

        
        */


      } else if (amode == InputHandler::kCCCoherent) {
        f3DHist_CCCoh->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                            d3dml_bx->sum_TProt, weight);
        f3DHist_EAvail_CCCoh->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        f3DHist_EAvail_CCCoh_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);

      } else {
        f3DHist_Other->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                            d3dml_bx->sum_TProt, weight);
        f3DHist_EAvail_Other->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        f3DHist_EAvail_Other_pzbins_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        /*
        f3DHist_EAvail_Other_pzbins_uniform->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        f3DHist_EAvail_Other_pzbins_uniform_double->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
        d3dml_bx->E_Avail, weight);
        
        */
      }

      // ---- oscillation-window subsample (true Enu between kEnuOscLo and kEnuOscHi) ----
      if (d3dml_bx->E_nu >= kEnuOscLo && d3dml_bx->E_nu < kEnuOscHi) {
        int cat = 5; // Other
        if (amode == InputHandler::kCCQE) {
          cat = 1;
        } else if (amode == InputHandler::kCC2p2h) {
          cat = 2;
        } else if ((amode == InputHandler::kCC1piponp) ||
                   (amode == InputHandler::kCC1pi0onn) ||
                   (amode == InputHandler::kCC1piponn)) {
          cat = 3;
        } else if ((amode == InputHandler::kCCmultipi) ||
                   (amode == InputHandler::kCCDIS)) {
          cat = 4;
        } else if (amode == InputHandler::kCCCoherent) {
          cat = 6;
        }
        fOscWin[0]->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
                         d3dml_bx->E_Avail, weight);
        fOscWin[cat]->Fill(d3dml_bx->p_perp, d3dml_bx->p_para,
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
    else {highenergy_nu=highenergy_nu+1;
            }
  }




  //********************************************************************
  bool isSignal(FitEvent *event) {
    //********************************************************************
    //return SignalDef::isCC0pi_MINERvAPTPZ(event, 14, EnuMin, EnuMax);
    //return SignalDef::isCCincLowRecoil_MINERvA(event, EnuMin, EnuMax);
    // CC numu inclusive: no cut on muon angle or energy
    return SignalDef::isCCINC(event, 14, EnuMin, EnuMax);
  }

  
  
  void Create2DHistogram(std::unique_ptr<TH3D> const &h){ //
        h->Write();
        h->SetDirectory(nullptr);
        TH1::SetDefaultSumw2(true);
        std::cout<<"making 2D hist"<<std::endl;
        int number_of_pt_bins = h->GetXaxis()->GetNbins();
        int number_of_pz_bins = h->GetYaxis()->GetNbins();
        int number_of_EAvail_bins = h->GetZaxis()->GetNbins();
        TH3* ptslice = (TH3*)h->Clone();
        //TH2D *ptslice_2D = new TH2D("ptslice_2D", "ptslice_2D", , number_of_pz_bins,number_of_EAvail_bins );
        for (int x = 0; x < number_of_pt_bins - 1 ; ++x) {
            std::cout<<"making slice in pt"<< x+1 << std::endl;
            ptslice->GetXaxis()->SetRange(x+1,x+1); //take slice of width of a bin in pt
            auto ptslice_2D= std::unique_ptr<TH2>(static_cast<TH2*>(ptslice->Project3D("yz"))); 

            std::stringstream ss;
            ss <<  "twodprojection_ptbin"<< x+1 <<  h->GetName()  ;

            std::stringstream title;
            title <<  "p_t"<< h->GetXaxis()->GetBinLowEdge(x + 1) << " < p_t < "
            << h->GetXaxis()->GetBinUpEdge(x + 1) << " [GeV/c], " << "p_z"<< 
            h->GetYaxis()->GetBinLowEdge(x + 1) << " < p_z < "<< h->GetYaxis()->GetBinUpEdge(x + 1) << " [GeV/c], "
            <<h->GetZaxis()->GetBinLowEdge(x + 1) << " < E_avail < "<< h->GetYaxis()->GetBinUpEdge(x + 1) << " [GeV/c], " ;
            
            std::string s = ss.str();
            std::cout<< "ss for 2D" << s <<std::endl;

            ptslice_2D->SetName(ss.str().c_str());
            ptslice_2D->SetTitle(title.str().c_str());
            ptslice_2D->GetYaxis()->SetTitle("pz (GeV) ");
            ptslice_2D->GetZaxis()->SetTitle("Available Energy (GeV) ");
            
            ptslice_2D->Write();
            ptslice_2D->SetDirectory(nullptr);
          
      }
      
      
  }
  

  void SplitWrite3D(std::unique_ptr<TH3D> const &h) {

  std::cout<<"number of Enu < 8GeV"<< lowenergy_nu << std::endl;
  std::cout<<"number of Enu > 8GeV"<< highenergy_nu << std::endl;
    h->Write();
    h->SetDirectory(nullptr);
    TH1::SetDefaultSumw2(true);
    std::ofstream outputFile1DALL("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/binerrors/all1Dbinerrors.txt"); // create a new output file or overwrite an existing one
       
    auto ptpzproj = std::unique_ptr<TH2>(static_cast<TH2*>(h->Project3D("yx")));
    std::string hist_2Dname = ptpzproj->GetName();
    std::ofstream outputFile(("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/binerrors/2Dbinerrors" +hist_2Dname+ ".txt").c_str()); // create a new output file or overwrite an existing one
    outputFile << " Bin errors for 2D projection histogram " << h->GetName() << "      \n" ;

      int number_of_pt_bins = h->GetXaxis()->GetNbins();
      int number_of_pz_bins = h->GetYaxis()->GetNbins();
      int number_of_EAvail_bins = h->GetZaxis()->GetNbins();

      std::cout << "The 3D histogram " << h->GetName() << "has" << number_of_pt_bins << "bins in pt (x axis) and " << number_of_pz_bins << "Pz bins" <<"and" << h->GetZaxis()->GetNbins()<< "energy bins"<<std::endl;
      for (int x = 0; x < h->GetXaxis()->GetNbins(); ++x) {
  
        for (int y = 0; y < h->GetYaxis()->GetNbins(); ++y) {
          TH1 *proj =
              h->ProjectionZ((std::string(h->GetName()) + "_x" +
                              std::to_string(x) + "_y" + std::to_string(y))
                                .c_str(),
                                x+1, x+1 , y+1 , y+1 , "e");
    
        
        std::stringstream ss;
        ss << h->GetXaxis()->GetBinLowEdge(x + 1) << " < p_t < "
           << h->GetXaxis()->GetBinUpEdge(x + 1) << " [GeV/c], "
           << h->GetYaxis()->GetBinLowEdge(y + 1) << " < p_z < "
           << h->GetYaxis()->GetBinUpEdge(y + 1) << " [GeV/c]";

        double proj_cell_area = (h->GetXaxis()->GetBinUpEdge(x + 1) - h->GetXaxis()->GetBinLowEdge(x + 1)) *
                                (h->GetYaxis()->GetBinUpEdge(y + 1) - h->GetYaxis()->GetBinLowEdge(y + 1));

        
        double x_bin_width_pt = h->GetXaxis()->GetBinUpEdge(x + 1) - h->GetXaxis()->GetBinLowEdge(x + 1);
        double y_bin_width_pz = h->GetYaxis()->GetBinUpEdge(y + 1) - h->GetYaxis()->GetBinLowEdge(y + 1);
      
        //proj->Scale(1.0); 
        proj->Scale(1.0/proj_cell_area,"WIDTH"); 
        //proj->Scale(fScaleFactor/proj_cell_area,"WIDTH"); 
       // proj->Scale(proj_cell_area/fScaleFactor , "WIDTH"); //, "WIDTH"    fScaleFactor /
        //std::cout<< "proj_cell_area  =  " << proj_cell_area  <<std::endl;
        proj->SetTitle(ss.str().c_str());
        proj->GetXaxis()->SetTitle("#Sigma T_{p} [GeV] or Available Energy ");

        std::string hist_1Dname = proj->GetName();
        std::ofstream outputFile1D(("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/binerrors/1Dbinerrors" + hist_1Dname + ".txt").c_str()); // create a new output file or overwrite an existing one
        
        outputFile1D << " Bin errors for 1D projection histogram " << proj->GetName() << "      \n" ;
        int binsin1Dhist=proj->GetNbinsX();

        if (outputFile1D.is_open()) { // check if the file was opened successfully
          for(int n =1 ; n< binsin1Dhist+1 ; n++){

          double sqrt_bincontent1D = std::sqrt(proj->GetBinContent(n));
          double standard_error_bincontent1D = 1.0/sqrt_bincontent1D;
          outputFile1D << "Hist_name " << proj->GetName()  << "  "  << " Histogram_noofbins " << binsin1Dhist << "  \n";
          outputFile1D << "Hist_binwidth_x" << proj->GetXaxis()->GetBinWidth(n)<< "  ";
          //outputFile1D << "Histogram_binwidth_y" << proj->GetYaxis()->GetBinWidth(n)<< "  ";
          outputFile1D << "Bin_number " << proj->GetBin(n) << "  ";
          outputFile1D << "proj_cell_area " << proj_cell_area << "  ";
          outputFile1D << "x_bin_width_pt " << x_bin_width_pt << "  ";
          outputFile1D << "y_bin_width_pz " << y_bin_width_pz << "  ";

          outputFile1D << "Bin_Content " << proj->GetBinContent(n) << "  " ;
          outputFile1D << "Sqrt_Bin_Content " << sqrt_bincontent1D << "  " ;
          outputFile1D << "Bin_Error_ROOT " << proj->GetBinError(n) << "  \n ";
          
          outputFile1DALL << "Histogram_name " << proj->GetName() << "  "  << " Histogram_noofbins " << binsin1Dhist << "  \n";
          outputFile1DALL << "Histogram_binwidth " <<proj->GetXaxis()->GetBinWidth(n)<< "  ";
          outputFile1DALL << "Bin_number " << proj->GetBin(n) << "  ";
          outputFile1DALL << "proj_cell_area " << proj_cell_area << "  ";
          outputFile1DALL << "x_bin_edge " << x_bin_width_pt << "  ";
          outputFile1DALL << "y_bin_edge " << y_bin_width_pz << "  ";

          outputFile1DALL << "Bin_Content " << proj->GetBinContent(n) << "  " ;
          outputFile1DALL << "Sqrt_Bin_Content " << sqrt_bincontent1D << "  " ;
          outputFile1DALL << "Bin_Error_ROOT " << proj->GetBinError(n) << "  \n";
         

          }
        }
        double sqrt_bincontent = std::sqrt(ptpzproj->GetBinContent(x + 1, y + 1));
        double standard_error_bincontent = 1.0/sqrt_bincontent;
        if (outputFile.is_open()) { // check if the file was opened successfully
          outputFile << " Bin number " << ptpzproj->GetBin(x + 1, y + 1) << "  ";
          outputFile << " Bin Content = " << ptpzproj->GetBinContent(x + 1, y + 1) << "  " ;
          outputFile << " Sqrt Bin Content = " << sqrt_bincontent << "  " ;
          outputFile << " Bin Error = " << ptpzproj->GetBinError(x + 1, y + 1) << "  ";
          outputFile << " Standard error (1/sqrt(n)) = " << standard_error_bincontent << "      \n";
        }
      
        outputFile1D.close();
        proj->Write();
        proj->SetDirectory(nullptr);
        delete proj;
      }

       
    }
    outputFile.close(); // close the file when done
    outputFile1DALL.close();
    ptpzproj->Scale(fScaleFactor,"WIDTH");  //get  root to divide out the bin width
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
    SplitWrite3D(f3DHist_CCCoh);

    SplitWrite3D(f3DHist_EAvail);
    SplitWrite3D(f3DHist_EAvail_CCQE);
    SplitWrite3D(f3DHist_EAvail_CC2p2h);
    SplitWrite3D(f3DHist_EAvail_CC1p1pi);
    SplitWrite3D(f3DHist_EAvail_CCDIS);
    SplitWrite3D(f3DHist_EAvail_Other);
    SplitWrite3D(f3DHist_EAvail_CCCoh);

    
    

    //SplitWrite3D(f3DHist_EAvail_pzbins_double);
    /*
    SplitWrite3D(f3DHist_EAvail_pzbins_uniform);
    SplitWrite3D(f3DHist_EAvail_pzbins_uniform_double);
    SplitWrite3D(f3DHist_EAvail_CCQE_pzbins_uniform);
    SplitWrite3D(f3DHist_EAvail_CC2p2h_pzbins_uniform);
    SplitWrite3D(f3DHist_EAvail_CC1p1pi_pzbins_uniform);
    SplitWrite3D(f3DHist_EAvail_CCDIS_pzbins_uniform);
    SplitWrite3D(f3DHist_EAvail_Other_pzbins_uniform);


    SplitWrite3D(f3DHist_EAvail_CCQE_pzbins_uniform_double);
    SplitWrite3D(f3DHist_EAvail_CC2p2h_pzbins_uniform_double);
    SplitWrite3D(f3DHist_EAvail_CC1p1pi_pzbins_uniform_double);
    SplitWrite3D(f3DHist_EAvail_CCDIS_pzbins_uniform_double);
    SplitWrite3D(f3DHist_EAvail_Other_pzbins_uniform_double);
    */
    SplitWrite3D(f3DHist_EAvail_pzbins_double);
    SplitWrite3D(f3DHist_EAvail_CCQE_pzbins_double);
    SplitWrite3D(f3DHist_EAvail_CC2p2h_pzbins_double);
    SplitWrite3D(f3DHist_EAvail_CC1p1pi_pzbins_double);
    SplitWrite3D(f3DHist_EAvail_CCDIS_pzbins_double);
    SplitWrite3D(f3DHist_EAvail_Other_pzbins_double);
    SplitWrite3D(f3DHist_EAvail_CCCoh_pzbins_double);
    
   

    SplitWrite3D(f3DHist_EAvail_lowW);
    SplitWrite3D(f3DHist_EAvail_midW);
    SplitWrite3D(f3DHist_EAvail_highW);

    // Oscillation-window subsample: writes the 3D hists, the (pT,pz)
    // projections and the per-slice E_avail histograms
    // (e.g. DUNE_3D_MINERvALike_EAvail_OscWin_CCQE_x3_y10)
    for (auto &h : fOscWin) SplitWrite3D(h);



    Create2DHistogram(f3DHist_EAvail);
    Create2DHistogram(f3DHist_EAvail_CCQE);
    Create2DHistogram(f3DHist_EAvail_CC2p2h);
    Create2DHistogram(f3DHist_EAvail_CC1p1pi);
    Create2DHistogram(f3DHist_EAvail_CCDIS);
    Create2DHistogram(f3DHist_EAvail_Other);
    Create2DHistogram(f3DHist_EAvail_CCCoh);
    
    f3DHist_EAvail->Write();
    f3DHist_EAvail_CCQE->Write();
    f3DHist_EAvail_CC2p2h->Write();
    f3DHist_EAvail_CC1p1pi->Write();
    f3DHist_EAvail_CCDIS->Write();
    f3DHist_EAvail_Other->Write();
    f3DHist_EAvail_CCCoh->Write();
    f1DHist_EAvail_true_nuenergy->Write();

    f3DHist_EAvail_pzbins_double->Write();
    f3DHist_EAvail_CCQE_pzbins_double->Write();
    f3DHist_EAvail_CC2p2h_pzbins_double->Write();
    f3DHist_EAvail_CC1p1pi_pzbins_double->Write();
    f3DHist_EAvail_CCDIS_pzbins_double->Write();
    f3DHist_EAvail_Other_pzbins_double->Write();
    f3DHist_EAvail_CCCoh_pzbins_double->Write();

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
    f3DHist_CCCoh->Reset();

    f3DHist_EAvail->Reset();
    f3DHist_EAvail_CCQE->Reset();
    f3DHist_EAvail_CC2p2h->Reset();
    f3DHist_EAvail_CC1p1pi->Reset();
    f3DHist_EAvail_CCDIS->Reset();
    f3DHist_EAvail_Other->Reset();
    f3DHist_EAvail_CCCoh->Reset();

    f1DHist_EAvail_true_nuenergy->Reset();

    /*
    f3DHist_EAvail_pzbins_uniform->Reset();
    f3DHist_EAvail_CCQE_pzbins_uniform->Reset();
    f3DHist_EAvail_CC2p2h_pzbins_uniform->Reset();
    f3DHist_EAvail_CC1p1pi_pzbins_uniform->Reset();
    f3DHist_EAvail_CCDIS_pzbins_uniform->Reset();
    f3DHist_EAvail_Other_pzbins_uniform->Reset();

    f3DHist_EAvail_pzbins_uniform_double->Reset();
    f3DHist_EAvail_CCQE_pzbins_uniform_double->Reset();
    f3DHist_EAvail_CC2p2h_pzbins_uniform_double->Reset();
    f3DHist_EAvail_CC1p1pi_pzbins_uniform_double->Reset();
    f3DHist_EAvail_CCDIS_pzbins_uniform_double->Reset();
    f3DHist_EAvail_Other_pzbins_uniform_double->Reset();

    f3DHist_EAvail_pzbins_uniform_double->Reset();
    */

    f3DHist_EAvail_pzbins_double->Reset();
    f3DHist_EAvail_CCQE_pzbins_double->Reset();
    f3DHist_EAvail_CC2p2h_pzbins_double->Reset();
    f3DHist_EAvail_CC1p1pi_pzbins_double->Reset();
    f3DHist_EAvail_CCDIS_pzbins_double->Reset();
    f3DHist_EAvail_Other_pzbins_double->Reset();
    f3DHist_EAvail_CCCoh_pzbins_double->Reset();
    
    f3DHist_EAvail_lowW->Reset();
    f3DHist_EAvail_midW->Reset();
    f3DHist_EAvail_highW->Reset();

    for (auto &h : fOscWin) h->Reset();

  }
  
};