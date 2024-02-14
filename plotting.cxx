#include "TRatioPlot.h"
//some constants
double NDFHC_IntFlux = 0.0010489263;
  double NDRHC_IntFlux =  0.00094423594;
double NucleonTonneScale = 6.02831;

//double_t cross_section_toeventrate =  NDFHC_IntFlux * 1E-4 *  NucleonTonneScale * 1.1E12;
          // Scale by Ar atoms * 1 year POT * 1E-38
        
//double_t cross_section_toeventrate = 8.85e49; //conversion factor for cross-section to event rate;
double_t cross_section_toeventrate = 1.610845e52; // 8.85e49^2/3 to account for 2D plane of 3D detector

// useful resources for color choices
//  https://personal.sron.nl/~pault/
int cols[] = {TColor::GetColor("#000000"), TColor::GetColor("#0077bb"),
              TColor::GetColor("#ee3377"), TColor::GetColor("#aa4499"),
               TColor::GetColor("#ee99aa"), TColor::GetColor("#009988")};

double fontsize = 0.05;
std::string componets[4] =  {"Total", "CCQE", "CCQE (b1-sigma ,b2,b3,b4)","CCQE (b1+sigma ,b2,b3,b4)"};

std::string bins[11] = {"[0 < P_{t} < 0.075]",	"[0.075 < P_{t} < 0.15]",	"[0.15  < P_{t} < 0.25]",	"[0.25 < P_{t} < 0.325]",
                        "[0.325 < P_{t} <0 .4]",	"[0.4 < P_{t} < 0.475]",  "[0.465 < P_{t} < 0.55]", "[0.55 <  P_{t} < 0.7]" ,	
                        "[0.7 < P_{t} < 0.85]" , "[0.85 < P_{t} < 1.0]"	, " [1.0 < P_{t} < 2.5]"};

std::string bins_pz[11] = { "[1.5 < P_{z} < 3.5]", "[3.5 <  P_{z} < 4.5]" ,	
                        "[4.5 < P_{z} < 7]" , "[7.0 < P_{z} < 8.0]"	, " [8.0 < P_{z} < 10]", " [10.0 < P_{z} < 20]"};

// [LOG Minmzr]:- |-> MicroBooNE_CC1Mu1p_XSec_1DDeltaPT_nu : 9.84183/13 [LOG
// Minmzr]:- |-> MicroBooNE_CC1Mu1p_XSec_1DMuonMomentum_nu         : 15.7307/10
// [LOG Minmzr]:- |-> MicroBooNE_CC1Mu1p_XSec_1DECal_nu : 14.1948/9
 

std::map<std::string, std::string> data_release_mc_legends = {
    {"DUNE_3D_MINERvALike_x0_y1", "DUNE"},
    {"DUNE_3D_MINERvALike_CCQE_x0_y1", "CCQE "},
    {"DUNE_3D_MINERvALike_CC2p2h_x0_y1",
     "CC2p2h"},
    {"DUNE_3D_MINERvALike_CC1pi_x0_y1",
     "CC1pi"},
     {"DUNE_3D_MINERvALike_CCDIS_x0_y1",
     "CCDIS"},
     {"DUNE_3D_MINERvALike_CCOther_x0_y1",
     "CCOther"},

     {"DUNE_3D_MINERvALike_x0_y1", "DUNE"},
    {"DUNE_3D_MINERvALike_CCQE_x0_y1", "CCQE "},
    {"DUNE_3D_MINERvALike_CC2p2h_x0_y1",
     "CC2p2h"},
    {"DUNE_3D_MINERvALike_CC1pi_x0_y1",
     "CC1pi"},
     {"DUNE_3D_MINERvALike_CCDIS_x0_y1",
     "CCDIS"},
     {"DUNE_3D_MINERvALike_CCOther_x0_y1",
     "CCOther"},


     {"DUNE_3D_MINERvALike_x0_y1", "DUNE"},
    {"DUNE_3D_MINERvALike_CCQE_x0_y1", "CCQE "},
    {"DUNE_3D_MINERvALike_CC2p2h_x0_y1",
     "CC2p2h"},
    {"DUNE_3D_MINERvALike_CC1pi_x0_y1",
     "CC1pi"},
     {"DUNE_3D_MINERvALike_CCDIS_x0_y1",
     "CCDIS"},
     {"DUNE_3D_MINERvALike_CCOther_x0_y1",
     "CCOther"},

/*
    {"MicroBooNE_CC1Mu1p_XSec_1DECal_nu_MC",
     "GENIE AR23i, #Chi^{2}_{CV} = 9.8/9"},
    {"Overlay9NuWro_ECalPlot", "NEUT, #Chi^{2}_{CV} = 10.2/9"},
    {"GiBUU_ECalPlot", "GiBUU, #Chi^{2}_{CV} = 6.6/9"},
    {"NEUT_ECalPlot", "NuWro, #Chi^{2}_{CV} = 9.9/9"},

    {"MicroBooNE_CC1Mu1p_XSec_1DMuonMomentum_nu_MC",
     "GENIE AR23i, #Chi^{2}_{CV} = 15.7/10"},
    {"Overlay9NuWro_MuonMomentumPlot", "NEUT"},
    {"GiBUU_MuonMomentumPlot", "GiBUU"},
    {"NEUT_MuonMomentumPlot", "NuWro"},
    */
};

std::map<std::string, std::pair<std::string, std::string>> data_release_titles =
    {
        {"DUNE_3D_MINERvALike_x0_y1",
         {"#frac{d#sigma}{d#delta#it{p}_{T}} [10^{-38} cm^{2} GeV/#it{c} /Ar] ",
          "Visible Energy [GeV/#it{c}] "}},
};


void DrawData(TH1 *data, TLegend *leg, std::string opt = ""){

  data->SetMarkerColor(kBlack);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(2);
  data->SetLineColor(kBlack);
  data->SetLineWidth(3);

  //data->GetYaxis()->SetRangeUser(0, data_release_maxima[data->GetName()]);
  data->GetYaxis()->SetLabelSize(1.1*fontsize);
  data->GetYaxis()->SetTitleSize(1.1* fontsize);
  data->GetYaxis()->SetTitleOffset(1.2);
  data->GetYaxis()->SetTitleFont(132);
  data->GetYaxis()->SetLabelFont(132);
  data->GetYaxis()->SetNdivisions(505);
  data->GetYaxis()->SetTitle(
      data_release_titles[data->GetName()].first.c_str());

  data->GetXaxis()->SetLabelSize(1.1* fontsize);
  data->GetXaxis()->SetTitleSize(1.1*fontsize);
  data->GetXaxis()->SetLabelFont(132);
  data->GetXaxis()->SetTitleFont(132);
  data->GetXaxis()->SetTitleOffset(1);
  data->GetXaxis()->SetNdivisions(505);
  data->GetXaxis()->SetTitle(
      data_release_titles[data->GetName()].second.c_str());

  if (!opt.size()) {
    leg->AddEntry(data, "", "lp");
  }
  data->Draw(opt.c_str());
}

void DrawStatisticalUncertainty(TH1 *mc, TLegend *leg , std::string name) {
  mc->SetLineWidth(0);
  mc->SetLineColorAlpha(cols[0], 0);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(cols[0], 0);
  mc->SetFillColorAlpha(cols[0], 0.1);
  mc->SetStats(0);
  leg->AddEntry(mc, name.c_str() , "f");
  mc->Draw("ehist");
}


void DrawMCBand(TH1 *mc, TLegend *leg , std::string name) {
  mc->SetLineWidth(0);
  mc->SetLineColorAlpha(cols[1], 0);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(cols[1], 0);
  mc->SetFillColorAlpha(cols[1], 0.5);
  mc->SetStats(0);
  leg->AddEntry(mc, name.c_str() , "f");
  mc->Draw("E2SAME");
}


/*
void DrawStatisticalUncertaintyBand(TH1 *mc, TLegend *leg) {
  mc->SetLineWidth(0);
  mc->SetLineColorAlpha(TColor::GetColor("#000000"), 0);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(TColor::GetColor("#000000"), 0);
  mc->SetFillColorAlpha(TColor::GetColor("#000000"), 0.5);

  leg->AddEntry(mc, "1 #sigma (stat)", "f");
  mc->Draw("E2SAME");
}
*/


int gcol = TColor::GetColor("#bbbbbb");
void DrawNonQE(TH1 *mc, TLegend *leg) {

  mc->SetLineWidth(0);
  mc->SetLineColorAlpha(gcol, 0);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(gcol, 0);
  mc->SetFillColorAlpha(gcol, 0.5);

  leg->AddEntry(mc, "Non-CCQE Contribution", "f");
  mc->Draw("HISTSAME");
}

void DrawMCCV(TH1 *mc, TLegend *leg, int col ,  std::string newhnistname) {
  mc->SetLineWidth(4);
  mc->SetLineColorAlpha(col, 1);
  mc->SetLineStyle(1);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(col, 0);
  mc->SetFillColorAlpha(col, 0);

  mc->SetStats(0);


  std::cout << mc->GetName() << " -> " << data_release_mc_legends[mc->GetName()]
            << std::endl;
  //mc->SetName("");
  //leg->AddEntry(mc, data_release_mc_legends[mc->GetName()].c_str(), "l");
  
  mc->SetName(newhnistname.c_str());
  
  mc->SetTitle(""); 
  //mc->GetYaxis()->SetTitle("cm^{2}/GeV^{3}/nucleon");
  mc->GetYaxis()->SetTitle("cm^{2}/GeV^{3}/POT");
  leg->AddEntry(mc, mc->GetName(), "l");
  
  mc->Draw("HISTSAME");
  

}

void plotting() {
  double event_rate_total;
  TFile *fdatar = TFile::Open("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/b1b2b3bb4all0.root");
  TFile *fdatar_plus1sigma = TFile::Open("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/b1plus.root");
  TFile *fdatar_minus1sigma = TFile::Open("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/b1minus.root");
  //events/dunesyst.root");
  //fdatar->cd("nominal_throw");
  //TFile *fthrows = TFile::Open("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/");
  //events/dunesyst.root");
     int i  = 0;
     int j=0;
  for (auto dist : {"_x0_y0;1", "_x1_y0;1", "_x2_y0;1", "_x3_y0;1", "_x4_y0;1", "_x5_y0;1", "_x6_y0;1", "_x7_y0;1","_x8_y0;1", "_x9_y0;1", "_x10_y0;1",
                    }){
        
        /*"_x0_y1;1", "_x1_y1;1", "_x2_y1;1", "_x3_y1;1", "_x4_y1;1", "_x5_y1;1", "_x6_y1;1", "_x7_y1;1","_x8_y1;1", "_x9_y1;1", "_x10_y1;1",
                    "_x0_y2;1", "_x1_y2;1", "_x2_y2;1", "_x3_y2;1", "_x4_y2;1", "_x5_y2;1", "_x6_y2;1", "_x7_y2;1","_x8_y2;1", "_x9_y2;1", "_x10_y2;1",
                    "_x0_y3;1", "_x1_y3;1", "_x2_y3;1", "_x3_y3;1", "_x4_y3;1", "_x5_y3;1", "_x6_y3;1", "_x7_y3;1","_x8_y3;1", "_x9_y3;1", "_x10_y3;1",
                    "_x0_y4;1", "_x1_y4;1", "_x2_y4;1", "_x3_y4;1", "_x4_y4;1", "_x5_y4;1", "_x6_y4;1", "_x7_y4;1","_x8_y4;1", "_x9_y4;1", "_x10_y4;1",
                    "_x0_y5;1", "_x1_y5;1", "_x2_y5;1", "_x3_y5;1", "_x4_y5;1", "_x5_y5;1", "_x6_y5;1", "_x7_y5;1","_x8_y5;1", "_x9_y5;1", "_x10_y5;1",*/
                     
    //, "", "_CCQE", "_CC2p2h", "_CC1pi", "_CCDIS", "_CCOther"
        //std::string hname = ("DUNE_3d_MINERvALike_" , dist ,"_x0_y1").c_str();
        //hname += hname + std::string(dist , "_x0_y1").c_str();
        //std::cout<<hname<<std::endl;
        std::string distribution = dist;
        std::string hname = "DUNE_3D_MINERvALike" + distribution ;
        //hname += dist.c_str() , "_x0_y1";
        //get hist
        
        auto fdata =
            fdatar->Get<TH1>(hname.c_str());

            if(!fdata){ std::cout << "Is the name: " << hname << " definitely correct? I couldn't find that hist in the file." << std::endl;}
            if(!fdata) { std::cout << "Failed to read hist <bla>" << std::endl; } 
        fdata->SetName((hname.c_str()));
        //fdata->Scale(cross_section_toeventrate , "width");
        std::cout<< " fdata->Integral(width) === " <<  fdata->Integral("width") <<std::endl;
        
        std::string hname2 = "DUNE_3D_MINERvALike_CCQE"+ distribution  + "Plot";
        if(!fdata){ std::cout << "Is the name: " << hname2 << " definitely correct? I couldn't find that hist in the file." << std::endl;}
        auto fdatar_mc_CCQE = fdatar->Get<TH1>(hname2.c_str());
        fdatar_mc_CCQE->SetName(hname2.c_str());


        std::string hname2_plus1sigma = "DUNE_3D_MINERvALike_CCQE"+ distribution  + "Plot";
        if(!fdatar_plus1sigma){ std::cout << "Is the name: " << hname2_plus1sigma << " definitely correct? I couldn't find that hist in the file." << std::endl;}
        auto fdatar_mc_CCQE_plus1sigma = fdatar_plus1sigma ->Get<TH1>(hname2_plus1sigma .c_str());
        fdatar_mc_CCQE_plus1sigma ->SetName(hname2_plus1sigma .c_str());

        std::string hname2_minus1sigma = "DUNE_3D_MINERvALike_CCQE"+ distribution  + "Plot";
        if(!fdatar_minus1sigma){ std::cout << "Is the name: " << hname2_minus1sigma << " definitely correct? I couldn't find that hist in the file." << std::endl;}
        auto fdatar_mc_CCQE_minus1sigma = fdatar_minus1sigma ->Get<TH1>(hname2_minus1sigma .c_str());
        fdatar_mc_CCQE_minus1sigma ->SetName(hname2_minus1sigma .c_str());
       
        //fdatar_mc_CCQE->Scale(cross_section_toeventrate , "width");

        std::string hname3 = "DUNE_3D_MINERvALike_CC2p2h" + distribution +  "Plot";
        if(!fdata){ std::cout << "Is the name: " << hname3 << " definitely correct? I couldn't find that hist in the file." << std::endl;}
        auto fdatar_mc_CC2p2h = fdatar->Get<TH1>(hname3.c_str());
        fdatar_mc_CC2p2h->SetName(hname3.c_str());
        //fdatar_mc_CC2p2h->Scale(cross_section_toeventrate  , "width");

        std::string hname4 = "DUNE_3D_MINERvALike_CC1pi"  + distribution + "Plot";
        if(!fdata){ std::cout << "Is the name: " << hname4 << " definitely correct? I couldn't find that hist in the file." << std::endl;}
        auto fdatar_mc_CC1pi = fdatar->Get<TH1>(hname4.c_str());
        fdatar_mc_CC1pi->SetName(hname4.c_str());
       // fdatar_mc_CC1pi->Scale(cross_section_toeventrate , "width");

        std::string hname5 = "DUNE_3D_MINERvALike_CCDIS" + distribution +  "Plot";
        if(!fdata){ std::cout << "Is the name: " << hname5 << " definitely correct? I couldn't find that hist in the file." << std::endl;}
        auto fdatar_mc_CCDIS = fdatar->Get<TH1>(hname5.c_str());
        fdatar_mc_CCDIS->SetName(hname5.c_str());
       // fdatar_mc_CCDIS->Scale(cross_section_toeventrate, "width");
        
        std::string hname6 = "DUNE_3D_MINERvALike_CCOther" +  distribution + "Plot";
        if(!fdata){ std::cout << "Is the name: " << hname6 << " definitely correct? I couldn't find that hist in the file." << std::endl;}
        auto fdatar_mc_CCOther = fdatar->Get<TH1>(hname6.c_str());
        fdatar_mc_CCOther->SetName(hname6.c_str());
        //fdatar_mc_CCOther->Scale(cross_section_toeventrate  , "width");

        /*
        std::string throw_hist = "error_bands/DUNE_3D_MINERvALike_CCQE" + distribution  + "_prof" ;
        auto fmc_throws = fthrows->Get<TH1>((throw_hist).c_str());
        std::cout << throw_hist << std::endl;
        fmc_throws->GetYaxis()->SetTitle(fmc_throws->GetYaxis()->GetTitle());
        //fmc_throws->Scale(cross_section_toeventrate   , "width");

        */
        
        /*
        auto fdatar_mc_GiBUU =
            fdatar->Get<TH1>((std::string("GiBUU_") + dist + "Plot").c_str());
        fdatar_mc_GiBUU->SetName((std::string("GiBUU_") + dist + "Plot").c_str());
        auto fdatar_mc_neut =
            fdatar->Get<TH1>((std::string("NEUT_") + dist + "Plot").c_str());
        fdatar_mc_neut->SetName((std::string("NEUT_") + dist + "Plot").c_str());
        */
    /*
        auto fmc_throws = fthrows->Get<TH1>(
            (std::string("error_bands/MicroBooNE_CC1Mu1p_XSec_1D") + dist +
            "_nu_MC")
                .c_str());
        fdata->GetYaxis()->SetTitle(fmc_throws->GetYaxis()->GetTitle());
        fmc_throws->Scale(1E38);
        auto fmc_cv = fthrows->Get<TH1>(
            (std::string("nominal_throw/MicroBooNE_CC1Mu1p_XSec_1D") + dist +
            "_nu_MC")
                .c_str());
        fmc_cv->Scale(1E38);
        auto fmc_cv_qe = fthrows->Get<TH1>(
            (std::string("nominal_throw/MicroBooNE_CC1Mu1p_XSec_1D") + dist +
            "_nu_MODES_CCQE;"+mc_name_cycle[fdata->GetName()])
                .c_str());
        fmc_cv_qe->Scale(1E38);

        auto fmc_cv_nonqe = dynamic_cast<TH1 *>(fmc_cv->Clone(
            (std::string("MicroBooNE_CC1Mu1p_XSec_1D") + dist + "_nu_nonQE")
                .c_str()));
        fmc_cv_nonqe->Add(fmc_cv_qe, -1.0);
        */
    TLegend *legendl = new TLegend(0.3, 0.60, 0.8, 0.85);
    legendl->SetTextFont(132);
    legendl->SetTextSize(fontsize * 0.95);
    legendl->SetBorderSize(0);
    legendl->SetFillStyle(0);
    legendl->SetNColumns(1);

    /*
    TPaveText *pt = new TPaveText(0.01,0.1,0.2,0.1);
    pt->SetTextSize(0.04);
    pt->SetFillColor(0);
    */
    //spt->SetTextAlign(12);

    TLatex latex;
    latex.SetNDC();
    
    latex.SetTextSize(0.05); 
    

    TCanvas c1("c1", "", 1200, 1200);

    c1.SetLeftMargin(0.1);
    c1.SetRightMargin(0.1);
    c1.SetTopMargin(0.1);
    c1.SetBottomMargin(0.1);

  
    TPad pad1("pad1", "", 0, 0.3, 1, 1);
    pad1.Draw();
    pad1.cd();
    /*
    auto ratio_plot_plus = new TRatioPlot(fdatar_mc_CCQE_plus1sigma, fdatar_mc_CCQE);
    auto ratio_plot_minus = new TRatioPlot(fdatar_mc_CCQE_minus1sigma, fdatar_mc_CCQE);
      
    ratio_plot_plus->Draw();
    ratio_plot_plus->GetLowYaxis()->SetNdivisions(505);

    ratio_plot_minus->Draw();
    ratio_plot_minus->GetLowYaxis()->SetNdivisions(505);
    */
    

    DrawStatisticalUncertainty(fdata, legendl, "Statistical Uncertainty");
    //DrawMCBand(fmc_throws, legendl, "z-Expansion Uncertainty");---------------------------------------
    
    /*
    DrawMCCV(fmc_cv, legendl, cols[0]);
    DrawMCBand(fmc_throws, legendl);
    // DrawNonQE(fmc_cv_nonqe, legendl);
    DrawMCCV(fdatar_mc_neut, legendl, cols[1]);
    DrawMCCV(fdatar_mc_nuwro, legendl, cols[2]);
    DrawMCCV(fdatar_mc_GiBUU, legendl, cols[3]);

    */
    //DrawStatisticalUncertaintyBand(fdata, legendl);
    DrawMCCV(fdata, legendl, cols[0] ,componets[0]);
    pad1.Update();
    
    DrawMCCV(fdatar_mc_CCQE, legendl, cols[1],componets[1]);
    pad1.Update();
    DrawMCCV(fdatar_mc_CCQE_plus1sigma, legendl, cols[2],componets[2]);
    pad1.Update();
    DrawMCCV(fdatar_mc_CCQE_minus1sigma, legendl, cols[3],componets[3]);
    pad1.Update();
    //DrawMCCV(fdatar_mc_CC1pi, legendl, cols[2],componets[2]);
   //DrawMCCV(fdatar_mc_CC2p2h, legendl, cols[3],componets[3]);
    

    //ratio_plot_plus->Draw();
    //ratio_plot_minus->Draw();
    //auto 3d_hist_total =  DUNE_3D_MINERvALike
    double_t integral =  fdata->Integral("width") ;
    std::cout<<"integral= " << integral<<std::endl;
    
    //DrawMCCV(fdatar_mc_CCDIS, legendl, cols[4],componets[4]);
    //DrawMCCV(fdatar_mc_CCOther, legendl, cols[5], componets[5]);
   
    //DrawData(fdata, legendl, "SAMES");
    std::cout<<bins[i].c_str()<<std::endl;
    
    legendl->Draw();
    //pt->AddText(bins[i].c_str());
    //pt->Draw();
    latex.DrawLatexNDC(0.14,0.94,"[1.5 < P_{z} < 3.5]");
    latex.DrawLatexNDC(0.55,0.94,bins[i].c_str());
    //latex.DrawLatexNDC(0.55,0.84,bins_pz[j].c_str());
    //latex.DrawLatexNDC(0.6,0.81,"N_{events} = 1.97 x 10^{14} y^{-1}");
    latex.SetTextSize(0.03); 
    latex.DrawLatexNDC(0.1,0.03,"#it{A Peake}");
    pad1.Update();
    c1.cd();
    TPad pad2("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2.Draw();
    pad2.cd(); //#pad2 becomes the current pad

    TH1D* ratioplus = (TH1D*)fdatar_mc_CCQE_plus1sigma ->Clone(); 
    ratioplus->Divide(fdatar_mc_CCQE);
    pad2.Update();
    TH1D* ratiominus = (TH1D*)fdatar_mc_CCQE_minus1sigma ->Clone(); 
    ratiominus->Divide(fdatar_mc_CCQE);
    ratioplus->Draw("p");
    ratiominus->Draw("p");
    ratioplus->SetTitle(""); 
    ratioplus->GetXaxis()->SetLabelSize(2* fontsize);
     ratioplus->GetXaxis()->SetTitleSize(2* fontsize); 
     ratioplus->GetYaxis()->SetLabelSize(2* fontsize); 
     ratioplus->GetYaxis()->SetTitleSize(2* fontsize); 
     ratioplus->GetYaxis()->SetTitle("Shift/CV"); 
     
     ratiominus->GetXaxis()->SetLabelSize(2* fontsize);
     ratiominus->GetXaxis()->SetTitleSize(2* fontsize); 
     ratiominus->GetYaxis()->SetLabelSize(2* fontsize); 
     ratiominus->GetYaxis()->SetTitleSize(2* fontsize); 
     ratiominus->GetYaxis()->SetTitle("Shift/CV"); 
    
      ratioplus->GetYaxis()->SetRangeUser(0.8,1.1);
      ratiominus->GetYaxis()->SetRangeUser(0.8,1.1);
  
    pad2.Update();
    c1.cd();
 
    c1.Update();
    c1.Draw();

    c1.SaveAs((std::string("/root/plots/varying_bparams/") + dist + "21M_b1_pm1sigma.pdf").c_str());


     //c1.cd();
      //c1.Update();
      //c1.Draw();

    //c1.Print((std::string("/root/plots/varying_bparams/") + dist + "21M_b1_pm1sigma.pdf").c_str());
    i = i+ 1;
    
    for(int i = 0; i < fdatar_mc_CCQE_plus1sigma->GetXaxis()->GetNbins(); ++i){ 
        std::cout << "bin " << i << " error = " << fdatar_mc_CCQE_plus1sigma->GetBinError(i+1) << std::endl;
         std::cout << "bin " << i << " content = " << fdatar_mc_CCQE_plus1sigma->GetBinContent(i+1) << std::endl;
        std::cout << "bin " << i << "  content / bin error = " << (fdatar_mc_CCQE_plus1sigma->GetBinContent(i+1)) /(fdatar_mc_CCQE_plus1sigma->GetBinError(i+1)) << std::endl;


    }
    for(int i = 0; i < fdatar_mc_CCQE->GetXaxis()->GetNbins(); ++i){ 
        std::cout << "bin " << i << " error = " << fdatar_mc_CCQE->GetBinError(i+1) << std::endl;
         std::cout << "bin " << i << " content = " << fdatar_mc_CCQE->GetBinContent(i+1) << std::endl;
        std::cout << "bin " << i << "  content / bin error = " << (fdatar_mc_CCQE->GetBinContent(i+1)) /(fdatar_mc_CCQE->GetBinError(i+1)) << std::endl;


    }
    
    //event_rate_total = event_rate_total * event_rate_perhist;

//* cross_section_toeventrate

   }
   //std::cout << "total event rate" << event_rate_total <<std::endl;

}

