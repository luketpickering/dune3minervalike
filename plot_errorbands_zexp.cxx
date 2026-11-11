#include "TRatioPlot.h"
#include "TGraphAsymmErrors.h"
//#define FMT_HEADER_ONLY
//#include "fmt-10.2.0/include/fmt/format.h"
#include <sstream>

using namespace std;

// ===================== global style =====================

void SetPublicationStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // Fonts: 43 = Helvetica (precision 3 = pixel size via SetTextSize)
  gStyle->SetTextFont(43);
  gStyle->SetLabelFont(43, "xyz");
  gStyle->SetTitleFont(43, "xyz");

  // Axis label / title sizes (in pixels, robust across canvas sizes)
  gStyle->SetLabelSize(18, "xyz");
  gStyle->SetTitleSize(20, "xyz");
  gStyle->SetTitleOffset(1.3, "y");
  gStyle->SetTitleOffset(1.1, "x");

  // Ticks on all four sides
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Clean pad
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadBottomMargin(0.13);

  // No border around the legend
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(16);
}

// ===================== colour palette =====================
// Colour-blind friendly (Paul Tol "Bright" set)
int cols[] = {
  TColor::GetColor("#000000"),  // 0 – black  (total)
  TColor::GetColor("#4477AA"),  // 1 – blue   (CCQE)
  TColor::GetColor("#EE6677"),  // 2 – red    (CC1p1pi)
  TColor::GetColor("#228833"),  // 3 – green  (CC2p2h)
  TColor::GetColor("#CCBB44"),  // 4 – yellow (CCDIS)
  TColor::GetColor("#66CCEE"),  // 5 – cyan   (Other)
  TColor::GetColor("#AA3377"),  // 6 – purple (z-exp PCA)
};

// ===================== drawing helpers =====================

// Shared axis formatting applied once to the "leading" histogram
void FormatAxes(TH1* h, const char* xtitle, const char* ytitle) {
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->GetXaxis()->SetTitleFont(43);
  h->GetXaxis()->SetLabelFont(43);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetLabelFont(43);
  h->GetXaxis()->SetTitleSize(20);
  h->GetXaxis()->SetLabelSize(18);
  h->GetYaxis()->SetTitleSize(20);
  h->GetYaxis()->SetLabelSize(18);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.5);
}

void DrawTotal(TH1* h, TLegend* leg, double ymax) {
  h->SetLineWidth(2);
  h->SetLineColor(cols[0]);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
  h->SetMarkerColor(cols[0]);
  h->GetYaxis()->SetRangeUser(0, 1.25 * ymax);
  // 2.0 lands exactly on a real EAvail_bins edge (..., 1.6, 2.0, 4.0, 8.0);
  // anything between 2.0 and 4.0 falls inside that single wide bin, so ROOT
  // would round outward and show the whole bin up to 4.0 regardless.
  h->GetXaxis()->SetRangeUser(0, 2.0);
  leg->AddEntry(h, "Total (stat. unc.)", "lep");
  h->Draw("E1");
}

void DrawComponent(TH1* h, TLegend* leg, int col,
                   const char* label) {
  h->SetLineWidth(2);
  h->SetLineStyle(1);
  h->SetLineColor(col);
  h->SetFillStyle(0);
  leg->AddEntry(h, label, "l");
  h->Draw("HIST SAME");
}

// Draws a filled band running from (hCentral - hSpread) to (hCentral + hSpread)
// at each bin, centered exactly on hCentral's own value -- NOT on inherited
// statistical errors, which is what an older two-histogram E2 approach
// incorrectly did, causing the band to float off the CCQE line.
TGraphAsymmErrors* DrawPCABand(TH1* hCentral, TH1* hSpread, int col, TLegend* leg) {
  int n = hCentral->GetNbinsX();
  // Two points per bin (left edge, right edge) with the same content, so the
  // band steps flat across each bin instead of drawing a diagonal line between
  // bin centers, which is what made it look pointy/jagged.
  auto* g = new TGraphAsymmErrors(2 * n);
  for (int b = 1; b <= n; ++b) {
    double xlo  = hCentral->GetBinLowEdge(b);
    double xhi  = hCentral->GetBinLowEdge(b) + hCentral->GetBinWidth(b);
    double y    = hCentral->GetBinContent(b);
    double err  = hSpread->GetBinContent(b);
    g->SetPoint(2 * (b - 1),     xlo, y);
    g->SetPoint(2 * (b - 1) + 1, xhi, y);
    g->SetPointError(2 * (b - 1),     0, 0, err, err);
    g->SetPointError(2 * (b - 1) + 1, 0, 0, err, err);
  }
  g->SetFillColorAlpha(col, 0.30);
  g->SetLineWidth(0);
  leg->AddEntry(g, "CCQE z-exp PCA (1#sigma)", "f");
  g->Draw("3 SAME");
  return g;
}

// Draw a DUNE / MINERvA style preliminary watermark
void DrawWatermark(const char* expt = "DUNE Simulation",
                   const char* prelim = "Preliminary") {
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(73);   // bold italic Helvetica
  latex.SetTextSize(18);
  latex.SetTextColor(kGray + 1);
  latex.DrawLatex(0.14, 0.885, expt);

  latex.SetTextFont(43);
  latex.SetTextSize(16);
  latex.SetTextColor(kGray + 1);
  latex.DrawLatex(0.14 + 0.20, 0.887, prelim);
}

// Draw the kinematic slice label in the top-right corner
void DrawSliceLabel(double pT_lo, double pT_hi,
                    double pZ_lo, double pZ_hi) {
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(15);
  latex.SetTextAlign(31);   // right-aligned

  char s_buf[256];
  snprintf(s_buf, sizeof(s_buf),
      "%.2f < #it{p}_{T} < %.2f GeV/c,  %.2f < #it{p}_{Z} < %.2f GeV/c",
      pT_lo, pT_hi, pZ_lo, pZ_hi);
  std::string s = s_buf;
  latex.DrawLatex(0.94, 0.92, s.c_str());
}

// ===================== main =====================
// NOTE: function renamed to match the file name (plot_errorbands_zexp.cxx)
// so that `root -l -b -q plot_errorbands_zexp.cxx` actually invokes it.
void plot_errorbands_zexp() {

  SetPublicationStyle();

  // ---------- files ----------
  TFile* fCV = TFile::Open("/home/abipeake/scratch/dune3minervalike/DUNE_3D_MINERvALike.out_throw_nominal_new.root");
  TFile* fb1 = TFile::Open("/home/abipeake/scratch/dune3minervalike/DUNE_3D_MINERvALike.out_throw_b1_new.root");
  TFile* fb2 = TFile::Open("/home/abipeake/scratch/dune3minervalike/DUNE_3D_MINERvALike.out_throw_b2_new.root");
  TFile* fb3 = TFile::Open("/home/abipeake/scratch/dune3minervalike/DUNE_3D_MINERvALike.out_throw_b3_new.root");
  TFile* fb4 = TFile::Open("/home/abipeake/scratch/dune3minervalike/DUNE_3D_MINERvALike.out_throw_b4_new.root");

  auto h3D = fCV->Get<TH3>("DUNE_3D_MINERvALike_EAvail");
  const int npt = h3D->GetXaxis()->GetNbins();
  const int npz = h3D->GetYaxis()->GetNbins();

  // ---------- output PDF ----------
  const std::string outpdf =
    "Highstatzepnew_pub.pdf";

  // Use a square-ish canvas; 600×500 prints cleanly at journal width
  TCanvas c1("c1", "", 700, 560);
  c1.Print((outpdf + "(").c_str());

  // ---------- loops ----------
  for (int iy = 0; iy < npz; ++iy) {
    for (int ix = 0; ix < npt; ++ix) {

      char dist_buf[32];
      // NOTE: must match SplitWrite3D's std::to_string(x)/(y) naming exactly,
      // which does NOT zero-pad (x0, x1, ..., x9, x10, ...) -- using %02d here
      // previously caused every slice with ix<10 or iy<10 to be silently
      // skipped as "missing", since the zero-padded name never matched.
      snprintf(dist_buf, sizeof(dist_buf), "_x%d_y%d;1", ix, iy);
      std::string dist = dist_buf;

      // --- fetch histograms ---
      auto hTot   = fCV->Get<TH1>(("DUNE_3D_MINERvALike_EAvail"         + dist).c_str());
      auto hCCQE  = fCV->Get<TH1>(("DUNE_3D_MINERvALike_EAvail_CCQE"    + dist).c_str());
      auto hCC1pi = fCV->Get<TH1>(("DUNE_3D_MINERvALike_EAvail_CC1p1pi" + dist).c_str());
      auto hCC2p2h= fCV->Get<TH1>(("DUNE_3D_MINERvALike_EAvail_CC2p2h"  + dist).c_str());
      auto hDIS   = fCV->Get<TH1>(("DUNE_3D_MINERvALike_EAvail_CCDIS"   + dist).c_str());

      auto hCCQE_b1 = fb1->Get<TH1>(("DUNE_3D_MINERvALike_EAvail_CCQE" + dist).c_str());
      auto hCCQE_b2 = fb2->Get<TH1>(("DUNE_3D_MINERvALike_EAvail_CCQE" + dist).c_str());
      auto hCCQE_b3 = fb3->Get<TH1>(("DUNE_3D_MINERvALike_EAvail_CCQE" + dist).c_str());
      auto hCCQE_b4 = fb4->Get<TH1>(("DUNE_3D_MINERvALike_EAvail_CCQE" + dist).c_str());

      if (!hTot || !hCCQE || !hCC1pi || !hCC2p2h || !hDIS ||
          !hCCQE_b1 || !hCCQE_b2 || !hCCQE_b3 || !hCCQE_b4) {
        std::cerr << "Skipping missing slice: " << dist << "\n";
        continue;
      }

      if (hTot->Integral() <= 0) {
        std::cerr << "Skipping empty slice: " << dist << "\n";
        continue;
      }

      // --- PCA uncertainty band (quadrature sum of throws) ---
      auto hpca = (TH1D*) hCCQE->Clone("hpca");
      hpca->Reset();

      for (TH1* hb : {hCCQE_b1, hCCQE_b2, hCCQE_b3, hCCQE_b4}) {
        auto diff = (TH1D*) hCCQE->Clone("tmp");
        diff->Add(hb, -1);
        for (int b = 1; b <= diff->GetNbinsX(); ++b) {
          double d = diff->GetBinContent(b);
          hpca->AddBinContent(b, d * d);
        }
        delete diff;
      }
      for (int b = 1; b <= hpca->GetNbinsX(); ++b)
        hpca->SetBinContent(b, std::sqrt(hpca->GetBinContent(b)));

      // --- draw ---
      c1.Clear();
      c1.cd();

      double ymax = hTot->GetMaximum();

      // Legend: top-right, compact
      TLegend leg(0.55, 0.60, 0.93, 0.88);
      leg.SetNColumns(1);

      DrawTotal(hTot, &leg, ymax);
      FormatAxes(hTot,
                 "E^{Had.}_{Avail.} [GeV]",
                 "d#sigma / dE^{Had.}_{Avail.} [cm^{2}/GeV/nucleon]");

      auto pcaGraph = DrawPCABand(hCCQE, hpca, cols[1], &leg);
      DrawComponent(hCCQE,   &leg, cols[1], "CCQE");
      DrawComponent(hCC1pi,  &leg, cols[2], "CC1p1#pi");
      DrawComponent(hCC2p2h, &leg, cols[3], "CC2p2h");
      DrawComponent(hDIS,    &leg, cols[4], "CCDIS");

      leg.Draw();

      // Slice boundaries
      double pT_lo = h3D->GetXaxis()->GetBinLowEdge(ix + 1);
      double pT_hi = h3D->GetXaxis()->GetBinUpEdge (ix + 1);
      double pZ_lo = h3D->GetYaxis()->GetBinLowEdge(iy + 1);
      double pZ_hi = h3D->GetYaxis()->GetBinUpEdge (iy + 1);

      DrawWatermark("DUNE Simulation", "Preliminary");
      DrawSliceLabel(pT_lo, pT_hi, pZ_lo, pZ_hi);

      c1.Update();
      c1.Print(outpdf.c_str());

      delete hpca;
      delete pcaGraph;

    }  // ix
  }  // iy

  c1.Print((outpdf + ")").c_str());
}