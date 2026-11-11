#include "TRatioPlot.h"
#include "TGraphAsymmErrors.h"
//#define FMT_HEADER_ONLY
//#include "fmt-10.2.0/include/fmt/format.h"
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

// =====================================================================
//  Outputs (all written by one run of plot_errorbands_zexp_slices()):
//
//   <stem>_fixedPZ_scanPT.pdf      one slice per page: pZ bin fixed, step
//                                  through all pT, then next pZ (old order)
//   <stem>_fixedPT_scanPZ.pdf      one slice per page: pT bin fixed, step
//                                  through all pZ, then next pT
//   <stem>_grid_onePTperPage.pdf   one page per pT bin, all its pZ slices
//                                  side by side on that page
//   <stem>_grid_onePZperPage.pdf   one page per pZ bin, all its pT slices
//                                  side by side on that page
//   <stem>_withOscWin_*.pdf        the four PDFs above again, with the
//                                  Enu-window overlay (grey in-window total,
//                                  dashed in-window components)
//   <stem>_oscwindow_fixedPT_scanPZ.pdf
//   <stem>_oscwindow_grid_onePTperPage.pdf
//   <stem>_oscwindow_grid_onePZperPage.pdf
//                                  the Enu-window subsample on its own: its
//                                  total and whichever in-window components
//                                  exist (kEnuWinLo < Enu < kEnuWinHi)
//
//  The coherent component and the Enu-window histograms must be written by
//  the NUISANCE sample; if they are absent the script says so once and
//  draws everything else as before.
//
//  Memory note: nothing read from a TFile is ever drawn directly. Each pad
//  draws private copies (PadCopy) that the pad owns and deletes when it is
//  cleared, so no object is shared between pads or canvases.
// =====================================================================

// ===================== configuration =====================

const std::string kDir     = "/home/abipeake/scratch/dune3minervalike/";
const std::string kPrefix  = "DUNE_3D_MINERvALike_EAvail";
const std::string kOutStem = "Highstatzepnew_pub";

// Neutrino-energy window around the DUNE first oscillation maximum
// (1300 km, dm2 ~ 2.5e-3 eV^2: first max ~2.6 GeV, minimum below it ~1.3 GeV).
// Window histograms carry kWinSuffix after the component name, e.g.
//   DUNE_3D_MINERvALike_EAvail_CCQE_OscWin_x3_y5
//   DUNE_3D_MINERvALike_EAvail_OscWin_x3_y5          (window total)
// The Enu cut itself is applied when the sample is filled; the two numbers
// below are only used for labels, so keep them in sync with the sample.
// kDrawWindow = false switches the Enu-window overlay (grey in-window total,
// dashed in-window components) off everywhere. main() writes the
// full-sample PDFs once without the overlay and once with it.
const bool        kDrawWindow = true;
const std::string kWinSuffix  = "_OscWin";
const double      kEnuWinLo   = 1.3;  // GeV
const double      kEnuWinHi   = 2.6;  // GeV

// x-axis range of every plot (E_avail, GeV)
const double kXMin = 0.0;
const double kXMax = 2.0;

const char* kXTitle ="E^{Had.}_{Avail.} [GeV]";
const char* kYTitle = "d#sigma / dE^{Had.}_{Avail.} [cm^{2}/GeV/nucleon]";

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
  TColor::GetColor("#66CCEE"),  // 5 – cyan   (CC coherent)
  TColor::GetColor("#AA3377"),  // 6 – purple (z-exp PCA)
};

// ===================== components =====================
// Optional components are drawn when their histogram exists and left out
// otherwise, so the script still runs on the old output files.
struct Component {
  std::string suffix;
  std::string label;
  int         col;
  bool        optional;
};

const std::vector<Component> kComponents = {
  {"_CCQE",    "CCQE",        cols[1], false},
  {"_CC1p1pi", "CC1p1#pi",    cols[2], false},
  {"_CC2p2h",  "CC2p2h",      cols[3], false},
  {"_CCDIS",   "CCDIS",       cols[4], false},
  {"_CCCoh",   "CC Coh. #pi", cols[5], true },
};

// ===================== input handling =====================

struct Inputs {
  TFile*              cv  = nullptr;
  std::vector<TFile*> throws;
  TH3*                h3D = nullptr;
};

std::set<std::string> gWarned;
void WarnOnce(const std::string& msg) {
  if (gWarned.insert(msg).second) std::cerr << msg << "\n";
}

// NOTE: must match SplitWrite3D's std::to_string(x)/(y) naming exactly,
// which does NOT zero-pad (x0, x1, ..., x9, x10, ...) -- using %02d here
// previously caused every slice with ix<10 or iy<10 to be silently
// skipped as "missing", since the zero-padded name never matched.
//
// Each histogram is read from disk exactly ONCE and cached here. Several PDFs
// use the same slices, and re-reading "name;1" makes ROOT replace (delete)
// the copy it already holds in memory while it is still in use -- that was
// the cause of the segfaults. With TH1::AddDirectory(kFALSE) set in main,
// the histograms belong to this cache, not to the file. Missing histograms
// are cached as nullptr so they are not looked up again.
std::map<std::string, TH1*> gSliceCache;

TH1* GetSlice(TFile* f, const std::string& comp, int ix, int iy) {
  std::string name = kPrefix + comp + "_x" + std::to_string(ix) +
                     "_y" + std::to_string(iy) + ";1";
  std::string key = std::string(f->GetName()) + "|" + name;
  auto it = gSliceCache.find(key);
  if (it != gSliceCache.end()) return it->second;
  TH1* h = f->Get<TH1>(name.c_str());
  gSliceCache[key] = h;
  return h;
}

// A private copy for drawing: detached from the file and owned by the pad
// it is drawn in (deleted when that pad is cleared).
TH1* PadCopy(const TH1* h) {
  TH1* c = (TH1*) h->Clone();
  c->SetDirectory(nullptr);
  c->SetBit(kCanDelete);
  return c;
}

std::string SliceName(int ix, int iy) {
  return "x" + std::to_string(ix) + " y" + std::to_string(iy);
}

std::string RangeText(const char* var, TAxis* ax, int bin0) {
  char buf[128];
  snprintf(buf, sizeof(buf), "%.2f < #it{p}_{%s} < %.2f GeV/c",
           ax->GetBinLowEdge(bin0 + 1), var, ax->GetBinUpEdge(bin0 + 1));
  return buf;
}

std::string WindowText() {
  char buf[128];
  snprintf(buf, sizeof(buf), "%.1f < E_{#nu} < %.1f GeV", kEnuWinLo, kEnuWinHi);
  return buf;
}

// ===================== drawing helpers =====================

// Shared axis formatting applied once to the frame histogram
// xLabels = false hides the x labels and title (grid pads that have another
// pad directly below them share that pad's x axis).
void FormatAxes(TH1* h, const char* xtitle, const char* ytitle,
                bool compact = false, bool xLabels = true) {
  h->GetXaxis()->SetTitle(xLabels ? xtitle : "");
  h->GetYaxis()->SetTitle(ytitle);
  h->GetXaxis()->SetTitleFont(43);
  h->GetXaxis()->SetLabelFont(43);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetLabelFont(43);
  h->GetXaxis()->SetTitleSize(compact ? 18 : 20);
  h->GetXaxis()->SetLabelSize(!xLabels ? 0 : (compact ? 16 : 18));
  h->GetYaxis()->SetTitleSize(compact ? 18 : 20);
  h->GetYaxis()->SetLabelSize(compact ? 16 : 18);
  h->GetXaxis()->SetTitleOffset(compact ? 1.2 : 1.1);
  h->GetYaxis()->SetTitleOffset(1.5);
  if (compact) {
    // Small panels: fewer numbers (0, 0.5, 1, 1.5, 2) so they stay readable
    h->GetXaxis()->SetNdivisions(504);
    h->GetYaxis()->SetNdivisions(505);
  }
}

void StyleTotal(TH1* h) {
  h->SetLineWidth(2);
  h->SetLineStyle(1);
  h->SetLineColor(cols[0]);
  h->SetFillStyle(0);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
  h->SetMarkerColor(cols[0]);
}

void StyleComponent(TH1* h, int col, int lineStyle) {
  h->SetLineWidth(2);
  h->SetLineStyle(lineStyle);   // 1 = full sample, 2 = in Enu window
  h->SetLineColor(col);
  h->SetFillStyle(0);
}

TLine gWinKey;  // dummy object, only used as the legend key for dashed lines

// Builds a filled band running from (hCentral - hSpread) to
// (hCentral + hSpread) at each bin, centered exactly on hCentral's own
// value -- NOT on inherited statistical errors, which is what an older
// two-histogram E2 approach incorrectly did, causing the band to float off
// the CCQE line.
TGraphAsymmErrors* MakePCABand(const TH1* hCentral, const TH1* hSpread, int col) {
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
  g->SetBit(kCanDelete);   // owned by the pad it is drawn in
  return g;
}

// Draw a DUNE / MINERvA style preliminary watermark
void DrawWatermark(const char* expt = "",
                   const char* prelim = "",
                   double x = 0.14, double y = 0.885, double dx = 0.20) {
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(73);   // bold italic Helvetica
  latex.SetTextSize(18);
  latex.SetTextColor(kGray + 1);
  latex.DrawLatex(x, y, expt);

  latex.SetTextFont(43);
  latex.SetTextSize(16);
  latex.SetTextColor(kGray + 1);
  latex.DrawLatex(x + dx, y + 0.002, prelim);
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

// ===================== slice content =====================
// Draws the (ix, iy) slice into the current pad and fills `leg`.
//   sfx = ""         : full sample (optionally with the Enu-window overlay)
//   sfx = kWinSuffix : the Enu-window subsample: its total plus whichever
//                      in-window components exist
//   compact          : small-pad version for the grid pages (no y title,
//                      no dashed per-component window lines)
//   xLabels          : draw x labels/title (false for grid pads with a pad
//                      directly below)
// Returns false, having drawn nothing, if the slice is missing or empty.
bool DrawSliceContent(const Inputs& in, int ix, int iy,
                      const std::string& sfx, bool overlayWindow,
                      bool compact, TLegend* leg, bool xLabels = true) {

  const bool winOnly = !sfx.empty();

  // ---------- fetch (file-owned objects: read only, never drawn) ----------
  // The full sample needs its total. The window subsample draws its own
  // total when the file has one, plus whichever window components exist.
  TH1* fTot = GetSlice(in.cv, sfx, ix, iy);
  if (!winOnly) {
    if (!fTot) {
      WarnOnce("Skipping missing slice: " + SliceName(ix, iy));
      return false;
    }
    if (fTot->Integral() <= 0) {
      WarnOnce("Skipping empty slice: " + SliceName(ix, iy));
      return false;
    }
  } else if (!fTot) {
    WarnOnce("Note: no " + kPrefix + kWinSuffix +
             " (window total) histograms found.");
  }

  std::vector<TH1*>             fComps;
  std::vector<const Component*> compDefs;
  TH1* fCCQE = nullptr;
  for (const auto& comp : kComponents) {
    TH1* h = GetSlice(in.cv, comp.suffix + sfx, ix, iy);
    if (!h) {
      if (comp.optional || winOnly) {
        WarnOnce("Note: no " + kPrefix + comp.suffix + sfx +
                 " histograms found; drawing without that component.");
        continue;
      }
      WarnOnce("Skipping missing slice: " + SliceName(ix, iy) +
               " (no " + comp.suffix + ")");
      return false;
    }
    fComps.push_back(h);
    compDefs.push_back(&comp);
    if (comp.suffix == "_CCQE") fCCQE = h;
  }

  // Window subsample: skip the slice if neither the total nor any component
  // has events in it
  double ymax = fTot ? fTot->GetMaximum() : 0.0;
  if (winOnly) {
    double sum = fTot ? fTot->Integral() : 0.0;
    for (TH1* h : fComps) {
      sum += h->Integral();
      ymax = std::max(ymax, h->GetMaximum());
    }
    if (sum <= 0) {
      WarnOnce("Skipping empty slice: " + SliceName(ix, iy) + " (window sample)");
      return false;
    }
  }

  // PCA throws: required for the full sample, optional for the window sample
  std::vector<TH1*> fThrows;
  if (fCCQE) {
    for (TFile* f : in.throws) {
      TH1* h = GetSlice(f, "_CCQE" + sfx, ix, iy);
      if (!h) {
        if (winOnly) {
          WarnOnce("Note: PCA throw files have no _CCQE" + sfx +
                   " histograms; window plots drawn without the PCA band.");
          fThrows.clear();
          break;
        }
        WarnOnce("Skipping missing slice: " + SliceName(ix, iy) +
                 " (missing PCA throw)");
        return false;
      }
      fThrows.push_back(h);
    }
  }

  // Enu-window overlay on the full sample (optional)
  TH1*              fWinTot = nullptr;
  std::vector<TH1*> fWinComps;
  std::vector<int>  winCols;
  if (kDrawWindow && overlayWindow && !winOnly) {
    fWinTot = GetSlice(in.cv, kWinSuffix, ix, iy);
    if (fWinTot) {
      if (!compact) {
        for (size_t k = 0; k < compDefs.size(); ++k) {
          TH1* h = GetSlice(in.cv, compDefs[k]->suffix + kWinSuffix, ix, iy);
          if (h) { fWinComps.push_back(h); winCols.push_back(compDefs[k]->col); }
        }
      }
    } else {
      WarnOnce("Note: no " + kPrefix + kWinSuffix +
               " histograms found; drawing without the Enu-window overlay.");
    }
  }

  // ---------- PCA band (quadrature sum of throws) ----------
  TGraphAsymmErrors* band = nullptr;
  if (fCCQE && !fThrows.empty()) {
    TH1* hpca = (TH1*) fCCQE->Clone();
    hpca->SetDirectory(nullptr);
    hpca->Reset();
    for (TH1* hb : fThrows) {
      for (int b = 1; b <= hpca->GetNbinsX(); ++b) {
        double d = fCCQE->GetBinContent(b) - hb->GetBinContent(b);
        hpca->AddBinContent(b, d * d);
      }
    }
    for (int b = 1; b <= hpca->GetNbinsX(); ++b)
      hpca->SetBinContent(b, std::sqrt(hpca->GetBinContent(b)));
    band = MakePCABand(fCCQE, hpca, cols[1]);   // band copies the values
    delete hpca;
  }

  // ---------- draw private copies, back to front ----------
  // Empty frame first: carries the axes, so the total can be drawn last
  // (on top) without being added to the pad twice.
  // DrawFrame gives an axis running exactly from kXMin to kXMax, whatever
  // the histogram binning is (SetRangeUser can only snap to bin edges, so
  // it could not force the axis to start at 0). Histograms drawn on top are
  // clipped to this range. The frame is owned by the pad.
  TH1F* frame = gPad->DrawFrame(kXMin, 0.0, kXMax, 1.25 * ymax);
  FormatAxes(frame, kXTitle, compact ? "" : kYTitle, compact, xLabels);

  TH1* win = nullptr;
  if (fWinTot) {
    win = PadCopy(fWinTot);
    win->SetLineColor(kGray);
    win->SetLineStyle(1);
    win->SetLineWidth(1);
    win->SetFillColorAlpha(kGray, 0.6);
    win->SetFillStyle(1001);
    win->Draw("HIST SAME");
  }

  if (band) band->Draw("3 SAME");

  std::vector<TH1*> comps;
  for (size_t k = 0; k < fComps.size(); ++k) {
    TH1* h = PadCopy(fComps[k]);
    StyleComponent(h, compDefs[k]->col, 1);
    h->Draw("HIST SAME");
    comps.push_back(h);
  }

  for (size_t k = 0; k < fWinComps.size(); ++k) {
    TH1* h = PadCopy(fWinComps[k]);
    StyleComponent(h, winCols[k], 2);
    h->Draw("HIST SAME");
  }

  TH1* tot = nullptr;
  if (fTot) {
    tot = PadCopy(fTot);
    StyleTotal(tot);
    tot->Draw("E1 SAME");
  }

  gPad->RedrawAxis();   // keep ticks visible over the filled histogram

  // ---------- legend, in reading order ----------
  if (winOnly && !compact) leg->SetHeader(WindowText().c_str());
  if (tot) leg->AddEntry(tot, "Total (stat. unc.)", "lep");
  if (win)  leg->AddEntry(win, ("Total, " + WindowText()).c_str(), "f");
  if (band) leg->AddEntry(band, "CCQE z-exp PCA (1#sigma)", "f");
  for (size_t k = 0; k < comps.size(); ++k)
    leg->AddEntry(comps[k], compDefs[k]->label.c_str(), "l");
  if (!fWinComps.empty()) {
    gWinKey.SetLineStyle(2);
    gWinKey.SetLineWidth(2);
    gWinKey.SetLineColor(kGray + 2);
    leg->AddEntry(&gWinKey, "Dashed: same, in E_{#nu} window", "l");
  }
  return true;
}

// ===================== single-slice pages =====================

bool DrawSlicePage(TCanvas& c, const Inputs& in, int ix, int iy,
                   const std::string& outpdf,
                   const std::string& sfx, bool overlayWindow) {
  c.Clear();
  c.cd();

  // Legend: top-right; height set from the number of entries
  auto* leg = new TLegend(0.55, 0.60, 0.93, 0.88);
  leg->SetNColumns(1);

  if (!DrawSliceContent(in, ix, iy, sfx, overlayWindow, false, leg)) {
    delete leg;   // never drawn, so not owned by the canvas
    return false;
  }

  leg->SetY1NDC(0.88 - 0.04 * leg->GetNRows());
  leg->SetBit(kCanDelete);
  leg->Draw();

  double pT_lo = in.h3D->GetXaxis()->GetBinLowEdge(ix + 1);
  double pT_hi = in.h3D->GetXaxis()->GetBinUpEdge (ix + 1);
  double pZ_lo = in.h3D->GetYaxis()->GetBinLowEdge(iy + 1);
  double pZ_hi = in.h3D->GetYaxis()->GetBinUpEdge (iy + 1);

  DrawWatermark("", "");
  DrawSliceLabel(pT_lo, pT_hi, pZ_lo, pZ_hi);

  c.Update();
  c.Print(outpdf.c_str());
  return true;
}

// pTOuter = true : fixed pT bin, step through all pZ, then the next pT bin
// pTOuter = false: fixed pZ bin, step through all pT, then the next pZ bin
void MakePDF(TCanvas& c, const Inputs& in, const std::string& outpdf,
             bool pTOuter, const std::string& sfx = "",
             bool overlayWindow = true) {
  const int npt = in.h3D->GetXaxis()->GetNbins();
  const int npz = in.h3D->GetYaxis()->GetNbins();
  const int nOuter = pTOuter ? npt : npz;
  const int nInner = pTOuter ? npz : npt;

  // "[" / "]" open and close the multi-page PDF without printing a page,
  // so there is no blank first page and no repeated last page.
  c.Print((outpdf + "[").c_str());

  int nPages = 0;
  for (int o = 0; o < nOuter; ++o) {
    for (int i = 0; i < nInner; ++i) {
      int ix = pTOuter ? o : i;   // pT bin
      int iy = pTOuter ? i : o;   // pZ bin
      if (DrawSlicePage(c, in, ix, iy, outpdf, sfx, overlayWindow)) ++nPages;
    }
  }

  c.Print((outpdf + "]").c_str());
  c.Clear();   // leave the canvas empty for the next PDF
  std::cout << "Wrote " << nPages << " pages to " << outpdf << "\n";
}

// ===================== grid pages =====================
// pTPerPage = true : pages for each pT bin, its pZ slices in a grid
// pTPerPage = false: pages for each pZ bin, its pT slices in a grid
//
// Landscape pages with kGridCols x kGridRows slices each; a bin with more
// non-empty slices than that continues on the next page ("page 1 of 2").
// Empty slices are left out. Layout, top to bottom: title strip, legend
// strip, grid.
//
// The pads are created ONCE and reused: each page only clears their
// contents and moves them into place. The canvas itself is never cleared,
// because that deletes the sub-pads (the cause of the earlier crash).
const int kGridCols = 5;
const int kGridRows = 3;

// True if the slice has anything to draw (empty slices are skipped)
bool SliceHasData(const Inputs& in, int ix, int iy, const std::string& sfx) {
  TH1* tot = GetSlice(in.cv, sfx, ix, iy);
  if (sfx.empty()) return tot && tot->Integral() > 0;
  double sum = tot ? tot->Integral() : 0.0;
  for (const auto& comp : kComponents) {
    TH1* h = GetSlice(in.cv, comp.suffix + sfx, ix, iy);
    if (h) sum += h->Integral();
  }
  return sum > 0;
}

void MakeGridPDF(const Inputs& in, const std::string& outpdf,
                 bool pTPerPage, const std::string& sfx = "",
                 bool overlayWindow = true) {
  TAxis* axPT = in.h3D->GetXaxis();
  TAxis* axPZ = in.h3D->GetYaxis();
  const int nFixed  = pTPerPage ? axPT->GetNbins() : axPZ->GetNbins();
  const int nSlices = pTPerPage ? axPZ->GetNbins() : axPT->GetNbins();
  const int nc = kGridCols, nr = kGridRows, perPage = nc * nr;

  const double titleB = 0.94;    // bottom of the title strip
  const double gridT  = 0.885;   // bottom of the legend strip = top of grid

  // Pad margins, as fractions of each pad. Only pads with no pad directly
  // below them carry x labels, so the others get a thin bottom margin. Those
  // pads are made taller to compensate, so every plot area is the same size.
  const double mLeft = 0.15, mRight = 0.02, mTop = 0.075;
  const double mBotIn = 0.025, mBotOut = 0.14;
  const double hIn  = 1.0 - mTop - mBotIn;    // plot fraction, inner pads
  const double hOut = 1.0 - mTop - mBotOut;   // plot fraction, x-labelled pads
  const double plotH = gridT / ((nr - 1) / hIn + 1.0 / hOut);

  // Unique canvas name per call: creating a canvas with an existing name
  // makes ROOT delete the old one.
  static int nGrid = 0;
  auto* cg = new TCanvas(Form("cg_%d", nGrid++), "", 2200, 1500);   // landscape
  cg->cd();

  auto* titlePad = new TPad("titlepad", "", 0.0, titleB, 1.0, 1.0);
  titlePad->Draw();
  cg->cd();
  auto* legPad = new TPad("legpad", "", 0.0, gridT, 1.0, titleB);
  legPad->Draw();

  std::vector<TPad*> pads;
  for (int k = 0; k < perPage; ++k) {
    cg->cd();
    auto* pad = new TPad(Form("gridpad_%d", k), "", 0.0, 0.0, 1.0, 1.0);
    pad->SetLeftMargin(mLeft);
    pad->SetRightMargin(mRight);
    pad->SetTopMargin(mTop);
    pad->SetFillStyle(4000);   // transparent, so unused slots hide nothing
    pad->Draw();
    pads.push_back(pad);
  }

  cg->Print((outpdf + "[").c_str());

  int nWritten = 0;
  for (int f = 0; f < nFixed; ++f) {

    // Non-empty slices for this fixed bin
    std::vector<int> slices;
    for (int s = 0; s < nSlices; ++s) {
      const int ix = pTPerPage ? f : s;   // pT bin
      const int iy = pTPerPage ? s : f;   // pZ bin
      if (SliceHasData(in, ix, iy, sfx)) slices.push_back(s);
    }
    if (slices.empty()) continue;
    const int nSub = ((int) slices.size() + perPage - 1) / perPage;

    for (int sub = 0; sub < nSub; ++sub) {
      const int first = sub * perPage;
      const int nOn   = std::min(perPage, (int) slices.size() - first);

      // Legend first: its entries point at objects in the grid pads
      legPad->Clear();
      titlePad->Clear();

      TLegend* shared = nullptr;

      for (int k = 0; k < perPage; ++k) {
        TPad* pad = pads[k];
        pad->Clear();

        const int  col = k % nc, row = k / nc;
        const bool used = (k < nOn);
        const bool xLab = used && (k + nc >= nOn);   // nothing directly below
        const double x1 = (double) col / nc, x2 = (col + 1.0) / nc;
        const double y2 = gridT - row * plotH / hIn;
        const double y1 = y2 - (xLab ? plotH / hOut : plotH / hIn);
        pad->SetPad(x1, std::max(0.0, y1), x2, y2);
        pad->SetBottomMargin(xLab ? mBotOut : mBotIn);
        pad->Modified();
        if (!used) continue;

        pad->cd();
        const int s  = slices[first + k];
        const int ix = pTPerPage ? f : s;   // pT bin
        const int iy = pTPerPage ? s : f;   // pZ bin
        std::string padTxt = pTPerPage ? RangeText("Z", axPZ, iy)
                                       : RangeText("T", axPT, ix);

        auto* leg = new TLegend(0.0, 0.0, 1.0, 1.0);
        if (DrawSliceContent(in, ix, iy, sfx, overlayWindow, true, leg, xLab)) {
          if (!shared) shared = leg; else delete leg;   // spare legends never drawn
        } else {
          delete leg;
          TLatex miss;
          miss.SetNDC();
          miss.SetTextFont(43);
          miss.SetTextSize(16);
          miss.SetTextAlign(22);
          miss.SetTextColor(kGray + 1);
          miss.DrawLatex(0.5, 0.5, "histograms missing");
        }

        // Slice label in the top margin, right-aligned
        TLatex lab;
        lab.SetNDC();
        lab.SetTextFont(43);
        lab.SetTextSize(17);
        lab.SetTextAlign(31);
        lab.DrawLatex(1.0 - mRight, 1.0 - mTop + 0.015, padTxt.c_str());
      }

      if (!shared) continue;   // nothing could be drawn: no page

      // Title strip: the fixed bin (+ window, + page count) and watermark,
      // with the axis quantities on a second line
      titlePad->cd();
      TLatex title;
      title.SetNDC();
      title.SetTextFont(43);
      title.SetTextSize(30);
      std::string ttl = pTPerPage ? RangeText("T", axPT, f)
                                  : RangeText("Z", axPZ, f);
      if (!sfx.empty()) ttl += ",   " + WindowText();
      if (nSub > 1) ttl += Form("   (page %d of %d)", sub + 1, nSub);
      title.DrawLatex(0.015, 0.50, ttl.c_str());

      TLatex axisNote;
      axisNote.SetNDC();
      axisNote.SetTextFont(43);
      axisNote.SetTextSize(17);
      axisNote.SetTextColor(kGray + 2);
      axisNote.DrawLatex(0.015, 0.10,
          "Each panel: d#sigma / dE^{Had.}_{Avail.} [cm^{2}/GeV/nucleon]"
          "  vs  E^{Had.}_{Avail.} [GeV]");
      DrawWatermark("DUNE Simulation", "Preliminary", 0.80, 0.50, 0.075);
      titlePad->Modified();

      // Legend strip: one or two rows across the page
      legPad->cd();
      const int nEnt = shared->GetListOfPrimitives()->GetSize();
      shared->SetNColumns(nEnt <= 4 ? std::max(1, nEnt) : (nEnt + 1) / 2);
      shared->SetTextSize(19);
      shared->SetMargin(0.22);
      shared->SetX1NDC(0.01);
      shared->SetX2NDC(0.99);
      shared->SetY1NDC(0.05);
      shared->SetY2NDC(0.95);
      shared->SetBit(kCanDelete);
      shared->Draw();
      legPad->Modified();

      cg->cd();
      cg->Modified();
      cg->Update();
      cg->Print(outpdf.c_str());
      ++nWritten;
    }
  }

  cg->Print((outpdf + "]").c_str());
  std::cout << "Wrote " << nWritten << " pages to " << outpdf << "\n";
  // cg is left for ROOT to clean up at exit, rather than deleted here,
  // so no pad is ever destroyed while its contents are still referenced.
}

// ===================== diagnostics =====================
// Prints the names of the Enu-window histograms the file holds for the first
// slice, so you can check they match what this script looks for
// (<prefix><component>_OscWin_x0_y0).
void ListWindowHists(TFile* f) {
  std::cout << "Enu-window histograms in " << f->GetName()
            << " for slice x0_y0:\n";
  int n = 0;
  TIter next(f->GetListOfKeys());
  while (TKey* key = (TKey*) next()) {
    std::string nm = key->GetName();
    const std::string tail = "_x0_y0";
    if (nm.find(kWinSuffix) != std::string::npos &&
        nm.size() >= tail.size() &&
        nm.compare(nm.size() - tail.size(), tail.size(), tail) == 0) {
      std::cout << "  " << nm << "\n";
      ++n;
    }
  }
  if (n == 0)
    std::cout << "  (none -- the window PDFs will be empty)\n";
}

// ===================== main =====================
// NOTE: function name matches the file name (plot_errorbands_zexp_slices.cxx)
// so that `root -l -b -q plot_errorbands_zexp_slices.cxx` actually invokes it.
void plot_errorbands_zexp_slices() {

  SetPublicationStyle();

  // Histograms read from the files (and any clones) are NOT attached to a
  // directory, so ROOT never deletes them behind the script's back.
  // Must be set before any file is opened or histogram read.
  TH1::AddDirectory(kFALSE);

  // ---------- files ----------
  Inputs in;
  in.cv = TFile::Open((kDir + "DUNE_3D_MINERvALike.out_throw_nominal_new.root").c_str());
  for (const char* b : {"b1", "b2", "b3", "b4"}) {
    in.throws.push_back(
      TFile::Open((kDir + "DUNE_3D_MINERvALike.out_throw_" + b + "_new.root").c_str()));
  }
  if (!in.cv) { std::cerr << "Could not open nominal file\n"; return; }
  for (TFile* f : in.throws)
    if (!f) { std::cerr << "Could not open a PCA throw file\n"; return; }

  in.h3D = in.cv->Get<TH3>(kPrefix.c_str());
  if (!in.h3D) { std::cerr << "Could not find " << kPrefix << "\n"; return; }

  ListWindowHists(in.cv);

  // Use a square-ish canvas; 700x560 prints cleanly at journal width
  auto* c1 = new TCanvas("c1", "", 700, 560);

  // Big PDFs, one slice per page, both orderings
  // Full sample WITHOUT the Enu-window overlay
  MakePDF(*c1, in, kOutStem + "_fixedPZ_scanPT.pdf", false, "", false);
  MakePDF(*c1, in, kOutStem + "_fixedPT_scanPZ.pdf", true,  "", false);

  // Overview pages: a whole row of slices on one page, both ways round
  MakeGridPDF(in, kOutStem + "_grid_onePTperPage.pdf", true,  "", false);
  MakeGridPDF(in, kOutStem + "_grid_onePZperPage.pdf", false, "", false);

  // The same full-sample PDFs WITH the Enu-window overlay
  MakePDF(*c1, in, kOutStem + "_withOscWin_fixedPZ_scanPT.pdf", false, "", true);
  MakePDF(*c1, in, kOutStem + "_withOscWin_fixedPT_scanPZ.pdf", true,  "", true);
  MakeGridPDF(in, kOutStem + "_withOscWin_grid_onePTperPage.pdf", true,  "", true);
  MakeGridPDF(in, kOutStem + "_withOscWin_grid_onePZperPage.pdf", false, "", true);

  // Enu-window subsample on its own (first osc. max to the minimum below it)
  MakePDF(*c1, in, kOutStem + "_oscwindow_fixedPT_scanPZ.pdf", true, kWinSuffix);
  MakeGridPDF(in, kOutStem + "_oscwindow_grid_onePTperPage.pdf", true,  kWinSuffix);
  MakeGridPDF(in, kOutStem + "_oscwindow_grid_onePZperPage.pdf", false, kWinSuffix);
}