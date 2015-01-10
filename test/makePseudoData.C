#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TString.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TRandom3.h"


const double alpha = 1 - 0.6827;

void makePseudoData()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.07);

  TRandom3 rand(37);

  // input file
  TFile *file = new TFile("test/histos_FullStats_JEC.root");

  // get event counts
  TH1F *eventCount = (TH1F*)file->Get("dijetAna/hEventCount");
  // get dijet mass distribution
  TH1F *dijetMass = (TH1F*)file->Get("dijetAna/hCorrDijetMass");
  
  float xsMC = 2.435e+09; // total xs in pb
  float intLumi = 100.;   // integrated lumi in /pb
  
  // scale simulated dijet mass distibution to 100 /pb
  dijetMass->Scale((xsMC*intLumi)/eventCount->GetBinContent(2));

  // differential counts
  TH1D *diffCounts = (TH1D*)dijetMass->Clone("diffCounts");
  diffCounts->Reset();
  
  // pseudo-data
  TH1D *hDijetMassPseudo = (TH1D*)dijetMass->Clone("hDijetMassPseudo");
  hDijetMassPseudo->Reset();

  for(int i=1; i<=dijetMass->GetNbinsX(); ++i)
  {
    double width = dijetMass->GetBinWidth(i);
    double n = rand.Poisson(dijetMass->GetBinContent(i));
    double err = TMath::Sqrt(n);
    //double l = 0.5*TMath::ChisquareQuantile(alpha/2,2*n);
    //double h = 0.5*TMath::ChisquareQuantile(1-alpha/2,2*(n+1));
    //double err = (h-l)/2;
    
    diffCounts->SetBinContent(i, n/width);
    diffCounts->SetBinError(i, err/width);
    
    hDijetMassPseudo->SetBinContent(i, n);
    hDijetMassPseudo->SetBinError(i, err);
  }

  diffCounts->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  diffCounts->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");
  diffCounts->SetMarkerStyle(20);
  diffCounts->SetMarkerSize(0.8);
  diffCounts->SetTitleOffset(1.4,"Y");

  // declare a fit function and fit
  //TF1 *fit = new TF1("fit", "[0]*pow(1-x/13000.0,[1])/pow(x/13000.,[2]+[3]*log(x/13000.))",1118,6099);
  //fit->SetParameter(0,2.94554e-04);
  //fit->SetParameter(1,5.57678e+00);
  //fit->SetParameter(2,7.09806e+00);
  //fit->SetParameter(3,3.48634e-01);

  //gStyle->SetOptFit(1111);

  //fit->SetLineWidth(2);
  //fit->SetLineColor(kRed);
  //std::cout << "*********************************************************" << std::endl;
  //TFitResultPtr s = diffCounts->Fit("fit","SRLI");
  //TString status_default = gMinuit->fCstatu.Data();
  // Results of the fit
  //std::cout << "*********************************************************" << std::endl;
  //Double_t chi_fit = fit->GetChisquare();
  //Double_t ndf_fit = fit->GetNDF();
  //std::cout << "Chi2/ndf: " << chi_fit << "/" << ndf_fit << " = " << chi_fit/ndf_fit << std::endl;
  //std::cout << "Status: "<<status_default<<std::endl;
  //std::cout << "*********************************************************" << std::endl;

  // Print fit results
  //s->Print("V");
  
  //TCanvas *c = new TCanvas("c", "",1000,800);
  //c->cd();

  //diffCounts->Draw();

  //TLegend *legend = new TLegend(.7,.5,.85,.6);
  //legend->SetBorderSize(0);
  //legend->SetFillColor(0);
  //legend->SetFillStyle(0);
  //legend->AddEntry(diffCounts, "QCD MC","lp");
  //legend->AddEntry(fit, "Fit","l");
  //legend->Draw();
  
  //c->SetLogy();
  //c->SaveAs("diffCounts.eps");

  // save pseudo-data
  hDijetMassPseudo->SaveAs("hDijetMassPseudo_JEC.root");

}
