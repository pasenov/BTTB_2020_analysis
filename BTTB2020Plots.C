#include <memory>
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <TMath.h>

#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH1F.h"

//#include "tdrstyle.C"
//#include "CMS_lumi.C"
#include "TH1.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TPad.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include "TPaveStats.h"

#define extendedRootTree true

void BTTB2020Plots()
{
  //double sc = 26814/17539;
  float sc = 1.529;

  TString cmsText     = "CMS";
  float cmsTextFont   = 61;  // default is helvetic-bold

  bool writeExtraText = false;
  TString extraText   = "Preliminary";
  TString extraText1   = "Data";
  TString extraText2   = "Simulation";
  float extraTextFont = 52;
  // default is helvetica-italics

  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  float lumiTextSize     = 0.6;  
  float extraTextSize    = 0.6;
  float lumiTextOffset   = 0.2;
  float cmsTextSize      = 0.75;
  float cmsTextOffset    = 0.1;  // only used in outOfFrame version

  float relPosX    = 0.045;
  float relPosY    = 0.035;
  float relExtraDY = 1.2;

  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;

  TString lumi_13TeV = "20.1 fb^{-1}";
  TString lumi_8TeV  = "19.7 fb^{-1}";
  TString lumi_7TeV  = "5.1 fb^{-1}";
  TString lumi_sqrtS = "";

  int alignX_=1;
  int alignY_=1;

  bool drawLogo      = false;

  int iPeriod=3;
  int iPosX=10;
  //setTDRStyle();

  TLatex latex;
				
  int n_ = 2;

  float x1_l = 0.92;
  float y1_l = 0.60;

  float dx_l = 0.30;
  float dy_l = 0.18;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;

  float ar_l = dy_l/dx_l;
		
  float x_l[1];
  float ex_l[1];
  float y_l[1];
  float ey_l[1];
		
  //    float gap_ = 0.09/ar_l;
  float gap_ = 1./(n_+1);
		
  float bwx_ = 0.12;
  float bwy_ = gap_/1.5;
		
  x_l[0] = 1.2*bwx_;
  //    y_l[0] = 1-(1-0.10)/ar_l;
  y_l[0] = 1-gap_;
  ex_l[0] = 0;
  ey_l[0] = 0.04/ar_l;

  int W = 800;
  int H = 600;
  //
  int H_ref = 600; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  float l = 0.1;
  float t = 0.05;
  float r = 0.4;

  TFile *f = new TFile("proton_GOOD.root");
  //TFile *f1 = new TFile("proton1.root");
  gROOT->Reset();

  gStyle->Reset();
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  gStyle->SetEndErrorSize(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetOptDate(0);
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.15);

  /*gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.02); */

  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.25);
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);

  /* // Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);

  gStyle->SetHatchesLineWidth(5);
  gStyle->SetHatchesSpacing(0.05);

  gStyle->cd(); */

  //All histograms
  
  //Cluster occupancy
  TH2F *hClOcPix2Top = (TH2F*)f->Get("14");
  TH2F *hClOcPix2Bot = (TH2F*)f->Get("13");

  TCanvas *cClOcPix2Top = new TCanvas("cClOcPix2Top", "Cluster occupancy for top module of pixel layer (simulation)", W, H);
  hClOcPix2Top->GetXaxis()->SetLabelFont(42);
  hClOcPix2Top->GetXaxis()->SetLabelSize(0.035);
  hClOcPix2Top->GetXaxis()->SetTitleSize(0.035);
  hClOcPix2Top->GetYaxis()->SetLabelFont(42);
  hClOcPix2Top->GetYaxis()->SetLabelSize(0.035);
  hClOcPix2Top->GetYaxis()->SetTitleSize(0.035);
  hClOcPix2Top->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cClOcPix2Top->SetFillColor(0);
  cClOcPix2Top->SetBorderMode(0);
  cClOcPix2Top->SetFrameFillStyle(0);
  cClOcPix2Top->SetFrameBorderMode(0);
  cClOcPix2Top->SetLeftMargin( L/W );
  cClOcPix2Top->SetRightMargin( R/W );
  cClOcPix2Top->SetTopMargin( T/H );
  cClOcPix2Top->SetBottomMargin( B/H );
  cClOcPix2Top->SetTickx(0);
  cClOcPix2Top->SetTicky(0);
  hClOcPix2Top->GetXaxis()->SetTitle("#font[42]{col}");
  hClOcPix2Top->GetYaxis()->SetTitle("#font[42]{row}");
  hClOcPix2Top->Draw();
  float cClOcPix2TopX = 5.;   //l +  relPosX*(1-l-r);
  float cClOcPix2TopY = 155.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cClOcPix2TopX,cClOcPix2TopY,extraText2);
  cClOcPix2Top->Print(Form("Cluster_occupancy_top_module_pixel_layer_2_simulation.eps"),"eps");
  cClOcPix2Top->Print(Form("Cluster_occupancy_top_module_pixel_layer_2_simulation.jpg"),"jpg");

  TCanvas *cClOcPix2Bot = new TCanvas("cClOcPix2Bot", "Cluster occupancy for bottom module of pixel layer (simulation)", W, H);
  hClOcPix2Bot->GetXaxis()->SetLabelFont(42);
  hClOcPix2Bot->GetXaxis()->SetLabelSize(0.035);
  hClOcPix2Bot->GetXaxis()->SetTitleSize(0.035);
  hClOcPix2Bot->GetYaxis()->SetLabelFont(42);
  hClOcPix2Bot->GetYaxis()->SetLabelSize(0.035);
  hClOcPix2Bot->GetYaxis()->SetTitleSize(0.035);
  hClOcPix2Bot->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cClOcPix2Top->SetFillColor(0);
  cClOcPix2Top->SetBorderMode(0);
  cClOcPix2Top->SetFrameFillStyle(0);
  cClOcPix2Top->SetFrameBorderMode(0);
  cClOcPix2Top->SetLeftMargin( L/W );
  cClOcPix2Top->SetRightMargin( R/W );
  cClOcPix2Top->SetTopMargin( T/H );
  cClOcPix2Top->SetBottomMargin( B/H );
  cClOcPix2Top->SetTickx(0);
  cClOcPix2Top->SetTicky(0);
  hClOcPix2Bot->GetXaxis()->SetTitle("#font[42]{col}");
  hClOcPix2Bot->GetYaxis()->SetTitle("#font[42]{row}");
  hClOcPix2Bot->Draw();
  float cClOcPix2BotX = 5.;   //l +  relPosX*(1-l-r);
  float cClOcPix2BotY = 155.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cClOcPix2BotX,cClOcPix2BotY,extraText2);
  cClOcPix2Bot->Print(Form("Cluster_occupancy_bottom_module_pixel_layer_2_simulation.eps"),"eps");
  cClOcPix2Bot->Print(Form("Cluster_occupancy_bottom_module_pixel_layer_2_simulation.jpg"),"jpg");

  //Cluster size
  TH1F *hClSize2SS1 = (TH1F*)f->Get("241");
  TH1F *hClSize2SS2 = (TH1F*)f->Get("242");
  TH1F *hClSizePix2 = (TH1F*)f->Get("240");

  //Time difference between beam arrival at scintillators
  TH1F *hT1 = (TH1F*)f->Get("235");
  TH1F *hT2 = (TH1F*)f->Get("236");
  TH1F *hDt = (TH1F*)f->Get("237");

  TCanvas *cTime = new TCanvas("cTime", "Time of beam arrival scintillators (simulation)", W, H);  
  float cTimeX = 0.3;   //l +  relPosX*(1-l-r);
  float cTimeY = 20000.; //1-t+lumiTextOffset*t;
  cTime->SetFillColor(0);
  cTime->SetBorderMode(0);
  cTime->SetFrameFillStyle(0);
  cTime->SetFrameBorderMode(0);
  cTime->SetLeftMargin( L/W );
  cTime->SetRightMargin( R/W );
  cTime->SetTopMargin( T/H );
  cTime->SetBottomMargin( B/H );
  cTime->SetTickx(0);
  cTime->SetTicky(0);
  hT1->SetLineColor(kBlue);
  hT1->GetXaxis()->SetLabelFont(42);
  hT1->GetXaxis()->SetLabelSize(0.035);
  hT1->GetXaxis()->SetTitleSize(0.035);
  hT1->GetYaxis()->SetLabelFont(42);
  hT1->GetYaxis()->SetLabelSize(0.035);
  hT1->GetYaxis()->SetTitleSize(0.035);
  hT1->GetYaxis()->SetTitleOffset(1.4);
  hT1->GetXaxis()->SetRangeUser(0,3);
  //hT1->GetYaxis()->SetRangeUser(0,1800);
  hT1->GetXaxis()->SetTitle("#font[42]{Time (ns)}");
  hT1->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  hT1->Draw("hist");
  hT2->SetLineColor(kGreen);
  hT2->GetXaxis()->SetLabelFont(42);
  hT2->GetXaxis()->SetLabelSize(0.035);
  hT2->GetXaxis()->SetTitleSize(0.035);
  hT2->GetYaxis()->SetLabelFont(42);
  hT2->GetYaxis()->SetLabelSize(0.035);
  hT2->GetYaxis()->SetTitleSize(0.035);
  hT2->GetYaxis()->SetTitleOffset(1.4);
  hT2->GetXaxis()->SetRangeUser(0,3);
  //hT2->GetYaxis()->SetRangeUser(0,1800);
  hT2->GetXaxis()->SetTitle("#font[42]{Time (ns)}");
  hT2->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  hT2->Draw("hist same");
  hDt->SetLineColor(kRed);
  hDt->GetXaxis()->SetLabelFont(42);
  hDt->GetXaxis()->SetLabelSize(0.035);
  hDt->GetXaxis()->SetTitleSize(0.035);
  hDt->GetYaxis()->SetLabelFont(42);
  hDt->GetYaxis()->SetLabelSize(0.035);
  hDt->GetYaxis()->SetTitleSize(0.035);
  hDt->GetYaxis()->SetTitleOffset(1.4);
  hDt->GetXaxis()->SetRangeUser(0,3);
  //hDt->GetYaxis()->SetRangeUser(0,1800);
  hDt->GetXaxis()->SetTitle("#font[42]{Time (ns)}");
  hDt->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  hDt->Draw("hist same");
  TLegend *legendTime = new TLegend(x0_l-0.12, y0_l+0.35, x1_l-0.32, y1_l);
  //legendTime->SetHeader("Legend"); 
  legendTime->SetTextFont(42);
  legendTime->SetTextAngle(0);
  legendTime->SetTextColor(kBlack);
  legendTime->SetTextSize(0.035);
  legendTime->SetTextAlign(12); 
  legendTime->AddEntry(hT1, "Time of arrival at first scintillator (t_{1})", "l");
  legendTime->AddEntry(hT2, "Time of arrival at second scintillator (t_{2})", "l");
  legendTime->AddEntry(hDt, "t_{2} - t_{1}", "l");
  legendTime->SetBorderSize(0);
  legendTime->Draw();
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cTimeX,cTimeY,extraText2);
  gStyle->SetLabelFont(42, "XYZ");
  cTime->Update();
  cTime->Write();
  cTime->Print(Form("Time_diff.eps"),"eps");
  cTime->Print(Form("Time_diff.jpg"),"jpg");

  //Kinetic energies at entrances and exits of volumes
  TH1F *hEkinAlex = (TH1F*)f->Get("220");
  TH1F *hEkin2SS1en = (TH1F*)f->Get("225");
  TH1F *hEkin2SS1ex = (TH1F*)f->Get("226");
  TH1F *hEkin2SS2en = (TH1F*)f->Get("227");
  TH1F *hEkin2SS2ex = (TH1F*)f->Get("228");
  TH1F *hEkinScint1en = (TH1F*)f->Get("221");
  TH1F *hEkinScint1ex = (TH1F*)f->Get("222");
  TH1F *hEkinScint2en = (TH1F*)f->Get("231");
  TH1F *hEkinScint2ex = (TH1F*)f->Get("232");
  TH1F *hEkinPixel2en = (TH1F*)f->Get("229");
  TH1F *hEkinPixel2ex = (TH1F*)f->Get("230");


  //Kinetic energy vs. energy loss in Scintillator 2
  TH2F *hEkinEloss2 = (TH2F*)f->Get("9");

  TCanvas *cEkinEloss2 = new TCanvas("cEkinEloss2", "Kinetic energy vs. energy deposition of primary protons in Scintillator 2 (simulation)", W, H);
  hEkinEloss2->GetXaxis()->SetLabelFont(42);
  hEkinEloss2->GetXaxis()->SetLabelSize(0.035);
  hEkinEloss2->GetXaxis()->SetTitleSize(0.035);
  hEkinEloss2->GetYaxis()->SetLabelFont(42);
  hEkinEloss2->GetYaxis()->SetLabelSize(0.035);
  hEkinEloss2->GetYaxis()->SetTitleSize(0.035);
  hEkinEloss2->GetYaxis()->SetTitleOffset(1.4);
  cEkinEloss2->SetFillColor(0);
  cEkinEloss2->SetBorderMode(0);
  cEkinEloss2->SetFrameFillStyle(0);
  cEkinEloss2->SetFrameBorderMode(0);
  cEkinEloss2->SetLeftMargin( L/W );
  cEkinEloss2->SetRightMargin( R/W );
  cEkinEloss2->SetTopMargin( T/H );
  cEkinEloss2->SetBottomMargin( B/H );
  cEkinEloss2->SetTickx(0);
  cEkinEloss2->SetTicky(0);
  hEkinEloss2->GetXaxis()->SetTitle("#font[42]{Energy deposition (MeV)}");
  hEkinEloss2->GetYaxis()->SetTitle("#font[42]{Kinetic energy (MeV)}");
  hEkinEloss2->GetYaxis()->SetRangeUser(0,25.2);
  hEkinEloss2->Draw();
  float cEkinEloss2X = 0.5;   //l +  relPosX*(1-l-r);
  float cEkinEloss2Y = 24.5; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cEkinEloss2X,cEkinEloss2Y,extraText2);
  cEkinEloss2->Print(Form("Energy_dep_ekin_Scint_2.eps"),"eps");
  cEkinEloss2->Print(Form("Energy_dep_ekin_Scint_2.jpg"),"jpg");

  //Kinetic energy vs. z
  //TH2F *hEkinZ1 = (TH2F*)f1->Get("10");
  TH2F *hEkinZ  = (TH2F*)f->Get("10");

  TCanvas *cEkinZ = new TCanvas("cEkinZ", "Kinetic energy vs. z (simulation)", W, H);
  hEkinZ->GetXaxis()->SetLabelFont(42);
  hEkinZ->GetXaxis()->SetLabelSize(0.035);
  hEkinZ->GetXaxis()->SetTitleSize(0.035);
  hEkinZ->GetYaxis()->SetLabelFont(42);
  hEkinZ->GetYaxis()->SetLabelSize(0.035);
  hEkinZ->GetYaxis()->SetTitleSize(0.035);
  hEkinZ->GetYaxis()->SetTitleOffset(1.4);
  cEkinZ->SetFillColor(0);
  cEkinZ->SetBorderMode(0);
  cEkinZ->SetFrameFillStyle(0);
  cEkinZ->SetFrameBorderMode(0);
  cEkinZ->SetLeftMargin( L/W );
  cEkinZ->SetRightMargin( R/W );
  cEkinZ->SetTopMargin( T/H );
  cEkinZ->SetBottomMargin( B/H );
  cEkinZ->SetTickx(0);
  cEkinZ->SetTicky(0);
  hEkinZ->GetXaxis()->SetTitle("#font[42]{z (cm)}");
  hEkinZ->GetYaxis()->SetTitle("#font[42]{Kinetic energy (MeV)}");
  hEkinZ->GetYaxis()->SetRangeUser(0,28.0);
  hEkinZ->Draw();
  float cEkinZX = 7.5;   //l +  relPosX*(1-l-r);
  float cEkinZY = 25.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cEkinZX,cEkinZY,extraText2);
  cEkinZ->Print(Form("Ekin_Z.eps"),"eps");
  cEkinZ->Print(Form("Ekin_Z.jpg"),"jpg");

  //Energy depositions in volumes
  TH1F *hEnDepAir = (TH1F*)f->Get("214");
  TH1F *hEnDepAl = (TH1F*)f->Get("215");
  TH1F *hEnDep2SS1 = (TH1F*)f->Get("216");
  TH1F *hEnDep2SS2 = (TH1F*)f->Get("217");
  TH1F *hEnDepPixSen = (TH1F*)f->Get("219");
  TH1F *hEnDepPixROC = (TH1F*)f->Get("234");
  TH1F *hEnDepScint1 = (TH1F*)f->Get("212");
  TH1F *hEnDepScint2 = (TH1F*)f->Get("213");

  TCanvas *cEnDepAl = new TCanvas("cEnDepAl", "Energy deposition in Al plate (simulation)", W, H);
  hEnDepAl->GetXaxis()->SetLabelFont(42);
  hEnDepAl->GetXaxis()->SetLabelSize(0.035);
  hEnDepAl->GetXaxis()->SetTitleSize(0.035);
  hEnDepAl->GetYaxis()->SetLabelFont(42);
  hEnDepAl->GetYaxis()->SetLabelSize(0.035);
  hEnDepAl->GetYaxis()->SetTitleSize(0.035);
  hEnDepAl->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cEnDepAl->SetFillColor(0);
  cEnDepAl->SetBorderMode(0);
  cEnDepAl->SetFrameFillStyle(0);
  cEnDepAl->SetFrameBorderMode(0);
  cEnDepAl->SetLeftMargin( L/W );
  cEnDepAl->SetRightMargin( R/W );
  cEnDepAl->SetTopMargin( T/H );
  cEnDepAl->SetBottomMargin( B/H );
  cEnDepAl->SetTickx(0);
  cEnDepAl->SetTicky(0);
  gPad->SetLogy();
  hEnDepAl->GetXaxis()->SetTitle("#font[42]{Energy (MeV)}");
  hEnDepAl->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  //hEnDepAl->GetYaxis()->SetRangeUser(0,28.0);
  hEnDepAl->Draw("HIST");
  float cEnDepAlX = 0.5;   //l +  relPosX*(1-l-r);
  float cEnDepAlY = 10000.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cEnDepAlX,cEnDepAlY,extraText2);
  cEnDepAl->Print(Form("EnDep_Al.eps"),"eps");
  cEnDepAl->Print(Form("EnDep_Al.jpg"),"jpg");

  TCanvas *cEnDepScint1 = new TCanvas("cEnDepScint1", "Energy deposition in first scintillator (simulation)", W, H);
  hEnDepScint1->GetXaxis()->SetLabelFont(42);
  hEnDepScint1->GetXaxis()->SetLabelSize(0.035);
  hEnDepScint1->GetXaxis()->SetTitleSize(0.035);
  hEnDepScint1->GetYaxis()->SetLabelFont(42);
  hEnDepScint1->GetYaxis()->SetLabelSize(0.035);
  hEnDepScint1->GetYaxis()->SetTitleSize(0.035);
  hEnDepScint1->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cEnDepScint1->SetFillColor(0);
  cEnDepScint1->SetBorderMode(0);
  cEnDepScint1->SetFrameFillStyle(0);
  cEnDepScint1->SetFrameBorderMode(0);
  cEnDepScint1->SetLeftMargin( L/W );
  cEnDepScint1->SetRightMargin( R/W );
  cEnDepScint1->SetTopMargin( T/H );
  cEnDepScint1->SetBottomMargin( B/H );
  cEnDepScint1->SetTickx(0);
  cEnDepScint1->SetTicky(0);
  gPad->SetLogy();
  hEnDepScint1->GetXaxis()->SetTitle("#font[42]{Energy (MeV)}");
  hEnDepScint1->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  //hEnDepAl->GetYaxis()->SetRangeUser(0,28.0);
  hEnDepScint1->Draw("HIST");
  float cEnDepScint1X = 0.5;   //l +  relPosX*(1-l-r);
  float cEnDepScint1Y = 10000.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cEnDepScint1X,cEnDepScint1Y,extraText2);
  cEnDepScint1->Print(Form("EnDep_Scint1.eps"),"eps");
  cEnDepScint1->Print(Form("EnDep_Scint1.jpg"),"jpg");

  TCanvas *cEnDepScint2 = new TCanvas("cEnDepScint2", "Energy deposition in second scintillator (simulation)", W, H);
  hEnDepScint2->GetXaxis()->SetLabelFont(42);
  hEnDepScint2->GetXaxis()->SetLabelSize(0.035);
  hEnDepScint2->GetXaxis()->SetTitleSize(0.035);
  hEnDepScint2->GetYaxis()->SetLabelFont(42);
  hEnDepScint2->GetYaxis()->SetLabelSize(0.035);
  hEnDepScint2->GetYaxis()->SetTitleSize(0.035);
  hEnDepScint2->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cEnDepScint2->SetFillColor(0);
  cEnDepScint2->SetBorderMode(0);
  cEnDepScint2->SetFrameFillStyle(0);
  cEnDepScint2->SetFrameBorderMode(0);
  cEnDepScint2->SetLeftMargin( L/W );
  cEnDepScint2->SetRightMargin( R/W );
  cEnDepScint2->SetTopMargin( T/H );
  cEnDepScint2->SetBottomMargin( B/H );
  cEnDepScint2->SetTickx(0);
  cEnDepScint2->SetTicky(0);
  gPad->SetLogy();
  hEnDepScint2->GetXaxis()->SetTitle("#font[42]{Energy (MeV)}");
  hEnDepScint2->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  //hEnDepAl->GetYaxis()->SetRangeUser(0,28.0);
  hEnDepScint2->Draw("HIST");
  float cEnDepScint2X = 0.5;   //l +  relPosX*(1-l-r);
  float cEnDepScint2Y = 1100.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cEnDepScint2X,cEnDepScint2Y,extraText2);
  cEnDepScint2->Print(Form("EnDep_Scint2.eps"),"eps");
  cEnDepScint2->Print(Form("EnDep_Scint2.jpg"),"jpg");

  TCanvas *cEnDep2SS1 = new TCanvas("cEnDep2SS1", "Energy deposition in first 2S sensor (simulation)", W, H);
  hEnDep2SS1->GetXaxis()->SetLabelFont(42);
  hEnDep2SS1->GetXaxis()->SetLabelSize(0.035);
  hEnDep2SS1->GetXaxis()->SetTitleSize(0.035);
  hEnDep2SS1->GetYaxis()->SetLabelFont(42);
  hEnDep2SS1->GetYaxis()->SetLabelSize(0.035);
  hEnDep2SS1->GetYaxis()->SetTitleSize(0.035);
  hEnDep2SS1->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cEnDep2SS1->SetFillColor(0);
  cEnDep2SS1->SetBorderMode(0);
  cEnDep2SS1->SetFrameFillStyle(0);
  cEnDep2SS1->SetFrameBorderMode(0);
  cEnDep2SS1->SetLeftMargin( L/W );
  cEnDep2SS1->SetRightMargin( R/W );
  cEnDep2SS1->SetTopMargin( T/H );
  cEnDep2SS1->SetBottomMargin( B/H );
  cEnDep2SS1->SetTickx(0);
  cEnDep2SS1->SetTicky(0);
  gPad->SetLogy();
  hEnDep2SS1->GetXaxis()->SetTitle("#font[42]{Energy (MeV)}");
  hEnDep2SS1->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  //hEnDepAl->GetYaxis()->SetRangeUser(0,28.0);
  hEnDep2SS1->Draw("HIST");
  float cEnDep2SS1X = 0.5;   //l +  relPosX*(1-l-r);
  float cEnDep2SS1Y = 10000.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cEnDep2SS1X,cEnDep2SS1Y,extraText2);
  cEnDep2SS1->Print(Form("EnDep_2SS1.eps"),"eps");
  cEnDep2SS1->Print(Form("EnDep_2SS1.jpg"),"jpg");

  TCanvas *cEnDep2SS2 = new TCanvas("cEnDep2SS2", "Energy deposition in second 2S sensor (simulation)", W, H);
  hEnDep2SS2->GetXaxis()->SetLabelFont(42);
  hEnDep2SS2->GetXaxis()->SetLabelSize(0.035);
  hEnDep2SS2->GetXaxis()->SetTitleSize(0.035);
  hEnDep2SS2->GetYaxis()->SetLabelFont(42);
  hEnDep2SS2->GetYaxis()->SetLabelSize(0.035);
  hEnDep2SS2->GetYaxis()->SetTitleSize(0.035);
  hEnDep2SS2->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cEnDep2SS2->SetFillColor(0);
  cEnDep2SS2->SetBorderMode(0);
  cEnDep2SS2->SetFrameFillStyle(0);
  cEnDep2SS2->SetFrameBorderMode(0);
  cEnDep2SS2->SetLeftMargin( L/W );
  cEnDep2SS2->SetRightMargin( R/W );
  cEnDep2SS2->SetTopMargin( T/H );
  cEnDep2SS2->SetBottomMargin( B/H );
  cEnDep2SS2->SetTickx(0);
  cEnDep2SS2->SetTicky(0);
  gPad->SetLogy();
  hEnDep2SS2->GetXaxis()->SetTitle("#font[42]{Energy (MeV)}");
  hEnDep2SS2->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  //hEnDepAl->GetYaxis()->SetRangeUser(0,28.0);
  hEnDep2SS2->Draw("HIST");
  float cEnDep2SS2X = 0.5;   //l +  relPosX*(1-l-r);
  float cEnDep2SS2Y = 30000.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cEnDepScint2X,cEnDepScint2Y,extraText2);
  cEnDep2SS2->Print(Form("EnDep_2SS2.eps"),"eps");
  cEnDep2SS2->Print(Form("EnDep_2SS2.jpg"),"jpg");

  TCanvas *cEnDepPixSen = new TCanvas("cEnDepPixSen", "Energy deposition in sensor of pixel module (simulation)", W, H);
  hEnDepPixSen->GetXaxis()->SetLabelFont(42);
  hEnDepPixSen->GetXaxis()->SetLabelSize(0.035);
  hEnDepPixSen->GetXaxis()->SetTitleSize(0.035);
  hEnDepPixSen->GetYaxis()->SetLabelFont(42);
  hEnDepPixSen->GetYaxis()->SetLabelSize(0.035);
  hEnDepPixSen->GetYaxis()->SetTitleSize(0.035);
  hEnDepPixSen->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cEnDepPixSen->SetFillColor(0);
  cEnDepPixSen->SetBorderMode(0);
  cEnDepPixSen->SetFrameFillStyle(0);
  cEnDepPixSen->SetFrameBorderMode(0);
  cEnDepPixSen->SetLeftMargin( L/W );
  cEnDepPixSen->SetRightMargin( R/W );
  cEnDepPixSen->SetTopMargin( T/H );
  cEnDepPixSen->SetBottomMargin( B/H );
  cEnDepPixSen->SetTickx(0);
  cEnDepPixSen->SetTicky(0);
  gPad->SetLogy();
  hEnDepPixSen->GetXaxis()->SetTitle("#font[42]{Energy (MeV)}");
  hEnDepPixSen->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  //hEnDepAl->GetYaxis()->SetRangeUser(0,28.0);
  hEnDepPixSen->Draw("HIST");
  float cEnDepPixSenX = 0.5;   //l +  relPosX*(1-l-r);
  float cEnDepPixSenY = 10000.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cEnDepPixSenX,cEnDepPixSenY,extraText2);
  cEnDepPixSen->Print(Form("EnDep_Pix_Sen.eps"),"eps");
  cEnDepPixSen->Print(Form("EnDep_Pix_Sen.jpg"),"jpg");

  TCanvas *cEnDepPixROC = new TCanvas("cEnDepPixROC", "Energy deposition in ROC of pixel module (simulation)", W, H);
  hEnDepPixROC->GetXaxis()->SetLabelFont(42);
  hEnDepPixROC->GetXaxis()->SetLabelSize(0.035);
  hEnDepPixROC->GetXaxis()->SetTitleSize(0.035);
  hEnDepPixROC->GetYaxis()->SetLabelFont(42);
  hEnDepPixROC->GetYaxis()->SetLabelSize(0.035);
  hEnDepPixROC->GetYaxis()->SetTitleSize(0.035);
  hEnDepPixROC->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cEnDepPixROC->SetFillColor(0);
  cEnDepPixROC->SetBorderMode(0);
  cEnDepPixROC->SetFrameFillStyle(0);
  cEnDepPixROC->SetFrameBorderMode(0);
  cEnDepPixROC->SetLeftMargin( L/W );
  cEnDepPixROC->SetRightMargin( R/W );
  cEnDepPixROC->SetTopMargin( T/H );
  cEnDepPixROC->SetBottomMargin( B/H );
  cEnDepPixROC->SetTickx(0);
  cEnDepPixROC->SetTicky(0);
  gPad->SetLogy();
  hEnDepPixROC->GetXaxis()->SetTitle("#font[42]{Energy (MeV)}");
  hEnDepPixROC->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  //hEnDepAl->GetYaxis()->SetRangeUser(0,28.0);
  hEnDepPixROC->Draw("HIST");
  float cEnDepPixROCX = 0.5;   //l +  relPosX*(1-l-r);
  float cEnDepPixROCY = 10000.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cEnDepPixROCX,cEnDepPixROCY,extraText2);
  cEnDepPixROC->Print(Form("EnDep_Pix_ROC.eps"),"eps");
  cEnDepPixROC->Print(Form("EnDep_Pix_ROC.jpg"),"jpg");

  //Number of electrons produced in the 2S sensors per event
  TH1F *hNbe2S = (TH1F*)f->Get("238");

  TCanvas *cNbe2S = new TCanvas("cNbe2S", "Number of electrons created in 2S sensors (simulation)", W, H);
  hNbe2S->GetXaxis()->SetLabelFont(42);
  hNbe2S->GetXaxis()->SetLabelSize(0.035);
  hNbe2S->GetXaxis()->SetTitleSize(0.035);
  hNbe2S->GetYaxis()->SetLabelFont(42);
  hNbe2S->GetYaxis()->SetLabelSize(0.035);
  hNbe2S->GetYaxis()->SetTitleSize(0.035);
  hNbe2S->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cNbe2S->SetFillColor(0);
  cNbe2S->SetBorderMode(0);
  cNbe2S->SetFrameFillStyle(0);
  cNbe2S->SetFrameBorderMode(0);
  cNbe2S->SetLeftMargin( L/W );
  cNbe2S->SetRightMargin( R/W );
  cNbe2S->SetTopMargin( T/H );
  cNbe2S->SetBottomMargin( B/H );
  cNbe2S->SetTickx(0);
  cNbe2S->SetTicky(0);
  //gPad->SetLogy();
  hNbe2S->GetXaxis()->SetTitle("#font[42]{Number of e}");
  hNbe2S->GetYaxis()->SetTitle("#font[42]{Entries/bin}");
  //hEnDepAl->GetYaxis()->SetRangeUser(0,28.0);
  hNbe2S->Draw("HIST");
  float cNbe2SX = 30.;   //l +  relPosX*(1-l-r);
  float cNbe2SY = 460.; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cNbe2SX,cNbe2SY,extraText2);
  cNbe2S->Print(Form("Nbe2S.eps"),"eps");
  cNbe2S->Print(Form("Nbe2S.jpg"),"jpg");

  //Beam profile at exit of each volume
  TH2F *hxyAl = (TH2F*)f->Get("2");
  TH2F *hxyScint1 = (TH2F*)f->Get("3");
  TH2F *hxyScint2 = (TH2F*)f->Get("8");
  TH2F *hxy2S1 = (TH2F*)f->Get("5");
  TH2F *hxy2S2 = (TH2F*)f->Get("6");
  TH2F *hxyPix2 = (TH2F*)f->Get("7");

  TCanvas *cxyAl = new TCanvas("cxyAl", "Beam profile at the exit of the Al plate (simulation)", W, H);
  hxyAl->GetXaxis()->SetLabelFont(42);
  hxyAl->GetXaxis()->SetLabelSize(0.035);
  hxyAl->GetXaxis()->SetTitleSize(0.035);
  hxyAl->GetYaxis()->SetLabelFont(42);
  hxyAl->GetYaxis()->SetLabelSize(0.035);
  hxyAl->GetYaxis()->SetTitleSize(0.035);
  hxyAl->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cxyAl->SetFillColor(0);
  cxyAl->SetBorderMode(0);
  cxyAl->SetFrameFillStyle(0);
  cxyAl->SetFrameBorderMode(0);
  cxyAl->SetLeftMargin( L/W );
  cxyAl->SetRightMargin( R/W );
  cxyAl->SetTopMargin( T/H );
  cxyAl->SetBottomMargin( B/H );
  cxyAl->SetTickx(0);
  cxyAl->SetTicky(0);
  hxyAl->GetXaxis()->SetTitle("#font[42]{x (cm)}");
  hxyAl->GetYaxis()->SetTitle("#font[42]{y (cm)}");
  hxyAl->Draw();
  float cxyAlX = -2.4;   //l +  relPosX*(1-l-r);
  float cxyAlY = 2.4; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cxyAlX,cxyAlY,extraText2);
  cxyAl->Print(Form("xy_Ex_Al.eps"),"eps");
  cxyAl->Print(Form("xy_Ex_Al.jpg"),"jpg");

  TCanvas *cxyScint1 = new TCanvas("cxyScint1", "Beam profile at the exit of the first scintillator (simulation)", W, H);
  hxyScint1->GetXaxis()->SetLabelFont(42);
  hxyScint1->GetXaxis()->SetLabelSize(0.035);
  hxyScint1->GetXaxis()->SetTitleSize(0.035);
  hxyScint1->GetYaxis()->SetLabelFont(42);
  hxyScint1->GetYaxis()->SetLabelSize(0.035);
  hxyScint1->GetYaxis()->SetTitleSize(0.035);
  hxyScint1->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cxyScint1->SetFillColor(0);
  cxyScint1->SetBorderMode(0);
  cxyScint1->SetFrameFillStyle(0);
  cxyScint1->SetFrameBorderMode(0);
  cxyScint1->SetLeftMargin( L/W );
  cxyScint1->SetRightMargin( R/W );
  cxyScint1->SetTopMargin( T/H );
  cxyScint1->SetBottomMargin( B/H );
  cxyScint1->SetTickx(0);
  cxyScint1->SetTicky(0);
  hxyScint1->GetXaxis()->SetTitle("#font[42]{x (cm)}");
  hxyScint1->GetYaxis()->SetTitle("#font[42]{y (cm)}");
  hxyScint1->Draw();
  float cxyScint1X = -2.4;   //l +  relPosX*(1-l-r);
  float cxyScint1Y = 2.4; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cxyScint1X,cxyScint1Y,extraText2);
  cxyScint1->Print(Form("xy_Ex_Scint1.eps"),"eps");
  cxyScint1->Print(Form("xy_Ex_Scint1.jpg"),"jpg");

  TCanvas *cxyScint2 = new TCanvas("cxyScint2", "Beam profile at the exit of the second scintillator (simulation)", W, H);
  hxyScint2->GetXaxis()->SetLabelFont(42);
  hxyScint2->GetXaxis()->SetLabelSize(0.035);
  hxyScint2->GetXaxis()->SetTitleSize(0.035);
  hxyScint2->GetYaxis()->SetLabelFont(42);
  hxyScint2->GetYaxis()->SetLabelSize(0.035);
  hxyScint2->GetYaxis()->SetTitleSize(0.035);
  hxyScint2->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cxyScint2->SetFillColor(0);
  cxyScint2->SetBorderMode(0);
  cxyScint2->SetFrameFillStyle(0);
  cxyScint2->SetFrameBorderMode(0);
  cxyScint2->SetLeftMargin( L/W );
  cxyScint2->SetRightMargin( R/W );
  cxyScint2->SetTopMargin( T/H );
  cxyScint2->SetBottomMargin( B/H );
  cxyScint2->SetTickx(0);
  cxyScint2->SetTicky(0);
  hxyScint2->GetXaxis()->SetTitle("#font[42]{x (cm)}");
  hxyScint2->GetYaxis()->SetTitle("#font[42]{y (cm)}");
  hxyScint2->Draw();
  float cxyScint2X = -2.4;   //l +  relPosX*(1-l-r);
  float cxyScint2Y = 2.4; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cxyScint2X,cxyScint2Y,extraText2);
  cxyScint2->Print(Form("xy_Ex_Scint2.eps"),"eps");
  cxyScint2->Print(Form("xy_Ex_Scint2.jpg"),"jpg");

  TCanvas *cxy2S1 = new TCanvas("cxy2S1", "Beam profile at the exit of the first 2S sensor (simulation)", W, H);
  hxy2S1->GetXaxis()->SetLabelFont(42);
  hxy2S1->GetXaxis()->SetLabelSize(0.035);
  hxy2S1->GetXaxis()->SetTitleSize(0.035);
  hxy2S1->GetYaxis()->SetLabelFont(42);
  hxy2S1->GetYaxis()->SetLabelSize(0.035);
  hxy2S1->GetYaxis()->SetTitleSize(0.035);
  hxy2S1->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cxy2S1->SetFillColor(0);
  cxy2S1->SetBorderMode(0);
  cxy2S1->SetFrameFillStyle(0);
  cxy2S1->SetFrameBorderMode(0);
  cxy2S1->SetLeftMargin( L/W );
  cxy2S1->SetRightMargin( R/W );
  cxy2S1->SetTopMargin( T/H );
  cxy2S1->SetBottomMargin( B/H );
  cxy2S1->SetTickx(0);
  cxy2S1->SetTicky(0);
  hxy2S1->GetXaxis()->SetTitle("#font[42]{x (cm)}");
  hxy2S1->GetYaxis()->SetTitle("#font[42]{y (cm)}");
  hxy2S1->Draw();
  float cxy2S1X = -2.4;   //l +  relPosX*(1-l-r);
  float cxy2S1Y = 2.4; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cxy2S1X,cxy2S1Y,extraText2);
  cxy2S1->Print(Form("xy_Ex_2S1.eps"),"eps");
  cxy2S1->Print(Form("xy_Ex_2S1.jpg"),"jpg");

  TCanvas *cxy2S2 = new TCanvas("cxy2S2", "Beam profile at the exit of the second 2S sensor (simulation)", W, H);
  hxy2S2->GetXaxis()->SetLabelFont(42);
  hxy2S2->GetXaxis()->SetLabelSize(0.035);
  hxy2S2->GetXaxis()->SetTitleSize(0.035);
  hxy2S2->GetYaxis()->SetLabelFont(42);
  hxy2S2->GetYaxis()->SetLabelSize(0.035);
  hxy2S2->GetYaxis()->SetTitleSize(0.035);
  hxy2S2->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cxy2S2->SetFillColor(0);
  cxy2S2->SetBorderMode(0);
  cxy2S2->SetFrameFillStyle(0);
  cxy2S2->SetFrameBorderMode(0);
  cxy2S2->SetLeftMargin( L/W );
  cxy2S2->SetRightMargin( R/W );
  cxy2S2->SetTopMargin( T/H );
  cxy2S2->SetBottomMargin( B/H );
  cxy2S2->SetTickx(0);
  cxy2S2->SetTicky(0);
  hxy2S2->GetXaxis()->SetTitle("#font[42]{x (cm)}");
  hxy2S2->GetYaxis()->SetTitle("#font[42]{y (cm)}");
  hxy2S2->Draw();
  float cxy2S2X = -2.4;   //l +  relPosX*(1-l-r);
  float cxy2S2Y = 2.4; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cxy2S2X,cxy2S2Y,extraText2);
  cxy2S2->Print(Form("xy_Ex_2S2.eps"),"eps");
  cxy2S2->Print(Form("xy_Ex_2S2.jpg"),"jpg");

  TCanvas *cxyPix2 = new TCanvas("cxyPix2", "Beam profile at the exit of the pixel module (simulation)", W, H);
  hxyPix2->GetXaxis()->SetLabelFont(42);
  hxyPix2->GetXaxis()->SetLabelSize(0.035);
  hxyPix2->GetXaxis()->SetTitleSize(0.035);
  hxyPix2->GetYaxis()->SetLabelFont(42);
  hxyPix2->GetYaxis()->SetLabelSize(0.035);
  hxyPix2->GetYaxis()->SetTitleSize(0.035);
  hxyPix2->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat("mr");
  cxyPix2->SetFillColor(0);
  cxyPix2->SetBorderMode(0);
  cxyPix2->SetFrameFillStyle(0);
  cxyPix2->SetFrameBorderMode(0);
  cxyPix2->SetLeftMargin( L/W );
  cxyPix2->SetRightMargin( R/W );
  cxyPix2->SetTopMargin( T/H );
  cxyPix2->SetBottomMargin( B/H );
  cxyPix2->SetTickx(0);
  cxyPix2->SetTicky(0);
  hxyPix2->GetXaxis()->SetTitle("#font[42]{x (cm)}");
  hxyPix2->GetYaxis()->SetTitle("#font[42]{y (cm)}");
  hxyPix2->Draw();
  float cxyPix2X = -2.4;   //l +  relPosX*(1-l-r);
  float cxyPix2Y = 2.4; //1-t+lumiTextOffset*t;
  latex.SetTextSize(0.028);
  latex.SetTextAlign(12);
  latex.SetTextFont(52);
  latex.DrawLatex(cxyPix2X,cxyPix2Y,extraText2);
  cxyPix2->Print(Form("xy_Ex_Pix2.eps"),"eps");
  cxyPix2->Print(Form("xy_Ex_Pix2.jpg"),"jpg");

}
