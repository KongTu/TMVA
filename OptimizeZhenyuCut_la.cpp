#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector>
#include <iostream>
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include <sstream>

#include "TF1.h"
#include "TH1.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"

#include "RooFit.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHist.h"

using namespace std;

TSystem->Load("libRooFit");

void OptimizeZhenyuCut_la(){

TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_corrected/JetTripPb_data_June4_2014.root");
//TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb/JetTripPb_Data.root");
TTree* theTree = (TTree*)file->Get("ana/v0_Lambda");

using namespace RooFit;

TCanvas* s1 = new TCanvas("s1","test",600,800);

theTree->Draw("la_mass>>myHist(360,1.05,1.17)","la_pt > 9 && la_pt < 12 && TMath::Abs(la_dau1_dzos) > 1 & TMath::Abs(la_dau2_dzos)>1 & TMath::Abs(la_dau1_dxyos)>1 & TMath::Abs(la_dau2_dxyos)>1 & la_dau1_nhit > 3 & la_dau2_nhit >3 & la_agl > 0.999 & la_dlos > 5");
TH1F* myHist = (TH1F*) gDirectory->Get( "myHist" );

RooRealVar x("x","mass",1,1.2);
RooDataHist data("data","dataset",x,myHist);
RooPlot* xframe = x.frame(160);

data.plotOn(xframe, Name("data"));

RooRealVar mean("mean","mean",1.115,1.11,1.12);
RooRealVar sigma1("sigma1","sigma1",0.003,0.001,0.01);
RooRealVar sigma2("sigma2","sigma2",0.003,0.001,0.01);
RooRealVar sig1("sig1","signal1",10,0,10000000);
RooRealVar sig2("sig2","signal2",10,0,10000000);
RooRealVar a("a","a",0,-100000,100000);
RooRealVar b("b","b",0,-100000,100000);
RooRealVar cp("cp","cp",0,-100000,100000);
RooRealVar d("d","d",0,-100000,100000);

RooRealVar f("f","f",0,-100000,100000);
RooRealVar g("g","g",0,-100000,100000);
RooRealVar h("h","h",0,-100000,100000);
RooRealVar k("k","k",0,-100000,100000);

RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
RooRealVar polysig("polysig","polysig",10,0,10000000);
RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));


x.setRange("cut",1.085,1.165);

sum.fitTo(data,Range("cut"));
sum.fitTo(data,Range("cut"));
sum.fitTo(data,Range("cut"));
sum.fitTo(data,Range("cut"));
sum.fitTo(data,Range("cut"));
sum.fitTo(data,Range("cut"));

sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(0.5),LineColor(kRed));
sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(0.5),LineColor(kRed));

xframe->Draw();

double chi2  = xframe->chiSquare("sum","data");
double meanf  = mean.getVal();
double meanfe  = mean.getError();
double sigmaf1  = sigma1.getVal();
double sigmaf2  = sigma2.getVal();
double bkgf  = polysig.getVal();
double sigf1  = sig1.getVal();
double sigf2  = sig2.getVal();
double sigwf1  = sigf1 /(sigf1 +sigf2 );
double sigwf2  = sigf2 /(sigf1 +sigf2 );
double c1 = a.getVal();
double c2 = b.getVal();

double sigmaf  = sqrt(sigmaf1 **2*sigwf1  + sigmaf2 **2*sigwf2 );
double massmin  = meanf  - 2*sigmaf ;
double massmax  = meanf  + 2*sigmaf ;


int nmin  = myHist ->GetXaxis()->FindBin(massmin );
int nmax  = myHist ->GetXaxis()->FindBin(massmax );
int anmin  = myHist ->GetXaxis()->FindBin(1.0);
int anmax  = myHist ->GetXaxis()->FindBin(1.2);

double awyh1  = myHist ->Integral(anmin ,nmin );
double awyh2  = myHist ->Integral(nmax ,anmax );
double awyh  = awyh1  + awyh2 ;
double totyh  = myHist ->Integral(nmin ,nmax );

cout << "let's see the total integral value:" << totyh <<endl;

x.setRange("cut",massmin ,massmax );
RooAbsReal* ibkg  = poly.createIntegral(x,NormSet(x),Range("cut"));
RooAbsReal* isig1  = gaus1.createIntegral(x,NormSet(x),Range("cut"));
RooAbsReal* isig2  = gaus2.createIntegral(x,NormSet(x),Range("cut"));
double ibkgf  = ibkg ->getVal();
double bkgfe  = polysig.getError();
double isig1f  = isig1 ->getVal();
double isig2f  = isig2 ->getVal();

double bkgy  = ibkgf *bkgf ;
double bkgye  = ibkgf *bkgfe ;
double sigy1  = isig1f *sigf1 ;
double sigy2  = isig2f *sigf2 ;
double sigy  = sigy1  + sigy2 ;
double toty  = bkgy  + sigy ;

double abkgy  = (1-ibkgf )*bkgf ;
double asigy1  = (1-isig1f )*sigf1 ;
double asigy2  = (1-isig2f )*sigf2 ;
double asigy  = asigy1  + asigy2 ;
double awy  = abkgy  + asigy ;

double sigfrac  = sigy /toty ;
double bkgfrac  = bkgy /toty ;

double sigyh  = totyh  - bkgy ;
double sigfrach  = sigyh /totyh ;
double bkgfrach  = bkgy /totyh ;

double signif  = sigyh / sqrt( totyh );

cout << "let's see the signal yield: " << sigyh << endl;
cout << "let's see the background yield: " << bkgy << endl;
cout << "let's see the total yield: " << totyh <<endl;
cout << "let's see the signal significance: " << signif <<endl;


}