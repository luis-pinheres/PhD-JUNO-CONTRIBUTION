// Monte Carlo Toy Model: PMTS

#include <fstream> 
#include <string>
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <endian.h>
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"
#include "TF1.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "tgmath.h"
#include <arpa/inet.h>
#include "TStyle.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
 #include "TROOT.h"
 #include "TMath.h"
 #include "TRandom.h"
 #include "TRandom3.h"
 #include "TSystem.h"
 #include "TDirectory.h"
 #include "Math/QuantFuncMathCore.h"
 #include "TUUID.h"
 #include "TApplication.h"
#include "TText.h"
#include <stdlib.h>
#include <stdio.h>

using namespace std;

// Bellamy's article: Spectrum

//Parameters:
// par[0]=N: Number of Entry in the charge distribution
// par[1]=mu: Mean of the Poisson distribution
// par[2]=Q0: Mean of the Pedestal distribution
// par[3]= Q1: Gain
// par[4]=Sigma0: RMS of the Pedestal distribution
// par[5]= Sigma1 : RMS of the photoelectrons distribution
// Dx= binning
// npe= Number of photoelectrons 
Double_t fitfp(Double_t *x,Double_t *par){
 Int_t npe=20;
 Double_t fitvalp=0;
 Double_t Qn=0;
 Double_t Sn=0;
 Double_t poissonterm=0;
for (int n=0; n<npe; n++){
 Qn=par[2]+n*par[3];
 Sn=sqrt(pow(par[4],2)+n*pow(par[5],2));
 poissonterm=pow(par[1],n)/TMath::Factorial(n)*TMath::Exp(-par[1]);
 
 fitvalp += par[0]*poissonterm*TMath::Exp(-pow(x[0]-Qn,2)/(2*Sn*Sn))/(sqrt(2*TMath::Pi())*Sn);
}

return fitvalp;
}


// Leonida's article: Spectrum 

Double_t fitleon(Double_t *x,Double_t *par){
Single photoelectron response model
 Double_t fitl=0;
 Double_t gammaterm1=0;
 Double_t gammaterm2=0;

 Double_t gammaterm=0;
 Double_t exponentialterm=0;
 //Parameters:
 //par[0]= w
 //par[1]= alpha
 //par[2]= lambda
 //par[3]= theta
 // par[4]=N: Number of Entry in the charge distribution
 // par[5]=mu: Mean of the Poisson distribution
 // par[6]=Q0: Mean of the Pedestal distribution
 // par[7]= Q1: Gain
 // par[8]=Sigma0: RMS of the Pedestal distribution
 // par[9]= Sigma1 : RMS of the photoelectrons distribution
 // Dx= binning
 
 exponentialterm=par[0]*par[1]*TMath::Exp(-par[2]*x[0]);
 gammaterm1=par[2]*(1+par[3]);
 gammaterm2=pow(gammaterm1,par[3])/TMath::Factorial(par[3]);
 gammaterm= gammaterm1+ gammaterm2*pow(x[0],par[3]);
 fitl= exponentialterm + gammaterm*TMath::Exp(-gammaterm1*x[0]);

return fitl;
}

int main(int argc, char *argv[])
{
gROOT->Reset();
gROOT->SetStyle("Plain");

if (argc > 1) {
        cout << "argv[1] = " << argv[1] << endl;
    } else {
        cout << "No file name entered. Exiting..."<<endl;
        printf("Cannot Open file, Write again the filename \n");
        return -1;
    }

    TFile *f1= new TFile(argv[1]); open the file
    Bool_t stat= f1->IsOpen();
    if ( stat ) {
   cout<<"File is open:"<<endl;

//Parameters FOR MONTECARLO TOY
double Nentries=1.0e6; // N, DELTAX=1;
double mupoi= 0.8;  //mu
double meped=30.0;  //Q0
double gain=50.0;   //Q1
double sigped=5.0;  //RMS1
double sigpht=13.0; //RMS2

TFile *MonteCarloToy = new TFile("MonteCarloToy.pdf", "RECREATE");
cout<<"SELECT YOUR CHANNEL TO FIT"<<endl;

TH1F *h2 = (TH1F*)f1->Get(argv[2]);
 
  // Number of Entries, Mean, rms in the Charge Distribution
 
  double mean_charge=0;  
  double rms_charge=0;   
  mean_charge= h2->GetMean();
  rms_charge=h2->GetRMS();

  // Parameters
  double Ntotal=0; // NORMALIZATION FACTOR
  Ntotal=h2->Integral(0,350);  

  double Q0=0; // Mean Pedestal Distribution
  Q0=  h2->GetMaximumBin(); 
   
   // RMS Pedestal Distribution
   Double_t minimumBinContent = h2->GetMinimum(0);
   double bin1 = h2->FindFirstBinAbove(h2->GetMaximum()/2);
   double bin2 = h2->FindLastBinAbove(h2->GetMaximum()/2);
   double fwhm = h2->GetBinCenter(bin2) - h2->GetBinCenter(bin1);
   //S0=fwhm/(2*sqrt(2*TMath::Log(2)));
   cout<<"bin1, bin2, fwhm = "<<bin1<<"   "<<bin2<<"   "<<fwhm<<"   "<<endl;

  double Nped=0; // Number of Entries in the Pedestal Distribution
  Nped=h2->Integral(0,Q0+8);

  double mu=0;  // Mean number of Photoelectrons
  double Q1=10; // Gain
  double rn= Ntotal/Nped;
  double mus=0;
  mus= TMath::Log(rn);
  mu= (mean_charge-Q0)/Q1;
  double Npeds=0;
  Npeds=Ntotal*exp(-mu);
  //double Q1=0; // Gain
  //Q1= (mean_charge-Q0)/mu;
  double S1=0;  // RMS NPE Distribution
  S1=0.3*Q1;
  double fvari=0, fvari2=0;
  fvari=rms_charge*rms_charge;
  fvari2=1.09*Q1*Q1;
  double S0=0;
  S0=sqrt(fvari-(mu*fvari2));

  int binmax = h2->GetMaximumBin(); 
  double xy = h2->GetBinContent(binmax);

cout<<"--------------------------------------------"<<endl;
  cout<<"Mean given by the ChargeDist:"<<mean_charge<<endl;
  cout<<"RMS given by the  ChargeDist:"<<rms_charge<<endl;
  cout<<"Number of Entries           :"<<Ntotal<<endl<<endl;

  cout<<"*************PARAMETERS**************************************"<<endl;
  cout<<"Number of Entries Pedestal(mu /I)   :"<<Nped<<"	"<<Npeds<<endl;
  cout<<"Mean number of Photoelectrons:"<<mu<<"		"<<mus<<""<<endl;
  cout<<"Mean pedestal Distribution   :"<<Q0<<"          	"<<endl;
  cout<<"Gain                         :"<<Q1<<"          	"<<endl;
  cout<<"RMS pedestal Distribution    :"<<S0<<"          	"<<endl;
  cout<<"RMS NPE Distribution         :"<<S1<<"          	"<<endl;
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"X position and max    content:"<<binmax<<"             "<<xy<<endl;
  cout<<"-------------------------------------------------------"<<endl;

TF1 *funcp = new TF1 ("funcp", fitfp,0,350,6);
funcp->SetParameters(Ntotal, mus, Q0, Q1, S0, S1);
funcp->SetParNames("Event", "MeanPois", "MeanPed", "Gain", "RMSPED", "RMSPhot");
h2->Fit("funcp","","",0,150);

//Parameter leonida's paper
double theta=7, lambda=0,alpha=0, w=0;
lambda= 1/(rms_charge*(sqrt(8)));
alpha= 1/Q1;//mus/(mean_charge-Q0);
w=alpha*(lambda*Q1-1)/(lambda-alpha);
TF1 *funcleo = new TF1 ("funcleo", fitleon,0,350,4);
funcleo->SetParameters(w,alpha,lambda,theta);
funcleo->SetParNames("w", "#alpha", "#lambda", "#theta");
h2->SetLineColor(kGreen);
h2->SetLineWidth(2);
h2->Fit("funcleo","+","",0,150);

TCanvas* c1pos= new TCanvas("c1pos","Position",400,400);
 gROOT->SetStyle("Plain");
  gStyle->GetAttDate()->SetTextColor(1);
  gStyle->SetTextFont(132);
  gStyle->SetTitleFont(132,"XYZ");
  gStyle->SetStatH(0.12);
  gStyle->SetStatW(0.15);
  gStyle->SetOptFit(1111);
  gStyle->SetCanvasColor(2);
  gStyle->SetFrameBorderMode(2);
  gStyle->SetStatBorderSize(2);
  gStyle->SetOptStat(1111);

  h2->SetTitle("Monte Carlo Toy: Fitting Charge Distribution");
  h2->GetXaxis()->SetTitle("Charge[ADC] ");
  h2->GetYaxis()->SetTitle("Number of Entries ");
  h2->SetLineColor(kBlue);
  h2->SetLineWidth(2);
  h2->Draw();
  h2->Write();
  MonteCarloToy->cd();
  MonteCarloToy->Print();
  c1pos->SetLogy();
  c1pos->Print("ChargeFit.pdf","pdf");
  cout<<"***************************************************************"<<endl;
  double chi2=0;
  double degree=0;
  chi2 = funcp->GetChisquare();
  degree= funcp->GetNDF();
  cout<<"Chi Square        = "<<chi2<<endl;
  cout<<"Degree of freedom = "<<degree<<endl;
  cout<<" #chi /Ndf        = "<<chi2/degree<<endl;

  double Chi2=0;
  double Degree=0;
  Chi2 = funcleo->GetChisquare();
  Degree= funcleo->GetNDF();
  cout<<"Chi Square        = "<<Chi2<<endl;
  cout<<"Degree of freedom = "<<Degree<<endl;
  cout<<" #chi /Ndf        = "<<Chi2/Degree<<endl;

f1->Close();
                }

        return 0;
}

