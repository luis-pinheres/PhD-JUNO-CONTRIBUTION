//-----------------------------------------------------------------
//This Code estimate the Trigger Rate of 40K per PMT for 1e+8 events 
//------------------------------------------------------------------

// Using the Channel ID convention of 32 bit
//-----------------------------------------
// Sub-detector: 8 bits 
// not Used:     5 bits
//-----------------------------------------
// Layer:        3 bits -----> 3 Layers
// Column:       3 bits -----> 3 Columns
// Wall:         3 bits -----> 7 Walls
// PMT:          4 bits -----> 16 PMTs
// ----------------------------------------
// Strip:        6 bits
// ---------------------------------------

// PhD student Luis F Pineres (IPHC-CNRS)

#include <iostream>
#include <fstream>
#include <string>
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include <bitset>
#include <iostream>
#include <cmath>
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include<bits/stdc++.h>
#include <stdio.h>
#include <stdlib.h>
#include "TLegend.h"
#include "THStack.h"
using namespace std;

string itob(int bits, int n) {
    int c;
    char s[bits+1]; // +1 to append NULL character.

    s[bits] = '\0'; // The NULL character in a character array flags the end of the string, not appending it may cause problems.

    c = bits - 1; // If the length of a string is n, than the index of the last character of the string will be n - 1. Cause the index is 0 based not 1 based. Try yourself.

    do {
        if(n%2) s[c] = '1';
        else s[c] = '0';
        n /= 2;
        c--;
    } while (n>0);

    while(c > -1) {
        s[c] = '0';
        c--;
}

    return s;
}


//--------------------------------------------------------------------------------------
int main()
{

        
// Select and Open the File

TFile *f=new TFile("/sps/hep/juno/pineresr/SUBMIT/K/hadd/RateK40_hadd_user.root ");
//TFile *f=new TFile("/pbs/home/p/pineresr/sps_juno/jdeandre/TTonly/simprod/MediumProd_50cm_TT_40K_hadd_0_user.root"); 



f->ls();
	TTree *tr=(TTree*)f->Get("TTDigit");
	

//  Add a Branch 
  	int NChannels;
    	int Channel[1000];
    	tr->SetBranchAddress("NTouchedChannel",&NChannels);
    	tr->SetBranchAddress("TB_DMchannel",&Channel);
   
// Read the elements of TB_DMchannelC
        int nentries=tr->GetEntries();
        
         int lfchannel=0;
	 int lflayers;
	 int lfrow;
	 int lfcolumn;
	 int lfpmt;
 
         int lfchannels=0;
         int lflayerss;
         int lfrows;
         int lfcolumns;
         int lfpmts;
  
         double lfRate[3][3][7][16]={0};
         int lfpm[1000]={0};         
         double lfRat[1000]={0};

         
	 for (int n=0; n<=nentries;n++)
         { 
         tr->GetEntry(n);
           
          if(NChannels==0) continue;

          for(int ich=0; ich<NChannels; ich++){
           
             
             lfchannel=Channel[ich];        
             lflayers    = (lfchannel>>16) & 0x7;
	     lfrow       = (lfchannel>>10) & 0x7;
             lfcolumn    = (lfchannel>>13) & 0x7;
             lfpmt       = (lfchannel>>6)  & 0xF;
       
                                       
          bool lfpmte=true; 
          for(int j=0; j<ich; j++){      
             lfchannels= Channel[j];
             lflayerss    = (lfchannels>>16) & 0x7;
             lfrows       = (lfchannels>>10) & 0x7;
             lfcolumns    = (lfchannels>>13) & 0x7;
             lfpmts       = (lfchannels>>6)  & 0xF;
           
            if(lfpmts==lfpmt and lfcolumns==lfcolumn and lfrows==lfrow and lflayerss==lflayers) lfpmte=false;
                            }

            if(lfpmte==true)
            lfRate[lflayers][lfcolumn][lfrow][lfpmt]++;
           
                                  
                                               }                     
       
}

         
        
       
        double efficiency[3][3][7][16]={0};
        double TriggerRate[3][3][7][16]={0};
         
        // Define the paramater of rock
        double  mass_rock; //kg)
        double rock_depth = 0.5;// m
        double rock_density  = 2.6e3; //kg/m^3
        double cavern_height = 18.5;  //m 
        double cavern_side   = 48;    //m 
        double rock_volume;
        // volume of the rock:    
        rock_volume= rock_depth * (cavern_side*cavern_side + 4*cavern_side*cavern_height + 4*(cavern_side+cavern_height)*rock_depth + 4*rock_depth*rock_depth);
        // mass of the rock
        mass_rock= rock_density * rock_volume;

        double  element_rate;   //Bq (decays/s) 
//*********************** SELECT THE ISOTOPE *****************************************

        // Radioactivity inputs for  K40 
        double  element_activity=258.4; // Bq/kg/ppm
        double  element_abundance=5;   // ppm (parts per million) 
        double  element_chain=1;

       


        element_rate = element_chain*element_activity*element_abundance*mass_rock;
        
        TH1F *layer000  = new TH1F("layer000","Expected Rate from JUNO Simulation: From 40K dcy chain",100,0,40);//LAYER 000
        
        for (int i=0; i<3; i++){
         if(i==0){  //layers
           for(int j=0; j<3; j++){ //columns
             for(int k=0; k<7; k++){ //rows
               for(int l=0; l<16; l++){ //pmts
              efficiency[i][j][k][l]=lfRate[i][j][k][l]/nentries;
              TriggerRate[i][j][k][l]=element_rate*efficiency[i][j][k][l]/1000;
              if(TriggerRate[i][j][k][l]!=0){
              cout<<TriggerRate[i][j][k][l]<<endl;
              layer000->Fill(TriggerRate[i][j][k][l]);     
                                     } }
				     } 
				    }
                                   }
                                 }

          TH1F *layer001  = new TH1F("layer001","Middle Layer",100,0,40);//LAYER 001

        for (int i=0; i<3; i++){
         if(i==1){  //layers
           for(int j=0; j<3; j++){ //columns
             for(int k=0; k<7; k++){ //rows
               for(int l=0; l<16; l++){ //pmts

              efficiency[i][j][k][l]=lfRate[i][j][k][l]/nentries;
              TriggerRate[i][j][k][l]=element_rate*efficiency[i][j][k][l]/1000;
              layer001->Fill(TriggerRate[i][j][k][l]);
                                        }
                                       }
                                      }
                                     }
                                    }



         TH1F *layer010  = new TH1F("layer010","Top Layer",100,0,40);//LAYER 010

        for (int i=0; i<3; i++){
         if(i==2){  //layers
           for(int j=0; j<3; j++){ //columns
             for(int k=0; k<7; k++){ //rows
               for(int l=0; l<16; l++){ //pmts

              efficiency[i][j][k][l]=lfRate[i][j][k][l]/nentries;
              TriggerRate[i][j][k][l]=element_rate*efficiency[i][j][k][l]/1000;
              layer010->Fill(TriggerRate[i][j][k][l]);
                                        }
                                       }
                                      }
                                     }
                                    }




          TH1F *layer  = new TH1F("layer","All Layers",50,7000,20000);//LAYER 010

        for (int i=0; i<3; i++){
           for(int j=0; j<3; j++){ //columns
             for(int k=0; k<7; k++){ //rows
               for(int l=0; l<16; l++){ //pmts

              efficiency[i][j][k][l]=lfRate[i][j][k][l]/nentries;
              TriggerRate[i][j][k][l]=element_rate*efficiency[i][j][k][l];
              layer->Fill(TriggerRate[i][j][k][l]);
                                        }
                                       }
                                      }
                                     }
                                   


//---------------------------------------------------------------------------
     TCanvas* RatesK40= new TCanvas("RatesK40","Rate-K40",500,600);
     
     gStyle->SetCanvasColor(0);
     gStyle->SetFrameBorderMode(2);
     gStyle->SetStatBorderSize(2);
     gStyle->SetOptStat(0);
     
     layer000->SetLineWidth(2);
     layer000->SetLineStyle(1);
     layer000->SetMarkerStyle(25);
     layer000->SetLineColor(2);
     layer000->GetXaxis()->SetTitle("Hits Rate sent to Concentrator Card per Pmt (kHz)");
     layer000->GetXaxis()->SetRangeUser(7,20);
     layer000->GetYaxis()->SetRangeUser(0,50);
     layer000->Draw("hist E1");
   
     layer001->SetLineWidth(2);
     layer001->SetLineStyle(1);
     layer001->SetMarkerStyle(25);
     layer001->SetLineColor(1);
     layer001->GetXaxis()->SetRangeUser(7,20);
     layer001->GetYaxis()->SetRangeUser(0,50);
     layer001->Draw("SAMES hist E1");

     layer010->SetLineWidth(2);
     layer010->SetLineStyle(1);
     layer010->SetMarkerStyle(25);
     layer010->SetLineColor(9);
     layer010->GetXaxis()->SetRangeUser(7,20);
     layer010->GetYaxis()->SetRangeUser(0,50);
     layer010->Draw("SAMES hist E1");


     layer->SetLineWidth(2);
     layer->SetLineStyle(1);
     layer->SetMarkerStyle(25);
     layer->SetLineColor(30);
     layer->GetXaxis()->SetRangeUser(5000,45000);
     layer->GetYaxis()->SetRangeUser(0,80);
     //layer->Draw("SAMES");


     RatesK40->Print("RateK40.root","root");
    

//----------------------------------------------------------------------------------------------
     TCanvas* TriggerK40= new TCanvas("TriggerK40","TriggerRateK40",900,300);
     
     gStyle->SetCanvasColor(0);
     gStyle->SetFrameBorderMode(2);
     gStyle->SetStatBorderSize(2);
     gStyle->SetTitleFillColor(40);
     gStyle->SetOptStat(11111111);


     
     TriggerK40->Divide(1,3);

     TriggerK40->cd(1);
   
     layer000->SetLineWidth(2);
     layer000->SetLineStyle(1);
     layer000->SetMarkerStyle(25);
     layer000->SetLineColor(2);
     layer000->SetTitle("Layer 000 ");
     layer000->GetXaxis()->SetTitle("Rate/Pmt");
     layer000->GetXaxis()->SetRangeUser(8000,40000);
     layer000->GetYaxis()->SetRangeUser(0,50);
     layer000->Draw();
    
     TriggerK40->cd(2);
     
     layer001->SetLineWidth(2);
     layer001->SetLineStyle(1);
     layer001->SetMarkerStyle(25);
     layer001->SetLineColor(1);
     layer001->SetTitle("Layer 001 ");
     layer001->GetXaxis()->SetTitle("Rate/Pmt");
     layer001->GetXaxis()->SetRangeUser(8000,40000);
     layer001->GetYaxis()->SetRangeUser(0,50);
     layer001->Draw();
     
     TriggerK40->cd(3);
   
     layer010->SetLineWidth(2);
     layer010->SetLineStyle(1);
     layer010->SetMarkerStyle(25);
     layer010->SetLineColor(9);
     layer010->SetTitle("Layer 010");
     layer010->GetXaxis()->SetTitle("Rate/Pmt");
     layer010->GetXaxis()->SetRangeUser(8000,40000);
     layer010->GetYaxis()->SetRangeUser(0,50);
     layer010->Draw();
     
     TriggerK40->Print("TriggerRateK40_LAYER.root","root");
     

}
















