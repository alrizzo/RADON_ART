#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstring>
#include "TFile.h"
#include <TROOT.h>
#include <vector>
#include <deque>
#include <cmath>
#include <ctime>
#include <TChain.h>
#include <TTree.h>
#include <TF2.h>
#include <TF1.h>
#include <stdio.h>
#include <math.h>
#include <TROOT.h>
#include <TProof.h>
#include <TProofLog.h>
#include <TStopwatch.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm> // std::min_element
#include <iterator>  // std::begin, std::end

using namespace std;
using std::ios;

struct tree_struct_t{
//   int dect;
   int year;
   int month;
   int day;
   int hour;
   int min;
   double epoch;
   double f_d;
   double f_d_err;
   double f_dose_d;
   double f_dose_d_err;
   double f_dose_m;
   double f_dose_m_err;
   };


void CompAnalysis(){
// ----- tree variable definition -----
int year;
int month;
int day;
int hour;
int min;
double epoch;
double f_d;
double f_d_err;
double f_dose_d;
double f_dose_d_err;
double f_dose_m;
double f_dose_m_err;
// ----- write final histograms in file :
//TFile *f = new TFile("/Users/AlessandroR/Desktop/DoseAmbientaleCasaccia/ANALISI_GIUGNO_2020/FILES/root_data_histo/histo2analyze.root","RECREATE");
TFile *f = new TFile("histo2analyze.root","RECREATE");

// ----- Variable declaration -----
int dect;
int f_year;
int f_month;
int f_day;
int f_hour;
int f_min;
double f_time;   // numero di giorni dalla prima misura
double d;
double d_err;
double dose_h;
double dose_h_err;
double dose_d;
double dose_d_err;
double dose_m;
double dose_m_err;
double epoch_min=0;
double epoch_max=3000;
// ----- other variables -----
double f_time_old=5000;
double f_time_min=5000;
// ----- Boolean variables -----
Bool_t DEBUG=kTRUE;
// ----- calcolo dei mesi -----
int month_min = epoch_min/(60*60*24*30);
int month_max = epoch_max/(60*60*24*30);
int number_of_month = month_max-month_min;
// ----- Calcolo dei giorni -----
int day_min = 0;
int day_max = 3000;
int number_of_day = 3000;
//cout << " month_min " << epoch_min/(60*60*24*30) << " month max " << epoch_max/(60*60*24*30) << " number of months " <<(int)(epoch_max/(60*60*24*30) - epoch_min/(60*60*24*30)) << endl;
cout << " month_min " << month_min << " month max " << month_max << " number of months " << number_of_month << endl;

// ----- multithread enabling -----
int nthreads = 4;
ROOT::EnableThreadSafety(); 
// ----- Histogram Definition -----
// ----- TGraph Definition -----
vector<double> v1_time_m;
vector<double> v1_time_d;
vector<double> v1_dose_d;
vector<double> v1_dose_d_err;
vector<double> v1_dose_m;
vector<double> v1_dose_m_err;
int nbin_d=0;
int min_d;
int max_d;
int m_cont=0;
int mm;
int mm_old=0;
int dd;
int dd_old=0;
// ----- TChain Definition -----
TChain* T=new TChain("tf","tf");
//T->Add("/Users/AlessandroR/Desktop/DoseAmbientaleCasaccia/ANALISI_GIUGNO_2020/FILES/TotalTXTdata2analyze_copy/root_data/186_dal17al24NOV17.txt.root");
//T->Add("/Users/AlessandroR/Desktop/DoseAmbientaleCasaccia/ANALISI_GIUGNO_2020/FILES/TotalTXTdata2analyze_copy/root_data/*.root");
//T->Add("/Users/AlessandroR/Desktop/DoseAmbientaleCasaccia/DatiDaAnalizzare/TotalTXTdata2analyze_copy/root_data/*.root"); // --- chiama i dati ri-editi 26 Sett 2020
   // ------ Definizione della Tree da analizzare -----
   TTree *dt= new TTree("dt","dt");
   dt->SetMaxTreeSize(1024*1024*200); //200 Mb
   tree_struct_t tree_struct;
   dt->Branch("year",&year);
   dt->Branch("month",&month);
   dt->Branch("day",&day);
   dt->Branch("hour",&hour);
   dt->Branch("min",&min);
   dt->Branch("epoch",&epoch);
   dt->Branch("f_d",&f_d);
   dt->Branch("f_d_err",&f_d_err);
   dt->Branch("f_dose_d",&f_dose_d);
   dt->Branch("f_dose_d_err",&f_dose_d_err);
   dt->Branch("f_dose_m",&f_dose_m);
   dt->Branch("f_dose_m_err",&f_dose_m_err);


T->Add("/Users/AlessandroR/Desktop/DoseAmbientaleCasaccia/ANALISI_GIUGNO_2020/FILES/TotalTXTdata2analyze_copy/root_data/*.root"); // --- chiama i dati ri-editi 26 Sett 2020
int nentries;
nentries = T->GetEntries();
//nentries=10000;
   for(int i=0;i<nentries;i++){
   T->GetEntry(i+1);
   T->SetBranchAddress("f_year",&f_year);
   T->SetBranchAddress("f_month",&f_month);
   T->SetBranchAddress("f_day",&f_day);
   T->SetBranchAddress("f_hour",&f_hour);
   T->SetBranchAddress("f_min",&f_min);
   T->SetBranchAddress("f_time",&f_time); // ----- sostituito ad epoch -  f_time e' il numero di giorni dalla prima misura
   T->SetBranchAddress("d",&d);
   T->SetBranchAddress("d_err",&d_err);
   T->SetBranchAddress("dose_d",&dose_d);
   T->SetBranchAddress("dose_d_err",&dose_d_err);
   T->SetBranchAddress("dose_m",&dose_m);
   T->SetBranchAddress("dose_m_err",&dose_m_err);
   // ----- filling vectors for histos ----
   v1_dose_d.push_back(d);
   v1_dose_d_err.push_back(d_err);
   v1_time_d.push_back((int)f_time/60/60/24);
      if((f_month<13&&f_month>0)&&(f_year<20)){
      //if((f_month<13&&f_month>0)&&(f_year>2013&&f_year<2020)){
      year=f_year;
      tree_struct.year=year;
      month=f_month;
      tree_struct.month=month;
      day=f_day;
      tree_struct.day=day;
      hour=f_hour;
      tree_struct.hour=hour;
      min=f_min;
      tree_struct.min=min;
      epoch=f_time;          // ----- epoch meno l'offset
      tree_struct.epoch=epoch;
      f_d=d;
      tree_struct.f_d=f_d;
      f_d_err=d_err;
      tree_struct.f_d_err=f_d_err;
      f_dose_d=dose_d;
      tree_struct.f_dose_d=f_dose_d;
      f_dose_d_err=dose_d_err;
      tree_struct.f_dose_d_err=f_dose_d_err;
      f_dose_m=dose_m;
      tree_struct.f_dose_m=f_dose_m;
      f_dose_m_err=dose_m_err;
      tree_struct.f_dose_m_err=f_dose_m_err;
      dt->Fill();
      }
         if(f_d<0.1||f_d>0.2){
         cout << " epoch " << epoch << " f_d " << f_d << endl; 
         }
      
      // ----- Ciclo per trovare il minimo di time -----
      if(f_time/60/60/24<f_time_old){
      f_time_min=f_time/60/60/24;
      }
      

      if(DEBUG==kFALSE){
         if(d<1){ 
      cout << " -------------------------------------------------------- " << endl;
      cout << " f_year " << f_year << " f_month " << f_month << " f_day " << f_day << " f_hour " << f_hour << " f_min " << f_min << " f_time " << f_time << " d " << d << " d_err " << d_err << " dose_d " << dose_d << " dose_d_err " << dose_d_err << " dose_m " << dose_m << " dose_m_err " << dose_m_err << endl;
            
         }
      }

      if(i%10000==1){
      cout << f_day <<"/"<<f_month<<"/"<<f_year << endl;
      }
      mm=f_month;
      //`detector=dect;
      if(mm!=mm_old){
      m_cont++;
            if(d<1){ 
            v1_time_m.push_back((int)f_time/60/60/24);
            v1_dose_m.push_back(dose_m);
            v1_dose_m_err.push_back(dose_m_err);
            }
      mm_old=mm;
      }
   } // ---- chiude for su nentries
// ----- find min and max time for histo definition -----


   if(DEBUG==kFALSE){
   // ----- for daily dose histo -----  
   //cout << " daily dose - time min "  <<  v1_time_d.at(1) << " time_max " << v1_time_d.at((int)v1_time_d.size()-1) << " size " << v1_time_d.size()<< endl;
   cout << " daily dose - time min "  <<  f_time_min << " time_max " << v1_time_d.at((int)v1_time_d.size()-1) << " size " << v1_time_d.size()<< endl;
   // ----- for monthly dose histo ----- 
   cout << " monthly dose - time min " <<  v1_time_m.at(1) << " time_max " << v1_time_m.at((int)v1_time_m.size()-1) << " size " << v1_time_m.size()<< endl;
   }
// ----- Histograms definition  ----
TH1F *h1_dose_d = new TH1F("h1_dose_d","h1_dose_d",2355,0,2356);
TH1F *h1_dose_m = new TH1F("h1_dose_m","h1_dose_m",92,0,2356);
//TH1F *h1_dose_d = new TH1F("h1_dose_d","h1_dose_d",2355,0,2356);
// ----- Tree Filling -----
// ----- Histograms filling -----
   cout << "vector time size " << v1_time_d.size() <<  endl;
   for(int k=0;k<v1_time_d.size();k++){ 
   //cout << " v1_time " << v1_time_d.at(k) << " bin " << h1_dose_d->FindBin(v1_time_d.at(k)) << " dose " << v1_dose_d.at(k) << endl;
   h1_dose_d->SetBinContent(h1_dose_d->FindBin(v1_time_d.at(k)),v1_dose_d.at(k));
   }
   for(int k=0;k<v1_time_m.size();k++){ 
   cout << " v1_timei_m " << v1_time_m.at(k) << " bin " << h1_dose_m->FindBin(v1_time_m.at(k)) << " dose " << v1_dose_m.at(k) << endl;
   h1_dose_m->SetBinContent(h1_dose_m->FindBin(v1_time_m.at(k)),v1_dose_m.at(k));
   }

//TH1F *h1_dose_m = new TH1F("h1_dose_m","h1_dose_m",number_of_month,epoch_min,epoch_max);
   //if(dect==1){
   //   for(int i=0;i<v1_time_m.size();i++){
   //   h1_dose_m->SetBinContent(h1_dose_m->FindBin(v1_time_m.at(i)),v1_dose_m.at(i));
   //   h1_dose_m->SetBinError(h1_dose_m->FindBin(v1_time_m.at(i)),v1_dose_m_err.at(i));
   //   }
   //}

TCanvas *c1 = new TCanvas("c1","c1");   
gStyle->SetOptStat(0);
h1_dose_d->SetTitle("");
h1_dose_d->GetXaxis()->SetTitle("time [# day]");
h1_dose_d->GetYaxis()->SetTitle("dose rate (uSv/h)");
h1_dose_d->SetTitle("");
h1_dose_d->Draw();
//h1_dose_d->Draw("E SAME");
c1->Update();
TCanvas *c2 = new TCanvas("c2","c2");   
gStyle->SetOptStat(0);
h1_dose_m->SetTitle("");
h1_dose_m->GetXaxis()->SetTitle("epoch time [s]");
h1_dose_m->GetYaxis()->SetTitle("dose rate (nSv/h)");
h1_dose_m->SetTitle("");
h1_dose_m->Draw("P");
//h1_dose_d->Draw("E SAME");
c2->Update();
//c2->WaitPrimitive();
//c2->WaitPrimitive();
//c2->WaitPrimitive();
//c2->WaitPrimitive();
f->Write();
f->Close();

} // ------ chiude il void principale
