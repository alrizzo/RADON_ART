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
#include <TGraph.h> 
#include <TVector.h> 
#include <TGraphErrors.h> 
#include <TProofLog.h>
#include <TStopwatch.h>
#include <sstream> 
#include <stdlib.h> 
#include <stdio.h>
#include <TParameter.h>
#include <TObject.h>
using namespace std;
using std::ios;
struct tree_struct_t{
int f_year;
int f_month;
int f_day;
int f_hour;
int f_min;
double f_epoch;
double d;
double d_err;
double dose_d;
double dose_d_err;
double dose_m;
double dose_m_err;
double d_dev;
double d_theo;
   };
void RSanalysis(){
// ----- definizione variabili per leggere la tree d'ingresso ------
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
// ----- definizione variabili per scrivere la tree finale -----
int f_year;
int f_month;
int f_day;
int f_hour;
int f_min;
double f_epoch;
double d;
double d_err;
double dose_d;
double dose_d_err;
double dose_m;
double dose_m_err;
double d_dev;
double d_theo;
// ----- file 2 save
// ----- file 2 read
TFile *f = new TFile("histo2analyze.root");
// ----- Boolean Variable Definition -----
Bool_t PRINT=kFALSE;
Bool_t FIXslope=kFALSE;
// ----- Variable Definition -----
double par4=pow(10,-12);
//double par4=0;
double per_min = 0.5;
double per_max = 1.5;
//TH1F * h1 = new TH1F("Magliano Dei Marsi - monthly","Magliano dei Marsi- monthly",232, 1000000000, 1600000000);
//TH1F * h1_d = new TH1F("Magliano Dei Marsi - daily","Magliano dei Marsi - daily",6944, 1000000000, 1600000000);
TH1F * h3_m = new TH1F("Reuter-Stokes - monthly","Reuter-Stokes - monthly",92,0,2356);
TH1F * h3_d = new TH1F("Reuter-Stokes - daily","Reuter-Stokes - daily",2355,0,2356);
TH1F * h3_residuals = new TH1F("h3_residuals","h3_residuals",100, -0.02, 0.02);
TH1F * h3_dev = new TH1F("Reuter-Stokes - daily","Reuter-Stokes - daily",2355, 0, 2356);
//h3->GetXaxis()->SetRangeUser(1240000000,1550000000);
h3_m = (TH1F*)f->Get("h1_dose_m");
h3_d = (TH1F*)f->Get("h1_dose_d");
//double min_d=h3->GetXaxis()->GetXmin()/60/60/24;
//double max_d=h3->GetXaxis()->GetXmax()/60/60/24;
//TH1F * h3_m = new TH1F("Reuter-Stokes - monthly","Reuter-Stokes - monthly",232, min_d, max_d);
//TH1F * h3_m_shift = new TH1F("Reuter-Stokes - monthly shift","Reuter-Stokes - monthly shift",232, 0, max_d-min_d);
//cout << " min " << min_d << " max " << max_d << endl;
//h3_m->GetXaxis()->SetRangeUser(16000,17900);

//h3_m->GetYaxis()->SetRangeUser(100,140);
//TH1F * h3_day = new TH1F("Reuter-Stokes - daily","Reuter-Stokes - daily",6944, min_d, max_d);
//TH1F * h3_day_shift = new TH1F("Reuter-Stokes - daily shift","Reuter-Stokes - daily shift",6944, 0, max_d-min_d);
//TH1F * h3_dev = new TH1F("Reuter-Stokes - daily","Reuter-Stokes - daily",6944, min_d, max_d);
//TH1F * h3_dev_shift = new TH1F("Reuter-Stokes - daily shift","Reuter-Stokes - daily shift",6944, 0, max_d-min_d);
//h3_day->GetXaxis()->SetRangeUser(16000,18000);
//h3_day->GetYaxis()->SetRangeUser(100,140);
//h3_dev->GetXaxis()->SetRangeUser(16000,18000);
//h3_dev->GetYaxis()->SetRangeUser(100,140);
//TH1F * h3_residuals = new TH1F("Residuals","residuals",100, -10, 10);
////h3_residuals->GetXaxis()->SetRangeUser(16000,18000);
////h3_residuals->GetYaxis()->SetRangeUser(100,140);
//   for(int u=0;u<232;u++){
//   h3_m->SetBinContent(u,h3->GetBinContent(u));
//   h3_m->SetBinError(u,h3->GetBinError(u));
//   h3_m_shift->SetBinContent(u-149,h3_m->GetBinContent(u));
//   h3_m_shift->SetBinError(u-149,h3_m->GetBinError(u));
//   //cout << " bin u " << u << " valore h3 " << h3->GetBinContent(u) << " valore h3_m " << h3_m->GetBinContent(u) << " valore h3_m_shift " << h3_m_shift->GetBinContent(u) << endl;
//   }
//   for(int u=0;u<6944;u++){
//   h3_day->SetBinContent(u,h3_d->GetBinContent(u));
//   h3_day->SetBinError(u,h3_d->GetBinError(u));
//   h3_day_shift->SetBinContent(u-4446,h3_d->GetBinContent(u));
//   h3_day_shift->SetBinError(u-4446,h3_d->GetBinError(u));
//   cout << " bin u " << u << " valore h3_day " << h3_day->GetBinContent(u) <<  " valore h3_day_shift " << h3_day_shift->GetBinContent(u) << endl;
//   }
//h3_m_shift->GetXaxis()->SetRangeUser(0,1900);
//h3_m_shift->GetYaxis()->SetRangeUser(100,140);
//h3_day_shift->GetXaxis()->SetRangeUser(0,1900);
//h3_day_shift->GetYaxis()->SetRangeUser(100,140);




// ------ Histogram Preview ------
TCanvas *a1 = new TCanvas("a1","a1");
a1->Divide(1,2);
a1->cd(1);
a1->cd(1)->SetGridx();
a1->cd(1)->SetGridy();
//h3_m->SetMarkerSize(10);
h3_m->SetTitle("Monthly dose rate");
h3_m->GetXaxis()->SetTitle("time [# day]");
h3_m->GetYaxis()->SetTitle("dose rate (nSv/h)");
h3_m->SetMarkerStyle(28);
h3_m->SetMarkerColor(4);
h3_m->Draw("p");
a1->cd(1)->Update();
a1->cd(2);
a1->cd(2)->SetGridx();
a1->cd(2)->SetGridy();
h3_d->SetMarkerColor(1);
h3_d->SetMarkerSize(2);
h3_d->SetMarkerStyle(7);
h3_d->SetTitle("Daily dose rate");
h3_d->GetXaxis()->SetTitle("time [# day]");
h3_d->GetYaxis()->SetTitle("dose rate (nSv/h)");
h3_d->Draw("p");
h3_m->Draw("p SAME");
a1->cd(2)->Update();
a1->Print("Report.pdf(");

//a1->WaitPrimitive();
// ------ Prefit Step -----

TF1  *f3 = new TF1("f3","[0]+[1]*sin([2]*x+[3])+[4]*x",0,1850);
TF1  *d3 = new TF1("d3","[0]+[1]*sin([2]*x+[3])+[4]*x",0,1850);
//TF1  *f3 = new TF1("f3","[0]+[1]*cos([2]*x+[3])",16000,18000);
f3->SetLineWidth(1);
f3->SetLineStyle(4);
f3->SetParName(0,"Offset");
f3->SetParameter(0,0.115);
f3->SetParName(1,"Amplitude");
f3->SetParameter(1,0.01);
f3->SetParName(2,"Angular Frequency");
//f3->SetParameter(2,0.0178);
f3->SetParameter(2,0.0172);
f3->SetParName(3,"Phase");
f3->SetParameter(3,365);
//f3->SetParameter(3,0);
f3->SetParName(4,"Slope");
f3->SetParameter(4,pow(10,-3));
TCanvas *prova = new TCanvas("prova","prova");
prova->Divide(1,2);
prova->cd(1);
h3_m->SetMarkerStyle(28);
h3_m->SetMarkerColor(4);
h3_m->Draw("p");
prova->cd(2);
h3_d->Draw("p");
TCanvas *a2 = new TCanvas("a2","a2");
h3_m->Draw();
   if(PRINT==kTRUE){
   f3->SetLineColor(6);
   f3->Draw("SAME"); 
   a2->Update();
   a2->WaitPrimitive();
head:
   cout << " ----- Parameter to change ----- " << endl;
   cout << " 0 -  Offset - current value:" << f3->GetParameter(0)<< endl; 
   cout << " 1 -  Amplitude - current value:" << f3->GetParameter(1)<< endl; 
   cout << " 2 -  Ang. Freq. - current value:" << f3->GetParameter(2)<< endl; 
   cout << " 3 -  Phase - current value:" << f3->GetParameter(3)<< endl; 
   cout << " 4 -  Slope - current value:" << f3->GetParameter(4)<< endl;
   cout << " 5 - Nothing " << endl;
   int ans1;
   double value1;
   int ans2;
   cout << " Which parameter do you want to change ? " << endl;
   cin >> ans1;
      if((ans1==0)||(ans1==1)||(ans1==2)||(ans1==3)||(ans1==4)){ 
      cout << " Please insert the new value for parameter " << ans1 << endl;
      cin >> value1;
      f3->SetParameter(ans1,value1);
      f3->Draw("SAME"); 
      a2->Update();
      cout << " Is it ok ? (1 - yes / 0 - no ) " << endl;
      cin >> ans2;
         if(ans2==0){
         goto head;
         }
      cout << " Pre-print step is finished " << endl;
      cout << " Parameters are free to vary within:  " << per_min << " - " << per_max <<  endl;
      }
      if(ans1==5) goto tail;
      else{
      goto head;
      }
   }
tail:
   //for(int k=0;k<5;k++){
   //f3->SetParLimits(k,per_min*f3->GetParameter(k),per_max*f3->GetParameter(k));
   //}
   //f3->SetParLimits(4,0.99*f3->GetParameter(4),1.01*f3->GetParameter(4));
      if(FIXslope==kTRUE){
      f3->FixParameter(4,par4);
      }
   h3_m->Fit(f3,"WWRE0+");
   cout << " Reduced Chi Square " << f3->GetChisquare()/f3->GetNDF() << endl;
   vector<double> v_offset;
   vector<double> v_amplitude;
   vector<double> v_ang_freq;
   vector<double> v_phase;
   vector<double> v_slope;
   vector<double> v_offset_err;
   vector<double> v_amplitude_err;
   vector<double> v_ang_freq_err;
   vector<double> v_phase_err;
   vector<double> v_slope_err;
   double a_offset=0;
   double a_offset_old=0;
   double a_amplitude=0;
   double a_amplitude_old=0;
   double a_ang_freq=0;
   double a_ang_freq_old=0;
   double a_phase=0;
   double a_phase_old=0;
   double a_slope=0;
   double a_slope_old=0;
   int fit_offset=0;   
   int fit_amplitude=0;   
   int fit_ang_freq=0;   
   int fit_phase=0;   
   int fit_slope=0;
   int fit_tot=0;  
   double delta=pow(10,-5); 
   // ----- Iterative part ------
      for(int i=0;i<120;i++){
      f3->FixParameter(0,f3->GetParameter(0));
      f3->FixParameter(1,f3->GetParameter(1));
      f3->FixParameter(3,f3->GetParameter(3));
      f3->FixParameter(4,f3->GetParameter(4));
      f3->ReleaseParameter(2);
      f3->SetLineColor(i);
      h3_m->Fit(f3,"WWRE0+");
      a2->Update();
      f3->FixParameter(0,f3->GetParameter(0));
      f3->FixParameter(1,f3->GetParameter(1));
      f3->FixParameter(2,f3->GetParameter(2));
      f3->FixParameter(4,f3->GetParameter(4));
      f3->ReleaseParameter(3);
      f3->SetLineColor(i);
      h3_m->Fit(f3,"WWRE0+");
      a2->Update();
      f3->FixParameter(2,f3->GetParameter(2));
      f3->FixParameter(3,f3->GetParameter(3));
      f3->ReleaseParameter(0);
      f3->ReleaseParameter(1);
      f3->ReleaseParameter(4);
         if(FIXslope==kTRUE){
         f3->FixParameter(4,par4);
         }
      f3->SetLineColor(i);
      h3_m->Fit(f3,"WWRE0+");
       
      
 
      v_offset.push_back(f3->GetParameter(0));
      v_amplitude.push_back(f3->GetParameter(1));
      v_ang_freq.push_back(f3->GetParameter(2));
      v_phase.push_back(f3->GetParameter(3));
      v_slope.push_back(f3->GetParameter(4));
      v_offset_err.push_back(f3->GetParError(0));
      v_amplitude_err.push_back(f3->GetParError(1));
      v_ang_freq_err.push_back(f3->GetParError(2));
      v_phase_err.push_back(f3->GetParError(3));
      v_slope_err.push_back(f3->GetParError(4));
      a2->Update();
    

      a_offset=f3->GetParameter(0);
         if(fabs(a_offset-a_offset_old)/a_offset<delta){
         fit_offset=1;         
         }
      a_offset_old=a_offset;
      a_amplitude=f3->GetParameter(1);
         if(fabs(a_amplitude-a_amplitude_old)/a_amplitude<delta){
         fit_amplitude=1;         
         }
      a_amplitude_old=a_amplitude;
      a_ang_freq=f3->GetParameter(2);
         if(fabs(a_ang_freq-a_ang_freq_old)/a_ang_freq<delta){
         fit_ang_freq=1;         
         }
      a_ang_freq_old=a_ang_freq;
      a_phase=f3->GetParameter(3);
         if(fabs(a_phase-a_phase_old)/a_phase<delta){
         fit_phase=1;         
         }
      a_phase_old=a_phase;
      a_slope=f3->GetParameter(4);
         if(fabs(a_slope-a_slope_old)/a_slope<delta){
         fit_slope=1;         
         }
      a_slope_old=a_slope;
      fit_tot=fit_offset+fit_amplitude+fit_ang_freq+fit_phase+fit_slope;
         if(fit_tot==5){
         goto exit;
         }
      cout << " step " << i << " Offset " << a_offset << "-" << fit_offset << "- Amp " << a_amplitude << "-" << fit_amplitude << "- Ang. Freq " << a_ang_freq << "-"<< fit_ang_freq << "- phase " << a_phase <<"-" << fit_phase << "- slope " << a_slope << "-"<< fit_slope << endl;    
      }
exit:
      f3->ReleaseParameter(0);
      f3->ReleaseParameter(1);
      f3->ReleaseParameter(2);
      f3->ReleaseParameter(3);
      f3->ReleaseParameter(4);
         if(FIXslope==kTRUE){
         f3->FixParameter(4,par4);
         }
         //for(int k=0;k<5;k++){
         //f3->SetParLimits(k,0.10*f3->GetParameter(k),1.90*f3->GetParameter(k));
         //}
      //for(int k=0;k<5;k++){
      //f3->SetParLimits(k,per_min*f3->GetParameter(k),per_max*f3->GetParameter(k));
      //}
      //f3->FixParameter(4,0);
      cout << " ---------------------------------------- " << endl;
      gStyle->SetOptFit(11111111);
      h3_m->Fit(f3,"WWER+");
      v_offset.push_back(f3->GetParameter(0));
      v_amplitude.push_back(f3->GetParameter(1));
      v_ang_freq.push_back(f3->GetParameter(2));
      v_phase.push_back(f3->GetParameter(3));
      v_slope.push_back(f3->GetParameter(4));
      v_offset_err.push_back(f3->GetParError(0));
      v_amplitude_err.push_back(f3->GetParError(1));
      v_ang_freq_err.push_back(f3->GetParError(2));
      v_phase_err.push_back(f3->GetParError(3));
      v_slope_err.push_back(f3->GetParError(4));
      h3_m->Draw();
      f3->SetLineColor(4);
      d3->SetLineColor(4);
      f3->Draw("SAME");
      d3->SetParameter(0,f3->GetParameter(0));
      d3->SetParameter(1,f3->GetParameter(1));
      d3->SetParameter(2,f3->GetParameter(2));
      d3->SetParameter(3,f3->GetParameter(3));
      d3->SetParameter(4,f3->GetParameter(4));
      ofstream myfile;
      myfile.open("FitResults.txt",std::ios::out | std::ios::app);
      myfile <<  " RS " ; 
         for(int g=0;g<5;g++){
         myfile << f3->GetParameter(g) << " " ;
         }
         myfile <<"\n" ;
         myfile.close();
      cout << " Final Reduced Chi Square " << f3->GetChisquare()/f3->GetNDF() << endl;
      // ----- End of Iterative Part

      TCanvas *a3 = new TCanvas("a3","a3"); 
      a3->Divide(1,2);
      a3->cd(1);
      a3->cd(1)->SetGridx();
      a3->cd(1)->SetGridy();
      h3_m->Draw("p");
      d3->SetNpx(10000);
      d3->Draw("SAME");
      a3->cd(2);     
      a3->cd(2)->SetGridx();
      a3->cd(2)->SetGridy();
      h3_d->Draw("p");
      d3->Draw("SAME");
      a3->Update();
      //d3->Draw("SAME");
      //h3_d->SetRangeUser
      //f3->Draw("SAME");
      //f3->Draw("SAME");
      //f3->Draw("SAME");
      a3->Update();
      a3->Print("Report.pdf");
      // ----- Residuals Evaluation -----
      double res;
      double sigma_res;
      double off;
         for(int i=0;i<2356;i++){ 
         //for(int i=0;i<1800;i++){ 
            if(h3_d->GetBinContent(i)>0){
            cout << " bin " << i << " value = " << h3_d->GetBinContent(i) << " center " << h3_d->GetBinCenter(i) << " function " << f3->Eval(h3_d->GetBinCenter(i))<< endl;
            res = h3_d->GetBinContent(i)-(f3->Eval(h3_d->GetBinCenter(i))); 
            h3_residuals->Fill(res);
            }
         }
      TCanvas *a4 = new TCanvas("a4","a4");
      gStyle->SetOptStat("kKsSiourRmMen");
      h3_residuals->SetTitle("Model Residuals");
      h3_residuals->GetYaxis()->SetTitle("Counts");
      h3_residuals->GetXaxis()->SetTitle("dose rate (nSv/h)");
      h3_residuals->Draw();
      h3_residuals->Fit("gaus");
      sigma_res= h3_residuals->GetFunction("gaus")->GetParameter(2);
      off= h3_residuals->GetFunction("gaus")->GetParameter(1);
      //off= -0.0009;
      //cout << " sigma " <<  h3_residuals->GetFunction("gaus")->GetParameter(2) << endl;
         for(int k=0;k<2356;k++){
            if(h3_d->GetBinContent(k)>0){
               if(fabs(h3_d->GetBinContent(k)-f3->Eval(h3_d->GetBinCenter(k))-off)>1.96*sigma_res){
               h3_dev->SetBinContent(k,h3_d->GetBinContent(k));
               // ----- prova conversione -----
               
               std::time_t result = (h3_d->GetBinCenter(k)+4446)*24*60*60;
               
               //std::cout << std::asctime(std::localtime(&result)) << " epoch " << h3_day->GetBinCenter(k) << endl;
               cout <<  std::asctime(std::localtime(&result)) << " epoch " << h3_d->GetBinCenter(k) << " dose rate " <<   h3_d->GetBinContent(k) << " expected dose rate " << f3->Eval(h3_d->GetBinCenter(k))-off << " offset " << off << " sigma_res " << sigma_res << " diff " << fabs(h3_d->GetBinContent(k)-f3->Eval(h3_d->GetBinCenter(k))-off) << " 1.96 sigma " << 1.96*sigma_res <<   endl;
               //cout << " punto numero "<< k << " epoch " << min_d << " date " <<  endl;
               // ------ ----- ----- ----- ----
               //cout << " day number " << h3_day->GetBinCenter(k) << endl; 
               }
            }        
         }
      // ----- scrittura del file in uscita con tree -----
      TFile *final_result = new TFile("RS_final_result.root","RECREATE");
      TChain* T=new TChain("dt","dt");
      T->Add("histo2analyze.root");
      int nentries;
      TTree *data= new TTree("data","data");
      data->SetMaxTreeSize(1024*1024*200); //200 Mb
      tree_struct_t tree_struct;
      data->Branch("f_year",&f_year);
      data->Branch("f_month",&f_month);
      data->Branch("f_day",&f_day);
      data->Branch("f_hour",&f_hour);
      data->Branch("f_min",&f_min);
      data->Branch("f_epoch",&f_epoch);
      data->Branch("d",&d);
      data->Branch("d_err",&d_err);
      data->Branch("dose_d",&dose_d);
      data->Branch("dose_d_err",&dose_d_err);
      data->Branch("dose_m",&dose_m);
      data->Branch("dose_m_err",&dose_m_err);
      data->Branch("d_dev",&d_dev);
      data->Branch("d_theo",&d_theo);
      nentries=T->GetEntries();
      //nentries=10000;
      cout << " NENTRIES " << nentries << endl;
         for(int r=0;r<nentries;r++){
         T->GetEntry(r);
         T->SetBranchAddress("year",&year);
         T->SetBranchAddress("month",&month);
         T->SetBranchAddress("day",&day);
         T->SetBranchAddress("hour",&hour);
         T->SetBranchAddress("min",&min);
         T->SetBranchAddress("epoch",&epoch); 
         T->SetBranchAddress("f_d",&f_d);
         T->SetBranchAddress("f_d_err",&f_d_err);
         T->SetBranchAddress("f_dose_d",&f_dose_d);
         T->SetBranchAddress("f_dose_d_err",&f_dose_d_err);
         T->SetBranchAddress("f_dose_m",&f_dose_m);
         T->SetBranchAddress("f_dose_m_err",&f_dose_m_err);
         
  
         if(r%1000==0) cout << " event " << r << endl;
   
            if((year>13||year<20)&&(month>0&&month<13)){
            f_year=year;
            tree_struct.f_year=f_year;
            f_month=month;
            tree_struct.f_month=f_month;
            f_day=day;
            tree_struct.f_day=f_day;
            f_hour=hour;
            tree_struct.f_hour=f_hour;
            f_min=min;
            tree_struct.f_min=f_min;
            f_epoch=epoch;          // ----- epoch meno l'offset
            tree_struct.f_epoch=f_epoch;
            d=f_d;
            tree_struct.d=d;
            d_err=f_d_err;
            tree_struct.d_err=d_err;
            dose_d=f_dose_d;
            tree_struct.dose_d=dose_d;
            dose_d_err=f_dose_d_err;
            tree_struct.dose_d_err=dose_d_err;
            dose_m=f_dose_m;
            tree_struct.dose_m=dose_m;
            dose_m_err=f_dose_m_err;
            tree_struct.dose_m_err=dose_m_err;
            d_dev=h3_dev->GetBinContent(f_epoch/60/60/24);         
            tree_struct.d_dev=d_dev;
            d_theo=f3->Eval(f_epoch/60/60/24);         
            tree_struct.d_theo=d_theo;
         cout << " estrapolazione deviazione f_epoch " << f_epoch <<" dose rate "<< h3_dev->GetBinContent(f_epoch) << " per entry " << r << endl;  
         cout << " estrapolazione della funzione " << d3->Eval(f_epoch) << endl;
 
            data->Fill();
            } 
         }

      a4->Print("Report.pdf");
      TCanvas *a5 = new TCanvas("a5","a5");
      h3_dev->SetMarkerStyle(34);
      h3_dev->SetMarkerColor(2);
      h3_dev->GetYaxis()->SetRangeUser(0.08,0.14);
      h3_dev->Draw("p");
      h3_dev->Draw("p SAME");
      h3_d->Draw("p SAME");
      d3->Draw("SAME");
      a5->Update();
      a5->Print("Report.pdf)");

       

      TCanvas *a6 = new TCanvas("a6","a6");
      a6->Divide(1,5);
      int np = v_amplitude.size();
      double arr_offset[np];
      double arr_amplitude[np];
      double arr_ang_freq[np];
      double arr_phase[np];
      double arr_slope[np];
      double arr_offset_err[np];
      double arr_amplitude_err[np];
      double arr_ang_freq_err[np];
      double arr_phase_err[np];
      double arr_slope_err[np];
      double step[np];
      double step_err[np];
         for(int k=0;k<np;k++){
         arr_offset[k]=v_offset.at(k);
         arr_amplitude[k]=v_amplitude.at(k);
         arr_ang_freq[k]=v_ang_freq.at(k);
         arr_phase[k]=v_phase.at(k);
         arr_slope[k]=v_slope.at(k);
         arr_offset_err[k]=v_offset_err.at(k);
         arr_amplitude_err[k]=v_amplitude_err.at(k);
         arr_ang_freq_err[k]=v_ang_freq_err.at(k);
         arr_phase_err[k]=v_phase_err.at(k);
         arr_slope_err[k]=v_slope_err.at(k);
         step[k]=k;
         step_err[k]=0;
         }
      TGraphErrors *g_offset=new TGraphErrors(np,step,arr_offset,step_err,arr_offset_err);
      TGraphErrors *g_amplitude=new TGraphErrors(np,step,arr_amplitude,step_err,arr_amplitude_err);
      TGraphErrors *g_ang_freq=new TGraphErrors(np,step,arr_ang_freq,step_err,arr_ang_freq_err);
      TGraphErrors *g_phase=new TGraphErrors(np,step,arr_phase,step_err,arr_phase_err);
      TGraphErrors *g_slope=new TGraphErrors(np,step,arr_slope,step_err,arr_slope_err);

      a6->cd(1);
      g_offset->SetTitle("offset parameter convergence");
      g_offset->GetXaxis()->SetTitle("step number");
      g_offset->GetYaxis()->SetTitle("Offset parameter value");
      //g_offset->GetYaxis()->SetRangeUser(70,150);
      g_offset->SetMarkerColor(4);
      g_offset->SetMarkerStyle(21);
      g_offset->Draw("ALP");
      a6->cd(2);
      g_amplitude->SetTitle("amplitude parameter convergence");
      g_amplitude->GetXaxis()->SetTitle("step number");
      g_amplitude->GetYaxis()->SetTitle("amplitude parameter value");
      //g_amplitude->GetYaxis()->SetRangeUser(2,5);
      g_amplitude->SetMarkerColor(4);
      g_amplitude->SetMarkerStyle(21);
      g_amplitude->Draw("ALP");
      a6->cd(3);
      g_ang_freq->SetTitle("angular frequency parameter convergence");
      g_ang_freq->GetXaxis()->SetTitle("step number");
      g_ang_freq->GetYaxis()->SetTitle("angular frequency parameter value");
      //g_ang_freq->GetYaxis()->SetRangeUser(2,5);
      g_ang_freq->SetMarkerColor(4);
      g_ang_freq->SetMarkerStyle(21);
      g_ang_freq->Draw("ALP");
      a6->cd(4);
      g_phase->SetTitle("phase parameter convergence");
      g_phase->GetXaxis()->SetTitle("step number");
      g_phase->GetYaxis()->SetTitle("phase parameter value");
      //g_phase->GetYaxis()->SetRangeUser(2,5);
      g_phase->SetMarkerColor(4);
      g_phase->SetMarkerStyle(21);
      g_phase->Draw("ALP");
      a6->cd(5);
      g_slope->SetTitle("slope parameter convergence");
      g_slope->GetXaxis()->SetTitle("step number");
      g_slope->GetYaxis()->SetTitle("slope parameter value");
      //g_slope->GetYaxis()->SetRangeUser(0,0.003);
      g_slope->SetMarkerColor(4);
      g_slope->SetMarkerStyle(21);
      g_slope->Draw("ALP");

      final_result->Write();
      final_result->Close();
}// ----- chiude void principale
