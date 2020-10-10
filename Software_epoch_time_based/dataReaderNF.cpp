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
using namespace std;
using std::ios;

struct tree_structure_t{
//   int dect;
   int f_year;
   int f_month;
   int f_day;
   int f_hour;
   int f_min;
   double f_time;
   double d;
   double d_err;
   double dose_d;
   double dose_d_err;
   double dose_m;
   double dose_m_err;
   };



void dataReaderNF(){       // ---- VERSIONE GIUGNO 2020 per i nuovi dati ------
// ----- Tree Variables declaration -----
int f_year;
int f_month;
int f_day;
int f_hour;
int f_min;
double f_time;
double d;
double d_err;
double dose_d;
double dose_d_err;
double dose_m;
double dose_m_err;
// ----- general variable declaration -----
double time_offset=1384124400;
int file_num;
int month;
int day;
int year;
int hour;
int min; 
int sec;
double time;
// ----- Variables for hourly dose calculation -----
vector<double> v_dose;
int hh;
int hh_old=-1;
double sum_sigma_h;
double sigma_h;
// ----- Variables for daily dose calculation -----
int dd;
int dd_old=-1;
vector<double> v_dose_d;
double sum_dose_d;
double av_dose_d;
double sum_sigma_d;
double sigma_d;
double day_numb=1855; // contatore dei giorni
double day_numb_old; // var appoggio per contatore giorni 
// ----- Variables for monthly dose calculation -----
int mm;
int mm_old=0;
vector<double> v_dose_m;
double sum_dose_m;
double av_dose_m;
double sum_sigma_m;
double sigma_m;
// ----- Boolean variables declaration -----
Bool_t DEBUG=kTRUE;
// ----- Detector -----
//int dect;
//dect=1; // Magliano Dei Marsi
//dect=2; // Castel Del Monte
// ----- file to open and to save -----
ifstream list("file_list.txt");
string str1;
string ext=".root";
   while (std::getline(list,str1)){
   file_num++;
   cout << "lettura file list " << str1.c_str() << endl;
   unsigned found2 = str1.find_last_of("/\\");
   std::cout << " path: " << str1.substr(0,found2) << '\n';
   std::cout << " file: " << str1.substr(found2+1) << '\n';
   //ifstream file2read(str1.c_str());
   FILE *file2read;
   file2read = fopen(str1.c_str(), "r");
   string rootFile ="/root_data/";
   string path2save=str1.substr(0,found2)+rootFile;
   string fw= path2save+str1.substr(found2+1)+ext;
   cout << "*********************************************" << endl;
   cout << "Written file in: " << fw.c_str() << endl;
   cout << "*********************************************" << endl;
   TFile *f = new TFile(fw.c_str(),"RECREATE");
   // ----- Tree definition -----
   TTree *tf= new TTree("tf","tf");
   tf->SetMaxTreeSize(1024*1024*200); //200 Mb
   tree_structure_t tree_struct;
   //tf->Branch("dect",&dect);
   tf->Branch("f_year",&f_year);
   tf->Branch("f_month",&f_month);
   tf->Branch("f_day",&f_day);
   tf->Branch("f_hour",&f_hour);
   tf->Branch("f_min",&f_min);
   tf->Branch("f_time",&f_time);
   tf->Branch("d",&d);
   tf->Branch("d_err",&d_err);
   tf->Branch("dose_d",&dose_d);
   tf->Branch("dose_d_err",&dose_d_err);
   tf->Branch("dose_m",&dose_m);
   tf->Branch("dose_m_err",&dose_m_err);
      if (file2read != NULL) {
      std::string line;
         //while(fscanf(file2read, "%2u/%2u/%2u %d:%d  %lf",&day, &month, &year, &hour, &min, &d) != EOF) {   // per i dati dal 230 in poi
         while(fscanf(file2read, "%2u/%2u/%2u %2d:%2d  %lf",&month, &day, &year, &hour, &min, &d) != EOF) { //per i dati fino al 229
         int sec=0;
         if(day!=day_numb_old) day_numb++;
         day_numb_old = day;
             if(month<0||month>12){
             char pippo;
             cin>> pippo;

             }
             if(DEBUG==kTRUE){
             cout << year<<"/"<<month<<"/"<<day<<"  "<< hour<<":"<<min<<":"<<sec<<" dose " << d << endl;
             cout << " -------------------------------------------------- " << endl;
             cout << " day_numb " << day_numb << " file_num " << file_num <<endl;
             }
          f_year=year;
          f_month=month;
          f_day=day;
          f_hour=hour;
          f_min=min;
          // ------ calculate epoch time ------
	  //std::tm t = {0};
          struct tm t = {0};
	  char timestring[60];
          //sprintf(timestring,"%d-%d-%d %d:%d:%dZ",year,month,day,hour,min,sec);
          //sprintf(timestring,"20%d-%d-%d",year,month,day);
          sprintf(timestring,"20%d-%d-%d %d:%d",year,month,day,hour,min);   // funziona ma con inseriti anche minuti e ore
          //sprintf(timestring,"20%d-%d-%d",year,month,day);
	  ////std::istringstream ss("2010-11-04T23:23:01Z");
          cout.precision(17);
             if(DEBUG==kTRUE){
             cout << " timestring " << timestring << endl;
             cout << " year " << year << " month " << month << " day " << day << endl;
             }
	  std::istringstream ss(timestring);
	     if(ss >> std::get_time(&t, "%Y-%m-%d %H:%M")){
                if(DEBUG==kTRUE){
	        cout << "GUARDA QUI " << std::put_time(&t, "%c") << "-----\n" << std::mktime(&t) << "----\n";
                }
             time=std::mktime(&t);
	     }
	     else{
	     std::cout << "Parse failed\n";
             char pippo;
	     cin >>pippo;
	     }
             if(DEBUG==kTRUE){
             cout << " ------- ------ ------ ------- " << endl;
	     //cout << " registered time " <<((long int)time-time_offset)/60/60/24 <<  endl;
	     cout << " registered time " <<((long int)time-time_offset) <<  endl; // -- print dell'epoch time in secondi
             cout << " ------- ------ ------ ------- " << endl;
             }
          //f_time=(time-time_offset)/60/60/24;
          f_time=(time-time_offset);
          //   epoch=time; // ----- load on the tree variable
          //// ----- end epoch time calculation -----
          d_err=0.05*d;
          tree_struct.d=d;
          tree_struct.d_err=d_err;
          //// ----- Calculation of hourly dose  ------ 
          //hh=hour;
          //   if(hh==hh_old){
          //      if(d>0&&d<1000){
          //      v_dose.push_back(d);
          //      hh_old=hh;
          //      }
          //   }
          //   if(hh!=hh_old){
          //   // ----- Sigma hourly dose calculation -----
          //      for(int i=0;i<v_dose.size();i++){
          //      sum_sigma_h+=(v_dose.at(i) - dose_h)*(v_dose.at(i) - dose_h);
          //      }
          //      if(v_dose.size()>0){
          //      sigma_h=sqrt(sum_sigma_h/(v_dose.size()-1));
          //      }
          //   // ----- end of sigma_h calculation -----
          //   dose_h_err=sigma_h/sqrt(v_dose.size());
          //   tree_struct.dose_h_err=dose_h_err;
          //   hh_old=hh;
          //   sum_sigma_h=0;
          //   sigma_h=0;
          //   v_dose.erase(v_dose.begin(), v_dose.begin()+v_dose.size());
          //      if(d>0&&d<1000){
          //      v_dose.push_back(d);
          //      }
          //   }
          // ----- Calculation of daily dose  ------ 
          dd=day;
             if(dd==dd_old){
                if(d>0&&d<1000){
                v_dose_d.push_back(d);
                dd_old=dd;
                }
             }
             if(dd!=dd_old){
                for(int g=0;g<v_dose_d.size();g++){
                 sum_dose_d+=v_dose_d.at(g);
                }
             av_dose_d=sum_dose_d/v_dose_d.size();  
                if(v_dose_d.size()>0){
                dose_d=av_dose_d;
    		}
             tree_struct.dose_d=dose_d;  
             // ----- Sigma daily dose calculation -----
                for(int i=0;i<v_dose_d.size();i++){
                sum_sigma_d+=(v_dose_d.at(i) - av_dose_d)*(v_dose_d.at(i) - av_dose_d);
                }
                if(v_dose_d.size()>0){
                sigma_d=sqrt(sum_sigma_d/(v_dose_d.size()-1));
                dose_d_err=sigma_d/sqrt(v_dose_d.size());
                tree_struct.dose_d_err=dose_d_err;
                }
             dd_old=dd;
             sum_dose_d=0;
             av_dose_d=0;
             sum_sigma_d=0;
             sigma_d=0;
             v_dose_d.erase(v_dose_d.begin(), v_dose_d.begin()+v_dose_d.size()); 
                if(d>0&&d<1000){
                v_dose_d.push_back(d); 
                }
             }
          // ----- ----- ----- ----- ---- ----- ----
          // ----- Calculation of monthly dose  ------ 
          mm=month;
             if(mm==mm_old){
                if(d>0&&d<1000){
                v_dose_m.push_back(d);
                mm_old=mm;
                }
             }
             if(mm!=mm_old){
                for(int g=0;g<v_dose_m.size();g++){
                sum_dose_m += v_dose_m.at(g);
                }
             av_dose_m=sum_dose_m/v_dose_m.size();  
                if(v_dose_m.size()>0){
                dose_m=av_dose_m;
                }
	     tree_struct.dose_m=dose_m;  
             // ----- Sigma monthly dose calculation -----
                for(int i=0;i<v_dose_m.size();i++){
                sum_sigma_m+=(v_dose_m.at(i) - av_dose_m)*(v_dose_m.at(i) - av_dose_m);
                }
                if(v_dose_m.size()>0){
                sigma_m=sqrt(sum_sigma_m/(v_dose_m.size()-1));
                dose_m_err=sigma_m/sqrt(v_dose_m.size()); // errore sulla media e' std dev / sqrt(N)
                }
             tree_struct.dose_m_err=dose_m_err;
             //cout << " VERIFICA dose_d " << tree_struct.dose_d << " +/- " << dose_d_err << " epoch " << epoch << endl;         
             //cout << " VERIFICA dose_m " << tree_struct.dose_m << " +/- " << dose_m_err << " epoch " << epoch << endl;         
             mm_old=mm;
             sum_dose_m=0;
             av_dose_m=0;
             sum_sigma_m=0;
             sigma_m=0;
             v_dose_m.erase(v_dose_m.begin(), v_dose_m.begin()+v_dose_m.size()); 
                if(d>0&&d<1000){
                v_dose_m.push_back(d); 
                }
             }
          // ----- ----- ----- ----- ---- ----- ----
          //tree_struct.dect=dect;
          tree_struct.f_year=f_year;
          tree_struct.f_month=f_month;
          tree_struct.f_day=f_day;
          tree_struct.f_hour=f_hour;
          tree_struct.f_min=f_min;
          tree_struct.f_time=f_time;
          tf->Fill();
          }  // chiude il while sul singolo file
          //file2read.close(); // chiude il singolo file di dati letto
	  fclose(file2read);
      } // chiude il while su file_list
   cout << " uscito dal while " << endl;
   f->Write();
   f->Close();



   } // chiude il while su file_list




} // chiude il main void
