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

struct tree_structure_t{
   int f_year;
   int f_month;
   int f_day;
   int f_hour;
   int f_min;
   int f_sec;
   int f_time;
   double f_epoch;
   double f_lat;
   double f_lon;
   double f_depth;
   double f_M;
   double f_dist;
   double f_alfa_deg;
   };


void earthquakeDB(){
// ------ definizione variabili per lettura del file ------
long int EventID;
char  date[30];
char  date2[60];
char  time[30];
double lat;
double lon;
double depth;
char Author[30];
char CContributor[30];
char ContributorID[30];
char MagType[30];
double Magnitude;
char MagAuthor[30];
char EventLocationName[30];
int month;
int day;
int year;
int hour;
int min;
int sec;
double distance;
// ----- Definizione variabili per tree ----
int f_year;
int f_month;
int f_day;
int f_hour;
int f_min;
int f_sec;
int f_time;
double f_epoch;
double f_lat;
double f_lon;
double f_depth;
double f_M;
double f_dist;
double f_alfa_deg;
double tempo;
double epoch;
int R=6371000; // ----- earth radius in meters
double pi=3.141592;
double fi2;
double Dfi, Dlambda;
double a; 
double lat0=42.045679;
double fi1=lat0*pi/180;
double lon0=12.298524;
double c1;
double alfa_rad;
double alfa_deg;
// ----- Bool -----
Bool_t DEBUG=kTRUE;
int epoch_offset=1384124400;
// ----- Apertura file ------
ifstream file2read("earth1.txt");
//FILE *file2read;
//file2read = fopen("earth1.txt", "r");
TFile *f = new TFile("earth1.root","RECREATE");
// ------ tree definition -----
TTree *tf= new TTree("tf","tf");
tf->SetMaxTreeSize(1024*1024*200); //200 Mb
tree_structure_t tree_struct;
tf->Branch("f_year",&f_year);
tf->Branch("f_month",&f_month);
tf->Branch("f_day",&f_day);
tf->Branch("f_hour",&f_hour);
tf->Branch("f_min",&f_min);
tf->Branch("f_epoch",&f_epoch);
tf->Branch("f_lat",&f_lat);
tf->Branch("f_lon",&f_lon);
tf->Branch("f_depth",&f_depth);
tf->Branch("f_M",&f_M);
tf->Branch("f_dist",&f_dist);
tf->Branch("f_alfa_deg",&f_alfa_deg);
   //if (file2read != NULL) {
   for(int k=0;k<26351;k++){
   //for(int k=0;k<1000;k++){
   file2read >> EventID >> date >> time >> lat >> lon >> depth >> Author >> CContributor >> Magnitude ;
   cout << "EventID " << EventID << " date string " << date << " time string " << time << " latitude " << lat << " longitude " << lon << " depth " << depth << " Author " << Author << " CContributor "<< CContributor  <<  " magnitude " << Magnitude << endl;
   // ----- Estrazione data e ora dalla stringa ----
   sscanf(date,"%4u-%2u-%2u",&year,&month,&day);    
   cout << " date - year " << year << " month " << month << " day " << day <<endl; 
   sscanf(time,"%2u:%2u:%2u",&hour,&min,&sec);    
   cout << " time - hour " << hour << " min " << min << " sec " << sec <<endl;
   // ------ calculate epoch time ------
   struct tm t = {0};
   cout << " date 1 " << date << endl;
   strcat(date," ");   // per concatenare le due stringhe e leggere anche l'ora
   strcat(date,time);
   cout << " date 2 " << date << endl;
   std::istringstream ss(date);
      if(ss >> std::get_time(&t, "%Y-%m-%d %H:%M:%S")){
         if(DEBUG==kTRUE){
         cout << "GUARDA QUI " << std::put_time(&t, "%c") << "-----\n" << std::mktime(&t) << "----\n";
         }
      tempo=std::mktime(&t);
      }
      else{
      std::cout << "Parse failed\n";
      }
      if(DEBUG==kTRUE){
      cout << " ------- ------ ------ ------- " << endl;
      cout << " registered time " <<(long int)tempo <<  endl;
      cout << " ------- ------ ------ ------- " << endl;
      }
      //char pippo;
      //cin >>pippo;
   // ----- end epoch time calculation -----
   // ----- Distance calculation (Haversine formula)----
   fi2= lat*pi/180;
   Dfi=(lat-lat0)*pi/180;
   Dlambda=(lon-lon0)*pi/180; 
   a=pow(sin(Dfi/2),2)+cos(fi1)*cos(fi2)*pow(sin(Dlambda/2),2);
   distance=R*2*atan2(sqrt(a),sqrt(1-a));   

   //cout << " a " << a << " distance " << distance <<endl;
   // ----- calcolo dell'angolo rispetto il nord -----
      if((lat>=lat0)&&(lon>=lon0)){     // ------ primo quadrante
      c1=fabs(lat-lat0);
      alfa_rad=acos(c1/sqrt(pow((lon-lon0),2)+pow((lat-lat0),2)));
      alfa_deg=alfa_rad*180/pi;
         if(DEBUG==kTRUE){
         //cout << " primo quadrante alfa_rad " << alfa_rad << " deg " << alfa_deg << endl; 
         }
      }
      if((lat<lat0)&&(lon>=lon0)){     // ------ secondo quadrante
      c1=fabs(lat-lat0);
      alfa_rad=pi-acos(c1/sqrt(pow((lon-lon0),2)+pow((lat-lat0),2)));
      alfa_deg=alfa_rad*180/pi;
         if(DEBUG==kTRUE){
         //cout << " secondo quadrante alfa_rad " << alfa_rad << " deg " << alfa_deg << endl; 
         }
      }
      if((lat<lat0)&&(lon<lon0)){     // ------ terzo quadrante
      c1=fabs(lat-lat0);
      alfa_rad=pi+acos(c1/sqrt(pow((lon-lon0),2)+pow((lat-lat0),2)));
      alfa_deg=alfa_rad*180/pi;
         if(DEBUG==kTRUE){
         cout << " terzo quadrante alfa_rad " << alfa_rad << " deg " << alfa_deg << endl; 
         }
      }
      if((lat>lat0)&&(lon<lon0)){     // ------ quarto quadrante
      c1=fabs(lat-lat0);
      alfa_rad=3/2*pi+(pi-acos(c1/sqrt(pow((lon-lon0),2)+pow((lat-lat0),2))));
      alfa_deg=alfa_rad*180/pi;
      cout << " QUARTO QUADRANTE " << alfa_deg << endl;
         if(DEBUG==kTRUE){
         //cout << " quarto quadrante alfa_rad " << alfa_rad << " deg " << alfa_deg << endl; 
         }
      }


   // ----- Filling the tree -----
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
   f_sec=sec;
   tree_struct.f_sec=f_sec;
   //f_epoch=(tempo-epoch_offset)/60/60/24;   
   f_epoch=(tempo-epoch_offset)/60/60/24;   
   tree_struct.f_epoch=f_epoch;
   f_lat=lat;
   tree_struct.f_lat=f_lat;
   f_lon=lon;
   tree_struct.f_lon=f_lon;
   f_depth=depth;
   tree_struct.f_depth=f_depth;
   f_M=Magnitude;
   tree_struct.f_M=f_M;
   f_dist=distance;
   tree_struct.f_dist=f_dist;
   f_alfa_deg=alfa_deg;
   tree_struct.f_alfa_deg=f_alfa_deg;
   tf->Fill();
      if(year>2021||EventID<10000||day==0){
      //char pippo;
      //cin >> pippo;
      }
   } // chiude il for sul numero delle righe       
f->Write();
f->Close();
} // chiude il void principale
