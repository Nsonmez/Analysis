#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <time.h>

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "DelphesClasses.h"
#include "Math/GenVector/VectorUtil.h"

using namespace std;

typedef struct
{
    double Pt_lepton;
    double MET;
    double Pt_btag;
    double Mt_lv;
    double Mi_lvb;
    double N_jets;
    double R_lb;
    double Phi_lb;
    double Eta_lb;
    double alpha_t;
} MyEVENT;

typedef struct
{
    //	double PVz;
    //	double PVx;
    //	double PVy;
    double cross;
    double eff;
} SAMPLE;

typedef struct
{
    //	double PVz;
    //	double PVx;
    //	double PVy;
    double HT;
    double MET;
    double METEta;
    double METPhi;
} EVENT_VAR;



typedef struct
{
    double PT;
    double Eta;
    double Phi;
    double Mass;
    double DeltaEta;
    double DeltaPhi;
    int Charge;
    int BTag;
    int TauTag;
} JET;


typedef struct
{
    double	PT;
    double  Eta;
    double  Phi;
    int Charge;
    int ID;
    //double  Particle;
} LEPTON;




double lepton_invariant_massv2(Jet *bjet, Electron *ptlep, MissingET *miss, double MH2)
{
    TLorentzVector particle_lep;
    TLorentzVector particle_bjet;
    
    particle_lep.SetPtEtaPhiM(ptlep->PT, ptlep->Eta, ptlep->Phi, 0.000);
    particle_bjet.SetPtEtaPhiM(bjet->PT, bjet->Eta,  bjet->Phi,  0.005);
    
    double plx=particle_lep.Px();
    double ply=particle_lep.Py();
    double plz=particle_lep.Pz();
    double pl =particle_lep.P();
    
    double misspx=(miss->MET)*cos(miss->Phi);
    double misspy=(miss->MET)*sin(miss->Phi);
    
    double dummy=(double)(plx*misspx+ply*misspy)/(pl)+(double)MH2/(2*pl);
    double a=pow(((double)plz/pl),2)-1;
    double b=2*(dummy)*(plz/pl);
    double c=pow(dummy,2)-pow((miss->MET),2);
    double res_sqrt=(b*b-4*a*c);
    
    if (res_sqrt<0) return 10.0;
    
    double result1=(-b+sqrt(res_sqrt))/(2*a);
    double result2=(-b-sqrt(res_sqrt))/(2*a);
    
    TLorentzVector particle_miss1;
    TLorentzVector particle_miss2;
    particle_miss1.SetPtEtaPhiE( miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result1*result1) );
    particle_miss2.SetPtEtaPhiE( miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result2*result2) );
    double lept_mass1=(particle_lep + particle_miss1).M();
    double lept_mass2=(particle_lep + particle_miss2).M();
    
    if (miss->Eta * result1 > 0) return lept_mass1 ;
    else if (miss->Eta * result2 > 0) return lept_mass2 ;
    
    return 10.0;
}

double lepton_invariant_massv2(Jet *bjet, Muon *ptlep, MissingET *miss, double MH2)
{
    TLorentzVector particle_lep;
    TLorentzVector particle_bjet;
    
    particle_lep.SetPtEtaPhiM(ptlep->PT, ptlep->Eta, ptlep->Phi, 0.000);
    particle_bjet.SetPtEtaPhiM(bjet->PT, bjet->Eta,  bjet->Phi,  0.005);
    
    double plx=particle_lep.Px();
    double ply=particle_lep.Py();
    double plz=particle_lep.Pz();
    double pl =particle_lep.P();
    
    double misspx=(miss->MET)*cos(miss->Phi);
    double misspy=(miss->MET)*sin(miss->Phi);
    
    double dummy=(double)(plx*misspx+ply*misspy)/(pl)+(double)MH2/(2*pl);
    double a=pow(((double)plz/pl),2)-1;
    double b=2*(dummy)*(plz/pl);
    double c=pow(dummy,2)-pow((miss->MET),2);
    double res_sqrt=(b*b-4*a*c);
    
    if (res_sqrt<0) return 10.0;
    
    double result1=(-b+sqrt(res_sqrt))/(2*a);
    double result2=(-b-sqrt(res_sqrt))/(2*a);
    
    TLorentzVector particle_miss1;
    TLorentzVector particle_miss2;
    particle_miss1.SetPtEtaPhiE( miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result1*result1) );
    particle_miss2.SetPtEtaPhiE( miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result2*result2) );
    double lept_mass1=(particle_lep + particle_miss1).M();
    double lept_mass2=(particle_lep + particle_miss2).M();
    
    if (miss->Eta * result1 > 0) return lept_mass1 ;
    else if (miss->Eta * result2 > 0) return lept_mass2 ;
    
    return 10.0;
}

double lepton_invariant_mass(Muon *ptlep, MissingET *miss)
{
    double plx=(ptlep->PT)*cos(ptlep->Phi);
    double ply=(ptlep->PT)*sin(ptlep->Phi);
    double misspx=(miss->MET)*cos(miss->Phi);
    double misspy=(miss->MET)*sin(miss->Phi);
    
    double result=sqrt(2*(ptlep->PT)*(miss->MET)-2*(plx*misspx+ply*misspy ));
    return result;
}

double lepton_invariant_mass(Electron *ptlep, MissingET *miss)
{
    double plx=(ptlep->PT)*cos(ptlep->Phi);
    double ply=(ptlep->PT)*sin(ptlep->Phi);
    double misspx=(miss->MET)*cos(miss->Phi);
    double misspy=(miss->MET)*sin(miss->Phi);
    
    double result=sqrt(2*(ptlep->PT)*(miss->MET)-2*(plx*misspx+ply*misspy ));
    return result;
}

double top_invariant_mass(Jet *bjet, Electron *ptlep, MissingET *miss, double MH2)
{
    
    TLorentzVector particle_lep;
    TLorentzVector particle_bjet;
    
    particle_lep.SetPtEtaPhiM(ptlep->PT, ptlep->Eta, ptlep->Phi, 0.000);
    particle_bjet.SetPtEtaPhiM(bjet->PT, bjet->Eta,  bjet->Phi,  0.005);
    
    double plx=particle_lep.Px();
    double ply=particle_lep.Py();
    double plz=particle_lep.Pz();
    double pl =particle_lep.P();
    
    double misspx=(miss->MET)*cos(miss->Phi);
    double misspy=(miss->MET)*sin(miss->Phi);
    
    double dummy=(double)(plx*misspx+ply*misspy)/(pl)+(double)MH2/(2*pl);
    double a=pow(((double)plz/pl),2)-1;
    double b=2*(dummy)*(plz/pl);
    double c=pow(dummy,2)-pow((miss->MET),2);
    double res_sqrt=(b*b-4*a*c);
    
    if (res_sqrt<0) return 10.0;
    
    double result1=(-b+sqrt(res_sqrt))/(2*a);
    double result2=(-b-sqrt(res_sqrt))/(2*a);
    
    /*
     double result=(-b+sqrt(res_sqrt))/(2*a);
     TLorentzVector particle_miss;
     particle_miss.SetPxPyPzE(misspx, misspy, result, sqrt(misspx*misspx + misspy*misspy + result*result) );
     double top_mass=(particle_bjet + particle_lep + particle_miss).M();
     */
    
    TLorentzVector particle_miss1;
    TLorentzVector particle_miss2;
    particle_miss1.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result1*result1) );
    particle_miss2.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result2*result2) );
    double top_mass1=(particle_bjet + particle_lep + particle_miss1).M();
    double top_mass2=(particle_bjet + particle_lep + particle_miss2).M();
    
//    if (false)
//    {
//        cout.width(10);
//        cout << miss->Eta<< "\t";
//        cout.width(10);
//        cout << result1 << "\t";
//        cout.width(10);
//        cout << result2 << "\t";
//        cout.width(10);
//        cout << top_mass1 << "\t";
//        cout.width(10);
//        cout << top_mass2 << "\t";
//        cout << endl;
//    }
    
    
    //    if (top_mass1 < 0 && top_mass2 > 0 ) return top_mass2;
    //    else if (top_mass1 > 0 && top_mass2 < 0 ) return top_mass1;
    //    else if (top_mass1 > 0 && top_mass2 > 0 ) return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
    //    else if (top_mass1 < 0 && top_mass2 < 0 ) return 10;
    //    else return 5;
    
    //return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
    return min( abs(174. - abs(top_mass1)), abs(174. - abs(top_mass2)) )+174.;
    
}

double top_invariant_mass(Jet *bjet, Muon *ptlep, MissingET *miss, double MH2)
{
    
    TLorentzVector particle_lep;
    TLorentzVector particle_bjet;
    
    particle_lep.SetPtEtaPhiM(ptlep->PT, ptlep->Eta, ptlep->Phi, 0.105);
    particle_bjet.SetPtEtaPhiM(bjet->PT, bjet->Eta,  bjet->Phi,  0.005);
    
    double plx=particle_lep.Px();
    double ply=particle_lep.Py();
    double plz=particle_lep.Pz();
    double pl =particle_lep.P();
    
    double misspx=(miss->MET)*cos(miss->Phi);
    double misspy=(miss->MET)*sin(miss->Phi);
    
    double dummy=(double)(plx*misspx+ply*misspy)/(pl)+(double)MH2/(2*pl);
    double a=pow(((double)plz/pl),2)-1;
    double b=2*(dummy)*(plz/pl);
    double c=pow(dummy,2)-pow((miss->MET),2);
    double res_sqrt=(b*b-4*a*c);
    
    if (res_sqrt<0) return 10.0;
    
    double result1=(-b+sqrt(res_sqrt))/(2*a);
    double result2=(-b-sqrt(res_sqrt))/(2*a);
    
    
    /*
     double result=(-b+sqrt(res_sqrt))/(2*a);
     TLorentzVector particle_miss;
     particle_miss.SetPxPyPzE(misspx, misspy, result, sqrt(misspx*misspx + misspy*misspy + result*result) );
     double top_mass=(particle_bjet + particle_lep + particle_miss).M();
     */
    
    TLorentzVector particle_miss1;
    TLorentzVector particle_miss2;
    particle_miss1.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result1*result1) );
    particle_miss2.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result2*result2) );
    double top_mass1=(particle_bjet + particle_lep + particle_miss1).M();
    double top_mass2=(particle_bjet + particle_lep + particle_miss2).M();
    
//    if (false)
//    {
//        cout.width(10);
//        cout << miss->Eta<< "\t";
//        cout.width(10);
//        cout << result1 << "\t";
//        cout.width(10);
//        cout << result2 << "\t";
//        cout.width(10);
//        cout << top_mass1 << "\t";
//        cout.width(10);
//        cout << top_mass2 << "\t";
//        cout << endl;
//    }
    
    
    //    if (top_mass1 < 0 && top_mass2 > 0 ) return top_mass2;
    //    else if (top_mass1 > 0 && top_mass2 < 0 ) return top_mass1;
    //    else if (top_mass1 > 0 && top_mass2 > 0 ) return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
    //    else if (top_mass1 < 0 && top_mass2 < 0 ) return 10;
    //    else return 5;
    
   return min( abs(174. - abs(top_mass1)), abs(174. - abs(top_mass2)) )+174.;
    
}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
double missing_energy_pz(Electron *ptlep, MissingET *miss, double MH2)
{
    double plx=(ptlep->PT)*cos(ptlep->Phi);
    double ply=(ptlep->PT)*sin(ptlep->Phi);
    double plz=(ptlep->PT)*sinh(ptlep->Eta);
    double pl  =(ptlep->PT)*cosh(ptlep->Eta); //massless lepton energy
    double misspx=(miss->MET)*cos(miss->Phi);
    double misspy=(miss->MET)*sin(miss->Phi);
    
    double dummy=(double)(plx*misspx+ply*misspy)/(pl)+(double)MH2/(2*pl);
    double a=pow(((double)plz/pl),2)-1;
    double b=2*(dummy)*(plz/pl);
    double c=pow(dummy,2)-pow((miss->MET),2);
    double res_sqrt=(b*b-4*a*c);
    
    //cout << a << "\t" << b << "\t" << c << "\t"  << res_sqrt << endl;
    if (res_sqrt<0) return 10.0;
    else
    {
        double result=(-b+sqrt(res_sqrt))/(2*a);
        return result;
    }
}

double missing_energy_pz(Muon *ptlep, MissingET *miss, double MH2)
{
    double plx=(ptlep->PT)*cos(ptlep->Phi);
    double ply=(ptlep->PT)*sin(ptlep->Phi);
    double plz=(ptlep->PT)*sinh(ptlep->Eta);
    double pl=(ptlep->PT)*cosh(ptlep->Eta);
    double misspx=(miss->MET)*cos(miss->Phi);
    double misspy=(miss->MET)*sin(miss->Phi);
    
    double dummy=(double)(plx*misspx+ply*misspy)/(pl)+(double)MH2/(2*pl);
    double a=pow(((double)plz/pl),2)-1;
    double b=2*(dummy)*(plz/pl);
    double c=pow(dummy,2)-pow((miss->MET),2);
    double res_sqrt=(b*b-4*a*c);
    
    //cout << a << "\t" << b << "\t" << c << "\t"  << res_sqrt << endl;
    if (res_sqrt<0) return 10.0;
    else
    {
        double result=(-b+sqrt(res_sqrt))/(2*a);
        return result;
    }
    
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// taken from skimming_eventsv4


//double lepton_invariant_mass(Muon *ptlep, MissingET *miss)
//{
//	double plx=(ptlep->PT)*cos(ptlep->Phi);
//	double ply=(ptlep->PT)*sin(ptlep->Phi);
//	double misspx=(miss->MET)*cos(miss->Phi);
//	double misspy=(miss->MET)*sin(miss->Phi);
//
//	double result=sqrt(2*(ptlep->PT)*(miss->MET)-2*(plx*misspx+ply*misspy ));
//	return result;
//}
//
//double lepton_invariant_mass(Electron *ptlep, MissingET *miss)
//{
//	double plx=(ptlep->PT)*cos(ptlep->Phi);
//	double ply=(ptlep->PT)*sin(ptlep->Phi);
//	double misspx=(miss->MET)*cos(miss->Phi);
//	double misspy=(miss->MET)*sin(miss->Phi);
//
//	double result=sqrt(2*(ptlep->PT)*(miss->MET)-2*(plx*misspx+ply*misspy ));
//	return result;
//}
//
//double top_invariant_mass(Jet *bjet, Electron *ptlep, MissingET *miss, double MH2)
//{
//
//	TLorentzVector particle_lep;
//	TLorentzVector particle_bjet;
//
//	particle_lep.SetPtEtaPhiM(ptlep->PT, ptlep->Eta, ptlep->Phi, 0.000);
//	particle_bjet.SetPtEtaPhiM(bjet->PT, bjet->Eta,  bjet->Phi,  0.005);
//
//	double plx=particle_lep.Px();
//	double ply=particle_lep.Py();
//	double plz=particle_lep.Pz();
//	double pl =particle_lep.P();
//
//	double misspx=(miss->MET)*cos(miss->Phi);
//	double misspy=(miss->MET)*sin(miss->Phi);
//
//	double dummy=(double)(plx*misspx+ply*misspy)/(pl)+(double)MH2/(2*pl);
//	double a=pow(((double)plz/pl),2)-1;
//	double b=2*(dummy)*(plz/pl);
//	double c=pow(dummy,2)-pow((miss->MET),2);
//	double res_sqrt=(b*b-4*a*c);
//
//	if (res_sqrt<0) return 10.0;
//
//	double result1=(-b+sqrt(res_sqrt))/(2*a);
//	double result2=(-b-sqrt(res_sqrt))/(2*a);
//
//	/*
//	 double result=(-b+sqrt(res_sqrt))/(2*a);
//	 TLorentzVector particle_miss;
//	 particle_miss.SetPxPyPzE(misspx, misspy, result, sqrt(misspx*misspx + misspy*misspy + result*result) );
//	 double top_mass=(particle_bjet + particle_lep + particle_miss).M();
//	 */
//
//	TLorentzVector particle_miss1;
//	TLorentzVector particle_miss2;
//	particle_miss1.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result1*result1) );
//	particle_miss2.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result2*result2) );
//	double top_mass1=(particle_bjet + particle_lep + particle_miss1).M();
//	double top_mass2=(particle_bjet + particle_lep + particle_miss2).M();
//
//
//	if (false)
//	{
//		cout.width(10);
//		cout << miss->Eta<< "\t";
//		cout.width(10);
//		cout << result1 << "\t";
//		cout.width(10);
//		cout << result2 << "\t";
//		cout.width(10);
//		cout << top_mass1 << "\t";
//		cout.width(10);
//		cout << top_mass2 << "\t";
//		cout << endl;
//	}
//
//
//	//    if (top_mass1 < 0 && top_mass2 > 0 ) return top_mass2;
//	//    else if (top_mass1 > 0 && top_mass2 < 0 ) return top_mass1;
//	//    else if (top_mass1 > 0 && top_mass2 > 0 ) return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
//	//    else if (top_mass1 < 0 && top_mass2 < 0 ) return 10;
//	//    else return 5;
//
//	return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
//
//}
//
//double top_invariant_mass(Jet *bjet, Muon *ptlep, MissingET *miss, double MH2)
//{
//
//	TLorentzVector particle_lep;
//	TLorentzVector particle_bjet;
//
//	particle_lep.SetPtEtaPhiM(ptlep->PT, ptlep->Eta, ptlep->Phi, 0.000);
//	particle_bjet.SetPtEtaPhiM(bjet->PT, bjet->Eta,  bjet->Phi,  0.005);
//
//	double plx=particle_lep.Px();
//	double ply=particle_lep.Py();
//	double plz=particle_lep.Pz();
//	double pl =particle_lep.P();
//
//	double misspx=(miss->MET)*cos(miss->Phi);
//	double misspy=(miss->MET)*sin(miss->Phi);
//
//	double dummy=(double)(plx*misspx+ply*misspy)/(pl)+(double)MH2/(2*pl);
//	double a=pow(((double)plz/pl),2)-1;
//	double b=2*(dummy)*(plz/pl);
//	double c=pow(dummy,2)-pow((miss->MET),2);
//	double res_sqrt=(b*b-4*a*c);
//
//	if (res_sqrt<0) return 10.0;
//
//	double result1=(-b+sqrt(res_sqrt))/(2*a);
//	double result2=(-b-sqrt(res_sqrt))/(2*a);
//
//
//	/*
//	 double result=(-b+sqrt(res_sqrt))/(2*a);
//	 TLorentzVector particle_miss;
//	 particle_miss.SetPxPyPzE(misspx, misspy, result, sqrt(misspx*misspx + misspy*misspy + result*result) );
//	 double top_mass=(particle_bjet + particle_lep + particle_miss).M();
//	 */
//
//	TLorentzVector particle_miss1;
//	TLorentzVector particle_miss2;
//	particle_miss1.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result1*result1) );
//	particle_miss2.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result2*result2) );
//	double top_mass1=(particle_bjet + particle_lep + particle_miss1).M();
//	double top_mass2=(particle_bjet + particle_lep + particle_miss2).M();
//
//	if (false)
//	{
//		cout.width(10);
//		cout << miss->Eta<< "\t";
//		cout.width(10);
//		cout << result1 << "\t";
//		cout.width(10);
//		cout << result2 << "\t";
//		cout.width(10);
//		cout << top_mass1 << "\t";
//		cout.width(10);
//		cout << top_mass2 << "\t";
//		cout << endl;
//	}
//
//
//	//    if (top_mass1 < 0 && top_mass2 > 0 ) return top_mass2;
//	//    else if (top_mass1 > 0 && top_mass2 < 0 ) return top_mass1;
//	//    else if (top_mass1 > 0 && top_mass2 > 0 ) return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
//	//    else if (top_mass1 < 0 && top_mass2 < 0 ) return 10;
//	//    else return 5;
//
//	return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
//
//}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// skimming_eventv3

//double lepton_invariant_mass(Muon *ptlep, MissingET *miss)
//{
//	double plx=(ptlep->PT)*cos(ptlep->Phi);
//	double ply=(ptlep->PT)*sin(ptlep->Phi);
//	double misspx=(miss->MET)*cos(miss->Phi);
//	double misspy=(miss->MET)*sin(miss->Phi);
//
//	double result=sqrt(2*(ptlep->PT)*(miss->MET)-2*(plx*misspx+ply*misspy ));
//	return result;
//}
//
//double lepton_invariant_mass(Electron *ptlep, MissingET *miss)
//{
//	double plx=(ptlep->PT)*cos(ptlep->Phi);
//	double ply=(ptlep->PT)*sin(ptlep->Phi);
//	double misspx=(miss->MET)*cos(miss->Phi);
//	double misspy=(miss->MET)*sin(miss->Phi);
//
//	double result=sqrt(2*(ptlep->PT)*(miss->MET)-2*(plx*misspx+ply*misspy ));
//	return result;
//}
//
////Calculates the top invariant mass using ELECTRON
////and missing energy and higgs mass
//double top_invariant_mass(Electron *ptlep, MissingET *miss)
//{
//	double plx=(ptlep->PT)*cos(ptlep->Phi);
//	double ply=(ptlep->PT)*sin(ptlep->Phi);
//	double plz=(ptlep->PT)*sinh(ptlep->Eta);
//	double pl  =(ptlep->PT)*cosh(ptlep->Eta); //massless lepton energy
//	double misspx=(miss->MET)*cos(miss->Phi);
//	double misspy=(miss->MET)*sin(miss->Phi);
//
//	double dummy=(double)(plx*misspx+ply*misspy)/(pl)+(double)MH2/(2*pl);
//	double a=pow(((double)plz/pl),2)-1;
//	double b=2*(dummy)*(plz/pl);
//	double c=pow(dummy,2)-pow((miss->MET),2);
//	double res_sqrt=(b*b-4*a*c);
//
//	//cout << a << "\t" << b << "\t" << c << "\t"  << res_sqrt << endl;
//	if (res_sqrt<0) return 10.0;
//	else
//	{
//		double result1=(-b- sqrt(res_sqrt))/(2*a);
//		double result2=(-b+sqrt(res_sqrt))/(2*a);
//		if (result1>0&&result2<0) return result1;
//		else if (result1<0&&result2>0) return result2;
//		else if (result1>0&&result2>0) return min(result1,result2);
//		return 5;
//	}
//}
//
//double top_invariant_mass(Muon *ptlep, MissingET *miss)
//{
//	double plx=(ptlep->PT)*cos(ptlep->Phi);
//	double ply=(ptlep->PT)*sin(ptlep->Phi);
//	double plz=(ptlep->PT)*sinh(ptlep->Eta);
//	double pl=(ptlep->PT)*cosh(ptlep->Eta);
//	double misspx=(miss->MET)*cos(miss->Phi);
//	double misspy=(miss->MET)*sin(miss->Phi);
//
//	double dummy=(double)(plx*misspx+ply*misspy)/(pl)+(double)MH2/(2*pl);
//	double a=pow(((double)plz/pl),2)-1;
//	double b=2*(dummy)*(plz/pl);
//	double c=pow(dummy,2)-pow((miss->MET),2);
//	double res_sqrt=(b*b-4*a*c);
//
//	//cout << a << "\t" << b << "\t" << c << "\t"  << res_sqrt << endl;
//	if (res_sqrt<0) return 10.0;
//	else
//	{
//		double result1=(-b- sqrt(res_sqrt))/(2*a);
//		double result2=(-b+sqrt(res_sqrt))/(2*a);
//		if (result1>0&&result2<0) return result1;
//		else if (result1<0&&result2>0) return result2;
//		else if (result1>0&&result2>0) return min(result1,result2);
//		return 5;
//	}
//}



