//
//  analyze_events.h
//  
//
//  Created by Nasuf Sonmez on 8/4/15.
//
//

#ifndef _analyze_events_h
#define _analyze_events_h


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


//Matched cross section
const double nevents[] =
{
    
    //mhc100 icin degerler
    1000000,					// singletop_s
    1000000,	// singletop_s
    1000000,	// singletop_s
    1000000,	// singletop_s
    
    
    1000000,		//back_singletop_s
    1099999,		//back_singletop_t2
    1000000,		//back_singletop_t3
    1100000,		//back_singletop_tw
    
    845786,					// ttjets_lept
    853378,					// ttjets_semilept
    852616,					// ttjets_hadr
    
    30400000,		// wmp0jets
    3775793,		// wmp1jets
    2142360,    	// wmp2jets
    //5800.3387351386527,   // wmp3jets
    
    3006604,		// wc0jets
    1636475,		// wc1jets
    1247244,    	// wc2jets
    //1622.8438854931915,   // wc3jets
    
    1000000,		// wbb0jets
    1164558,		// wbb1jets
    664725,    	// wbb2jets
    //114.9			    	// wbb3jets
};



//void progress(int entries,  int i );
//
//
//
//void progress(int entries,  int i )
//{
//	double progress = 10.0*i/(1.0*entries);
//	int k = TMath::FloorNint(progress);
//    if (k > decade) {	std::cout << 10*k << "%\t"; gSystem->Exec("date '+%H:%M:%S'"); std::cout << std::endl;}
//	decade = k;
//}

double lepton_invariant_massv2(JET *bjet, LEPTON *ptlep, EVENT_VAR *miss, double MH2)
{
    TLorentzVector particle_lep;
    TLorentzVector particle_bjet;
    
    particle_lep.SetPtEtaPhiM(ptlep->PT, ptlep->Eta, ptlep->Phi, 0.000);
    particle_bjet.SetPtEtaPhiM(bjet->PT, bjet->Eta,  bjet->Phi,  0.005);
    
    double plx=particle_lep.Px();
    double ply=particle_lep.Py();
    double plz=particle_lep.Pz();
    double pl =particle_lep.P();
    
    double misspx=(miss->MET)*cos(miss->METPhi);
    double misspy=(miss->MET)*sin(miss->METPhi);
    
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
    particle_miss1.SetPtEtaPhiE( miss->MET, miss->METEta, miss->METPhi, sqrt( miss->MET*miss->MET + result1*result1) );
    particle_miss2.SetPtEtaPhiE( miss->MET, miss->METEta, miss->METPhi, sqrt( miss->MET*miss->MET + result2*result2) );
    double lept_mass1=(particle_lep + particle_miss1).M();
    double lept_mass2=(particle_lep + particle_miss2).M();
    
    if (miss->METEta * result1 > 0) return lept_mass1 ;
    else if (miss->METEta * result2 > 0) return lept_mass2 ;
    
    return 10.0;
}

double lepton_invariant_mass(LEPTON *ptlep, EVENT_VAR *miss)
{
    double plx=(ptlep->PT)*cos(ptlep->Phi);
    double ply=(ptlep->PT)*sin(ptlep->Phi);
    double misspx=(miss->MET)*cos(miss->METPhi);
    double misspy=(miss->MET)*sin(miss->METPhi);
    
    double result=sqrt(2*(ptlep->PT)*(miss->MET)-2*(plx*misspx+ply*misspy ));
    return result;
}


double top_invariant_mass(JET *bjet, LEPTON *ptlep, EVENT_VAR *miss, double MH2)
{
    
    TLorentzVector particle_lep;
    TLorentzVector particle_bjet;
    
    particle_lep.SetPtEtaPhiM(ptlep->PT, ptlep->Eta, ptlep->Phi, 0.000);
    particle_bjet.SetPtEtaPhiM(bjet->PT, bjet->Eta,  bjet->Phi,  0.005);
    
    double plx=particle_lep.Px();
    double ply=particle_lep.Py();
    double plz=particle_lep.Pz();
    double pl =particle_lep.P();
    
    double misspx=(miss->MET)*cos(miss->METPhi);
    double misspy=(miss->MET)*sin(miss->METPhi);
    
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
    particle_miss1.SetPtEtaPhiE(miss->MET, miss->METEta, miss->METPhi, sqrt( miss->MET*miss->MET + result1*result1) );
    particle_miss2.SetPtEtaPhiE(miss->MET, miss->METEta, miss->METPhi, sqrt( miss->MET*miss->MET + result2*result2) );
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
    //    else if (top_mass1 > 0 && top_mass2 > 0 ) return min( fabs(174.-top_mass1), fabs(174.-top_mass2) )+174.;
    //    else if (top_mass1 < 0 && top_mass2 < 0 ) return 10;
    //    else return 5;
    
    //return min( fabs(174.-top_mass1), fabs(174.-top_mass2) )+174.;
    return fmin( fabs(174. - fabs(top_mass1)), fabs(174. - fabs(top_mass2)) )+174.;
    
}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
double missing_energy_pz(LEPTON *ptlep, EVENT_VAR *miss, double MH2)
{
    double plx=(ptlep->PT)*cos(ptlep->Phi);
    double ply=(ptlep->PT)*sin(ptlep->Phi);
    double plz=(ptlep->PT)*sinh(ptlep->Eta);
    double pl  =(ptlep->PT)*cosh(ptlep->Eta); //massless lepton energy
    double misspx=(miss->MET)*cos(miss->METPhi);
    double misspy=(miss->MET)*sin(miss->METPhi);
    
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





#endif
