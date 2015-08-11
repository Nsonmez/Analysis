#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <time.h>

#include <TSystem.h>
#include <TChain.h>
#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TEventList.h>
#include "TClonesArray.h"
#include <TH1D.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TMath.h>
#include <TProfile.h>
#include <TString.h>
#include "TVectorD.h"
#include "TLorentzVector.h"

#include "Utilities.h"
#include "analyze_events.h"

using namespace std;


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

int main(int argc, char*argv[])
{
    //----------------------------------------------------------------------------------------
    /////////////////////////////////////////////  Read Input Parameters  ////////////////////
    //----------------------------------------------------------------------------------------
    if (argc < 2) {printf("******Error in input parameters\n");return 1;}
	

    CommandLine c1;
    c1.parse(argc,argv);
    gROOT->ProcessLine("#include >vector<");

    string InputFileName	= c1.getValue<string>	("InputFileName");
    string OutputFileName	= c1.getValue<string>	("OutputFileName");
    //string OutputFileTag	= c1.getValue<string>	("OutputFileTag");
    int NENTRIES			= c1.getValue<int>		("NEntries");
    double cross			= c1.getValue<double>	("CrossSection");
    double efficiency		= c1.getValue<double>	("Efficiency");
	bool DATA				= c1.getValue<bool>		("DATA");
	double MH2				= c1.getValue<double>	("MH2");
	bool SaveData			= c1.getValue<bool>		("SaveData");
	
	int NJET_MAX=6;
	
    if ( !c1.check() ) return 0;
    c1.print();				// Printing the options
    
    char OutputFileTag[100];
    sprintf(OutputFileTag,"mhc%.0f",sqrt(MH2));
    
    string outputfile = "plots/"+OutputFileName + "_" + OutputFileTag + ".root";
    TFile *outf = new TFile(outputfile.c_str(),"RECREATE");
    //TFile *outf = new TFile(OutputFileName.c_str(),"RECREATE");
    
    cout << "________________________________________________________________\n";
    cout << "\n";
    cout << "time start  " << endl;
    gSystem->Exec("date '+%H:%M:%S'");
    
    //----------------------------------------------------------------------------------------
    ///////////////////////////////////////////////  Histogram Output  ///////////////////////
    //----------------------------------------------------------------------------------------
    // Book histograms
    
    
    TH1 *histJet_pt[3];
    TH1 *histJet_eta[3];
    TH1 *histJet_phi[3];
	char hist_name[100];
	char hist_title[100];
	
    //----------------------------------------------------------------------------------------
    TDirectory *jetdir= outf->mkdir("jets");
    jetdir->cd();
    
    for (int i=0; i<3; i++)
    {
        sprintf(hist_name,"jet_pt%i",i);
        histJet_pt[i] = new TH1F(hist_name, "jet P_{T}", 1000, 0.0, 1000.0);
        
        sprintf(hist_name,"jet_eta%i",i);
        histJet_eta[i] = new TH1F(hist_name, "jet eta", 100, -5.0, 5.0);
        
        sprintf(hist_name,"jet_phi%i",i);
        histJet_phi[i] = new TH1F(hist_name, "jet phi", 100, -5.0, 5.0);
    }
    
    //TH1 *jet_size 		 = new TH1F("jet_size",     "Number of Jets", 10, 0, 10.0);
    TH1 *histJet_btag 	 = new TH1F("histJet_btag", "Number of B-tagged jets", 10, 0, 10.0);
    TH1 *histJet_btag_pt = new TH1F("histJet_btag_pt", "PT of B-tagged jets", 100, 0.0, 500.0);
    TH1 *histJet_btag_eta = new TH1F("histJet_btag_eta", "ETA of B-tagged jets", 100, -5.0, +5.0);

	TH1 *histJet_btag_pt_afterbjet1cut = new TH1F("histJet_btag_pt_afterbjet1cut", "PT of B-tagged jets", 100, 0.0, 500.0);
	TH1 *histJet_btag_eta_afterbjet1cut = new TH1F("histJet_btag_eta_afterbjet1cut", "ETA of B-tagged jets", 100, -5.0, +5.0);

	
	int btag[12]={35,40,45,50,55,60,65,70,80,90,100,110};
	TH1F *histJet_btag_ptbigger[12];
	for (int i=0; i<12; i++)
	{
		sprintf(hist_name,"histJet_btag%i",btag[i]);
		sprintf(hist_title,"Number of B-tagged jets for PT > %i",btag[i]);
		histJet_btag_ptbigger[i] = new TH1F(hist_name, hist_title, 10, 0.0, 10.0);
	}

	TH1 *jet_size_cut8   = new TH1F("jet_size_cut8", "Number of Jets with 30<PT GeV", 10, 0, 10.0);
    TH1 *jet_size_cut9   = new TH1F("jet_size_cut9", "Number of Jets with 15<PT<30 GeV", 10, 0, 10.0);
    TH1 *hist_before_jet_eta = new TH1F("hist_before_jet_eta", "Jet Eta Before cut 9", 100, -5.0, 5.0);
    TH1 *hist_alpha_t	 = new TH1F("hist_alpha_t", "Apha_T PT_2/M12", 100, 0, 5.0);
    
    /* TH1 *jet_size_cut8_35	 = new TH1F("jet_size_cut8_35", "jet size PT > 35", 10, 0, 10);
     TH1 *jet_size_cut8_40	 = new TH1F("jet_size_cut8_40", "jet size PT > 40", 10, 0, 10);
     TH1 *jet_size_cut8_45	 = new TH1F("jet_size_cut8_45", "jet size PT > 45", 10, 0, 10);
     TH1 *jet_size_cut8_50	 = new TH1F("jet_size_cut8_50", "jet size PT > 50", 10, 0, 10);
     
     TH1 *jet_size_cut9_35	 = new TH1F("jet_size_cut9_35", "jet size PT > 35", 10, 0, 10);
     TH1 *jet_size_cut9_40	 = new TH1F("jet_size_cut9_40", "jet size PT > 40", 10, 0, 10);
     TH1 *jet_size_cut9_45	 = new TH1F("jet_size_cut9_45", "jet size PT > 45", 10, 0, 10);
     TH1 *jet_size_cut9_50	 = new TH1F("jet_size_cut9_50", "jet size PT > 50", 10, 0, 10);
     */
    //----------------------------------------------------------------------------------------
    TDirectory *lepdir= outf->mkdir("lepton");
    lepdir->cd();
    
    
    //TH2 *hist_gen_lepton2D		= new TH2I("hist_gen_lepton2D"	, "Number of GEN Leptons ;GEN;SIM", 10, 0, 10.0,10,0,10);
    
    //TH2 *hist_gen_elec_SIM1GEN0	= new TH2F("hist_gen_elec_SIM1GEN0"	, "Number of Elec : GEN=0  and SIM=1 ;PT;ETA", 100, 0, 100,100,-5,+5);
    //TH2 *hist_gen_elec_SIM0GEN1	= new TH2F("hist_gen_elec_SIM0GEN1"	, "Number of Elec : GEN=1  and SIM=0 ;PT;ETA", 100, 0, 100,100,-5,+5);
    
    //TH2 *hist_gen_muon_SIM1GEN0	= new TH2F("hist_gen_muon_SIM1GEN0"	, "Number of Muon : GEN=0  and SIM=1 ;PT;ETA", 100, 0, 100,100,-5,+5);
    //TH2 *hist_gen_muon_SIM0GEN1	= new TH2F("hist_gen_muon_SIM0GEN1"	, "Number of Muon : GEN=1  and SIM=0 ;PT;ETA", 100, 0, 100,100,-5,+5);
    
    //TH1 *hist_gen_elec_pt		= new TH1F("hist_gen_elec_pt"	, "Elec Pt " , 100, 0, 250);
    //TH1 *hist_gen_elec_phi		= new TH1F("hist_gen_elec_phi"	, "Elec Phi ", 100, -5, 5);
    //TH1 *hist_gen_elec_eta		= new TH1F("hist_gen_elec_eta"	, "Elec Eta ", 100, -5, 5);
    
    //TH1 *hist_gen_muon_pt		= new TH1F("hist_gen_muon_pt"  , "Muon Pt " , 100, 0, 250);
    //TH1 *hist_gen_muon_phi		= new TH1F("hist_gen_muon_phi" , "Muon Phi ", 100, -5, 5);
    //TH1 *hist_gen_muon_eta		= new TH1F("hist_gen_muon_eta" , "Muon Eta ", 100, -5, 5);
    
    TH1 *histElec_pt            = new TH1F("elec_pt1"	 , "1st elec P_{T}", 100, 0.0, 500.0);
    TH1 *histElec_phi           = new TH1F("elec_pt1_phi", "1st elec Phi  ", 100, -5.0, 5.0);
    TH1 *histElec_eta           = new TH1F("elec_pt1_eta", "1st elec Eta  ", 100, -5.0, 5.0);
    
    TH1 *histMuon_pt        = new TH1F("muon_pt1"     , "1st mu P_{T}", 100, 0.0, 500.0);
    TH1 *histMuon_phi       = new TH1F("muon_pt1_phi" , "1st mu Phi  ", 100, -5.0, 5.0);
    TH1 *histMuon_eta       = new TH1F("muon_pt1_eta" , "1st mu Eta  ", 100, -5.0, 5.0);
    
    TH1 *histLepton_pt 	= new TH1F("lepton_pt"  , "lepton P_{T} ", 100, 0.0, 500.0);
    TH1 *histLepton_eta	= new TH1F("lepton_eta" , "lepton Eta   ", 100, -5.0, 5.0);
    TH1 *histLepton_phi	= new TH1F("lepton_phi" , "lepton  Phi  ", 100, -5.0, 5.0);
    
    
    //TH1 *hist_gen_elec_deltaR	= new TH1F("hist_gen_elec_deltaR" , "DeltaR Elec Eta ", 1000, 0, 0.1);
    //TH1 *hist_gen_muon_deltaR	= new TH1F("hist_gen_muon_deltaR" , "DeltaR Muon Eta ", 1000, 0, 0.1);
    
    //TH1 *hist_gen_lepton             = new TH1F("hist_gen_lepton"            , "Number of GEN Leptons ;GEN;SIM", 10, 0, 10.0);
    //TH1 *hist_lepton_before_trig     = new TH1F("numb_lepton_before_trig"    , "Number of SIM Leptons PT>20/30(elec/muon)", 10, 0, 10.0);
    //TH1 *hist_lepton_pass_trig       = new TH1F("numb_lepton_pass_trig"      , "Number of SIM Leptons PT>20/30(elec/muon)", 10, 0, 10.0);
    //TH1 *hist_lepton_before_hardphot = new TH1F("hist_lepton_before_hardphot", "Number of SIM Leptons |eta|<2.5 and PT>20/30(elec/muon)", 10, 0, 10.0);
    //TH1 *hist_lepton_before_10gev      = new TH1F("numb_lepton_before_10gev" , "Number of SIM Leptons and PT>10(elec/muon)", 10, 0, 10.0);
    
    TH1 *lepton_invmass = new TH1F("lepton_invmass", "lepton inv mass", 100, 0, 500);
    
    TH2 *bjet_lepton_delta_eta 	= new TH2F("bjet_lepton_delta_eta"   , "Delta Eta between lepton and btagjet ", 100, 0, 250, 100, 0,   +5);
    TH2 *bjet_lepton_delta_phi 	= new TH2F("bjet_lepton_delta_phi"   , "Delta Phi between lepton and btagjet ", 100, 0, 250, 100, 0, +6.3);
    TH2 *bjet_lepton_deltaR     = new TH2F("bjet_lepton_deltaR"      , "Delta R between lepton and btagjet   ", 100, 0, 250, 100, 0,  +10);
    TH1 *hist_qlep_etajet       = new TH1F("hist_qlep_etajet"      , " lepton charge x jet eta  ", 100, -3, 3);
    TH2 *bjet_qeta				= new TH2F("bjet_qeta"      ,"Bjet Eta #cdot  lepton charge   ", 100, 0, 250, 100, 0,  +10);

	sprintf(hist_name,"top inv mass assume mhc=%f.2",sqrt(MH2));
    TH1 *hist_top_inv_mass_mhc 	= new TH1F("top_inv_mass_mhc"   , hist_name, 100, 0, 500);

	sprintf(hist_name,"miss_pz assume mhc=%f.2",sqrt(MH2));
    TH1 *hist_miss_pz_mhc 	= new TH1F("miss_pz_mhc"   ,hist_name, 100, 0, 500);
	
    //----------------------------------------------------------------------------------------
    
    TDirectory *metdir= outf->mkdir("met");
    metdir->cd();
    
    TH1 *histMET_et  = new TH1F("histMET_et" , "MET",  500,  0.0, 500.0);
    TH1 *histMET_eta = new TH1F("histMET_eta", "MET Eta", 100, -5.0, 5.0);
    TH1 *histMET_phi = new TH1F("histMET_phi", "histMET_phi Phi", 100, -5.0, 5.0);
    
    //---------------------------------------------------------------------------------------------------------
    //////////////////////////////////   BRANCHES AND CHAIN OVER INPUT FILES //////////////////////////////////
    //---------------------------------------------------------------------------------------------------------
	

	string inputfile = "../filtered_events_v3/"+InputFileName + ".root";

	TFile *inf  = new TFile(inputfile.c_str());
    TTree *data = (TTree*)inf->Get("DATA");
	
	
    TBranch *b_event = data->GetBranch("eventVar");
    EVENT_VAR event;
    b_event->SetAddress(&event);


    TBranch *b_jet[10];
    JET jet[10];
    char branch_name[10];
    
    for ( int a=0; a<NJET_MAX; a++ )
    {
        sprintf(branch_name,"jet%iVar",a+1);
        b_jet[a]= (TBranch*)data->GetBranch(branch_name);
        b_jet[a]->SetAddress( &jet[a] );
    }

    TBranch *b_lep = (TBranch*)data->GetBranch("mu1Var");
    LEPTON lep;
    b_lep->SetAddress( &lep );
    
    Long64_t numberOfEntries = b_event->GetEntries();

    cout << "Set Branch Addresses" << endl;
    
    //----------------------------------------------------------------------------------------
    /////////////////////////////////////  LOOP  Over the EVENTS  ////////////////////////////
    //----------------------------------------------------------------------------------------

    unsigned entries = 0;
    
    if (NENTRIES==-1) entries=numberOfEntries;
    else entries=NENTRIES;
    
    cout << "Reading TREE: " << numberOfEntries << " events available, \n";
    cout << "\t\t" << NENTRIES << " of them will be analyzed." << endl;
    
    //	int status[500];
    //	int pid[500];
    //	int mother[500];
    
    int event_counter0_after_trigger=0;
    int event_counter1_after_hardphoton=0;
    int event_counter2_after_10gev=0;
    int event_counter3_after_leptonpt55=0;
    int event_counter4_after_met=0;
    int event_counter4_after_alpha_t    =0;
    int event_counter5_after_btag=0;
    int event_counter5a_after_deltaEta=0;
    int event_counter6_after_topinvmass=0;
    int event_counter7_after_leptinvmass=0;
    int event_counter8_after_onejet=0;
    int event_counter9_after_onejet=0;
    int event_counter10_after_jeteta=0;
    int final_counter=0;

    //int part_counter=0;
	//int signal_counter=0;
    int nevents=0;
    int counter_saved=0;
	int decade = 0;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//___________________________________________________________________________________________________________
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Loop over all events
    for(Long64_t entry = 0; entry < entries; ++entry)
    {
        nevents++;
		//progress( entries, entry );
		double progress = 10.0*entry/(1.0*entries);
		int k = TMath::FloorNint(progress);
		if (k > decade) {   cout << 10*k << "%\t"; gSystem->Exec("date '+%H:%M:%S'"); cout << endl;	}
		decade = k;

        // Load selected branches with data from specified event
        for ( unsigned int a =0; a<NJET_MAX; a++ )	b_jet[a]->GetEvent(entry);

		b_event->GetEvent(entry);
        b_lep->GetEvent(entry);
        ////////////////////////////////////////////////////////////////////////////////////
        // 0. TRIGGER EMULATION
        ////////////////////////////////////////////////////////////////////////////////////
		/*
		1. Analysis over the lepton cut for PT > MUONPT_CUT or PT > ELECTRONPT_CUT
         in the central region then if there is a soft lepton get rid of this event.
		2. All THESE events are just for single lepton trigger passed ones.
		*/
		
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////

        if(lep.ID==11 && lep.PT < 25) continue;
        if(lep.ID==13 && lep.PT < 35) continue;

        if(lep.ID==11)
        {
            histElec_pt->Fill(lep.PT);
            histElec_eta->Fill(lep.Eta);
            histElec_phi->Fill(lep.Phi);
        }
        else
            if(lep.ID==13)
            {
                histMuon_pt->Fill(lep.PT);
                histMuon_eta->Fill(lep.Eta);
                histMuon_phi->Fill(lep.Phi);
            }

        histLepton_pt->Fill(lep.PT);
        histLepton_eta->Fill(lep.Eta);
        histLepton_phi->Fill(lep.Phi);

        //////////////////////////////////////////////////////////
        // KINEMATICAL CUT on Leptons PT
        // filter the event if there is lepton PT>55
        //////////////////////////////////////////////////////////

		//if ( lep.PT  < 50 && lep.PT > 30 ) continue;
		if ( lep.PT > 45 ) continue;

        event_counter3_after_leptonpt55++;

        ////////////////////////////////////////////////////////////////////
        //  MET cut
        ////////////////////////////////////////////////////////////////////

        histMET_et->Fill(event.MET);
        histMET_eta->Fill(event.METEta);
        histMET_phi->Fill(event.METPhi);
		
		// double MET_CUT=40;
		double MET_CUT=55;
        if (event.MET < MET_CUT ) continue;
        event_counter4_after_met++;

        ////////////////////////////////////////////////////////////////////
        //	ALPHA_T FILTER
        ////////////////////////////////////////////////////////////////////
        
        double alpha_t=0;
        double jet_inv_mass12=sqrt(2*jet[0].PT*jet[1].PT*(cosh(jet[0].Eta-jet[1].Eta)-cos(jet[0].Phi-jet[1].Phi)));
        alpha_t=(double)(jet[1].PT)/jet_inv_mass12;
        hist_alpha_t->Fill(alpha_t);
        
        //if( alpha_t > 0.4 ) continue;
        event_counter4_after_alpha_t++;

        ////////////////////////////////////////////////////////////////////
        //	BTAG FILTER
        ////////////////////////////////////////////////////////////////////

        int counter_btag=0;
		for (int j=0; j<12; j++)
		{
			counter_btag=0;
			// filter bjet tagged events
			for (int i=0; i<NJET_MAX; i++)
			{
				if( jet[i].BTag==1 && jet[i].PT > 30 && jet[i].PT < btag[j] ) counter_btag++;
			}
			histJet_btag_ptbigger[j]->Fill(counter_btag);
        }
        
        
        
        
        JET this_bjet;
        counter_btag=0;

		for (int i=0; i<NJET_MAX; i++)
		{
			if( jet[i].PT < 100 && jet[i].BTag==1 )
			{
				histJet_btag_pt->Fill(jet[i].PT);
				histJet_btag_eta->Fill(jet[i].Eta);
				counter_btag++;
				this_bjet=jet[i];
			}
		}
		histJet_btag->Fill(counter_btag);

		if ( counter_btag != 1 ) continue;

		histJet_btag_pt_afterbjet1cut->Fill(this_bjet.PT);
		histJet_btag_eta_afterbjet1cut->Fill(this_bjet.Eta);


		if ( this_bjet.PT <40 ) continue;
        event_counter5_after_btag++;

        
        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        //counter_saved++;
        //if ( SaveData && (double)cross*200*efficiency < counter_saved ) break;
        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        

		//////////////////////////////////////////////////////////////////////////////////////
		// BJET-LEPTON ETA DISTRIBUTION
		//////////////////////////////////////////////////////////////////////////////////////

		double deltaR2=pow(abs( (this_bjet.Eta)-(lep.Eta) ),2)+pow(abs( (this_bjet.Phi)-(lep.Phi) ),2);
        double deltaEta=abs( (this_bjet.Eta)-(lep.Eta) );
		bjet_lepton_delta_eta->Fill( lep.PT, deltaEta );
        bjet_lepton_delta_phi->Fill( lep.PT, abs( (this_bjet.Phi)-(lep.Phi) ));
        bjet_lepton_deltaR->Fill( lep.PT, sqrt(deltaR2) );

		hist_qlep_etajet->Fill(lep.Charge * jet[0].Eta);
        bjet_qeta->Fill(this_bjet.PT, lep.Charge * this_bjet.Eta);
    
        if ( deltaEta > 1.5 ) continue;
        event_counter5a_after_deltaEta++;
        
        //////////////////////////////////////////////////////////////////////////////////////
        // TOP INVARIANT MASS
        //////////////////////////////////////////////////////////////////////////////////////

        double miss_pz_mhc=missing_energy_pz(&lep,  &event, MH2);
			
        double top_inv_mass_mhc=top_invariant_mass(&this_bjet, &lep, &event, MH2);
			
		hist_miss_pz_mhc->Fill(miss_pz_mhc);
			
        hist_top_inv_mass_mhc->Fill(top_inv_mass_mhc);
		
        //if ( top_inv_mass_mhc < 230  ) continue;
        event_counter6_after_topinvmass++;

		//////////////////////////////////////////////////////////////////////////////////////
		// LEPTON INVARIANT MASS
		//////////////////////////////////////////////////////////////////////////////////////

		double leptoninvariantmass=lepton_invariant_mass(&lep, &event);
		lepton_invmass->Fill(leptoninvariantmass);

        if ( leptoninvariantmass > 65  ) continue;
        event_counter7_after_leptinvmass++;
        
		//////////////////////////////////////////////////////////////////////////////////////
        //	JET BRANCH
		//////////////////////////////////////////////////////////////////////////////////////
        int numb_hardjet=0;
        
        // eliminate hard jets
        JET this_hardjet;
        for (int i=0; i< NJET_MAX; i++)
        {
            if ( (jet[i].PT)>35 && abs(jet[i].Eta)< 4.9 )
            {
                this_hardjet=jet[i];
                numb_hardjet++;
            }
                
            if (i>2)continue;
            histJet_pt[i]->Fill(jet[i].PT);
            histJet_eta[i]->Fill(jet[i].Eta);
            histJet_phi[i]->Fill(jet[i].Phi);
        }
        jet_size_cut8->Fill(numb_hardjet);
        
        if (numb_hardjet !=1) continue;
        //if (numb_hardjet == 0) continue;
        event_counter8_after_onejet++;
        
        // eliminate soft jets
        JET this_softjet;
        int numb_softjet=0;
        
        for (int i=0; i< NJET_MAX; i++)
        {
            if ( (jet[i].PT)>15 && (jet[i].PT)<35 && abs(jet[i].Eta)< 4.9 )
            {
                numb_softjet++;
                this_softjet=jet[i];
            }
        }
		
        jet_size_cut9->Fill(numb_softjet);
        
		//if (numb_softjet!=1) continue;
        event_counter9_after_onejet++;
        
        // eliminate jets with central pseudorapidity distributions
        hist_before_jet_eta->Fill(this_hardjet.Eta);
		//if ( abs(this_hardjet.Eta) < 2.5 ) continue;
        event_counter10_after_jeteta++;
        
        final_counter++;
		
    }	//end of event loop
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TVectorD *crossdata = (TVectorD*)inf->Get("cross");
    TVectorD *effdata   = (TVectorD*)inf->Get("eff");
    
	//double test=(*crossdata)[0];

	outf->cd();
	TVectorD cross1(1);
	cross1[0] = cross;
	cross1.Write("cross");
	
	TVectorD eff1(12);
	eff1[0]  =  (*effdata)[0];
	eff1[1]  =  (double)event_counter3_after_leptonpt55/nevents;
	eff1[2]  =  (double)event_counter4_after_met/nevents;
	eff1[3]  =  (double)event_counter4_after_alpha_t/nevents;
	eff1[4]  =  (double)event_counter5_after_btag/nevents;
	eff1[5]  =  (double)event_counter5a_after_deltaEta/nevents;
	eff1[6]  =  (double)event_counter6_after_topinvmass/nevents;
	eff1[7]  =  (double)event_counter7_after_leptinvmass/nevents;
	eff1[8]  =  (double)event_counter8_after_onejet/nevents;
	eff1[9]  =  (double)event_counter9_after_onejet/nevents;
	eff1[10] =  (double)event_counter10_after_jeteta/nevents;
	eff1[11] =  (double)final_counter/nevents;
	eff1.Write("eff");
	
    //---------------------------------------------------------------------------------------------------------
    ////////////////////////////////  Saving the histograms into output  //////////////////////////////////////////
    //---------------------------------------------------------------------------------------------------------
    // end of event loop
    
    cout << "nevents : " << nevents << endl;
    cout << "_____________________________________________________________________" << endl;
    cout << " 0  event after trigger     : "    << event_counter0_after_trigger     << "\t" << (double)event_counter0_after_trigger/nevents		<< "\t" << (double)event_counter0_after_trigger/nevents*cross*1000 << endl;
    cout << " 1  event after hard photon : "    << event_counter1_after_hardphoton  << "\t" << (double)event_counter1_after_hardphoton/nevents  << "\t" << (double)event_counter1_after_hardphoton/nevents*cross*1000 << endl;
    cout << " 2  event after 10gev       : "    << event_counter2_after_10gev       << "\t" << (double)event_counter2_after_10gev/nevents       << "\t" << (double)event_counter2_after_10gev/nevents*cross*1000 << endl;
    cout << " 3  event after Lep(PT)>55  : "    << event_counter3_after_leptonpt55  << "\t" << (double)event_counter3_after_leptonpt55/nevents  << "\t" << (double)event_counter3_after_leptonpt55/nevents*cross*1000 << endl;
    cout << " 4  event after met > 50    : "    << event_counter4_after_met         << "\t" << (double)event_counter4_after_met/nevents         << "\t" << (double)event_counter4_after_met/nevents*cross*1000 << endl;
    cout << " 4  event after alp_t > 0.5 : "    << event_counter4_after_alpha_t     << "\t" << (double)event_counter4_after_alpha_t/nevents     << "\t" << (double)event_counter4_after_alpha_t/nevents*cross*1000 << endl;
    cout << " 5  event after btag = 1    : "    << event_counter5_after_btag        << "\t" << (double)event_counter5_after_btag/nevents        << "\t" << (double)event_counter5_after_btag/nevents*cross*1000 << endl;
    cout << " 5a event after delta Eta   : "    << event_counter5a_after_deltaEta   << "\t" << (double)event_counter5a_after_deltaEta/nevents   << "\t" << (double)event_counter5a_after_deltaEta/nevents*cross*1000 << endl;
    cout << " 6  event after topinvmass  : "    << event_counter6_after_topinvmass  << "\t" << (double)event_counter6_after_topinvmass/nevents  << "\t" << (double)event_counter6_after_topinvmass/nevents*cross*1000 << endl;
    cout << " 7  event after leptinvmass : "    << event_counter7_after_leptinvmass << "\t" << (double)event_counter7_after_leptinvmass/nevents << "\t" << (double)event_counter7_after_leptinvmass/nevents*cross*1000 << endl;
    cout << " 8  event after numb_jet =1 : "    << event_counter8_after_onejet      << "\t" << (double)event_counter8_after_onejet/nevents      << "\t" << (double)event_counter8_after_onejet/nevents*cross*1000 << endl;
    cout << " 9  event after numb_jet =1 : "    << event_counter9_after_onejet      << "\t" << (double)event_counter9_after_onejet/nevents      << "\t" << (double)event_counter9_after_onejet/nevents*cross*1000 << endl;
    cout << " 10 event after jet_eta>2.5 : "    << event_counter10_after_jeteta     << "\t" << (double)event_counter10_after_jeteta/nevents     << "\t" << (double)event_counter10_after_jeteta/nevents*cross*1000 << endl;
    cout << " 10 final counter           : "    << final_counter     << "\t" << (double)final_counter/nevents     << endl;

    cout << "_____________________________________________________________________" << endl;
	cout << " ntuple create eff : "<< (*effdata)[0];
	cout << "    analyze event eff : " << (double)final_counter/nevents;
	cout << "    final eff " << (*effdata)[0]*final_counter/nevents << endl;
    cout << " number of expct events L=1fb^-1  : "    << "\t" << (double)final_counter/nevents*(*effdata)[0]*cross*1000 ;
    cout << "\t" << InputFileName.c_str() << endl;
    
    
    //myfile.close();
    outf->Write();
    outf->Close();
    
    //end of main loop
}


