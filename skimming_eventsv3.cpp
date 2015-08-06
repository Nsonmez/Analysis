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
#include "TLorentzVector.h"
#include "TCanvas.h"
#include <TMath.h>
#include <TProfile.h>
#include "TStyle.h"
#include <time.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include "TLatex.h"
#include "TLegend.h"
#include "Utilities.h"
#include "TVectorD.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "DelphesClasses.h"

#include "skimming_events.h"

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
    //string JetAlgo		= c1.getValue<string>	("JetAlgo");
    //vector<string> HLTbit	= c1.getVector<string> 	("HLTbit","");
    int NENTRIES			= c1.getValue<int>		("NEntries");
    bool DATA				= c1.getValue<bool>		("DATA");
    bool SaveData			= c1.getValue<bool>		("SaveData");
    double MH2				= c1.getValue<double>	("MH2");
    double cross			= c1.getValue<double>	("CrossSection");
    double efficiency		= c1.getValue<double>	("Efficiency");
    
    double MUONPT_CUT=23.;
    double ELECPT_CUT=33.;
    double MET_CUT=40.;
	
	if (!c1.check()) return 0;
	c1.print(); // Printing the options
	
    char OutputFileTag[100];
    sprintf(OutputFileTag,"mhc%.0f",sqrt(MH2));
    
    string outputfile = OutputFileName + "_" + OutputFileTag + ".root";
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
	
	TH1 *jet_size 		 = new TH1F("jet_size",     "Number of Jets", 10, 0, 10.0);
	TH1 *histJet_btag 	 = new TH1F("histJet_btag", "Number of B-tagged jets", 10, 0, 10.0);
    TH1 *histJet_btag_pt = new TH1F("histJet_btag_pt", "PT of B-tagged jets", 100, 0.0, 500.0);
    TH1 *histJet_btag_eta = new TH1F("histJet_btag_eta", "ETA of B-tagged jets", 100, -5.0, +5.0);
    TH1 *histJet_btag_pt_aftercut = new TH1F("histJet_btag_pt_aftercut", "PT of B-tagged jets after the cut", 100, 0.0, 500.0);
    TH1 *histJet_btag_eta_aftercut = new TH1F("histJet_btag_eta_aftercut", "ETA of B-tagged jets after the cut", 100, -5.0, +5.0);
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
	

    TH2 *hist_gen_lepton2D		= new TH2I("hist_gen_lepton2D"	, "Number of GEN Leptons ;GEN;SIM", 10, 0, 10.0,10,0,10);
    
    TH2 *hist_gen_elec_SIM1GEN0	= new TH2F("hist_gen_elec_SIM1GEN0"	, "Number of Elec : GEN=0  and SIM=1 ;PT;ETA", 100, 0, 100,100,-5,+5);
    TH2 *hist_gen_elec_SIM0GEN1	= new TH2F("hist_gen_elec_SIM0GEN1"	, "Number of Elec : GEN=1  and SIM=0 ;PT;ETA", 100, 0, 100,100,-5,+5);
    
    TH2 *hist_gen_muon_SIM1GEN0	= new TH2F("hist_gen_muon_SIM1GEN0"	, "Number of Muon : GEN=0  and SIM=1 ;PT;ETA", 100, 0, 100,100,-5,+5);
    TH2 *hist_gen_muon_SIM0GEN1	= new TH2F("hist_gen_muon_SIM0GEN1"	, "Number of Muon : GEN=1  and SIM=0 ;PT;ETA", 100, 0, 100,100,-5,+5);
    
    
    TH1 *hist_gen_elec_pt		= new TH1F("hist_gen_elec_pt"	, "Elec Pt " , 100, 0, 250);
    TH1 *hist_gen_elec_phi		= new TH1F("hist_gen_elec_phi"	, "Elec Phi ", 100, -5, 5);
    TH1 *hist_gen_elec_eta		= new TH1F("hist_gen_elec_eta"	, "Elec Eta ", 100, -5, 5);
 
    TH1 *hist_gen_muon_pt		= new TH1F("hist_gen_muon_pt"  , "Muon Pt " , 100, 0, 250);
    TH1 *hist_gen_muon_phi		= new TH1F("hist_gen_muon_phi" , "Muon Phi ", 100, -5, 5);
    TH1 *hist_gen_muon_eta		= new TH1F("hist_gen_muon_eta" , "Muon Eta ", 100, -5, 5);

    TH1 *histElec_pt            = new TH1F("elec_pt1"	 , "1st elec P_{T}", 100, 0.0, 500.0);
    TH1 *histElec_phi           = new TH1F("elec_pt1_phi", "1st elec Phi  ", 100, -5.0, 5.0);
    TH1 *histElec_eta           = new TH1F("elec_pt1_eta", "1st elec Eta  ", 100, -5.0, 5.0);
    
    TH1 *histMuon_pt        = new TH1F("muon_pt1"     , "1st mu P_{T}", 100, 0.0, 500.0);
    TH1 *histMuon_phi       = new TH1F("muon_pt1_phi" , "1st mu Phi  ", 100, -5.0, 5.0);
    TH1 *histMuon_eta       = new TH1F("muon_pt1_eta" , "1st mu Eta  ", 100, -5.0, 5.0);
    
    TH1 *histLepton_pt 	= new TH1F("lepton_pt"  , "lepton P_{T} ", 100, 0.0, 500.0);
    TH1 *histLepton_eta	= new TH1F("lepton_eta" , "lepton Eta   ", 100, -5.0, 5.0);
    TH1 *histLepton_phi	= new TH1F("lepton_phi" , "lepton  Phi  ", 100, -5.0, 5.0);
    
    
    TH1 *hist_gen_elec_deltaR	= new TH1F("hist_gen_elec_deltaR" , "DeltaR Elec Eta ", 1000, 0, 0.1);
    TH1 *hist_gen_muon_deltaR	= new TH1F("hist_gen_muon_deltaR" , "DeltaR Muon Eta ", 1000, 0, 0.1);

    TH1 *hist_gen_lepton             = new TH1F("hist_gen_lepton"            , "Number of GEN Leptons ;GEN;SIM", 10, 0, 10.0);
    TH1 *hist_lepton_before_trig     = new TH1F("numb_lepton_before_trig"    , "Number of SIM Leptons PT>20/30(elec/muon)", 10, 0, 10.0);
    TH1 *hist_lepton_pass_trig       = new TH1F("numb_lepton_pass_trig"      , "Number of SIM Leptons PT>20/30(elec/muon)", 10, 0, 10.0);
    TH1 *hist_lepton_before_hardphot = new TH1F("hist_lepton_before_hardphot", "Number of SIM Leptons |eta|<2.5 and PT>20/30(elec/muon)", 10, 0, 10.0);
	TH1 *hist_lepton_before_10gev      = new TH1F("numb_lepton_before_10gev" , "Number of SIM Leptons and PT>10(elec/muon)", 10, 0, 10.0);

	TH1 *lepton_invmass = new TH1F("lepton_invmass", "lepton inv mass", 100, 0, 500);
	
	TH2 *bjet_lepton_delta_eta 	= new TH2F("bjet_lepton_delta_eta"   , "Delta Eta between lepton and btagjet ", 100, 0, 250, 100, 0,   +5);
	TH2 *bjet_lepton_delta_phi 	= new TH2F("bjet_lepton_delta_phi"   , "Delta Phi between lepton and btagjet ", 100, 0, 250, 100, 0, +6.3);
	TH2 *bjet_lepton_deltaR     = new TH2F("bjet_lepton_deltaR"      , "Delta R between lepton and btagjet   ", 100, 0, 250, 100, 0,  +10);
    TH1 *hist_qlep_etajet       = new TH1F("hist_qlep_etajet"      , " lepton charge x jet eta  ", 100, -3, 3);
	TH2 *bjet_qeta				= new TH2F("bjet_qeta"      ,"Bjet Eta #cdot  lepton charge   ", 100, 0, 250, 100, 0,  +10);

    TH1 *hist_top_inv_mass_mhc80 	= new TH1F("top_inv_mass_mhc80"   , "top inv mass assume mhc=80  ", 100, 0, 500);
    TH1 *hist_top_inv_mass_mhc100 	= new TH1F("top_inv_mass_mhc100"   , "top inv mass assume mhc=100  ", 100, 0, 500);
    TH1 *hist_top_inv_mass_mhc130 	= new TH1F("top_inv_mass_mhc130"   , "top inv mass  assume mhc=130 ", 100, 0, 500);
    
    TH1 *hist_miss_pz_mhc80 	= new TH1F("miss_pz_mhc80"   , "miss_pz assume mhc=80   ", 100, 0, 500);
    TH1 *hist_miss_pz_mhc100 	= new TH1F("miss_pz_mhc100"  , "miss_pz assume mhc=100  ", 100, 0, 500);
    TH1 *hist_miss_pz_mhc130 	= new TH1F("miss_pz_mhc130"  , "miss_pz assume mhc=130  ", 100, 0, 500);
    
	//----------------------------------------------------------------------------------------

    TDirectory *metdir= outf->mkdir("met");
	metdir->cd();
	
	TH1 *histMET_et  = new TH1F("histMET_et" , "MET",  500,  0.0, 500.0);
	TH1 *histMET_eta = new TH1F("histMET_eta", "MET Eta", 100, -5.0, 5.0);
	TH1 *histMET_phi = new TH1F("histMET_phi", "histMET_phi Phi", 100, -5.0, 5.0);
	
	//----------------------------------------------------------------------------------------
	////////////////////////////////////////  My Data STructure  /////////////////////////////
	//----------------------------------------------------------------------------------------
	
	TH1F::SetDefaultSumw2(true);
	
	gSystem->Load("libExRootAnalysis");
	gSystem->Load("libDelphes");
	
	//string OutputFileName="skimming_delphes.root";
	//TFile *outf = new TFile(OutputFileName.c_str(),"RECREATE");
	
	//ofstream myfile;
	//myfile.open ("skimming_mass70_1.txt");
	
	// Create chain of root trees
	TChain chain("Delphes");
	//chain.Add(InputFileName.c_str());
	char filename[1000];
	FILE *input;
	input = fopen(InputFileName.c_str(),"r");
	if (input != NULL)
	{
		// lets read each line and get the filename from it
		while (fscanf(input,"%s\n",filename) != EOF)
		{
			printf("%s\n",filename);
			chain.Add(filename);
		}
	}
	
	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	
	// Get pointers to branches used in this analysis
	TClonesArray *branchJet 	= treeReader->UseBranch("Jet");
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon 	= treeReader->UseBranch("Muon");
	//TClonesArray *branchPhoton 	= treeReader->UseBranch("Photon");
	TClonesArray *branchMET 	= treeReader->UseBranch("MissingET");
	
	//----------------------------------------------------------------------------------------
	/////////////////////////////////////  LOOP  Over the EVENTS  ////////////////////////////
	//----------------------------------------------------------------------------------------
	
	
	cout << "Set Branch Addresses" << endl;
	
	int decade = 0;
	unsigned entries = 0;
	
	if (NENTRIES==-1) entries=numberOfEntries;
	else entries=NENTRIES;

	cout << "Reading TREE: " << numberOfEntries << " events available, \n";
	cout << "\t\t" << NENTRIES << " of them will be analyzed." << endl;
	
    GenParticle *part, *mother1, *this_elec_part=0, *this_muon_part=0;//, *status;
    Jet *jet[10];
    Electron *elec[4]; //, *elec2, *elec3;
    Muon *mu[4]; //, *mu2, *mu3;
    MissingET *met=0;
    GenParticle *partm, *partp, *motherm, *motherp, *mother2, *mother3;
	
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
	
	//int part_counter=0;
	int signal_counter=0;
	int nevents=0;
    int counter_saved=0;
	
	// Loop over all events
	for(Long64_t entry = 0; entry < entries; ++entry)
	{
        nevents++;
		double progress = 10.0*entry/(1.0*entries);
		int k = TMath::FloorNint(progress);
		if (k > decade) {   cout << 10*k << "%\t"; gSystem->Exec("date '+%H:%M:%S'"); cout << endl;	}
		decade = k;
		
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		
		int part_counter_elec=0;
		int part_counter_muon=0;
		int genevent_size=0;
		
		genevent_size=branchParticle->GetEntries();
		
		bool the_signal=false;
		
        // counter for the number of events in the DATA
        // which has a tau decays leptonically
	
		if( genevent_size > 0)
		{
			for(int i=0; i<genevent_size; i++)
			{
				part = (GenParticle *) branchParticle->At(i);
				
				int motherPID=0;
				if (part->M1>1)
				{
					mother1 = (GenParticle *) branchParticle->At(part->M1);
					motherPID = mother1->PID;
				}
				else motherPID=0;
				
				//lepton status=1
				if ( (abs(part->PID)==11 && (part->Status)==1 && abs(motherPID)==15 ) ||
					(abs(part->PID)==13 && (part->Status)==1 && abs(motherPID)==15 ) )
				{
					partm   = (GenParticle *) branchParticle->At(i-1);
					partp   = (GenParticle *) branchParticle->At(i+1);
					motherm = (GenParticle *) branchParticle->At(partm->M1);
					mother2 = (GenParticle *) branchParticle->At(motherm->M1);
					mother3 = (GenParticle *) branchParticle->At(mother2->M1);
					motherp = (GenParticle *) branchParticle->At(partp->M1);
					
					if ( (partm->M1)==(partp->M1)
						&& abs(motherm->PID)==15
						&& abs(motherp->PID)==15
						&& abs(partm->PID)==16
						&& abs(partp->PID)==(abs(part->PID)+1)
						
						&& abs(motherm->PID)==15
						&& abs(mother2->PID)==15
						&& ( abs(mother3->PID)==15 || abs(mother3->PID)==37)
					) the_signal=true;
				}
				

                /*
				 if( (abs(part->PID)== 11 && (part->Status)==1 && abs(motherPID)==15 ) ||
				 (abs(part->PID)== 13 && (part->Status)==1 && abs(motherPID)==15 ) )
				 {
					the_signal=true;
					lept_counter++;
        		    cout.width(4);
					cout << i << "\t";
					cout.width(10);
                     cout << part->PID<< "\t";
                     cout.width(10);
                     cout << part->M1 << "\t";
                     cout.width(10);
                     cout << part->Status << "\t" << endl;
                     }
         			*/
				
				// count the number of leptons in the event
				if(   abs(part->PID)==11 && (part->Status)== 1
				   && abs(motherPID)==15 &&  part->PT > ELECPT_CUT
				   //&& abs(part->Eta) < 2.5
				   )
				{
					part_counter_elec++;
					hist_gen_elec_pt->Fill(part->PT);
					hist_gen_elec_eta->Fill(part->Eta);
					hist_gen_elec_phi->Fill(part->Phi);
					this_elec_part=part;
				}
				
				if(   abs(part->PID)==13 && (part->Status)== 1
				   && abs(motherPID)==15 &&  part->PT > MUONPT_CUT
				   //&& abs(part->Eta) < 2.5
				   )
				{
					part_counter_muon++;
					hist_gen_muon_pt ->Fill(part->PT);
					hist_gen_muon_eta->Fill(part->Eta);
					hist_gen_muon_phi->Fill(part->Phi);
					this_muon_part=part;
				}
			}
		}

        if ( the_signal == false && DATA == true ) continue;
        // if this an event where
        // higgs_charged > tau nu_tau
        // tau > leptonic decay
		
		if ( the_signal )
        {
            hist_gen_lepton->Fill(part_counter_muon+part_counter_elec);
            signal_counter++;
        }

		/*
		 if ( mysignal==true )
		 //if ( this_is_the_signal == true )
		 {
			if ( count_print < 10 )
			{
		 count_print++;
		 cout << "EVENT NUMBER   :   " << entry << endl;
		 for(int i=0; i<genevent_size; i++)
		 {
		 part = (GenParticle *) branchParticle->At(i);
		 
		 if (abs(part->PID)>10 && abs(part->PID)<17)
		 {
            cout.width(4);
            cout << i << "\t";
            cout.width(10);
            cout << part->PID<< "\t";
            cout.width(10);
            cout << part->M1 << "\t";
            cout.width(10);
            cout << part->Status << "\t" << endl;
            }
		 }
		 cout << endl;
		 }
		 }
		 */
    
        ////////////////////////////////////////////////////////////////////////////////////
		// 0. TRIGGER EMULATION
        ////////////////////////////////////////////////////////////////////////////////////

        Electron * this_elec=0;
        Muon     * this_muon=0;
		
		//int numb_muon=0;
		
		//filter lepton=1 events
		int numb_elec_pass_cuts=0;
		int numb_muon_pass_cuts=0;
		
        if( branchElectron->GetEntries() > 0)
        {
            for(int i=0; i<branchElectron->GetEntries(); i++)
            {
                elec[i] = (Electron *) branchElectron->At(i);
                if(  (elec[i]->PT)> 25 )
                {
                    numb_elec_pass_cuts++;
                }
            }
        }

        if( branchMuon->GetEntries() > 0)
        {
            for(int i=0; i<branchMuon->GetEntries(); i++)
            {
                mu[i] = (Muon *) branchMuon->At(i);
                if(  (mu[i]->PT) > 20 )
                {
                    numb_muon_pass_cuts++;
                }
            }
        }

        // comparison for generated Leptons and triggered leptons
        if (the_signal)
        {
            hist_gen_lepton2D->Fill(part_counter_elec,  numb_elec_pass_cuts);
            hist_gen_lepton2D->Fill(part_counter_muon+3,numb_muon_pass_cuts+3);
        }

        // number of lepton before the trigger
		hist_lepton_before_trig->Fill(numb_muon_pass_cuts + numb_elec_pass_cuts);
 
        // if the total number of lepton is higher than 0 fire for this event
        if ((numb_muon_pass_cuts + numb_elec_pass_cuts)==0) continue;
        event_counter0_after_trigger++;

        // number of leptons after the trigger cut
        hist_lepton_pass_trig->Fill(numb_muon_pass_cuts + numb_elec_pass_cuts);

        ////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////

        
        
        ////////////////////////////////////////////////////////////////////////////////////
        // 1a. HARD LEPTON IN Central REGION
        ////////////////////////////////////////////////////////////////////////////////////
        
        
        // will reuse these constants so set them to 0
        numb_muon_pass_cuts=0;
        numb_elec_pass_cuts=0;

        if( branchElectron->GetEntries() > 0)
        {
            for(int i=0; i<branchElectron->GetEntries(); i++)
            {
                elec[i] = (Electron *) branchElectron->At(i);
                if(  (elec[i]->PT)>ELECPT_CUT && abs(elec[i]->Eta)<2.5 )
                {
                    numb_elec_pass_cuts++;
                    this_elec=elec[i];
                }
            }
        }
        if( branchMuon->GetEntries() > 0)
        {
            for(int i=0; i<branchMuon->GetEntries(); i++)
            {
                mu[i] = (Muon *) branchMuon->At(i);
                if(  (mu[i]->PT)>MUONPT_CUT && abs(mu[i]->Eta)<2.5 )
                {
                    numb_muon_pass_cuts++;
                    this_muon=mu[i];
                }
            }
        }

        // GEN vs SIM comparison shopuld be done over the SIGNAL ONLY
        if (the_signal)
        {
            // Gen=1 and SIM=0
            if (part_counter_elec==1 && numb_elec_pass_cuts==0)	hist_gen_elec_SIM0GEN1->Fill(this_elec_part->PT,this_elec_part->Eta);
            if (part_counter_muon==1 && numb_muon_pass_cuts==0)	hist_gen_muon_SIM0GEN1->Fill(this_muon_part->PT,this_muon_part->Eta);
        
            // Gen=0 and SIM=1
            if (part_counter_elec==0 && numb_elec_pass_cuts==1)	hist_gen_elec_SIM1GEN0->Fill(this_elec->PT,this_elec->Eta);
            if (part_counter_muon==0 && numb_muon_pass_cuts==1)	hist_gen_muon_SIM1GEN0->Fill(this_muon->PT,this_muon->Eta);
       
            if (part_counter_elec==1 && numb_elec_pass_cuts==1)
            {
                TLorentzVector v1;
                TLorentzVector v2;
                v1.SetPxPyPzE(this_elec_part->Px,this_elec_part->Py,this_elec_part->Pz,this_elec_part->E);
                v2.SetPtEtaPhiM(this_elec->PT,this_elec->Eta,this_elec->Phi,0.000511);
                hist_gen_elec_deltaR->Fill( ROOT::Math::VectorUtil::DeltaR(v1,v2) );
                //double deltar= ROOT::Math::VectorUtil::DeltaR(v1,v2);
                //cout << " deltaR   : " << deltar << endl;
            }
        
            if (part_counter_muon==1  && numb_muon_pass_cuts==1)
            {
                TLorentzVector v1;
                TLorentzVector v2;
                v1.SetPxPyPzE(this_muon_part->Px,this_muon_part->Py,this_muon_part->Pz,this_muon_part->E);
                v2.SetPtEtaPhiM(this_muon->PT,this_muon->Eta,this_muon->Phi,0.100);
                hist_gen_muon_deltaR->Fill(ROOT::Math::VectorUtil::DeltaR(v1,v2));
            }
        }
        
        hist_lepton_before_hardphot->Fill((numb_muon_pass_cuts + numb_elec_pass_cuts));

        
        // if the total number of lepton is higher than 0 fire for this event
        if ((numb_muon_pass_cuts + numb_elec_pass_cuts)!=1) continue;
        event_counter1_after_hardphoton++;


        ////////////////////////////////////////////////////////////////////////////////////
        // 1b. SOFT LEPTON IN whole REGION
        // ELIMINATE EVENTS with SOFT LEPTONS (PT > 10 GeV)
        // if there is any soft lepton in the event in whole
        // eta region this event should be discarded
        ////////////////////////////////////////////////////////////////////////////////////
        
        // will reuse these constants so set them to 0
        numb_muon_pass_cuts=0;
        numb_elec_pass_cuts=0;

        
        if( branchElectron->GetEntries() > 0)
        {
            for(int i=0; i<branchElectron->GetEntries(); i++)
            {
                elec[i] = (Electron *) branchElectron->At(i);
                if(  (elec[i]->PT)>10 ) numb_elec_pass_cuts++;
            }
        }
        
        if( branchMuon->GetEntries() > 0)
        {
            for(int i=0; i<branchMuon->GetEntries(); i++)
            {
                mu[i] = (Muon *) branchMuon->At(i);
                if(  (mu[i]->PT)>10) numb_muon_pass_cuts++;
            }
        }
        
        hist_lepton_before_10gev->Fill(numb_muon_pass_cuts + numb_elec_pass_cuts);

        // filter the event if there is more than one lepton
        if ((numb_muon_pass_cuts + numb_elec_pass_cuts)!=1) continue;
        
        // number of leptons after the trigger cut
        event_counter2_after_10gev++;

		
        if(numb_elec_pass_cuts==1)
        {
            histElec_pt->Fill(this_elec->PT);
            histElec_eta->Fill(this_elec->Eta);
            histElec_phi->Fill(this_elec->Phi);
            histLepton_pt->Fill(this_elec->PT);
            histLepton_eta->Fill(this_elec->Eta);
            histLepton_phi->Fill(this_elec->Phi);
        }
        else
            if(numb_muon_pass_cuts ==1 )
            {
                histMuon_pt->Fill(this_muon->PT);
                histMuon_eta->Fill(this_muon->Eta);
                histMuon_phi->Fill(this_muon->Phi);
                histLepton_pt->Fill(this_muon->PT);
                histLepton_eta->Fill(this_muon->Eta);
                histLepton_phi->Fill(this_muon->Phi);
            }

        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        // KINEMATICAL CUT on Leptons PT
        // filter the event if there is lepton PT>55
        //////////////////////////////////////////////////////////
		
        if (numb_elec_pass_cuts==1 && this_elec->PT  < 50 && this_elec->PT > 30 ) continue;
        else
        if (numb_muon_pass_cuts==1 && this_muon->PT  < 50 && this_muon->PT > 30 ) continue;

        event_counter3_after_leptonpt55++;
		
        ////////////////////////////////////////////////////////////////////
        //  MET cut, if MET < 60 gev then discard this event
        ////////////////////////////////////////////////////////////////////
		
		double this_met=0;
		if(branchMET->GetEntries() > 0)
		{
			met = (MissingET *) branchMET->At(0);
			histMET_et->Fill(met->MET);
			histMET_eta->Fill(met->Eta);
			histMET_phi->Fill(met->Phi);
			this_met=met->MET;
		}
        
        if (this_met < MET_CUT ) continue;
        event_counter4_after_met++;

		////////////////////////////////////////////////////////////////////
		//	JET BRANCH
		////////////////////////////////////////////////////////////////////
		
        double  alpha_t=0;
        if(branchJet->GetEntries() > 1 )
        {
            jet[0] = (Jet*) branchJet->At(0);
            jet[1] = (Jet*) branchJet->At(1);
            
            double jet_inv_mass12=sqrt(2*jet[0]->PT*jet[1]->PT*(cosh(jet[0]->Eta-jet[1]->Eta)-cos(jet[0]->Phi-jet[1]->Phi)));
            alpha_t=(double)(jet[1]->PT)/jet_inv_mass12;
            hist_alpha_t->Fill(alpha_t);
        }
        
        // go to the end of the loop for alpha_t>0.5,
        // because we assume these are qcd events
        //if( alpha_t < 0.5 ) continue;
        event_counter4_after_alpha_t++;

        Jet *this_bjet=0;
        Jet *leadin_jet=0;
		int counter_btag=0;
		// filter bjet tagged events
		if(branchJet->GetEntries() > 0)
		{
            leadin_jet=(Jet*) branchJet->At(0);
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				jet[i] = (Jet*) branchJet->At(i);
                
				if (jet[i]->PT > 30 && jet[i]->BTag==1)
				{
                    histJet_btag_pt->Fill(jet[i]->PT);
                    histJet_btag_eta->Fill(jet[i]->Eta);

                    counter_btag++;
					this_bjet=jet[i];
				}
			}
			histJet_btag->Fill(counter_btag);
		}

		if ( counter_btag != 1 ) continue;
		event_counter5_after_btag++;

		histJet_btag_pt_aftercut->Fill(this_bjet->PT);
        histJet_btag_eta_aftercut->Fill(this_bjet->Eta);

        
        
        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        counter_saved++;
        if ( SaveData && (double)cross*200*efficiency < counter_saved ) break;
        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////

 
        
        double deltaR2=999.;
        double deltaEta=999.;
        if (numb_elec_pass_cuts==1)
		{
			deltaR2=pow(abs( (this_bjet->Eta)-(this_elec->Eta) ),2)+pow(abs( (this_bjet->Phi)-(this_elec->Phi) ),2);
            deltaEta=abs( (this_bjet->Eta)-(this_elec->Eta) );
			bjet_lepton_delta_eta->Fill( this_elec->PT, deltaEta );
			bjet_lepton_delta_phi->Fill( this_elec->PT, abs( (this_bjet->Phi)-(this_elec->Phi) ));
			bjet_lepton_deltaR->Fill( this_elec->PT, sqrt(deltaR2) );
		}
		else if (numb_muon_pass_cuts==1)
		{
			deltaR2=pow( (this_bjet->Eta)-(this_muon->Eta) ,2)+pow( (this_bjet->Phi)-(this_muon->Phi) ,2);
            deltaEta=abs( (this_bjet->Eta)-(this_muon->Eta) );
            bjet_lepton_delta_eta->Fill( this_muon->PT, deltaEta );
			bjet_lepton_delta_phi->Fill( this_muon->PT, abs( (this_bjet->Phi)-(this_muon->Phi) ));
			bjet_lepton_deltaR->Fill( this_muon->PT, sqrt(deltaR2) );
		}

        if ( numb_muon_pass_cuts == 1 )
        {
            hist_qlep_etajet->Fill(this_muon->Charge * leadin_jet->Eta);
            bjet_qeta->Fill(this_bjet->PT, this_muon->Charge * this_bjet->Eta);
        }
        else if ( numb_elec_pass_cuts == 1 )
        {
            hist_qlep_etajet->Fill(this_elec->Charge * leadin_jet->Eta);
            bjet_qeta->Fill(this_bjet->PT, this_elec->Charge * this_bjet->Eta);

        }
        
        if ( deltaEta > 1.5 ) continue;
        event_counter5a_after_deltaEta++;
        
        //////////////////////////////////////////////////////////////////////////////////////
        // TOP INVARIANT MASS
        //////////////////////////////////////////////////////////////////////////////////////
        
        double top_inv_mass_mhc80=0;
        double top_inv_mass_mhc100=0;
        double top_inv_mass_mhc130=0;
        
        double miss_pz_mhc80=0;
        double miss_pz_mhc100=0;
        double miss_pz_mhc130=0;

        if ( numb_muon_pass_cuts == 1 )
        {
            miss_pz_mhc80=missing_energy_pz(this_muon,  met, 6400);
            miss_pz_mhc100=missing_energy_pz(this_muon, met, 10000);
            miss_pz_mhc130=missing_energy_pz(this_muon, met, 16900);
            
            top_inv_mass_mhc80=top_invariant_mass(this_bjet,  this_muon, met,  6400);
            top_inv_mass_mhc100=top_invariant_mass(this_bjet, this_muon, met, 10000);
            top_inv_mass_mhc130=top_invariant_mass(this_bjet, this_muon, met, 16900);
            
            hist_miss_pz_mhc80->Fill(miss_pz_mhc80);
            hist_miss_pz_mhc100->Fill(miss_pz_mhc100);
            hist_miss_pz_mhc130->Fill(miss_pz_mhc130);
            
            hist_top_inv_mass_mhc80->Fill(top_inv_mass_mhc80);
            hist_top_inv_mass_mhc100->Fill(top_inv_mass_mhc100);
            hist_top_inv_mass_mhc130->Fill(top_inv_mass_mhc130);
            //top_invariant_massv2(this_muon, met,  6400);
        }
        else
        if ( numb_elec_pass_cuts == 1 )
        {
                miss_pz_mhc80=missing_energy_pz(this_elec,  met, 6400);
                miss_pz_mhc100=missing_energy_pz(this_elec, met, 10000);
                miss_pz_mhc130=missing_energy_pz(this_elec, met, 16900);
                
                top_inv_mass_mhc80=top_invariant_mass(this_bjet, this_elec,  met, 6400);
                top_inv_mass_mhc100=top_invariant_mass(this_bjet, this_elec, met, 10000);
                top_inv_mass_mhc130=top_invariant_mass(this_bjet, this_elec, met, 16900);
                
                hist_miss_pz_mhc80->Fill(miss_pz_mhc80);
                hist_miss_pz_mhc100->Fill(miss_pz_mhc100);
                hist_miss_pz_mhc130->Fill(miss_pz_mhc130);
                
                hist_top_inv_mass_mhc80->Fill(top_inv_mass_mhc80);
                hist_top_inv_mass_mhc100->Fill(top_inv_mass_mhc100);
                hist_top_inv_mass_mhc130->Fill(top_inv_mass_mhc130);
                //top_invariant_massv2(this_elec, met, 6400);
                
            }

        double leptoninvariantmass=0;
		if ( numb_muon_pass_cuts == 1 )
		{
			leptoninvariantmass=lepton_invariant_mass(this_muon, met);
			lepton_invmass->Fill(leptoninvariantmass);
		}
		else
        if ( numb_elec_pass_cuts == 1 )
        {
			leptoninvariantmass=lepton_invariant_mass(this_elec, met);
			lepton_invmass->Fill(leptoninvariantmass);
        }
		
		if ( top_inv_mass_mhc100 < 220  ) continue;
		event_counter6_after_topinvmass++;
		
		if ( leptoninvariantmass > 50  ) continue;
		event_counter7_after_leptinvmass++;
		
		
		////////////////////////////////////////////////////////////////////
		//	JET BRANCH
		////////////////////////////////////////////////////////////////////
		int numb_hardjet=0;
        
        // eliminate hard jets
        Jet *this_hardjet=0;
		if(branchJet->GetEntries() > 0)
		{
			jet_size->Fill( branchJet->GetEntries() );
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				jet[i] = (Jet*) branchJet->At(i);
                if ( (jet[i]->PT)>35 && abs(jet[i]->Eta)< 4.9 )
                {
                    this_hardjet=jet[i];
                    numb_hardjet++;
                }

                if (i>2)continue;
                histJet_pt[i]->Fill(jet[i]->PT);
                histJet_eta[i]->Fill(jet[i]->Eta);
                histJet_phi[i]->Fill(jet[i]->Phi);
            }
		}
        jet_size_cut8->Fill(numb_hardjet);

        if (numb_hardjet !=1) continue;
        //if (numb_hardjet == 0) continue;
		event_counter8_after_onejet++;

        // eliminate soft jets
        Jet *this_softjet=0;
        int numb_softjet=0;
        
		if( branchJet->GetEntries() > 0 )
		{
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				jet[i] = (Jet*) branchJet->At(i);
                if ( (jet[i]->PT)>15 && (jet[i]->PT)<35 && abs(jet[i]->Eta)< 4.9 )
                {
                    numb_softjet++;
                    this_softjet=jet[i];
                }
			}
		}
		
		jet_size_cut9->Fill(numb_softjet);

        if (numb_softjet!=1) continue;
		event_counter9_after_onejet++;

        // eliminate jets with central pseudorapidity distributions
		hist_before_jet_eta->Fill(this_hardjet->Eta);
		if ( abs(this_hardjet->Eta) < 2.5 ) continue;
		event_counter10_after_jeteta++;

	}	//end of event loop

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	//---------------------------------------------------------------------------------------------------------
	////////////////////////////////  Saving the histograms into output  //////////////////////////////////////////
	//---------------------------------------------------------------------------------------------------------
	// end of event loop
	
	cout << "nevents : " << nevents << endl;
	cout << "signal_counter :  " << signal_counter << "  %  " << (double)signal_counter/nevents << endl;
	cout << "_____________________________________________________________________" << endl;
    if (DATA) nevents=signal_counter;
	cout << " 0  event after trigger     : "    << event_counter0_after_trigger     << "\t" << (double)event_counter0_after_trigger/nevents     << endl;
	cout << " 1  event after hard photon : "    << event_counter1_after_hardphoton  << "\t" << (double)event_counter1_after_hardphoton/nevents  << endl;
	cout << " 2  event after 10gev       : "    << event_counter2_after_10gev       << "\t" << (double)event_counter2_after_10gev/nevents       << endl;
	cout << " 3  event after Lep(PT)>55  : "    << event_counter3_after_leptonpt55  << "\t" << (double)event_counter3_after_leptonpt55/nevents  << endl;
	cout << " 4  event after met > 50    : "    << event_counter4_after_met         << "\t" << (double)event_counter4_after_met/nevents         << endl;
	cout << " 4  event after alp_t > 0.5 : "    << event_counter4_after_alpha_t     << "\t" << (double)event_counter4_after_alpha_t/nevents     << endl;
    cout << " 5  event after btag = 1    : "    << event_counter5_after_btag        << "\t" << (double)event_counter5_after_btag/nevents        << endl;
    cout << " 5a event after delta Eta   : "    << event_counter5a_after_deltaEta   << "\t" << (double)event_counter5a_after_deltaEta/nevents   << endl;
	cout << " 6  event after topinvmass  : "    << event_counter6_after_topinvmass  << "\t" << (double)event_counter6_after_topinvmass/nevents  << endl;
	cout << " 7  event after leptinvmass : "    << event_counter7_after_leptinvmass << "\t" << (double)event_counter7_after_leptinvmass/nevents << endl;
	cout << " 8  event after numb_jet =1 : "    << event_counter8_after_onejet      << "\t" << (double)event_counter8_after_onejet/nevents      << endl;
	cout << " 9  event after numb_jet =1 : "    << event_counter9_after_onejet      << "\t" << (double)event_counter9_after_onejet/nevents      << endl;
	cout << " 10 event after jet_eta>2.5 : "    << event_counter10_after_jeteta     << "\t" << (double)event_counter10_after_jeteta/nevents     << endl;
	cout << "_____________________________________________________________________" << endl;
    cout << " number of expct events L=1fb^-1  : "    << "\t" << (double)event_counter10_after_jeteta/nevents*cross*1000 ;
    cout << "\t" << InputFileName.c_str() << endl;

	//myfile.close();
	outf->Write();
	outf->Close();
	
	//end of main loop
}


