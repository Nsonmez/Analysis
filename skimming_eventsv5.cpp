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
#include "Math/GenVector/VectorUtil.h"

#include "skimming_events.h"

using namespace std;

//----------------------------------------------------------------------------------------
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

	double MUONPT_CUT=23;
	double ELECPT_CUT=33;
	
	
	cout << "________________________________________________________________\n";
	cout << "This is a SIGNAL(true)/BACKGROUND(false) : " << DATA << endl;
	
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
	
	
	//TH1 *histJet_pt[10];
	//TH1 *histJet_eta[10];
	//TH1 *histJet_phi[10];
	char hist_name[100];
	
	//----------------------------------------------------------------------------------------
	TDirectory *histo= outf->mkdir("histograms");
	histo->cd();
	
	TH1 *hist_gen_lepton	= new TH1F("gen_lepton"	, "Number of GEN Leptons ;GEN;SIM", 10, 0, 10.0);
	
	TH1 *hist_gen_elec_pt	= new TH1F("hist_gen_elec_pt"	, "Elec Pt " , 100, 0, 250);
	TH1 *hist_gen_elec_eta	= new TH1F("hist_gen_elec_eta"	, "Elec Eta ", 100, -5, 5);
	
	TH1 *hist_gen_muon_pt	= new TH1F("hist_gen_muon_pt"  , "Muon Pt " , 100, 0, 250);
	TH1 *hist_gen_muon_eta	= new TH1F("hist_gen_muon_eta" , "Muon Eta ", 100, -5, 5);

    TH1 *hist_lept_inv_massv2	= new TH1F("hist_lept_inv_massv2" , "Lepton Inv Mass ", 100, 0, 500);

	TH1 *hist_alpha_t		= new TH1F("hist_alpha_t" , "Alpha_t ", 100, 0, 2);

	TH1F::SetDefaultSumw2(true);
	
	//----------------------------------------------------------------------------------------
	////////////////////////////////////////  My Data STructure  /////////////////////////////
	//----------------------------------------------------------------------------------------
	

	TTree *outputTree = new TTree("DATA","a tree of event info and jet info");
	outputTree->SetDirectory(outf);
	
	// new class for all my leaves
	MyEVENT my_event;

	outputTree->Branch("TreeS", &my_event,"Pt_lepton/D:MET:Pt_btag:Mt_lv:Mi_lvb:N_jets:R_lb:Phi_lb:Eta_lb:alpha_t");
	
	gSystem->Load("libExRootAnalysis");
	gSystem->Load("libDelphes");
	
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
	TClonesArray *branchJet 	 = treeReader->UseBranch("Jet");
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon 	 = treeReader->UseBranch("Muon");
	//TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchMET 	 = treeReader->UseBranch("MissingET");
	
	//----------------------------------------------------------------------------------------
	/////////////////////////////////////  LOOP  Over the EVENTS  ////////////////////////////
	//----------------------------------------------------------------------------------------
	
	
	cout << "Set Branch Addresses" << endl;
	
	int decade = 0;
	unsigned entries = 0;
	
	if (NENTRIES==-1) entries=numberOfEntries;
	else entries=NENTRIES;
	
	outf->cd();
	
	cout << "Reading TREE: " << numberOfEntries << " events available, \n";
	cout << "\t\t" << NENTRIES << " of them will be analyzed." << endl;
	
	GenParticle *part, *mother1, *this_elec_part=0, *this_muon_part=0;//, *status;
	Jet *jet[10];
	Electron *elec[4]; //, *elec2, *elec3;
	Muon *mu[4]; //, *mu2, *mu3;
	MissingET *met=0;
	GenParticle *partm, *partp, *motherm, *motherp, *mother2, *mother3;
	
	int event_counter1_after_trigger    =0;
	int event_counter2_after_10gev      =0;
	int event_counter4_after_alpha_t    =0;
	int event_counter5_after_btag       =0;
	
    /*
	int event_counter3_after_leptonpt55 =0;
	int event_counter6_after_topinvmass =0;
	int event_counter7_after_leptinvmass=0;
	int event_counter8_after_onejet     =0;
	int event_counter9_after_onejet     =0;
	int event_counter10_after_jeteta    =0;
	*/
    
	int part_counter=0;
	int signal_counter=0;
	int counter_saved=0;
	
	// Loop over all events
	for(Long64_t entry = 0; entry < entries; ++entry)
	{
		
		double progress = 10.0*entry/(1.0*entries);
		int k = TMath::FloorNint(progress);
		if (k > decade) {   cout << 10*k << "%\t"; gSystem->Exec("date '+%H:%M:%S'"); cout << endl;	}
		decade = k;
		
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		

        ////////////////////////////////////////////////////////////////////
        // RUN OVER THE SIGNAL DATA ONLY
        ////////////////////////////////////////////////////////////////////
        
        int part_counter_elec=0;
		int part_counter_muon=0;
		int genevent_size=0;
		
		genevent_size=branchParticle->GetEntries();
		
		bool the_signal=false;
		//bool mysignal  =false;
		
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

				if( abs(part->PID)== 11 && (part->Status)== 1
				   && abs(motherPID)==15 && part->PT > ELECPT_CUT
				   && abs(part->Eta) < 2.5 )
				{
					part_counter_elec++;
					hist_gen_elec_pt->Fill(part->PT);
					hist_gen_elec_eta->Fill(part->Eta);
					this_elec_part=part;
				}
				
				if( abs(part->PID)== 13 && (part->Status)== 1
				   && abs(motherPID)==15 && part->PT > MUONPT_CUT
				   && abs(part->Eta) < 2.5 )
				{
					part_counter_muon++;
					hist_gen_muon_pt ->Fill(part->PT);
					hist_gen_muon_eta->Fill(part->Eta);
					this_muon_part=part;
				}
			}
			hist_gen_lepton->Fill(part_counter);
			hist_gen_lepton->Fill(part_counter_muon+part_counter_elec+5);
		}
		
		if ( the_signal == false && DATA==true ) continue;
		
        // if this an event where
        // higgs_charged > tau nu_tau
        // tau > leptonic decay
        
        // counter for the number of events in the DATA
        // which has a tau decays leptonically
        if ( the_signal ) signal_counter++;


        
        
        
        
        
   
        //////////////////////////////////////////
        // TRIGGER EMULATION
        //////////////////////////////////////////
        
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

        // if the total number of muons is higher than 0 fire for this event
        if ( (numb_muon_pass_cuts + numb_elec_pass_cuts)==0 ) continue;
        event_counter1_after_trigger++;
        
        // will reuse these constants so set them to 0
        numb_muon_pass_cuts=0;
        numb_elec_pass_cuts=0;
        
        //////////////////////////////////////////
        // NUMBER OF LEPTONS PT > 10 GeV
        //////////////////////////////////////////
        
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

        //////////////////////////////////////////
        // FILTER THE EVENTS FOR LEPTON=1
        // filter the event if there is more than one lepton
        
        if ((numb_muon_pass_cuts + numb_elec_pass_cuts)!=1) continue;
        
        // number of leptons after the trigger cut
        event_counter2_after_10gev++;

        
        //////////////////////////////////////////
        // MY PRESELECTION CUTS
        //////////////////////////////////////////

		
		int numb_jet=0;
		int counter_btag=0;
		//int numb_muon=0;
		
		double this_met=0;
		if(branchMET->GetEntries() > 0)
		{
			met = (MissingET *) branchMET->At(0);
			this_met=met->MET;
			
		}

		//if (this_met < 50. ) continue;
		//event_counter4_after_met++;
		
		
		////////////////////////////////////////////////////////////////////
		//	JET BRANCH
		////////////////////////////////////////////////////////////////////
		Jet *this_bjet;

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
		
		int counter_btagfilter=0;
		// filter bjet tagged events
		if(branchJet->GetEntries() > 0)
		{
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				jet[i] = (Jet*) branchJet->At(i);
				if (jet[i]->BTag==1)
				{
					counter_btag++;
					if (jet[i]->PT < 100)
					{
						this_bjet=jet[i];
						counter_btagfilter++;
					}
				}
			}
		}

		if ( counter_btagfilter!=1 ) continue;
 		event_counter5_after_btag++;
		
		
		///////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////
		counter_saved++;
		
		if ( SaveData && (double)cross*20000*efficiency < counter_saved ) break;
		
		
		
	///////////////////////////////////////////////////////////////////////////
	///////////////////////// end of PRESELCTION
	///////////////////////////////////////////////////////////////////////////        
		my_event.MET = this_met;
		my_event.alpha_t=alpha_t;

        if(numb_elec_pass_cuts==1)	my_event.Pt_lepton=this_elec->PT;
		else if(numb_muon_pass_cuts ==1 ) my_event.Pt_lepton=this_muon->PT;


		if (numb_elec_pass_cuts==1)
		{
			double deltaR2=pow(abs( (this_bjet->Eta)-(this_elec->Eta) ),2)+pow(abs( (this_bjet->Phi)-(this_elec->Phi) ),2);
// 			bjet_lepton_delta_eta->Fill( this_elec->PT, abs( (this_jet->Eta)-(this_elec->Eta) ) );
// 			bjet_lepton_delta_phi->Fill( this_elec->PT, abs( (this_jet->Phi)-(this_elec->Phi) ));
// 			bjet_lepton_deltaR->Fill( this_elec->PT, sqrt(deltaR2) );
			my_event.Eta_lb=abs( (this_bjet->Eta)-(this_elec->Eta) );
			my_event.Phi_lb=abs( (this_bjet->Phi)-(this_elec->Phi) );
			my_event.R_lb=sqrt(deltaR2);
		}
		else if (numb_muon_pass_cuts==1)
		{
			double deltaR2=pow( (this_bjet->Eta)-(this_muon->Eta) ,2)+pow( (this_bjet->Phi)-(this_muon->Phi) ,2);
// 			bjet_lepton_delta_eta->Fill( this_muon->PT, abs( (this_jet->Eta)-(this_muon->Eta) ) );
// 			bjet_lepton_delta_phi->Fill( this_muon->PT, abs( (this_jet->Phi)-(this_muon->Phi) ));
// 			bjet_lepton_deltaR->Fill( this_muon->PT, sqrt(deltaR2) );
			my_event.Eta_lb=abs( (this_bjet->Eta)-(this_muon->Eta) );
			my_event.Phi_lb=abs( (this_bjet->Phi)-(this_muon->Phi) );
			my_event.R_lb=sqrt(deltaR2);
		}
		
		// TODO make a comparison plot for batg-jet
		// multiplicity and pt distribution before and after
		//Phi_lb=abs((this_jet->Phi)-(this_muon->Phi));
		//Teta_lb=abs( (this_jet->Eta)-(this_muon->Eta) );
		
		//////////////////////////////////////////////////////////////////////////////////////
		// TOP INVARIANT MASS and LEPTONIC TRANSVERSE MASS
		//////////////////////////////////////////////////////////////////////////////////////
		double topinvariantmass=0;
		double leptoninvariantmass=0;
		
		if ( numb_muon_pass_cuts == 1 )	// there is a elec in the event
		{
			// cout << "e  "<< top_invariant_mass(elec1, this_met);
			topinvariantmass=top_invariant_mass(this_bjet,this_muon, met,MH2);
            // top_invmass->Fill(topinvariantmass);
			my_event.Mi_lvb=topinvariantmass;

            // cout << "e  "<< lepton_invariant_mass(elec1, this_met);
            leptoninvariantmass=lepton_invariant_mass(this_muon, met);
            my_event.Mt_lv=leptoninvariantmass;
            
            hist_lept_inv_massv2->Fill( lepton_invariant_massv2(this_bjet,this_muon, met,MH2) );
		}
		else
		if ( numb_elec_pass_cuts == 1 )	// there is a muon in the event
		{
            // cout << "m  "<< top_invariant_mass(mu1, this_met);
			topinvariantmass=top_invariant_mass(this_bjet,this_elec, met,MH2);
            // top_invmass->Fill(topinvariantmass);
            my_event.Mi_lvb=topinvariantmass;
            
            // cout << "m  "<< lepton_invariant_mass(mu1, this_met);
            leptoninvariantmass=lepton_invariant_mass(this_elec, met);
            // lepton_invmass->Fill(leptoninvariantmass);
            my_event.Mt_lv=leptoninvariantmass;

            hist_lept_inv_massv2->Fill( lepton_invariant_massv2(this_bjet,this_elec, met,MH2) );
        }
		//cout<< endl;
		

		/*if ( topinvariantmass < 280  ) continue;
		event_counter6_after_topinvmass++;
		
		if ( topinvariantmass < 85  ) continue;
		event_counter7_after_leptinvmass++;*/
		
		//cout<< endl;
		
		////////////////////////////////////////////////////////////////////
		//	JET BRANCH
		////////////////////////////////////////////////////////////////////

		numb_jet=0;
		if(branchJet->GetEntries() > 0)
		{
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				jet[i] = (Jet*) branchJet->At(i);
				if ( (jet[i]->PT)>20 && abs(jet[i]->Eta)< 4.9 ) numb_jet++;
			}
		}
		
		my_event.Pt_btag=this_bjet->PT;
		
		my_event.N_jets=numb_jet;
		
		if (SaveData) outputTree->Fill();
		
		
//		cout.width(10);
//		cout << my_event.Pt_lepton;
//		cout.width(10);
//		cout << my_event.MET;
//		cout.width(10);
//		cout << my_event.Pt_btag;
//		cout.width(10);
//		cout << my_event.Mt_lv;
//		cout.width(10);
//		cout << my_event.Mi_lvb;
//		cout.width(10);
//		cout << my_event.N_jets;
//		cout.width(10);
//		cout << my_event.R_lb;
//		cout.width(10);
//		cout << my_event.Eta_lb;
//		cout.width(10);
//		cout << my_event.alpha_t;
//		cout << endl;
		
		
	}	//end of event loop
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	//---------------------------------------------------------------------------------------------------------
	////////////////////////////////  Saving the histograms into output  //////////////////////////////////////////
	//---------------------------------------------------------------------------------------------------------
	// end of event loop
	
	cout << " signal_counter :  " << signal_counter << "  %  " << (double)signal_counter/entries << endl;
    cout << "Number of Events Needed (for 1 fb^-1) : " << (double)cross*200*efficiency << endl;
	cout << "Total nevents : " << entries << endl;
	cout << "Counter_saved : " << counter_saved << endl;
	cout << "1  event after trigger       : " << event_counter1_after_trigger     << "\t" << (double)event_counter1_after_trigger/entries	<<  endl;
	cout << "2  event after 10gev         : " << event_counter2_after_10gev       << "\t" << (double)event_counter2_after_10gev/entries		<< endl;
	//cout << "3  event after Lep(PT)>55  : " << event_counter3_after_leptonpt55  << "\t" << (double)event_counter3_after_leptonpt55/entries << endl;
	cout << "4  event after alpha_t > 0.5 : " << event_counter4_after_alpha_t     << "\t" << (double)event_counter4_after_alpha_t/entries	<< endl;
	cout << "5  event after btag = 1      : " << event_counter5_after_btag        << "\t" << (double)event_counter5_after_btag/entries		<< endl;
	//cout << "6  event after topinvmass  : " << event_counter6_after_topinvmass  << "\t" << (double)event_counter6_after_topinvmass/entries << endl;
	//cout << "7  event after leptinvmass : " << event_counter7_after_leptinvmass << "\t" << (double)event_counter7_after_leptinvmass/entries << endl;
	//cout << "8  event after numb_jet =1 : " << event_counter8_after_onejet      << "\t" << (double)event_counter8_after_onejet/entries << endl;
	//cout << "9  event after numb_jet =1 : " << event_counter9_after_onejet      << "\t" << (double)event_counter9_after_onejet/entries << endl;
	//cout << "10 event after jet_eta<4.9 : " << event_counter10_after_jeteta     << "\t" << (double)event_counter10_after_jeteta/entries << endl;
	
	//myfile.close();
	outf->cd();
	outf = outputTree->GetCurrentFile(); //to get the pointer to the current file
	outf->Write();
	outf->Close();
	
	//end of main loop
}


