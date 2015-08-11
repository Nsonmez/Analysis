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
	int NENTRIES			= c1.getValue<int>		("NEntries");
	bool DATA				= c1.getValue<bool>		("DATA");
	bool SaveData			= c1.getValue<bool>		("SaveData");
	double MH2				= c1.getValue<double>	("MH2");
	double cross			= c1.getValue<double>	("CrossSection");
	double efficiency		= c1.getValue<double>	("Efficiency");

	double MUONPT_CUT=20;
	double ELECPT_CUT=30;
	unsigned MAX_NJETS =8;
	unsigned MAX_NLEPTONS=1;
	
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
	////////////////////////////////////////  My Data STructure  /////////////////////////////
	//----------------------------------------------------------------------------------------
	

	TTree *outputTree = new TTree("DATA","a tree of event info and jet info");
	//TTree::SetMaxTreeSize(100000000000LL);	//the maximum tree size
	outputTree->SetDirectory(outf);

	EVENT_VAR aux_event;
	JET aux_jet[MAX_NJETS];
	LEPTON aux_lepton[MAX_NLEPTONS];
	//SAMPLE aux_sample;
	char name[2000];
	

	//outputTree->Branch("sampleVar",&aux_sample,"cross/D:eff");

	
	//outputTree->Branch("eventVar",&aux_event,"PVz/D:PVx:PVy:HT:MET:METEta:METPhi:Weight");
	outputTree->Branch("eventVar",&aux_event,"HT/D:MET:METEta:METPhi");

	for(unsigned i=0;i<MAX_NJETS;i++)
	{
		sprintf(name,"jet%dVar",i+1);
		outputTree->Branch(name, &aux_jet[i],
					"PT/D:Eta:Phi:Mass:DeltaEta:DeltaPhi:Charge/I:BTag:TauTag");
	}

	
	for(unsigned i=0;i<MAX_NLEPTONS;i++)
	{
		sprintf(name,"mu%dVar",i+1);
		outputTree->Branch(name, &aux_lepton[i], "PT/D:Eta:Phi:Charge/I:ID");
	}
	
	cout<<"event tree is booked"<<endl;

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
	//TClonesArray *branchVertex = treeReader->UseBranch("Vertex");
	TClonesArray *branchMET 	 = treeReader->UseBranch("MissingET");
	TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");

	//----------------------------------------------------------------------------------------
	/////////////////////////////////////  LOOP  Over the EVENTS  ////////////////////////////
	//----------------------------------------------------------------------------------------
	
	cout << "Set Branch Addresses" << endl;
	
	int decade = 0;
	unsigned entries = 0;

	if (NENTRIES==-1) entries=numberOfEntries;
	else entries=NENTRIES;
	
	outf->cd();
	
	cout << "Reading TREE:  " << OutputFileName.c_str() << "   \t " << numberOfEntries << " events available, \n";
	cout << "\t\t" << NENTRIES << " of them will be analyzed." << endl;


	GenParticle *part, *mother1, *this_elec_part=0, *this_muon_part=0;//, *status;
	Jet *jet[10];
	Electron *elec[4];
	Electron *this_elec=0; //, *elec2, *elec3;
	Muon *mu[4];
	Muon *this_muon=0; //, *mu2, *mu3;
	GenParticle *partm, *partp, *motherm, *motherp, *mother2, *mother3;
	
    int event_counter0_after_trigger=0;
    int event_counter1_after_hardlepton=0;
    int event_counter2_after_10gev=0;
	int signal_counter=0;
	int counter_saved=0;
	
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

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
					this_elec_part=part;
				}
				
				if( abs(part->PID)== 13 && (part->Status)== 1
				   && abs(motherPID)==15 && part->PT > MUONPT_CUT
				   && abs(part->Eta) < 2.5 )
				{
					part_counter_muon++;
					this_muon_part=part;
				}
			}
		}
		
		if (the_signal==true)signal_counter++;
		
		// if this an event where
		// higgs_charged > tau nu_tau
		// tau > leptonic decay

		if ( DATA==true && the_signal == false ) continue;
		
		
        //////////////////////////////////////////
        // TRIGGER EMULATION
        //////////////////////////////////////////
		
        
        //filter lepton=1 events
        int numb_elec_pass_cuts=0;
        int numb_muon_pass_cuts=0;
        
        
        if( branchElectron->GetEntries() > 0)
        {
            for(int i=0; i<branchElectron->GetEntries(); i++)
            {
                elec[i] = (Electron *) branchElectron->At(i);
                if(  (elec[i]->PT)>25 ) numb_elec_pass_cuts++;
            }
        }
        
        if( branchMuon->GetEntries() > 0)
        {
            for(int i=0; i<branchMuon->GetEntries(); i++)
            {
                mu[i] = (Muon *) branchMuon->At(i);
                if(  (mu[i]->PT)>20) numb_muon_pass_cuts++;
            }
        }
        
        if ( (numb_muon_pass_cuts + numb_elec_pass_cuts)==0 ) continue;
        
        event_counter0_after_trigger++;

        //////////////////////////////////////////
        // HARD LEPTONs at THE CENTRAL REGION
        //////////////////////////////////////////

        
        if( branchElectron->GetEntries() > 0)
        {
            numb_elec_pass_cuts=0;
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
            numb_muon_pass_cuts=0;
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

        event_counter1_after_hardlepton++;

        // will reuse these constants so set them to 0
        int numb_softmuon_pass_cuts=0;
        int numb_softelec_pass_cuts=0;
        
        //////////////////////////////////////////
        // SECONDARY SOFT LEPTONS WITH PT > 10 GeV
        //////////////////////////////////////////
        
        if( branchElectron->GetEntries() > 0)
        {
            for(int i=0; i<branchElectron->GetEntries(); i++)
            {
                elec[i] = (Electron *) branchElectron->At(i);
                if( (elec[i]->PT)>10 ) numb_softelec_pass_cuts++;
            }
        }
        
        if( branchMuon->GetEntries() > 0)
        {
            for(int i=0; i<branchMuon->GetEntries(); i++)
            {
                mu[i] = (Muon *) branchMuon->At(i);
                if( (mu[i]->PT)>10) numb_softmuon_pass_cuts++;
            }
        }

        //////////////////////////////////////////
        // FILTER THE EVENTS FOR LEPTON=1
        // filter the event if there is more than one lepton
        
        if ((numb_softmuon_pass_cuts + numb_softelec_pass_cuts)!=1) continue;
		
        event_counter2_after_10gev++;
        
		if ( branchJet-> GetEntries() == 0 ) continue;

		EVENT_VAR event;

		if ( branchScalarHT-> GetEntries() > 0)
		{
			ScalarHT * scalarht = (ScalarHT *) branchScalarHT->At(0);
			event.HT= scalarht->HT;
		} else
			event.HT= 0;
		
		
		if ( branchMET-> GetEntries() > 0)
		{
			MissingET * met = (MissingET *) branchMET->At(0);
			event.MET=met->MET;
			event.METEta=met->Eta;
			event.METPhi=met->Phi;
		}
		else
		{
			event.MET=-999;
			event.METEta=-999;
			event.METPhi=-999;
		}
		
		vector<JET> all_jets;
	
		if(branchJet->GetEntries() > 0)
		{
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				JET jetvar;
				jet[i] = (Jet*) branchJet->At(i);
				if (jet[i]->PT<20) continue;
				jetvar.PT=jet[i]->PT;
				jetvar.Eta=jet[i]->Eta;
				jetvar.Phi=jet[i]->Phi;
				jetvar.Mass=jet[i]->Mass;
				jetvar.DeltaEta=jet[i]->DeltaEta;
				jetvar.DeltaPhi=jet[i]->DeltaPhi;
				jetvar.Charge=jet[i]->Charge;
				jetvar.BTag=jet[i]->BTag;
				jetvar.TauTag=jet[i]->TauTag;
				//jetvar.PartonPt=jet[i]->PartonPt();
				//jetvar.PdgId=jet[i]->PdgId();
				all_jets.push_back(jetvar);

			}
		}
		
		else
		{
			JET jetvar;
			jetvar.PT=-999;
			jetvar.Eta=-999;
			jetvar.Phi=-999;
			jetvar.Mass=-999;
			jetvar.DeltaEta=-999;
			jetvar.DeltaPhi=-999;
			jetvar.Charge=-999;
			jetvar.BTag=-999;
			jetvar.TauTag=-999;
			//jetvar.PartonPt=-999;
			//jetvar.PdgId=-999;
			all_jets.push_back(jetvar);
		}

		LEPTON lepton;

		if( numb_elec_pass_cuts > 0 )
		{
			lepton.PT=this_elec->PT;
			lepton.Eta=this_elec->Eta;
			lepton.Phi=this_elec->Phi;
			lepton.Charge=this_elec->Charge;
			//lepton.Particle=this_elec->Particle;
			lepton.ID=11;
		}
		
		else if( numb_muon_pass_cuts > 0 )
		{
			lepton.PT=this_muon->PT;
			lepton.Eta=this_muon->Eta;
			lepton.Phi=this_muon->Phi;
			lepton.Charge=this_muon->Charge;
			//lepton.Particle=this_muon->Particle;
			lepton.ID=13;
		}
		else
		{
			lepton.PT=-999;
			lepton.Eta=-999;
			lepton.Phi=-999;
			lepton.Charge=-999;
			//lepton.Particle=-999;
			lepton.ID=-999;
		}

		
		aux_event=event;
	
		//save calojets
		for(unsigned j=0; j < MAX_NJETS;j++)
		{
			if (all_jets.size() <= j) continue;
			aux_jet[j] = all_jets[j];
			
		}

		for(unsigned j=0; j < MAX_NLEPTONS;j++)
		{
			aux_lepton[j] = lepton;
		}

		if (SaveData) outputTree->Fill();
		
		counter_saved++;
		
	}	//end of event loop


	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
    TVectorD cross1(1);
    cross1[0] = cross;
    cross1.Write("cross");
    
    TVectorD eff1(1);
    eff1[0] = (double)counter_saved/entries;
    eff1.Write("eff");
    
    TVectorD neve(1);
    neve[0] = entries;
    neve.Write("neve");
    
	//---------------------------------------------------------------------------------------------------------
	////////////////////////////////  Saving the histograms into output  //////////////////////////////////////
	//---------------------------------------------------------------------------------------------------------
	// end of event loop
	cout << "_____________________________________________"<< endl;
	cout << "signal_counter : ";
	cout.width(10);
	cout << signal_counter << "  %  " << (double)signal_counter/entries << endl;
	cout << "counter_saved  : ";
 	cout.width(10);
	cout << counter_saved << "  %  " << (double)counter_saved/entries << endl;
    cout << "_____________________________________________________________________" << endl;
    cout << " 0  event after trigger     : "    << event_counter0_after_trigger     << "\t" << (double)event_counter0_after_trigger/entries		<< "\t" << (double)event_counter0_after_trigger/entries*cross*1000      << endl;
    cout << " 1  event after hard lepton : "    << event_counter1_after_hardlepton  << "\t" << (double)event_counter1_after_hardlepton/entries  << "\t" << (double)event_counter1_after_hardlepton/entries*cross*1000   << endl;
    cout << " 2  event after 10gev       : "    << event_counter2_after_10gev       << "\t" << (double)event_counter2_after_10gev/entries       << "\t" << (double)event_counter2_after_10gev/entries*cross*1000        << endl;

    cout << "_____________________________________________"<< endl;


	
	//myfile.close();
	outf->cd();
	outf = outputTree->GetCurrentFile(); //to get the pointer to the current file
	outf->Write();
	outf->Close();
	
	//end of main loop
}


