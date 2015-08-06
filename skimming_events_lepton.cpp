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
	int NENTRIES			= c1.getValue<int>	("NEntries");
	double MUONPT_CUT=24;
	double ELECPT_CUT=34;
	
	if (!c1.check()) return 0;
	c1.print(); // Printing the options

	//string outputfile = OutputFileName + "_" + OutputFileTag;
	//TFile *outf = new TFile(outputfile.c_str(),"RECREATE");
	TFile *outf = new TFile(OutputFileName.c_str(),"RECREATE");

	cout << "________________________________________________________________\n";
	cout << "\n";
	cout << "time start  " << endl;
	gSystem->Exec("date '+%H:%M:%S'");

//----------------------------------------------------------------------------------------
///////////////////////////////////////////////  Histogram Output  ///////////////////////
//----------------------------------------------------------------------------------------
// Book histograms

    //----------------------------------------------------------------------------------------
	TDirectory *lepdir= outf->mkdir("lepton");
	lepdir->cd();

	TH1 *hist_gen_lepton		= new TH1F("gen_lepton"			, "Number of GEN Leptons ;GEN;SIM", 10, 0, 10.0);
	
	TH1 *hist_gen_elec_pt		= new TH1F("hist_gen_elec_pt"	, "Elec Pt " , 100, 0, 250);
	TH1 *hist_gen_elec_phi		= new TH1F("hist_gen_elec_phi"	, "Elec Phi ", 100, -4, 4);
	TH1 *hist_gen_elec_eta		= new TH1F("hist_gen_elec_eta"	, "Elec Eta ", 100, -5, 5);

	TH1 *hist_gen_muon_pt		= new TH1F("hist_gen_muon_pt"  , "Muon Pt " , 100, 0, 250);
	TH1 *hist_gen_muon_phi		= new TH1F("hist_gen_muon_phi" , "Muon Phi ", 100, -4, 4);
	TH1 *hist_gen_muon_eta		= new TH1F("hist_gen_muon_eta" , "Muon Eta ", 100, -5, 5);


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
	
	GenParticle *part, *mother, *this_elec_part=0, *this_muon_part=0;//, *status;
	
	////////////////////////////////////////////////////////////////
	//	FILTER the events
	////////////////////////////////////////////////////////////////
    
	
    // activate the branches
//    chain.SetBranchStatus("Event",1);
//    chain.SetBranchStatus("Jet",1);
//	chain.SetBranchStatus("Particle",1);
//	chain.SetBranchStatus("Electron",1);
//	chain.SetBranchStatus("Muon",1);
//	chain.SetBranchStatus("MissingET",1);

	//output file
	TFile *newfile = new TFile("small.root","recreate");

	// I clone an EMPTY tree
	TTree* theclonetree = chain.CloneTree(0);

	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////

	// Loop over all events
	for(Long64_t entry = 0; entry < entries; ++entry)
	{
		double progress = 10.0*entry/(1.0*entries);
		int k = TMath::FloorNint(progress);
		if (k > decade)
		{
			cout << 10*k << "%\t";
			gSystem->Exec("date '+%H:%M:%S'");
			cout << endl;
		}
		decade = k;

  		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		
		int part_counter=0;
		int part_counter_elec=0;
		int part_counter_muon=0;
		int genevent_size=0;

        genevent_size=branchParticle->GetEntries();
		//cout << "particle size : " << genev_size << "   ";
		if( genevent_size > 0 )
		{
			for(int i=0; i<genevent_size; i++)
            {
				part = (GenParticle *) branchParticle->At(i);
				//cout << i << "\t" << part->PID<< "\t" << part->M1 << "\t" << part->Status << "\t" << endl;
				
				int motherPID=0;
				if (part->M1>1)
				{
					mother = (GenParticle *) branchParticle->At(part->M1);
					motherPID = mother->PID;
				}
				else motherPID=0;

				if( (abs(part->PID)== 11 && (part->Status)== 1 && abs(motherPID)==15 ) ||
					(abs(part->PID)== 13 && (part->Status)== 1 && abs(motherPID)==15 ) )
				{
					part_counter++;
					//cout << motherPID << "    ";
					//cout<<endl;
					//cout << i << "\t" << part->PID<< "\t" << part->M1 << "\t" << part->Status << "\t" << endl;
				}
				
				if( abs(part->PID)== 11 && (part->Status)== 1 && abs(motherPID)==15 && part->PT > ELECPT_CUT && abs(part->Eta) < 2.5 )
				{
					part_counter_elec++;
					hist_gen_elec_pt->Fill(part->PT);
					hist_gen_elec_eta->Fill(part->Eta);
					hist_gen_elec_phi->Fill(part->Phi);
					this_elec_part=part;
				}
				if( abs(part->PID)== 13 && (part->Status)== 1 && abs(motherPID)==15 && part->PT > MUONPT_CUT && abs(part->Eta) < 2.5 )
				{
					part_counter_muon++;
					hist_gen_muon_pt ->Fill(part->PT);
					hist_gen_muon_eta->Fill(part->Eta);
					hist_gen_muon_phi->Fill(part->Phi);
					this_muon_part=part;
				}
			}
			hist_gen_lepton->Fill(part_counter);
			hist_gen_lepton->Fill(part_counter_muon+part_counter_elec+5);
			
		}
        
		
		// fill the tree
        if (part_counter_muon+part_counter_elec > 0)
        {
            cout << "event : " << entry << "\t" << part_counter_muon+part_counter_elec << endl;
            theclonetree->Fill();
        }													//		//delete newfile;

	}	//end of event loop

	//---------------------------------------------------------------------------------------------------------
	////////////////////////////////  Saving the histograms into output  //////////////////////////////////////////	
	//---------------------------------------------------------------------------------------------------------
	// end of event loop

	//myfile.close();
	newfile->Write();
	newfile->Close();
	outf->Write();
	outf->Close(); 

	//end of main loop
}


