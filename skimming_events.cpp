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
#include <TH1F.h>
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

using namespace std;

double MH2=80.0*80.0;

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

//Calculates the top invariant mass using ELECTRON 
//and missing energy and higgs mass
double top_invariant_mass(Electron *ptlep, MissingET *miss)
{
	double plx=(ptlep->PT)*cos(ptlep->Phi);
	double ply=(ptlep->PT)*sin(ptlep->Phi);
	double plz=(ptlep->PT)*sinh(ptlep->Eta);
	double pl  =(ptlep->PT)*cosh(ptlep->Eta); //massless lepton energy

	double mh2=MH2;
	double mt=173.7;
	double misspx=(miss->MET)*cos(miss->Phi);
	double misspy=(miss->MET)*sin(miss->Phi);
	double a=pow(((double)plz/pl),2)-1;
	double b=2*((plx*misspx+ply*misspy)/(pl)+mh2/(2*pl))*(plz/pl);
	double c=((plx*misspx+ply*misspy)/(pl)+mh2/(2*pl))-pow((miss->MET),2);
	double res_sqrt=(b*b-4*a*c);

	//cout << a << "\t" << b << "\t" << c << "\t"  << res_sqrt << endl; 
	if (res_sqrt<0) return 10.0;
	else 
	{
		double result1=(-b- sqrt(res_sqrt))/(2*a);
		double result2=(-b+sqrt(res_sqrt))/(2*a);
		        if (result1>0&&result2<0) return result1;
		else if (result1<0&&result2>0) return result2;
		else if (result1>0&&result2>0) return min(result1-mt,result2-mt) +mt;
		return 0;
	}
}

double top_invariant_mass(Muon *ptlep, MissingET *miss)
{
	double plx=(ptlep->PT)*cos(ptlep->Phi);
	double ply=(ptlep->PT)*sin(ptlep->Phi);
	double plz=(ptlep->PT)*sinh(ptlep->Eta);
	double pl=(ptlep->PT)*cosh(ptlep->Eta);

	double mh2=MH2;
	double mt=173.7;
	double misspx=(miss->MET)*cos(miss->Phi);
	double misspy=(miss->MET)*sin(miss->Phi);
	double a=pow(((double)plz/pl),2)-1;
	double b=2*((plx*misspx+ply*misspy)/(pl)+mh2/(2*pl))*(plz/pl);
	double c=((plx*misspx+ply*misspy)/(pl)+mh2/(2*pl))-pow((miss->MET),2);
	double res_sqrt=(b*b-4*a*c);

	//cout << a << "\t" << b << "\t" << c << "\t"  << res_sqrt << endl; 
	if (res_sqrt<0) return 10.0;
	else 
	{

		double result1=(-b- sqrt(res_sqrt))/(2*a);
		double result2=(-b+sqrt(res_sqrt))/(2*a);
		        if (result1>0&&result2<0) return result1;
		else if (result1<0&&result2>0) return result2;
		else if (result1>0&&result2>0) return min(result1-mt,result2-mt) +mt;
		return 0;
	}
}


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


	TH1 *histJet_pt[10];
	TH1 *histJet_eta[10];
	TH1 *histJet_phi[10];
	char hist_name[100];
	TDirectory *jetdir= outf->mkdir("jets");
	jetdir->cd();

	for (int i=0; i<10; i++)
	{
		sprintf(hist_name,"jet_pt%i",i);
		histJet_pt[i] = new TH1F(hist_name, "jet P_{T}", 1000, 0.0, 1000.0);

		sprintf(hist_name,"jet_eta%i",i);
		histJet_eta[i] = new TH1F(hist_name, "jet eta", 100, -5.0, 5.0);

		sprintf(hist_name,"jet_phi%i",i);
		histJet_phi[i] = new TH1F(hist_name, "jet phi", 100, -5.0, 5.0);
	}
	
	//TH1 *histMass 	 = new TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0);
	TH1 *jet_size 		 = new TH1F("jet_size", "Number of Jets", 10, 0, 10.0);
	TH1 *histJet_btag 	 = new TH1F("histJet_btag", "Number of B-tagged jets", 10, 0, 10.0);
	TH1 *histJet_btag_pt = new TH1F("histJet_btag_pt", "PT of B-tagged jets", 100, 0.0, 500.0);

	TDirectory *lepdir= outf->mkdir("lepton");
	lepdir->cd();

 	TH1 *hist_gen_lepton  	   = new TH1F("gen_lepton"		, "Number of GEN Leptons ", 10, 0, 10.0);
    TH1 *hist_lepton_before_trig = new TH1F("numb_lepton_before_trig", "Number of SIM Leptons |eta|<2.5 and PT>20/30(elec/muon)", 10, 0, 10.0);
        TH1 *hist_lepton_pass_trig = new TH1F("numb_lepton_pass_trig"	 , "Number of SIM Leptons |eta|<2.5 and PT>20/30(elec/muon)", 10, 0, 10.0);
        TH1 *hist_lepton_pass_10gev = new TH1F("numb_lepton_pass_10gev"	 , "Number of SIM Leptons and PT>10(elec/muon)", 10, 0, 10.0);

	TH1 *histElec_pt 	= new TH1F("elec_pt1"	, "1st elec P_{T}", 100, 0.0, 500.0);
	TH1 *histElec_phi	= new TH1F("elec_pt1_phi", "1st elec Phi", 100, -5.0, 5.0);
	TH1 *histElec_eta	= new TH1F("elec_pt1_eta", "1st elec Eta", 100, -5.0, 5.0);

	TH1 *histMuon_pt 	= new TH1F("muon_pt1"	, "1st mu P_{T}", 100, 0.0, 500.0);  
	TH1 *histMuon_phi = new TH1F("muon_pt1_phi", "1st mu Phi", 	100, -5.0, 5.0);
	TH1 *histMuon_eta = new TH1F("muon_pt1_eta", "1st mu Eta", 100, -5.0, 5.0);
	
	TH1 *histLepton_pt 	= new TH1F("lepton_pt"	, " lepton P_{T}", 100, 0.0, 500.0);
	TH1 *histLepton_eta	= new TH1F("lepton_eta", "lepton Eta", 100, -5.0, 5.0);
	TH1 *histLepton_phi	= new TH1F("lepton_phi", "lepton  Phi", 100, -5.0, 5.0);	
	
	TH1 *lepton_invmass 	= new TH1F("lepton_invmass", "lepton inv mass", 100, 0, 500);    
	TH1 *top_invmass 	= new TH1F("top_invmass", "top inv mass", 100, 0, 500);    

	TDirectory *metdir= outf->mkdir("met");
	metdir->cd();

	TH1 *histMET_et  = new TH1F("histMET_et", "MET",  500,  0.0, 500.0);
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
	//unsigned int counter =0;
	
	if (NENTRIES==-1) entries=numberOfEntries;
	else entries=NENTRIES;
	
	outf->cd();
	TVectorD v(entries);
	v[0] = entries;
	v.Write("nevent");

	cout << "Reading TREE: " << numberOfEntries << " events available, \n";
	cout << "\t\t" << NENTRIES << " of them will be analyzed." << endl;
	
	GenParticle *part, *mother;//, *status;
	Jet *jet[10]; 
	Electron *elec[4]; //, *elec2, *elec3;
	Muon *mu[4]; //, *mu2, *mu3;
	MissingET *met;

	int counter_muon=0;
    int counter_elec=0;
    
    int event_counter_after_trigger=0;
    int event_counter_after_10gev=0;
    int event_counter_after_leptonpt55=0;
    int event_counter_after_met=0;
    int event_counter_after_btag=0;
    
	// Loop over all events
	for(Long64_t entry = 0; entry < entries; ++entry)
	{

		double progress = 10.0*entry/(1.0*entries);
		int k = TMath::FloorNint(progress);
		if (k > decade) {	cout << 10*k << "%\t"; gSystem->Exec("date '+%H:%M:%S'"); cout << endl;	}
		decade = k;

  		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		int part_counter=0;
		int genevent_size=0;
		genevent_size=branchParticle->GetEntries();
		//cout << "particle size : " << genev_size << "   ";

		if( genevent_size > 0)
		{
			for(int i=0; i<genevent_size; i++)
            {
				part = (GenParticle *) branchParticle->At(i);
				//cout << i << "\t" << part->PID<< "\t" << part->M1 << "\t" << part->M2 << "\t" << endl;
				int motherPID=0;				
				if (part->M1>0)
				{
					mother = (GenParticle *) branchParticle->At(part->M1);
					motherPID = mother->PID;
				}
				else motherPID=0;

				if( (abs(part->PID)== 11 && abs(motherPID)== 24) ||
					(abs(part->PID)== 13 && abs(motherPID)== 24) ||
					(abs(part->PID)== 15 && abs(motherPID)== 24) ) 
					{
						part_counter++; 
					//cout << motherPID << "    ";
				}
			}
			hist_gen_lepton->Fill(part_counter);
		}
		
		//////////////////////////////////////////
		// TRIGGER EMULATION
		//////////////////////////////////////////

		int numb_jet=0;
		int counter_btag=0;
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
				if(  (elec[i]->PT)>30 && abs(elec[i]->Eta)<2.5 ) 
				{
					numb_elec_pass_cuts++; 
					this_elec=elec[i];
				}
			}
			counter_elec+=numb_elec_pass_cuts;
		}

  		if( branchMuon->GetEntries() > 0)
        {
			for(int i=0; i<branchMuon->GetEntries(); i++)
         	{
	           	mu[i] = (Muon *) branchMuon->At(i);
            	if(  (mu[i]->PT)>20 && abs(mu[i]->Eta)<2.5 ) 
            	{
            		numb_muon_pass_cuts++;          
         			this_muon=mu[i];
				}
			}
			counter_muon+=numb_muon_pass_cuts;
		}

		// number of lepton before the trigger
		hist_lepton_before_trig->Fill((numb_muon_pass_cuts + numb_elec_pass_cuts));
		
		// if the total number of muons is higher than 0 fire for this event
        if ((numb_muon_pass_cuts + numb_elec_pass_cuts)==0) continue;
		event_counter_after_trigger++;

		// number of leptons after the trigger cut
		hist_lepton_pass_trig->Fill(numb_muon_pass_cuts + numb_elec_pass_cuts);	
		
		// will reuse these constants so set them to 0
		numb_muon_pass_cuts=0;
		numb_elec_pass_cuts=0;
		counter_elec=0;
		counter_muon=0;


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
			counter_elec+=numb_elec_pass_cuts;
		}
 
  		if( branchMuon->GetEntries() > 0)
        {
			for(int i=0; i<branchMuon->GetEntries(); i++)
         	{
	           	mu[i] = (Muon *) branchMuon->At(i);
            	if(  (mu[i]->PT)>10) numb_muon_pass_cuts++;                
			}
			counter_muon+=numb_muon_pass_cuts;
		}

		hist_lepton_pass_10gev->Fill(numb_muon_pass_cuts + numb_elec_pass_cuts);	

		//////////////////////////////////////////
		// FILTER THE EVENTS FOR LEPTON=1
		// filter the event if there is more than one lepton

        if ((numb_muon_pass_cuts + numb_elec_pass_cuts)!=1) continue;
		
		// number of leptons after the trigger cut
		event_counter_after_10gev++;


		//////////////////////////////////////////
		//////////////////////////////////////////
		//	ELECTRON AND MUON BRANCH
		//////////////////////////////////////////

		if(numb_elec_pass_cuts==1)
		{
      		//cout << this_elec->PT << endl;
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
      		//mu1 = (Muon *) branchMuon->At(this_muon);
      		//cout << mu1->PT << endl;
      		histMuon_pt->Fill(this_muon->PT);
      		histMuon_eta->Fill(this_muon->Eta);
      		histMuon_phi->Fill(this_muon->Phi); 	
			histLepton_pt->Fill(this_muon->PT);
			histLepton_eta->Fill(this_muon->Eta);
			histLepton_phi->Fill(this_muon->Phi);      		
		}


		// filter the event if there is lepton PT>55
		if ( numb_elec_pass_cuts==1 && this_elec->PT > 55 ) continue;
		if ( numb_muon_pass_cuts==1 && this_muon->PT > 55 ) continue;
		event_counter_after_leptonpt55++;
		
		//////////////////////////////////////////

// TODO
// make a comparison between the signal MET and background MET
// make a comparison and employ a cut on the background


		double this_met=0;
		if(branchMET->GetEntries() > 0)
	    { 
	    	met = (MissingET *) branchMET->At(0);
	    	histMET_et->Fill(met->MET);
	    	histMET_eta->Fill(met->Eta);
	    	histMET_phi->Fill(met->Phi);
			this_met=met->MET;
		}

	// TODO filter over this met ???
	//		if (this_met < ) continue;

	////////////////////////////////////////////////////////////////////	
	//	JET BRANCH
	////////////////////////////////////////////////////////////////////
	
	int counter_btagfilter=0;
	// filter bjet tagged events
 	if(branchJet->GetEntries() > 0)
	{
		for (int i=0; i< branchJet->GetEntries(); i++)
	   	{
  	  		jet[i] = (Jet*) branchJet->At(i);
	   		if (jet[i]->BTag==1) 
	   		{
	   			histJet_btag_pt->Fill(jet[i]->PT);	    			
	   			counter_btag++;
				if (jet[i]->PT < 75) counter_btagfilter++;
	   		}
		}
		histJet_btag->Fill(counter_btag);
	}


	// FILTER events having MET <50
	if (this_met <=50. ) continue;
	event_counter_after_met++;

	// or the number of b-tagged jets is not 1
	if ( counter_btagfilter!=1 ) continue;
    event_counter_after_btag++;            
	// TODO make a comparison plot for batg-jet 
	// multiplicity and pt distribution before and after

	//////////////////////////////////////////////////////////////////////////////////////
	// TOP INVARIANT MASS and LEPTONIC TRANSVERSE MASS
	//////////////////////////////////////////////////////////////////////////////////////
		

	if ( numb_muon_pass_cuts == 1 )	// there is a elec in the event
	{
		//cout << "e  "<< top_invariant_mass(elec1, this_met);
		top_invmass->Fill(top_invariant_mass(this_muon, met));
	}	
	else 
	if (numb_elec_pass_cuts == 1)	// there is a muon in the event
	{
		//cout << "m  "<< top_invariant_mass(mu1, this_met);
		top_invmass->Fill(top_invariant_mass(this_elec, met));
	}
	//cout<< endl;


	if ( numb_muon_pass_cuts == 1 )	// there is a elec in the event
	{
		//cout << "e  "<< lepton_invariant_mass(elec1, this_met);
		lepton_invmass->Fill(lepton_invariant_mass(this_muon, met));
	}
	else 
	if (numb_elec_pass_cuts == 1)	// there is a muon in the event
	{
		//cout << "m  "<< lepton_invariant_mass(mu1, this_met);
		lepton_invmass->Fill(lepton_invariant_mass(this_elec, met));
	}
	//cout<< endl;

	////////////////////////////////////////////////////////////////////	
	//	JET BRANCH
	////////////////////////////////////////////////////////////////////
  
	if(branchJet->GetEntries() > 0)
 	{
		jet_size->Fill( branchJet->GetEntries() );
		for (int i=0; i< branchJet->GetEntries(); i++)
	   	{
     		jet[i] = (Jet*) branchJet->At(i);
            if (i>9)continue;
			// Take first jet
      		histJet_pt[i]->Fill(jet[i]->PT);
      		histJet_eta[i]->Fill(jet[i]->Eta);
      		histJet_phi[i]->Fill(jet[i]->Phi);
      		if ( (jet[i]->PT)>30 && abs(jet[i]->Eta)< 4.9 )numb_jet++;
		}
	}
	if (numb_jet > 1)continue;



	}	//end of event loop

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//---------------------------------------------------------------------------------------------------------
	////////////////////////////////  Saving the histograms into output  //////////////////////////////////////////	
	//---------------------------------------------------------------------------------------------------------

	// end of event loop

	cout << "counter elec : " << counter_elec << endl;
	cout << "counter muon : " << counter_muon << endl;
	cout << "nevents : " << entries << endl;
	cout << "event after trigger : "  << event_counter_after_trigger << "\t" << (double)event_counter_after_trigger/entries <<  endl;
	cout << "event after 10gev : "    << event_counter_after_10gev << "\t" << (double)event_counter_after_10gev/entries << endl;
	cout << "event after Lepton(PT)>55 : "   << event_counter_after_leptonpt55 << "\t" << (double)event_counter_after_leptonpt55/entries << endl;
	cout << "event after met > 50 : " << event_counter_after_met << "\t" << (double)event_counter_after_met/entries << endl;
	cout << "event after btag = 1 : " << event_counter_after_btag << "\t" << (double)event_counter_after_btag/entries << endl;

	//myfile.close();
	outf->Write();	
	outf->Close(); 

	//end of main loop
}

