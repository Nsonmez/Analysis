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

/*
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

*/


/*
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
	
	TLorentzVector particle_miss1;
	TLorentzVector particle_miss2;
	particle_miss1.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result1*result1) );
	particle_miss2.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result2*result2) );
	double top_mass1=(particle_bjet + particle_lep + particle_miss1).M();
	double top_mass2=(particle_bjet + particle_lep + particle_miss2).M();
	
 
	if (false)
	{
		cout.width(10);
		cout << miss->Eta<< "\t";
		cout.width(10);
		cout << result1 << "\t";
		cout.width(10);
		cout << result2 << "\t";
		cout.width(10);
		cout << top_mass1 << "\t";
		cout.width(10);
		cout << top_mass2 << "\t";
		cout << endl;
	}
	
	
	//    if (top_mass1 < 0 && top_mass2 > 0 ) return top_mass2;
	//    else if (top_mass1 > 0 && top_mass2 < 0 ) return top_mass1;
	//    else if (top_mass1 > 0 && top_mass2 > 0 ) return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
	//    else if (top_mass1 < 0 && top_mass2 < 0 ) return 10;
	//    else return 5;
	
	return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
	
}

double top_invariant_mass(Jet *bjet, Muon *ptlep, MissingET *miss, double MH2)
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
	
	

//	 double result=(-b+sqrt(res_sqrt))/(2*a);
//	 TLorentzVector particle_miss;
//	 particle_miss.SetPxPyPzE(misspx, misspy, result, sqrt(misspx*misspx + misspy*misspy + result*result) );
//	 double top_mass=(particle_bjet + particle_lep + particle_miss).M();

	
	TLorentzVector particle_miss1;
	TLorentzVector particle_miss2;
	particle_miss1.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result1*result1) );
	particle_miss2.SetPtEtaPhiE(miss->MET, miss->Eta, miss->Phi, sqrt( miss->MET*miss->MET + result2*result2) );
	double top_mass1=(particle_bjet + particle_lep + particle_miss1).M();
	double top_mass2=(particle_bjet + particle_lep + particle_miss2).M();
	
	if (false)
	{
		cout.width(10);
		cout << miss->Eta<< "\t";
		cout.width(10);
		cout << result1 << "\t";
		cout.width(10);
		cout << result2 << "\t";
		cout.width(10);
		cout << top_mass1 << "\t";
		cout.width(10);
		cout << top_mass2 << "\t";
		cout << endl;
	}
	
	
	//    if (top_mass1 < 0 && top_mass2 > 0 ) return top_mass2;
	//    else if (top_mass1 > 0 && top_mass2 < 0 ) return top_mass1;
	//    else if (top_mass1 > 0 && top_mass2 > 0 ) return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
	//    else if (top_mass1 < 0 && top_mass2 < 0 ) return 10;
	//    else return 5;
	
	return min( abs(174.-top_mass1), abs(174.-top_mass2) )+174.;
	
}

double top_invariant_massv2(Muon *ptlep, MissingET *miss, double MH2)
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
	
	if (res_sqrt<0) return 10.0;
	
	double result1=(-b+sqrt(res_sqrt))/(2*a);
	double result2=(-b-sqrt(res_sqrt))/(2*a);
	
	if (true)
	{
		cout.width(10);
		cout << miss->Eta<< "\t";
		cout.width(10);
		cout << result1 << "\t";
		cout.width(10);
		cout << result2 << "\t";
		cout << endl;
	}
	return result1;
}

double top_invariant_massv2(Electron *ptlep, MissingET *miss, double MH2)
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
	
	if (res_sqrt<0) return 10.0;
	
	double result1=(-b+sqrt(res_sqrt))/(2*a);
	double result2=(-b-sqrt(res_sqrt))/(2*a);
	
	if (true)
	{
		cout.width(10);
		cout << miss->Eta<< "\t";
		cout.width(10);
		cout << result1 << "\t";
		cout.width(10);
		cout << result2 << "\t";
		cout << endl;
	}
	return result1;
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


*/
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
	bool DATA				= c1.getValue<bool>	("DATA");
	double noneed			= c1.getValue<double>("MH2");
	
	
	double MUONPT_CUT=24;
	double ELECPT_CUT=34;
	
	
	cout << "________________________________________________________________\n";
	cout << "This is a SIGNAL(true)/BACKGROUND(false) : " << DATA << endl;
	
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
	
	
	//TH1 *histJet_pt[10];
	//TH1 *histJet_eta[10];
	//TH1 *histJet_phi[10];
	char hist_name[100];
	
	//----------------------------------------------------------------------------------------
	TDirectory *jetdir= outf->mkdir("jets");
	jetdir->cd();
	
	for (int i=0; i<10; i++)
	{
		sprintf(hist_name,"jet_pt%i",i);
		//histJet_pt[i] = new TH1F(hist_name, "jet P_{T}", 1000, 0.0, 1000.0);
		
		sprintf(hist_name,"jet_eta%i",i);
		//histJet_eta[i] = new TH1F(hist_name, "jet eta", 100, -5.0, 5.0);
		
		sprintf(hist_name,"jet_phi%i",i);
		//histJet_phi[i] = new TH1F(hist_name, "jet phi", 100, -5.0, 5.0);
	}
	
	TH1 *jet_size 		 = new TH1F("jet_size",     "Number of Jets", 10, 0, 10.0);
	TH1 *histJet_btag 	 = new TH1F("histJet_btag", "Number of B-tagged jets", 10, 0, 10.0);
	TH1 *histJet_btag_pt = new TH1F("histJet_btag_pt", "PT of B-tagged jets", 100, 0.0, 500.0);
	TH1 *jet_size_cut8   = new TH1F("jet_size_cut8", "Number of Jets with 30<PT GeV", 10, 0, 10.0);
	TH1 *jet_size_cut9   = new TH1F("jet_size_cut9", "Number of Jets with 15<PT<30 GeV", 10, 0, 10.0);
	TH1 *hist_before_jet_eta = new TH1F("hist_before_jet_eta", "Jet Eta Before cut 9", 100, -5.0, 5.0);
	TH1 *hist_alpha_t	 = new TH1F("hist_alpha_t", "Apha_T PT_2/M12", 100, 0, 5.0);
	
	
	//----------------------------------------------------------------------------------------
	TDirectory *lepdir= outf->mkdir("lepton");
	lepdir->cd();
	
	TH1 *hist_gen_lepton		= new TH1F("gen_lepton"			, "Number of GEN Leptons ;GEN;SIM", 10, 0, 10.0);
	TH2 *hist_gen_lepton2D		= new TH2I("hist_gen_lepton2D"	, "Number of GEN Leptons ;GEN;SIM", 10, 0, 10.0,10,0,10);
	
	TH2 *hist_gen_elec_SIM1GEN0	= new TH2F("hist_gen_elec_SIM1GEN0"	, "Number of Elec : GEN=0  and SIM=1 ;PT;ETA", 100, 0, 100,100,-5,+5);
	TH2 *hist_gen_elec_SIM0GEN1	= new TH2F("hist_gen_elec_SIM0GEN1"	, "Number of Elec : GEN=1  and SIM=0 ;PT;ETA", 100, 0, 100,100,-5,+5);
	
	TH2 *hist_gen_muon_SIM1GEN0	= new TH2F("hist_gen_muon_SIM1GEN0"	, "Number of Muon : GEN=0  and SIM=1 ;PT;ETA", 100, 0, 100,100,-5,+5);
	TH2 *hist_gen_muon_SIM0GEN1	= new TH2F("hist_gen_muon_SIM0GEN1"	, "Number of Muon : GEN=1  and SIM=0 ;PT;ETA", 100, 0, 100,100,-5,+5);
	
	TH1 *hist_gen_elec_pt		= new TH1F("hist_gen_elec_pt"	, "Elec Pt " , 100, 0, 250);
	TH1 *hist_gen_elec_phi		= new TH1F("hist_gen_elec_phi"	, "Elec Phi ", 100, -4, 4);
	TH1 *hist_gen_elec_eta		= new TH1F("hist_gen_elec_eta"	, "Elec Eta ", 100, -5, 5);
	TH1 *hist_gen_elec_deltaR	= new TH1F("hist_gen_elec_deltaR" , "DeltaR Elec Eta ", 100, 0, 5);
	
	TH1 *hist_gen_muon_pt		= new TH1F("hist_gen_muon_pt"  , "Muon Pt " , 100, 0, 250);
	TH1 *hist_gen_muon_phi		= new TH1F("hist_gen_muon_phi" , "Muon Phi ", 100, -4, 4);
	TH1 *hist_gen_muon_eta		= new TH1F("hist_gen_muon_eta" , "Muon Eta ", 100, -5, 5);
	TH1 *hist_gen_muon_deltaR	= new TH1F("hist_gen_muon_deltaR" , "DeltaR Muon Eta ", 100, 0, 5);
	
	
	TH1 *hist_lepton_before_trig = new TH1F("numb_lepton_before_trig", "Number of SIM Leptons |eta|<2.5 and PT>20/30(elec/muon)", 10, 0, 10.0);
	TH1 *hist_lepton_pass_trig   = new TH1F("numb_lepton_pass_trig"	 , "Number of SIM Leptons |eta|<2.5 and PT>20/30(elec/muon)", 10, 0, 10.0);
	TH1 *hist_lepton_pass_10gev  = new TH1F("numb_lepton_pass_10gev" , "Number of SIM Leptons and PT>10(elec/muon)", 10, 0, 10.0);
	
	TH1 *histElec_pt 	= new TH1F("elec_pt1"	 , "1st elec P_{T}", 100, 0.0, 500.0);
	TH1 *histElec_phi	= new TH1F("elec_pt1_phi", "1st elec Phi  ", 100, -5.0, 5.0);
	TH1 *histElec_eta	= new TH1F("elec_pt1_eta", "1st elec Eta  ", 100, -5.0, 5.0);
	
	TH1 *histMuon_pt  = new TH1F("muon_pt1"     , "1st mu P_{T}", 100, 0.0, 500.0);
	TH1 *histMuon_phi = new TH1F("muon_pt1_phi" , "1st mu Phi  ", 100, -5.0, 5.0);
	TH1 *histMuon_eta = new TH1F("muon_pt1_eta" , "1st mu Eta  ", 100, -5.0, 5.0);
	
	TH1 *histLepton_pt 	= new TH1F("lepton_pt"  , "lepton P_{T} ", 100, 0.0, 500.0);
	TH1 *histLepton_eta	= new TH1F("lepton_eta" , "lepton Eta   ", 100, -5.0, 5.0);
	TH1 *histLepton_phi	= new TH1F("lepton_phi" , "lepton  Phi  ", 100, -5.0, 5.0);
	
	TH1 *lepton_invmass = new TH1F("lepton_invmass", "lepton inv mass", 100, 0, 500);
	//TH1 *top_invmass 	= new TH1F("top_invmass"   , "top inv mass   ", 100, 0, 500);
	
	TH1 *hist_top_inv_mass_mhc80 	= new TH1F("top_inv_mass_mhc80"   , "top inv mass assume mhc=80  ", 100, 0, 500);
	TH1 *hist_top_inv_mass_mhc100 	= new TH1F("top_inv_mass_mhc100"   , "top inv mass assume mhc=100  ", 100, 0, 500);
	TH1 *hist_top_inv_mass_mhc130 	= new TH1F("top_inv_mass_mhc130"   , "top inv mass  assume mhc=130 ", 100, 0, 500);
	
	TH1 *hist_miss_pz_mhc80 	= new TH1F("miss_pz_mhc80"   , "miss_pz assume mhc=80  ", 100, 0, 500);
	TH1 *hist_miss_pz_mhc100 	= new TH1F("miss_pz_mhc100"   , "miss_pz assume mhc=100  ", 100, 0, 500);
	TH1 *hist_miss_pz_mhc130 	= new TH1F("miss_pz_mhc130"   , "miss_pz assume mhc=130  ", 100, 0, 500);
	
	TH2 *bjet_lepton_delta_eta 	= new TH2F("bjet_lepton_delta_eta"   , "Delta Eta between lepton and btagjet ", 100, 0, 250, 100, -5,   +5);
	TH2 *bjet_lepton_delta_phi 	= new TH2F("bjet_lepton_delta_phi"   , "Delta Phi between lepton and btagjet ", 100, 0, 250, 100, -6.3, +6.3);
	TH2 *bjet_lepton_deltaR     = new TH2F("bjet_lepton_deltaR"      , "Delta R between lepton and btagjet   ", 100, 0, 250, 100, 0,  +20);
	
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
	//unsigned int counter =0;
	
	if (NENTRIES==-1) entries=numberOfEntries;
	else entries=NENTRIES;
	
	outf->cd();
	TVectorD v(entries);
	v[0] = entries;
	v.Write("nevent");
	
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
	
	int event_counter1_after_trigger    =0;
	int event_counter2_after_10gev      =0;
	int event_counter3_after_leptonpt55 =0;
	int event_counter4_after_met        =0;
	int event_counter5_after_btag       =0;
	int event_counter6_after_topinvmass =0;
	int event_counter7_after_leptinvmass=0;
	int event_counter8_after_onejet     =0;
	int event_counter9_after_onejet     =0;
	int event_counter10_after_jeteta    =0;
	
	int part_counter=0;
	int signal_counter=0;
	//int count_print=0;
	
	// Loop over all events
	for(Long64_t entry = 0; entry < 0.5*entries; ++entry)
	{
		
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
					hist_gen_elec_phi->Fill(part->Phi);
					this_elec_part=part;
				}
				
				if( abs(part->PID)== 13 && (part->Status)== 1
				   && abs(motherPID)==15 && part->PT > MUONPT_CUT
				   && abs(part->Eta) < 2.5 )
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
		
		if ( the_signal ) signal_counter++;
		
		if ( the_signal == false && DATA==true ) continue;
		
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
		
		hist_gen_lepton2D->Fill(part_counter_elec,  numb_elec_pass_cuts);
		hist_gen_lepton2D->Fill(part_counter_muon+3,numb_muon_pass_cuts+3);
		
		
		if (part_counter_elec==1 && numb_elec_pass_cuts==1)
		{
			TLorentzVector v1;
			TLorentzVector v2;
			v1.SetPxPyPzE(this_elec_part->Px,this_elec_part->Py,this_elec_part->Pz,this_elec_part->E);
			v2.SetPtEtaPhiM(this_elec->PT,this_elec->Eta,this_elec->Phi,0.00000511);
			hist_gen_elec_deltaR->Fill( ROOT::Math::VectorUtil::DeltaR(v1,v2) );
			//double deltar= ROOT::Math::VectorUtil::DeltaR(v1,v2);
			//cout << " deltaR   : " << deltar << endl;
		}
		
		
		if (part_counter_muon==1  && numb_muon_pass_cuts==1)
		{
			TLorentzVector v1;
			TLorentzVector v2;
			v1.SetPxPyPzE(this_muon_part->Px,this_muon_part->Py,this_muon_part->Pz,this_muon_part->E);
			v2.SetPtEtaPhiM(this_muon->PT,this_muon->Eta,this_muon->Phi,0.000100);
			hist_gen_muon_deltaR->Fill(ROOT::Math::VectorUtil::DeltaR(v1,v2));
		}
		
		// Gen=1 and SIM=0
		if (part_counter_elec==1 && numb_elec_pass_cuts==0)	hist_gen_elec_SIM0GEN1->Fill(this_elec_part->PT,this_elec_part->Eta);
		if (part_counter_muon==1 && numb_muon_pass_cuts==0)	hist_gen_muon_SIM0GEN1->Fill(this_muon_part->PT,this_muon_part->Eta);
		
		// Gen=0 and SIM=1
		if (part_counter_elec==0 && numb_elec_pass_cuts==1)	hist_gen_elec_SIM1GEN0->Fill(this_elec->PT,this_elec->Eta);
		if (part_counter_muon==0 && numb_muon_pass_cuts==1)	hist_gen_muon_SIM1GEN0->Fill(this_muon->PT,this_muon->Eta);
		
		if ( part_counter_elec == 1 && numb_elec_pass_cuts == 1 )	hist_lepton_before_trig->Fill( 5);
		if ( part_counter_muon == 1 && numb_muon_pass_cuts == 1 )	hist_lepton_before_trig->Fill( 6);
		
		
		// number of lepton before the trigger
		hist_lepton_before_trig->Fill(numb_muon_pass_cuts + numb_elec_pass_cuts);
		
		// if the total number of muons is higher than 0 fire for this event
		if ((numb_muon_pass_cuts + numb_elec_pass_cuts)==0) continue;
		event_counter1_after_trigger++;
		
		// number of leptons after the trigger cut
		hist_lepton_pass_trig->Fill(numb_muon_pass_cuts + numb_elec_pass_cuts);
		
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
		
		hist_lepton_pass_10gev->Fill(numb_muon_pass_cuts + numb_elec_pass_cuts);
		
		//////////////////////////////////////////
		// FILTER THE EVENTS FOR LEPTON=1
		// filter the event if there is more than one lepton
		
		if ((numb_muon_pass_cuts + numb_elec_pass_cuts)!=1) continue;
		
		// number of leptons after the trigger cut
		event_counter2_after_10gev++;
		
		
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
		event_counter3_after_leptonpt55++;
		
		//////////////////////////////////////////
		
		double this_met=0;
		if(branchMET->GetEntries() > 0)
		{
			met = (MissingET *) branchMET->At(0);
			histMET_et->Fill(met->MET);
			histMET_eta->Fill(met->Eta);
			histMET_phi->Fill(met->Phi);
			this_met=met->MET;
		}
		
		// FILTER events having MET <50
		if (this_met < 50. ) continue;
		event_counter4_after_met++;
		
		////////////////////////////////////////////////////////////////////
		//	JET BRANCH
		////////////////////////////////////////////////////////////////////
		Jet *this_jet=0;
		
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
					if (jet[i]->PT < 75){
						this_jet=jet[i];
						counter_btagfilter++;
					}
				}
			}
			histJet_btag->Fill(counter_btag);
		}
		
		if (branchJet->GetEntries() > 1)
		{
			jet[0] = (Jet*) branchJet->At(0);
			jet[1] = (Jet*) branchJet->At(1);
			
			double inv_jet_mass12=sqrt(2*jet[0]->PT*jet[1]->PT*(cosh(jet[0]->Eta-jet[1]->Eta)-cos(jet[0]->Phi-jet[1]->Phi)));
			double alpha_t=(double)(jet[1]->PT)/inv_jet_mass12;
			hist_alpha_t->Fill(alpha_t);
		}
		
		
		// or the number of b-tagged jets is not 1
		if ( counter_btagfilter!=1 ) continue;
		event_counter5_after_btag++;
		
		if (numb_elec_pass_cuts==1)
		{
			double deltaR2=pow(abs( (this_jet->Eta)-(this_elec->Eta) ),2)+pow(abs( (this_jet->Phi)-(this_elec->Phi) ),2);
			bjet_lepton_delta_eta->Fill( this_elec->PT, abs( (this_jet->Eta)-(this_elec->Eta) ) );
			bjet_lepton_delta_phi->Fill( this_elec->PT, abs( (this_jet->Phi)-(this_elec->Phi) ));
			bjet_lepton_deltaR->Fill( this_elec->PT, sqrt(deltaR2) );
		}
		else if (numb_muon_pass_cuts==1)
		{
			double deltaR2=pow( (this_jet->Eta)-(this_muon->Eta) ,2)+pow( (this_jet->Phi)-(this_muon->Phi) ,2);
			bjet_lepton_delta_eta->Fill( this_muon->PT, abs( (this_jet->Eta)-(this_muon->Eta) ) );
			bjet_lepton_delta_phi->Fill( this_muon->PT, abs( (this_jet->Phi)-(this_muon->Phi) ));
			bjet_lepton_deltaR->Fill( this_muon->PT, sqrt(deltaR2) );
		}
		
		// TODO make a comparison plot for batg-jet
		// multiplicity and pt distribution before and after
		
		//////////////////////////////////////////////////////////////////////////////////////
		// TOP INVARIANT MASS
		//////////////////////////////////////////////////////////////////////////////////////
		
		double top_inv_mass_mhc80=0;
		double top_inv_mass_mhc100=0;
		double top_inv_mass_mhc130=0;
		
		double miss_pz_mhc80=0;
		double miss_pz_mhc100=0;
		double miss_pz_mhc130=0;
		
		if ( numb_muon_pass_cuts == 1 )	// there is a elec in the event
		{
			miss_pz_mhc80=missing_energy_pz(this_muon,  met, 6400);
			miss_pz_mhc100=missing_energy_pz(this_muon, met, 10000);
			miss_pz_mhc130=missing_energy_pz(this_muon, met, 16900);
			
			top_inv_mass_mhc80=top_invariant_mass(this_jet, this_muon, met,  6400);
			top_inv_mass_mhc100=top_invariant_mass(this_jet, this_muon, met, 10000);
			top_inv_mass_mhc130=top_invariant_mass(this_jet, this_muon, met, 16900);
			
			hist_miss_pz_mhc80->Fill(miss_pz_mhc80);
			hist_miss_pz_mhc100->Fill(miss_pz_mhc100);
			hist_miss_pz_mhc130->Fill(miss_pz_mhc130);
			
			hist_top_inv_mass_mhc80->Fill(top_inv_mass_mhc80);
			hist_top_inv_mass_mhc100->Fill(top_inv_mass_mhc100);
			hist_top_inv_mass_mhc130->Fill(top_inv_mass_mhc130);
			//top_invariant_massv2(this_muon, met,  6400);
		}
		else
			if ( numb_elec_pass_cuts == 1 ) // there is a muon in the event
			{
				miss_pz_mhc80=missing_energy_pz(this_elec,  met, 6400);
				miss_pz_mhc100=missing_energy_pz(this_elec, met, 10000);
				miss_pz_mhc130=missing_energy_pz(this_elec, met, 16900);
				
				top_inv_mass_mhc80=top_invariant_mass(this_jet, this_elec,  met, 6400);
				top_inv_mass_mhc100=top_invariant_mass(this_jet, this_elec, met, 10000);
				top_inv_mass_mhc130=top_invariant_mass(this_jet, this_elec, met, 16900);
				
				hist_miss_pz_mhc80->Fill(miss_pz_mhc80);
				hist_miss_pz_mhc100->Fill(miss_pz_mhc100);
				hist_miss_pz_mhc130->Fill(miss_pz_mhc130);
				
				hist_top_inv_mass_mhc80->Fill(top_inv_mass_mhc80);
				hist_top_inv_mass_mhc100->Fill(top_inv_mass_mhc100);
				hist_top_inv_mass_mhc130->Fill(top_inv_mass_mhc130);
				//top_invariant_massv2(this_elec, met, 6400);

			}
		
		//if ( topinvariantmass < 120 || topinvariantmass > 280  ) continue;
		event_counter6_after_topinvmass++;
		
		
		////////////////////////////////////////////////////////////////////
		//	LEPTONIC TRANSVERSE MASS
		////////////////////////////////////////////////////////////////////
		
		double leptoninvariantmass=0;
		
		if ( numb_muon_pass_cuts == 1 )	// there is a elec in the event
		{
			leptoninvariantmass=lepton_invariant_mass(this_muon, met);
			lepton_invmass->Fill(leptoninvariantmass);
		}
		else
			if ( numb_elec_pass_cuts == 1 )	// there is a muon in the event
			{
				leptoninvariantmass=lepton_invariant_mass(this_elec, met);
				lepton_invmass->Fill(leptoninvariantmass);
			}
		
		if ( leptoninvariantmass > 65 && leptoninvariantmass < 115 ) continue;
		event_counter7_after_leptinvmass++;
		
		
		////////////////////////////////////////////////////////////////////
		//	JET BRANCH
		////////////////////////////////////////////////////////////////////
  
		numb_jet=0;
		if(branchJet->GetEntries() > 0)
		{
			jet_size->Fill( branchJet->GetEntries() );
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				jet[i] = (Jet*) branchJet->At(i);
				if ( (jet[i]->PT)>30 && abs(jet[i]->Eta)< 4.9 ) numb_jet++;
			}
		}
		jet_size_cut8->Fill(numb_jet);
		if (numb_jet != 1) continue;
		event_counter8_after_onejet++;
		
		numb_jet=0;
		if(branchJet->GetEntries() > 0)
		{
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				jet[i] = (Jet*) branchJet->At(i);
				if ( (jet[i]->PT)>15 && (jet[i]->PT)<30 && abs(jet[i]->Eta)< 4.9 ) numb_jet++;
			}
		}
		
		jet_size_cut9->Fill(numb_jet);
		if (numb_jet==1) continue;
		event_counter9_after_onejet++;
		
		
		double jet_eta=0.;
		if(branchJet->GetEntries() > 0)
		{
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				jet[i] = (Jet*) branchJet->At(i);
				hist_before_jet_eta->Fill(jet[i]->Eta);
				if ( (jet[i]->PT)>15 && abs(jet[i]->Eta) > 4.9 ) jet_eta++;
			}
		}
		
		if (jet_eta > 0) continue;
		event_counter10_after_jeteta++;
		
	}	//end of event loop
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	//---------------------------------------------------------------------------------------------------------
	////////////////////////////////  Saving the histograms into output  //////////////////////////////////////////
	//---------------------------------------------------------------------------------------------------------
	// end of event loop
	
	cout << " signal_counter :  " << signal_counter << endl;
	
	cout << "nevents : " << entries << endl;
	cout << "1  event after trigger     : " << event_counter1_after_trigger     << "\t" << (double)event_counter1_after_trigger/entries <<  endl;
	cout << "2  event after 10gev       : " << event_counter2_after_10gev       << "\t" << (double)event_counter2_after_10gev/entries << endl;
	cout << "3  event after Lep(PT)>55  : " << event_counter3_after_leptonpt55  << "\t" << (double)event_counter3_after_leptonpt55/entries << endl;
	cout << "4  event after met > 50    : " << event_counter4_after_met         << "\t" << (double)event_counter4_after_met/entries << endl;
	cout << "5  event after btag = 1    : " << event_counter5_after_btag        << "\t" << (double)event_counter5_after_btag/entries << endl;
	cout << "6  event after topinvmass  : " << event_counter6_after_topinvmass  << "\t" << (double)event_counter6_after_topinvmass/entries << endl;
	cout << "7  event after leptinvmass : " << event_counter7_after_leptinvmass << "\t" << (double)event_counter7_after_leptinvmass/entries << endl;
	cout << "8  event after numb_jet =1 : " << event_counter8_after_onejet      << "\t" << (double)event_counter8_after_onejet/entries << endl;
	cout << "9  event after numb_jet =1 : " << event_counter9_after_onejet      << "\t" << (double)event_counter9_after_onejet/entries << endl;
	cout << "10 event after jet_eta<4.9 : " << event_counter10_after_jeteta     << "\t" << (double)event_counter10_after_jeteta/entries << endl;
	
	//myfile.close();
	outf->Write();
	outf->Close();
	
	//end of main loop
}


