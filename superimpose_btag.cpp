//
//  superimpose_btag.cpp
//  
//
//  Created by Nasuf Sonmez on 8/9/15.
//
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <time.h>

#include <TSystem.h>
#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH1D.h>

#include <TF1.h>
#include <TH2F.h>

#include "TLorentzVector.h"
#include "TCanvas.h"
#include <TMath.h>
#include <TProfile.h>
#include "TStyle.h"
#include <TString.h>
#include "TLatex.h"
#include "TLegend.h"
#include "TVectorD.h"
#include "Utilities.h"

//#include "/Users/nsonmez/root/macros/tdrStyle.C"

using namespace std;

#include "superimpose.h"
//int main(int argc, char*argv[])

int main()
{
	gStyle->SetOptStat("n");
	
	//setTDRStyle();
	char h_name[1000];
	char hist_name[100];

	int mhc=100;
	
	
	//______________________________________________________________
	// read the root files
	TFile * rootfile[20];
	cout << "__________________________________________\n";
	for (int a=0;a<number_of_datasets;a++)
	{
		sprintf(h_name,"plots/plots_%s_mhc%i.root", data_name[a], mhc);
		cout << h_name << endl;
		rootfile[a] = new TFile( h_name, "READ" );
	}
	cout << "__________________________________________"<<endl;
	
	
	double scale[number_of_datasets];
	for (int a=0;a<number_of_datasets;a++)
	{
		double eff_lumi=(double)nevents[a]/(cross_sections[a]);
		scale[a]=1.0/eff_lumi;
	}
	
	//______________________________________________________________
	// get the lepton
	
	int btag[12]={30,35,40,45,50,55,60,65,70,80,90,100};
	TH1F *histJet_btag_ptbigger[number_of_datasets][12];
	TH1F *histJet_btag_eta_a[number_of_datasets];
	TH1F *histJet_btag_pt_a[number_of_datasets];
	
	for (int a=0;a<number_of_datasets;a++)
	{
		gen_lepton2D[a]				= (TH2F*)rootfile[a]->Get("lepton/hist_gen_lepton2D");
		
		gen_lepton[a]               = (TH1F*)rootfile[a]->Get("lepton/hist_gen_lepton");
		lepton_numb_before_trig[a]	= (TH1F*)rootfile[a]->Get("lepton/numb_lepton_before_trig");
		lepton_numb_pass_trig[a]	= (TH1F*)rootfile[a]->Get("lepton/numb_lepton_pass_trig");
		lepton_before_hardphot[a]	= (TH1F*)rootfile[a]->Get("lepton/hist_lepton_before_hardphot");
		lepton_numb_before_10gev[a]	= (TH1F*)rootfile[a]->Get("lepton/numb_lepton_before_10gev");
		
		
		lepton_pt[a]				= (TH1F*)rootfile[a]->Get("lepton/lepton_pt");
		lepton_eta[a]				= (TH1F*)rootfile[a]->Get("lepton/lepton_eta");
		lepton_phi[a]				= (TH1F*)rootfile[a]->Get("lepton/lepton_phi");
		
		met_pt[a]					= (TH1F*)rootfile[a]->Get("met/histMET_et");
		
		alpha_t[a]					= (TH1F*)rootfile[a]->Get("jets/hist_alpha_t");
		
		jet_btag[a]					= (TH1F*)rootfile[a]->Get("jets/histJet_btag");
		jet_btag_pt[a]				= (TH1F*)rootfile[a]->Get("jets/histJet_btag_pt");
		jet_btag_eta[a]				= (TH1F*)rootfile[a]->Get("jets/histJet_btag_eta");
		
		histJet_btag_eta_a[a]	    = (TH1F*)rootfile[a]->Get("jets/histJet_btag_eta_afterbjet1cut");
		histJet_btag_pt_a[a]	    = (TH1F*)rootfile[a]->Get("jets/histJet_btag_pt_afterbjet1cut");

		for (int i=0; i<12; i++)
		{
			sprintf(hist_name,"jets/histJet_btag%i",btag[i]);
			histJet_btag_ptbigger[a][i] = (TH1F*)rootfile[a]->Get(hist_name);
		}
	}
	cout << "Histos reading ..." << endl;
	//______________________________________________________________________
	
	
	//______________________________________________________________________

	bool comparison_for_btag=true;

	
	
	/////////////////////////////////////////////////////////////////////////
	//// COMPARISON FOR BTAG
	/////////////////////////////////////////////////////////////////////////
	
	if(comparison_for_btag)
	{
		TH1F * btag_numb_sig  = new TH1F("btag_numb_sig" ,"BTAG Jet Numb Distrib Sig;Number of BTAG Jets;Relative Occurence", 10, 0, 10);
		TH1F * btag_numb_back = new TH1F("btag_numb_back","BTAG Jet Numb Distrib Back;Number of BTAG Jets;Relative Occurence", 10, 0, 10);

		TCanvas *canvas7a = new TCanvas("name7a","title",2500,0,600,500);

		for (int i=0; i<12; i++)
		{
			cout << " --> comparison_for_btag" << endl;
		
			// analysis over the signal
			for (int a=0;a<4;a++)
			{
				histJet_btag_ptbigger[a][i]->Scale( scale[a] );
				btag_numb_sig->Add(histJet_btag_ptbigger[a][i]);
			}
			btag_numb_sig->Scale(1.0/btag_numb_sig->Integral());
			
			
			for (int a=4;a<number_of_datasets;a++)
			{
				histJet_btag_ptbigger[a][i]->Scale( scale[a] );
				btag_numb_back->Add(histJet_btag_ptbigger[a][i]);
			}
			btag_numb_back->Scale(1.0/btag_numb_back->Integral());

			double signal=btag_numb_sig->GetBinContent(2);
			double backgr=btag_numb_back->GetBinContent(2);

			cout << "pt  " << btag[i] << "  ratio " << (double)signal/backgr <<  "\t" << signal << endl;
			
			btag_numb_sig->GetYaxis()->SetRangeUser(0.00001,1);
			btag_numb_back->GetYaxis()->SetRangeUser(0.00001,1);
			btag_numb_back->GetXaxis()->SetRangeUser(0,5);
			btag_numb_sig->GetXaxis()->SetRangeUser(0,5);
			btag_numb_sig->SetLineColor(2);
			btag_numb_sig->SetLineWidth(3);
			btag_numb_back->SetLineWidth(3);
			//btag_numb_sig->Draw("HTEXT0");
			//btag_numb_back->Draw("same HTEXT0");
			btag_numb_sig->Draw("");
			btag_numb_back->Draw("same");
			btag_numb_sig->SetFillColor(2);
			btag_numb_sig->SetFillStyle(3002);
			canvas7a->BuildLegend();
			//canvas7a->SetLogy();
			canvas7a->SetLeftMargin(0.127517);
			canvas7a->SetRightMargin(0.02013423);
			btag_numb_back->SetStats(0);
			btag_numb_sig->SetStats(0);
			sprintf(h_name,"6btag_numb_mhc%i_pt%i.eps", mhc,btag[i]);
			canvas7a->Print(h_name);
		}
		
		
	}

	
	
	
	TCanvas *canvas7a = new TCanvas("name7a","title",2500,0,600,500);

	cout << " --> comparison_for_btag pt" << endl;
	
	TH1F * btag_sig_pt_a  = new TH1F("btag_sig_pt_a" ,"BTAG Jet Numb Distrib Sig;Number of BTAG Jets;Relative Occurence", 100, 0.0, 500.0);
	TH1F * btag_back_pt_a = new TH1F("btag_back_pt_a","BTAG Jet Numb Distrib Back;Number of BTAG Jets;Relative Occurence", 100, 0.0, 500.0);
	
	for (int a=0;a<4;a++)
	{
		histJet_btag_pt_a[a]->Scale( scale[a] );
		btag_sig_pt_a->Add(histJet_btag_pt_a[a]);
	}
	btag_sig_pt_a->Scale(1.0/btag_sig_pt_a->Integral());
	
	
	for (int a=4;a<number_of_datasets;a++)
	{
		histJet_btag_pt_a[a]->Scale( scale[a] );
		btag_back_pt_a->Add(histJet_btag_pt_a[a]);
	}
	btag_back_pt_a->Scale(1.0/btag_back_pt_a->Integral());
	
	
	//btag_sig_pt_a->GetYaxis()->SetRangeUser(0.00001,1);
	//btag_back_pt_a->GetYaxis()->SetRangeUser(0.00001,1);
	//btag_back_pt_a->GetXaxis()->SetRangeUser(0,5);
	//btag_sig_pt_a->GetXaxis()->SetRangeUser(0,5);
	btag_sig_pt_a->SetLineColor(2);
	btag_sig_pt_a->SetLineWidth(3);
	btag_back_pt_a->SetLineWidth(3);
	//btag_numb_sig->Draw("HTEXT0");
	//btag_numb_back->Draw("same HTEXT0");
	btag_sig_pt_a->Draw("");
	btag_back_pt_a->Draw("same");
	btag_sig_pt_a->SetFillColor(2);
	btag_sig_pt_a->SetFillStyle(3002);
	canvas7a->BuildLegend();
	//canvas7a->SetLogy();
	canvas7a->SetLeftMargin(0.127517);
	canvas7a->SetRightMargin(0.02013423);
	btag_back_pt_a->SetStats(0);
	btag_sig_pt_a->SetStats(0);
	sprintf(h_name,"6btag_pt_mhc%i.eps", mhc);
	canvas7a->Print(h_name);





}