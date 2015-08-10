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
	
    const int btag[12]={35,40,45,50,55,60,65,70,80,90,100,110};
	TH1F *histJet_btag_ptbigger[number_of_datasets][12];
	TH1F *histJet_btag_eta_a[number_of_datasets];
	TH1F *histJet_btag_pt_a[number_of_datasets];

	for (int a=0;a<number_of_datasets;a++)
	{
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
	
    for (int a=0;a<number_of_datasets;a++)
    {
        cout.width(3);
        cout << a;
        cout.width(20);
        cout << data_name[a];
        for (int i=0; i< 10;i++)
        {
            cout.width(10);
            cout << histJet_btag_ptbigger[a][4]->GetBinContent(i);
        }
        cout << endl;
        
        cout.width(23);
        cout << "";
        for (int i=0; i< 10;i++)
        {
            cout.width(10);
            cout << jet_btag[a]->GetBinContent(i);
        }
        cout << endl;
    }
    
    
	//______________________________________________________________________

	bool comparison_for_btag=true;

	
	
	/////////////////////////////////////////////////////////////////////////
	//// COMPARISON FOR BTAG
	/////////////////////////////////////////////////////////////////////////
	

	TCanvas *canvas7a = new TCanvas("name7a","title",2500,0,600,500);
    TH1F * btag_numb_sig  = new TH1F("btag_numb_sig" ,"BTAG Jet Numb Distrib Sig;Number of BTAG Jets;Relative Occurence", 10, 0, 10);
    TH1F * btag_numb_back = new TH1F("btag_numb_back","BTAG Jet Numb Distrib Back;Number of BTAG Jets;Relative Occurence", 10, 0, 10);

	for (int i=0; i<12; i++)
	{

        for (int j=0;j<10;j++)
        {
            btag_numb_sig->SetBinContent(j,0);
            btag_numb_back->SetBinContent(j,0);
        }
        cout << " --> comparison_for_btag" << endl;

        for (int a=0;a<4;a++)
		{
            if (histJet_btag_ptbigger[a][i]->GetEntries()==0) continue;
			histJet_btag_ptbigger[a][i]->Scale( scale[a] );
			btag_numb_sig->Add(histJet_btag_ptbigger[a][i]);
		}
		btag_numb_sig->Scale(1.0/btag_numb_sig->Integral());
		
		
		for (int a=4;a<number_of_datasets;a++)
        {
            if (histJet_btag_ptbigger[a][i]->GetEntries()==0) continue;
            histJet_btag_ptbigger[a][i]->Scale( scale[a] );
			btag_numb_back->Add(histJet_btag_ptbigger[a][i]);
        }
		btag_numb_back->Scale(1.0/btag_numb_back->Integral());

        //double signal=btag_numb_sig->GetBinContent(2);
        //double backgr=btag_numb_back->GetBinContent(2);
		//cout << "pt  " << btag[i] << "  ratio " << (double)signal/backgr <<  "\t" << signal << endl;
			
		btag_numb_sig->GetYaxis()->SetRangeUser(0.00001,1);
		btag_numb_back->GetYaxis()->SetRangeUser(0.00001,1);
		btag_numb_back->GetXaxis()->SetRangeUser(0,5);
		btag_numb_sig->GetXaxis()->SetRangeUser(0,5);
		btag_numb_sig->SetLineColor(2);
		btag_numb_sig->SetLineWidth(3);
		btag_numb_back->SetLineWidth(3);
		btag_numb_sig->Draw("HTEXT0");
		btag_numb_back->Draw("same HTEXT0");
		//btag_numb_sig->Draw("");
		//btag_numb_back->Draw("same");
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
		
	
	
    TCanvas *canvas7b = new TCanvas("name7b","title",2500,0,600,500);

    TH1F * btag_numb_sigv2  = new TH1F("btag_numb_sigv2" ,"BTAG Jet Numb Distrib Sig;Number of BTAG Jets;Relative Occurence", 10, 0, 10);
    TH1F * btag_numb_backv2 = new TH1F("btag_numb_backv2","BTAG Jet Numb Distrib Back;Number of BTAG Jets;Relative Occurence", 10, 0, 10);
    // analysis over the signal
    for (int a=0;a<4;a++)
    {
        if (jet_btag[a]->GetEntries()==0) continue;
        jet_btag[a]->Scale( scale[a] );
        btag_numb_sigv2->Add(jet_btag[a]);
    }
    btag_numb_sigv2->Scale(1.0/btag_numb_sigv2->Integral());
    
    for (int a=4;a<number_of_datasets;a++)
    {
        if (jet_btag[a]->GetEntries()==0) continue;
        jet_btag[a]->Scale( scale[a] );
        btag_numb_backv2->Add(jet_btag[a]);
    }
    btag_numb_backv2->Scale(1.0/btag_numb_backv2->Integral());
    
    btag_numb_sigv2->GetYaxis()->SetRangeUser(0.00001,1);
    btag_numb_backv2->GetYaxis()->SetRangeUser(0.00001,1);
    btag_numb_backv2->GetXaxis()->SetRangeUser(0,5);
    btag_numb_sigv2->GetXaxis()->SetRangeUser(0,5);
    btag_numb_sigv2->SetLineColor(2);
    btag_numb_sigv2->SetLineWidth(3);
    btag_numb_backv2->SetLineWidth(3);
    btag_numb_sigv2->Draw("HTEXT0");
    btag_numb_backv2->Draw("same HTEXT0");
    //btag_numb_sig->Draw("");
    //btag_numb_back->Draw("same");
    btag_numb_sigv2->SetFillColor(2);
    btag_numb_sigv2->SetFillStyle(3002);
    canvas7b->BuildLegend();
    //canvas7a->SetLogy();
    canvas7b->SetLeftMargin(0.127517);
    canvas7b->SetRightMargin(0.02013423);
    btag_numb_backv2->SetStats(0);
    btag_numb_sigv2->SetStats(0);
    sprintf(h_name,"6btag_numb_mhc%i.eps", mhc);
    canvas7b->Print(h_name);
    
    
    
    cout << " --> comparison_for_btag pt" << endl;

    
   	TCanvas *canvas7c = new TCanvas("name7c","title",2500,0,600,500);
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


	btag_sig_pt_a->GetYaxis()->SetRangeUser(0.00001,0.2);
	btag_back_pt_a->GetYaxis()->SetRangeUser(0.00001,0.2);
	btag_back_pt_a->GetXaxis()->SetRangeUser(0,300);
	btag_sig_pt_a->GetXaxis()->SetRangeUser(0,300);
	btag_sig_pt_a->SetLineColor(2);
	btag_sig_pt_a->SetLineWidth(3);
	btag_back_pt_a->SetLineWidth(3);
	//btag_numb_sig->Draw("HTEXT0");
	//btag_numb_back->Draw("same HTEXT0");
	btag_sig_pt_a->Draw("");
	btag_back_pt_a->Draw("same");
	btag_sig_pt_a->SetFillColor(2);
	btag_sig_pt_a->SetFillStyle(3002);
	canvas7c->BuildLegend();
	//canvas7c->SetLogy();
	canvas7c->SetLeftMargin(0.127517);
	canvas7c->SetRightMargin(0.02013423);
	btag_back_pt_a->SetStats(0);
	btag_sig_pt_a->SetStats(0);
	sprintf(h_name,"6btag_pt_mhc%i.C", mhc);
	canvas7c->Print(h_name);




}