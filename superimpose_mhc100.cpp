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


	//______________________________________________________________
	// get the lepton


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

		lepton_invmass[a]			= (TH1F*)rootfile[a]->Get("lepton/lepton_invmass");
		top_invmass[a]				= (TH1F*)rootfile[a]->Get("lepton/top_inv_mass_mhc");
		miss_pz_mhc[a]				= (TH1F*)rootfile[a]->Get("lepton/miss_pz_mhc");

		jet_pt0[a]					= (TH1F*)rootfile[a]->Get("jets/jet_pt0");
		jet_pt1[a]					= (TH1F*)rootfile[a]->Get("jets/jet_pt1");
		jet_pt2[a]					= (TH1F*)rootfile[a]->Get("jets/jet_pt2");

		jet_eta0[a]					= (TH1F*)rootfile[a]->Get("jets/jet_eta0");
		jet_eta1[a]					= (TH1F*)rootfile[a]->Get("jets/jet_eta1");
		jet_eta2[a]					= (TH1F*)rootfile[a]->Get("jets/jet_eta2");

		jet_size[a]					= (TH1F*)rootfile[a]->Get("jets/jet_size");
		jet_size_cut8[a]			= (TH1F*)rootfile[a]->Get("jets/jet_size_cut8");
		jet_size_cut9[a]			= (TH1F*)rootfile[a]->Get("jets/jet_size_cut9");
		
        jet_size_cut8_35[a]			= (TH1F*)rootfile[a]->Get("jets/jet_size_cut8_35");
        jet_size_cut8_40[a]			= (TH1F*)rootfile[a]->Get("jets/jet_size_cut8_40");
        jet_size_cut8_45[a]			= (TH1F*)rootfile[a]->Get("jets/jet_size_cut8_45");
        jet_size_cut8_50[a]			= (TH1F*)rootfile[a]->Get("jets/jet_size_cut8_50");

		jet_size_cut9_35[a]			= (TH1F*)rootfile[a]->Get("jets/jet_size_cut9_35");
		jet_size_cut9_40[a]			= (TH1F*)rootfile[a]->Get("jets/jet_size_cut9_40");
		jet_size_cut9_45[a]			= (TH1F*)rootfile[a]->Get("jets/jet_size_cut9_45");
		jet_size_cut9_50[a]			= (TH1F*)rootfile[a]->Get("jets/jet_size_cut9_50");

	
        bjet_lepton_dR[a]           =(TH2F*)rootfile[a]->Get("lepton/bjet_lepton_deltaR");
        bjet_lepton_dEta[a]         =(TH2F*)rootfile[a]->Get("lepton/bjet_lepton_delta_eta");
        bjet_lepton_dPhi[a]         =(TH2F*)rootfile[a]->Get("lepton/bjet_lepton_delta_phi");

        
	}
	cout << "Histos reading ..." << endl;
    //______________________________________________________________________

    // 2D GEN vs SIM leptons, basically acceptance for the signal
    if (false)
    {
        TCanvas *canvas2D = new TCanvas("canvas2D","acceptance",1901,0,900,800);
        canvas2D->Divide(2,2);

        TLatex *t = new TLatex();
        t->SetNDC();
        t->SetTextAlign(22);
        t->SetTextFont(63);
        t->SetTextSizePixels(22);

    
        for (int a=0;a<4;a++)
        {
            canvas2D->cd(a+1);
            gPad->SetLogz();
            gPad->SetGridx(1);
            gPad->SetGridy(1);

            gen_lepton2D[a]->GetXaxis()->SetRangeUser(0,7);
            gen_lepton2D[a]->GetYaxis()->SetRangeUser(0,7);
            gen_lepton2D[a]->SetMarkerSize(1.5);
            //gen_lepton2D[a]->Scale( 0.1 );
      
            sprintf(h_name,"Number of Gen vs SIM Leptons in %s;GEN;SIM", data_name[a]);
            gen_lepton2D[a]->SetTitle(h_name);
            gen_lepton2D[a]->Draw("COLZ TEXT");
            t->DrawLatex(0.25,0.5,"Electron");
            t->DrawLatex(0.60,0.82,"Muon");

        }
		
		sprintf(h_name,"lepton_acceptance_mhc%i.eps", mhc);
        canvas2D->Print(h_name);
    
    }
	//______________________________________________________________________

	
	bool lepton_acceptance=false;
	bool comparison_for_lepton=false;
	bool comparison_lepton_pt=true;
	bool comparison_for_met=true;
	bool comparison_for_alphat=true;
	bool comparison_for_btag=true;
	bool comparison_for_invmass=true;
	bool comparison_for_jetcut=true;
    bool comparison_for_jetsize=false;
    bool comparison_for_bjet_lepton=true;
	
	
	//______________________________________________________________________
	// LEPTON spectrum for background and the signal

    if (lepton_acceptance)
    {
		
 
        TH1F * gen_lepton_numb_sig          = new TH1F("lepton_gen_numb_sig"  ,"Gen level ;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
        TH1F * lepton_before_trig_sig       = new TH1F("lepton_before_trig_sig"  ,"Before trigg ;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
        TH1F * lepton_pass_trig_sig         = new TH1F("lepton_pass_trig_sig"  ,"Pass trigger ;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
        TH1F * lepton_before_hardphot_sig   = new TH1F("lepton_before_hardphot_sig"  ,"Hard photon cut ;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
        TH1F * lepton_before_10gev_sig      = new TH1F("lepton_before_10gev_sig"  ,"Before Soft photon cut ;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
        
        for (int a=0;a<4;a++)
        {
            double scale= (1.0/gen_lepton[a]->GetEntries())*(cross_sections[a])/total_cross_s;
            gen_lepton[a]->Scale( scale );
            gen_lepton_numb_sig->Add(gen_lepton[a]);
        }
        for (int a=0;a<4;a++)
        {
            double scale= (1.0/lepton_numb_before_trig[a]->GetEntries())*(cross_sections[a])/total_cross_s;
            lepton_numb_before_trig[a]->Scale( scale );
            lepton_before_trig_sig->Add(lepton_numb_before_trig[a]);
        }
        for (int a=0;a<4;a++)
        {
            double scale= (1.0/lepton_numb_pass_trig[a]->GetEntries())*(cross_sections[a])/total_cross_s;
            lepton_numb_pass_trig[a]->Scale( scale );
            lepton_pass_trig_sig->Add(lepton_numb_pass_trig[a]);
        }
        for (int a=0;a<4;a++)
        {
            double scale= (1.0/lepton_before_hardphot[a]->GetEntries())*(cross_sections[a])/total_cross_s;
            lepton_before_hardphot[a]->Scale( scale );
            lepton_before_hardphot_sig->Add(lepton_before_hardphot[a]);
        }
        for (int a=0;a<4;a++)
        {
            double scale= (1.0/lepton_numb_before_10gev[a]->GetEntries())*(cross_sections[a])/total_cross_s;
            lepton_numb_before_10gev[a]->Scale( scale );
            lepton_before_10gev_sig->Add(lepton_numb_before_10gev[a]);
        }

        TCanvas *canvas0b = new TCanvas("name0abcd","gen_lepton_numb_diff1",0,0,600,500);

		// number of leptons in gen level
        gen_lepton_numb_sig->Draw();
        gen_lepton_numb_sig->SetLineWidth(3);
        gen_lepton_numb_sig->SetLineColor(1);
		
		// number of leptons in sim level before trigger emulation
        lepton_before_trig_sig->Draw("same");
        lepton_before_trig_sig->SetLineWidth(3);
        lepton_before_trig_sig->SetLineColor(2);
        
        // number of leptons in sim level after trigger emulation
        lepton_pass_trig_sig->Draw("same");
        lepton_pass_trig_sig->SetLineWidth(3);
        lepton_pass_trig_sig->SetLineColor(4);

        canvas0b->BuildLegend();
		canvas0b->SetLogy();
		gen_lepton_numb_sig->SetStats(0);
		
		gen_lepton_numb_sig->GetXaxis()->SetRangeUser(0,5);
		lepton_before_trig_sig->GetXaxis()->SetRangeUser(0,5);
		
		canvas0b->SetLeftMargin(0.127517);
		canvas0b->SetRightMargin(0.02013423);
		sprintf(h_name,"lepton_numb_gen_mhc%i.eps", mhc);
		canvas0b->Print(h_name);
		
		
		
		TCanvas *canvas0c = new TCanvas("name0c","gen_lepton_numb_diff2",100,0,600,500);


		// number of leptons in sim level after hard photon cut
        lepton_before_hardphot_sig->Draw("");
        lepton_before_hardphot_sig->SetLineWidth(3);
        lepton_before_hardphot_sig->SetLineColor(2);
       
		// number of leptons in sim level after soft photon cut
        lepton_before_10gev_sig->Draw("same");
        lepton_before_10gev_sig->SetLineWidth(3);
        lepton_before_10gev_sig->SetLineColor(4);

        canvas0c->BuildLegend();
        canvas0c->SetLogy();
        lepton_pass_trig_sig->SetStats(0);
        
        lepton_pass_trig_sig->GetXaxis()->SetRangeUser(0,5);
        lepton_before_hardphot_sig->GetXaxis()->SetRangeUser(0,5);
        lepton_before_10gev_sig->GetXaxis()->SetRangeUser(0,5);

        canvas0c->SetLeftMargin(0.127517);
        canvas0c->SetRightMargin(0.02013423);

		sprintf(h_name,"lepton_numb_trig_mhc%i.eps", mhc);
		canvas0c->Print(h_name);
    }


    if (comparison_for_lepton)
    {
            
		cout << " --> comparison_for_lepton" << endl;
		
        //______________________________________________________________________
        // done
        
        TH1F * lepton_numb_before_trig_sig   = new TH1F("lepton_numb_before_trig_sig"  ,"sig ;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
        TH1F * lepton_numb_before_trig_back  = new TH1F("lepton_numb_before_trig_back" ,"back;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
		cout << "testing" << endl;

        for (int a=0;a<4;a++)
        {
            double scale= (1.0/lepton_numb_before_trig[a]->GetEntries())*(cross_sections[a])/total_cross_s;
            lepton_numb_before_trig[a]->Scale( scale );
            lepton_numb_before_trig_sig->Add(lepton_numb_before_trig[a]);
        }
		
        for (int a=4;a<number_of_datasets;a++)
        {
            if (lepton_numb_before_trig[a]->GetEntries()==0)continue;
            double scale= (1.0/lepton_numb_before_trig[a]->GetEntries())*(cross_sections[a])/total_cross_b;
            lepton_numb_before_trig[a]->Scale( scale );
            lepton_numb_before_trig_back->Add(lepton_numb_before_trig[a]);
        }

        TCanvas *canvas1a = new TCanvas("name1abc","beforetrigg",2500,0,600,500);
        lepton_numb_before_trig_sig->Draw();
        lepton_numb_before_trig_sig->SetLineColor(2);
        lepton_numb_before_trig_back->Draw("same");
        lepton_numb_before_trig_sig->SetLineWidth(3);
        lepton_numb_before_trig_back->SetLineWidth(3);
        lepton_numb_before_trig_sig->SetFillColor(2);
        lepton_numb_before_trig_sig->SetFillStyle(3002);

        canvas1a->BuildLegend();
        canvas1a->SetLogy();
        lepton_numb_before_trig_sig->SetStats(0);
        lepton_numb_before_trig_back->SetStats(0);
        lepton_numb_before_trig_sig->GetXaxis()->SetRangeUser(0,5);
        lepton_numb_before_trig_back->GetXaxis()->SetRangeUser(0,5);

        canvas1a->SetLeftMargin(0.127517);
        canvas1a->SetRightMargin(0.02013423);
		
		sprintf(h_name,"lepton_numb_before_trig_mhc%i.eps", mhc);
		canvas1a->Print(h_name);

        //______________________________________________________________________
        //
	
		/*
		TCanvas *canvas2a = new TCanvas("name2abc","after10GeV",2500,0,600,500);

        TH1F * signal_added_10gev_sig   = new TH1F("signal_added_10gev_sig"  ,"sig ;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
        TH1F * signal_added_10gev_back  = new TH1F("signal_added_10gev_back" ,"back;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);

        for (int a=0;a<4;a++)
        {
            double scale= (1.0/lepton_numb_before_10gev[a]->GetEntries())*(cross_sections[a])/total_cross_s;
            lepton_numb_before_10gev[a]->Scale( scale );
            signal_added_10gev_sig->Add(lepton_numb_before_10gev[a]);
        }
        for (int a=4;a<number_of_datasets;a++)
        {
            double scale= (1.0/lepton_numb_before_10gev[a]->GetEntries())*(cross_sections[a])/total_cross_b;
            lepton_numb_before_10gev[a]->Scale( scale );
            signal_added_10gev_back->Add(lepton_numb_before_10gev[a]);
        }
        signal_added_10gev_sig->Draw();
        signal_added_10gev_sig->SetLineColor(2);
        signal_added_10gev_sig->SetLineWidth(2);
        signal_added_10gev_back->SetLineWidth(2);
        signal_added_10gev_back->Draw("same");
        //canvas2a->BuildLegend(0.7055058,0.3357271,0.9475032,0.8617594);
        canvas2a->BuildLegend();
        //canvas2a->SetLogy();
        signal_added_10gev_sig->SetStats(0);
        signal_added_10gev_back->SetStats(0);
        canvas2a->SetLeftMargin(0.127517);
        canvas2a->SetRightMargin(0.02013423);
 
		sprintf(h_name,"signal_added_10gevmhc%i.eps", mhc);
		canvas2a->Print(h_name);
		 */
        //______________________________________________________________________

		/*
        TCanvas *canvas3a = new TCanvas("name3abc","afterTrigg",2500,0,600,500);

        TH1F * lepton_numb_after_trig_sig   = new TH1F("lepton_numb_after_trig_sig"  ,"sig ;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
        TH1F * lepton_numb_after_trig_back  = new TH1F("lepton_numb_after_trig_back" ,"back;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);

        for (int a=0;a<4;a++)
        {
            double scale= (1.0/lepton_numb_pass_trig[a]->GetEntries())*(cross_sections[a])/total_cross_s;
            lepton_numb_pass_trig[a]->Scale( scale );
            lepton_numb_after_trig_sig->Add(lepton_numb_pass_trig[a]);
        }
        for (int a=4;a<number_of_datasets;a++)
        {
            double scale= (1.0/lepton_numb_pass_trig[a]->GetEntries())*(cross_sections[a])/total_cross_b;
            lepton_numb_pass_trig[a]->Scale( scale );
            lepton_numb_after_trig_back->Add(lepton_numb_pass_trig[a]);
        }
        lepton_numb_after_trig_sig->Draw();
        lepton_numb_after_trig_sig->SetLineColor(2);
        lepton_numb_after_trig_back->Draw("same");
        lepton_numb_after_trig_sig->SetLineWidth(2);
        lepton_numb_after_trig_back->SetLineWidth(2);
        canvas3a->BuildLegend();
        //canvas3a->SetLogy();
        lepton_numb_after_trig_sig->SetStats(0);
        lepton_numb_after_trig_back->SetStats(0);
        canvas3a->SetLeftMargin(0.127517);
        canvas3a->SetRightMargin(0.02013423);

		sprintf(h_name,"lepton_numb_after_trig_mhc%i.eps", mhc);
		canvas3a->Print(h_name);
		*/
    }
	/////////////////////////////////////////////////////////////////////////
	//// COMPARISON FOR LEPTON PT
	/////////////////////////////////////////////////////////////////////////
 
	if (comparison_lepton_pt)
	{
		cout << " --> comparison_lepton_pt" << endl;
		
		TH1F * lepton_pt_sig  = new TH1F("lepton_pt_signal" ,"Lepton PT Distrib Sig;P_{T}^{l} GeV [GeV];Relative Occurence", 100, 0.0, 500.0);
		TH1F * lepton_pt_back = new TH1F("lepton_pt_back"	,"Lepton PT Distrib Back;P_{T}^{l} GeV [GeV];Relative Occurence", 100, 0.0, 500.0);
		// analysis over the signal
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/lepton_pt[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			lepton_pt[a]->Scale( scale );
			lepton_pt_sig->Add(lepton_pt[a]);
		}
		
		for (int a=4;a<number_of_datasets;a++)
		{
			if (lepton_pt[a]->GetEntries()==0)continue;
			double scale= (1.0/lepton_pt[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			lepton_pt[a]->Scale( scale );
			lepton_pt_back->Add(lepton_pt[a]);
		}
		TCanvas *canvas5a = new TCanvas("name5a","title",2500,0,600,500);
		lepton_pt_sig->GetYaxis()->SetRangeUser(0,0.30);
		lepton_pt_back->GetYaxis()->SetRangeUser(0,0.30);
		lepton_pt_sig->GetXaxis()->SetRangeUser(0,120);
		lepton_pt_back->GetXaxis()->SetRangeUser(0,120);
		lepton_pt_sig->SetLineColor(2);
		lepton_pt_sig->SetLineWidth(3);
		lepton_pt_back->SetLineWidth(3);
		lepton_pt_sig->Draw();
		lepton_pt_sig->SetFillColor(2);
		lepton_pt_sig->SetFillStyle(3002);
		lepton_pt_back->Draw("same");
		canvas5a->BuildLegend();
		//canvas5a->SetLogy();
		canvas5a->SetLeftMargin(0.127517);
		canvas5a->SetRightMargin(0.02013423);
		lepton_pt_sig->SetStats(0);
		lepton_pt_back->SetStats(0);
		
		sprintf(h_name,"3lepton_pt_signalback_mhc%i.eps", mhc);
		canvas5a->Print(h_name);
		
	}
	

/////////////////////////////////////////////////////////////////////////
// COMPARISON FOR MET 
/////////////////////////////////////////////////////////////////////////

    if (comparison_for_met)
    {
		
		cout << " --> comparison_for_met" << endl;

        TH1F * signal_met = new TH1F("signal_met" , "Missing Energy Sig;MET [GeV];Relative Occurence", 500, 0, 500);
        TH1F * back_met   = new TH1F("back_met"   , "Missing Energy Back;MET [GeV];Relative Occurence", 500, 0, 500);

        for (int a=0;a<4;a++)
        {
            double scale= (1.0/met_pt[a]->GetEntries())*(cross_sections[a])/total_cross_s;
            met_pt[a]->Scale( scale );
            signal_met->Add(met_pt[a]);
        }
        for (int a=4;a<number_of_datasets;a++)
        {
            if (met_pt[a]->GetEntries()==0)continue;
            double scale= (1.0/met_pt[a]->GetEntries())*(cross_sections[a])/total_cross_b;
            met_pt[a]->Scale( scale );
            back_met->Add(met_pt[a]);
        }
		
        TCanvas *canvas4a = new TCanvas("name4a","title",2500,0,600,500);
        signal_met->Rebin(5);
        back_met->Rebin(5);
        signal_met->SetLineColor(2);
        back_met->SetLineWidth(2);
        signal_met->GetYaxis()->SetRangeUser(0.00001,0.25);
        back_met->GetYaxis()->SetRangeUser(0.00001,0.25);
        signal_met->GetXaxis()->SetRangeUser(0.0,250);
        back_met->GetXaxis()->SetRangeUser(0.0,250);

        back_met->SetLineWidth(3);
        signal_met->SetLineWidth(3);
        signal_met->Draw();
		signal_met->SetFillColor(2);
		signal_met->SetFillStyle(3002);
		
        back_met->Draw("same");
        canvas4a->BuildLegend();
        //canvas4a->SetLogy();
        canvas4a->SetLeftMargin(0.127517);
        canvas4a->SetRightMargin(0.02013423);
        signal_met->SetStats(0);
        back_met->SetStats(0);
		sprintf(h_name,"4met_signalback_mhc%i.eps", mhc);
		canvas4a->Print(h_name);
	}

/////////////////////////////////////////////////////////////////////////
//// COMPARISON FOR ALPHA_T
/////////////////////////////////////////////////////////////////////////

	if (comparison_for_alphat)
	{
		
		cout << " --> comparison_for_alphat" << endl;

		TH1F * alpha_t_sig = new TH1F("alpha_t_sig" , "alpha_t Sig;alpha_t;Relative Occurence", 100, 0, 5);
		TH1F * alpha_t_back= new TH1F("alpha_t_back", "alpha_t Back;alpha_t;Relative Occurence", 100, 0, 5);
		
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/alpha_t[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			alpha_t[a]->Scale( scale );
			alpha_t_sig->Add(alpha_t[a]);
		}
		for (int a=4;a<number_of_datasets;a++)
		{
            if (alpha_t[a]->GetEntries()==0)continue;
			double scale= (1.0/alpha_t[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			alpha_t[a]->Scale( scale );
			alpha_t_back->Add(alpha_t[a]);
		}
		
		TCanvas *canvas4b = new TCanvas("name4b","title",2500,0,600,500);
//		alpha_t_sig->Rebin(5);
//		alpha_t_back->Rebin(5);
		alpha_t_sig->SetLineColor(2);
		alpha_t_back->SetLineWidth(2);
		alpha_t_sig->GetYaxis()->SetRangeUser(0.00001,0.15);
		alpha_t_back->GetYaxis()->SetRangeUser(0.00001,0.15);
		alpha_t_sig->GetXaxis()->SetRangeUser(0,2.5);
		alpha_t_back->GetXaxis()->SetRangeUser(0,2.5);
		alpha_t_back->SetLineWidth(3);
		alpha_t_sig->SetLineWidth(3);
		alpha_t_sig->Draw();
		alpha_t_sig->SetFillColor(2);
		alpha_t_sig->SetFillStyle(3002);
		
		alpha_t_back->Draw("same");
		canvas4b->BuildLegend();
		//canvas4b->SetLogy();
		canvas4b->SetLeftMargin(0.127517);
		canvas4b->SetRightMargin(0.02013423);
		alpha_t_sig->SetStats(0);
		alpha_t_back->SetStats(0);
		
		sprintf(h_name,"5alpha_t_mhc100%i.eps", mhc);
		canvas4b->Print(h_name);
	}

/////////////////////////////////////////////////////////////////////////
//// COMPARISON FOR BTAG
/////////////////////////////////////////////////////////////////////////

    if(comparison_for_btag)
    {
		
		cout << " --> comparison_for_btag" << endl;

		TH1F * btag_pt_sig  = new TH1F("btag_pt_sig"  ,"BTAG Jet PT Distrib Sig;P_{T}^{b-jet} [GeV];Relative Occurence", 100, 0, 500);
		TH1F * btag_pt_back = new TH1F("btag_pt_back" ,"BTAG Jet PT Distrib Back;P_{T}^{b-jet} [GeV];Relative Occurence", 100, 0, 500);
		// analysis over the signal
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/jet_btag_pt[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			jet_btag_pt[a]->Scale( scale );
			btag_pt_sig->Add(jet_btag_pt[a]);
		}
		for (int a=4;a<number_of_datasets;a++)
		{
            if (jet_btag_pt[a]->GetEntries()==0)continue;
            double scale= (1.0/jet_btag_pt[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			jet_btag_pt[a]->Scale( scale );
			btag_pt_back->Add(jet_btag_pt[a]);
		}
		TCanvas *canvas6a = new TCanvas("name6a","title",2500,0,600,500);
		btag_pt_sig->GetYaxis()->SetRangeUser(0.00001,0.15);
		btag_pt_back->GetYaxis()->SetRangeUser(0.00001,0.15);
        btag_pt_sig->GetXaxis()->SetRangeUser(0,250);
        btag_pt_back->GetXaxis()->SetRangeUser(0,250);
		btag_pt_sig->SetLineColor(2);
		btag_pt_sig->SetLineWidth(3);
		btag_pt_back->SetLineWidth(3);
		btag_pt_sig->Draw("");
		btag_pt_back->Draw("same");
		btag_pt_sig->SetFillColor(2);
		btag_pt_sig->SetFillStyle(3002);
		canvas6a->BuildLegend();
		//canvas6a->SetLogy();
		canvas6a->SetLeftMargin(0.127517);
		canvas6a->SetRightMargin(0.02013423);
		btag_pt_back->SetStats(0);
		btag_pt_sig->SetStats(0);

		sprintf(h_name,"6btag_pt_mhc%i.eps", mhc);
		canvas6a->Print(h_name);
	
		TH1F * btag_numb_sig  = new TH1F("btag_numb_sig" ,"BTAG Jet Numb Distrib Sig;Number of BTAG Jets;Relative Occurence", 10, 0, 10);
		TH1F * btag_numb_back = new TH1F("btag_numb_back","BTAG Jet Numb Distrib Back;Number of BTAG Jets;Relative Occurence", 10, 0, 10);
		// analysis over the signal
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/jet_btag[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			jet_btag[a]->Scale( scale );
			btag_numb_sig->Add(jet_btag[a]);
            //btag_numb_sig->Scale( 1.0/jet_btag[a]->Integral() );
		}
		for (int a=4;a<number_of_datasets;a++)
		{
            if (jet_btag[a]->GetEntries()==0)continue;
            double scale= (1.0/jet_btag[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			jet_btag[a]->Scale( scale );
			btag_numb_back->Add(jet_btag[a]);
            //btag_numb_back->Scale( 1.0/jet_btag[a]->Integral() );
		}

		TCanvas *canvas7a = new TCanvas("name7a","title",2500,0,600,500);
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
		sprintf(h_name,"6btag_numb_mhc%i.eps", mhc);
		canvas7a->Print(h_name);
	

		
		TH1F * btag_eta_sig  = new TH1F("btag_eta_sig" ,"BTAG Jet Eta Distrib Sig;ETA;Relative Occurence", 100, -5, +5);
		TH1F * btag_eta_back = new TH1F("btag_eta_back","BTAG Jet Eta Distrib Back;ETA;Relative Occurence", 100, -5, +5);

		// analysis over the signal
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/jet_btag_eta[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			jet_btag_eta[a]->Scale( scale );
			btag_eta_sig->Add(jet_btag_eta[a]);
		}
		for (int a=4;a<number_of_datasets;a++)
		{
            if (jet_btag_eta[a]->GetEntries()==0)continue;
            double scale= (1.0/jet_btag_eta[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			jet_btag_eta[a]->Scale( scale );
			btag_eta_back->Add(jet_btag_eta[a]);
		}
		
		TCanvas *canvas7b = new TCanvas("name7b","title",2500,0,600,500);
		btag_eta_sig->GetYaxis()->SetRangeUser(0.0,0.1);
		btag_eta_back->GetYaxis()->SetRangeUser(0.0,0.1);
		btag_eta_sig->Rebin(2);
		btag_eta_back->Rebin(2);
		btag_eta_sig->SetLineColor(2);
		btag_eta_sig->SetLineWidth(3);
		btag_eta_back->SetLineWidth(3);
		btag_eta_sig->Draw();
		btag_eta_back->Draw("same");
		btag_eta_sig->SetFillColor(2);
		btag_eta_sig->SetFillStyle(3002);
		canvas7b->BuildLegend();
		//canvas7a->SetLogy();
		canvas7b->SetLeftMargin(0.127517);
		canvas7b->SetRightMargin(0.02013423);
		btag_eta_back->SetStats(0);
		btag_eta_sig->SetStats(0);
		//btag_eta_back->GetXaxis()->SetRangeUser(0,5);
		//btag_eta_sig->GetXaxis() ->SetRangeUser(0,5);
		sprintf(h_name,"6btag_eta_mhc%i.eps", mhc);
		canvas7b->Print(h_name);

	}

	/////////////////////////////////////////////////////////////////////////
	//// LEPTON - BJET angle
	/////////////////////////////////////////////////////////////////////////
	
	if(comparison_for_bjet_lepton)
	{
		cout << " --> comparison_for_bjet_lepton" << endl;
		
		TCanvas *canvasBjet0 = new TCanvas("canvasBjet0","title",2500,0,1000,500);
		canvasBjet0->Divide(2);
		
		TH2F * bjet_lepton_dR_sig  = new TH2F("bjet_lepton_dR_sig" ,"Bjet lepton dR Sig;P_T^{lepton};Delta R", 100, 0, 250, 100, 0,  +10);
		TH2F * bjet_lepton_dR_back = new TH2F("bjet_lepton_dR_back","Bjet lepton dR Back;P_T^{lepton};Delta R", 100, 0, 250, 100, 0,  +10);
		// analysis over the signal
		
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/bjet_lepton_dR[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			bjet_lepton_dR[a]->Scale( scale );
			bjet_lepton_dR_sig->Add(bjet_lepton_dR[a]);
		}
		
		for (int a=4;a<number_of_datasets;a++)
		{
			if (bjet_lepton_dR[a]->GetEntries()==0)continue;
			double scale= (1.0/bjet_lepton_dR[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			bjet_lepton_dR[a]->Scale( scale );
			bjet_lepton_dR_back->Add(bjet_lepton_dR[a]);
		}
		canvasBjet0->cd(1);
		//        bjet_lepton_dR_sig ->GetYaxis()->SetRangeUser(0.00001,1);
		//        bjet_lepton_dR_sig ->GetXaxis()->SetRangeUser(0,6);
		//        bjet_lepton_dR_back->GetYaxis()->SetRangeUser(0.00001,1);
		//        bjet_lepton_dR_back->GetXaxis()->SetRangeUser(0,6);
		bjet_lepton_dR_sig->Draw("COLZ");
		gPad->SetLogz();
		canvasBjet0->cd(2);
		bjet_lepton_dR_back->Draw("COLZ");
		gPad->SetLogz();
		//bjet_lepton_dR_back->SetStats(0);
		//bjet_lepton_dR_sig->SetStats(0);
		
		//        canvasBjet0->SetLeftMargin(0.127517);
		//        canvasBjet0->SetRightMargin(0.02013423);
		sprintf(h_name,"7bjet_lepton_dR_mhc%i.C", mhc);
		canvasBjet0->Print(h_name);
		
		TH2F * bjet_lepton_dEta_sig  = new TH2F("bjet_lepton_dEta_sig" ,"Bjet lepton dEta Sig;P_T^{lepton};Delta Eta", 100, 0, 250, 100, 0,   +5);
		TH2F * bjet_lepton_dEta_back = new TH2F("bjet_lepton_dEta_back","Bjet lepton dEta Back;P_T^{lepton};Delta Eta", 100, 0, 250, 100, 0,   +5);
		
		
		// analysis over the signal
		
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/bjet_lepton_dEta[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			bjet_lepton_dEta[a]->Scale( scale );
			bjet_lepton_dEta_sig->Add(bjet_lepton_dEta[a]);
		}
		
		for (int a=4;a<number_of_datasets;a++)
		{
			if (bjet_lepton_dR[a]->GetEntries()==0)continue;
			double scale= (1.0/bjet_lepton_dEta[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			bjet_lepton_dEta[a]->Scale( scale );
			bjet_lepton_dEta_back->Add(bjet_lepton_dEta[a]);
		}
		
		canvasBjet0->cd(1);
		//        bjet_lepton_dR_sig ->GetYaxis()->SetRangeUser(0.00001,1);
		//        bjet_lepton_dR_sig ->GetXaxis()->SetRangeUser(0,6);
		//        bjet_lepton_dR_back->GetYaxis()->SetRangeUser(0.00001,1);
		//        bjet_lepton_dR_back->GetXaxis()->SetRangeUser(0,6);
		bjet_lepton_dEta_sig->Draw("COLZ");
		gPad->SetLogz();
		canvasBjet0->cd(2);
		bjet_lepton_dEta_back->Draw("COLZ");
		gPad->SetLogz();
		//bjet_lepton_dEta_back->SetStats(0);
		//bjet_lepton_dEta_sig->SetStats(0);
		
		//        canvasBjet0->SetLeftMargin(0.127517);
		//        canvasBjet0->SetRightMargin(0.02013423);
		sprintf(h_name,"7bjet_lepton_dEta_mhc%i.C", mhc);
		canvasBjet0->Print(h_name);
		
	}
	
/////////////////////////////////////////////////////////////////////////
//// LEPTON, TOP INVARIANT MASS, PZ
/////////////////////////////////////////////////////////////////////////
	
    if(comparison_for_invmass)
    {
		cout << " --> lepton_invmass_sig" << endl;

		TH1F * lepton_invmass_sig  = new TH1F("lepton_invmass_sig" ,
											  "Lepton Invariant Mass Sig;M_{T}^{l\nu};Relative Occurence", 100, 0, 500);
		TH1F * lepton_invmass_back  = new TH1F("lepton_invmass_back",
											   "Lepton Invariant Mass Back;M_{T}^{l\nu};Relative Occurence", 100, 0, 500);

		// analysis over the signal

		for (int a=0;a<4;a++)
		{
			if (lepton_invmass[a]->GetEntries()==0)continue;
			double scale= (1.0/lepton_invmass[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			lepton_invmass[a]->Scale( scale );
			lepton_invmass_sig->Add(lepton_invmass[a]);
		}

		for (int a=4;a<number_of_datasets;a++)
		{
            if (lepton_invmass[a]->GetEntries()==0)continue;
            double scale= (1.0/lepton_invmass[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			lepton_invmass[a]->Scale( scale );
			lepton_invmass_back->Add(lepton_invmass[a]);
		}

		TCanvas *canvas8a = new TCanvas("name8a","title",2500,0,600,500);
		lepton_invmass_sig ->GetYaxis()->SetRangeUser(0.00001,0.25);
		lepton_invmass_back->GetYaxis()->SetRangeUser(0.00001,0.25);
		lepton_invmass_sig ->GetXaxis()->SetRangeUser(0,250);
		lepton_invmass_back->GetXaxis()->SetRangeUser(0,250);
		lepton_invmass_sig->SetLineColor(2);
		lepton_invmass_sig->SetLineWidth(3);
		lepton_invmass_back->SetLineWidth(3);
		lepton_invmass_sig->Draw();
		lepton_invmass_sig->SetFillColor(2);
		lepton_invmass_sig->SetFillStyle(3002);
		lepton_invmass_back->Draw("same");
		canvas8a->BuildLegend();
		//canvas8a->SetLogy();
		canvas8a->SetLeftMargin(0.127517);
		canvas8a->SetRightMargin(0.02013423);
		lepton_invmass_back->SetStats(0);
		lepton_invmass_sig->SetStats(0);
		sprintf(h_name,"8lepton_invmass_mhc%i.eps", mhc);
		canvas8a->Print(h_name);
	
		TH1F * top_invmass_sig  = new TH1F("top_invmass_sig" ,"Top Invariant Mass Sig;M_{i}^{blv};Relative Occurence",  100, 0, 500);
		TH1F * top_invmass_back = new TH1F("top_invmass_back","Top Invariant Mass Back;M_{i}^{blv};Relative Occurence", 100, 0, 500);
		// analysis over the signal
		for (int a=0;a<4;a++)
		{
			if (top_invmass[a]->GetEntries()==0) continue;
			double scale= (1.0/top_invmass[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			top_invmass[a]->Scale( scale );
			top_invmass_sig->Add(top_invmass[a]);
		}

		for (int a=4;a<number_of_datasets;a++)
		{
            if (top_invmass[a]->GetEntries()==0)continue;
			double scale= (1.0/top_invmass[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			top_invmass[a]->Scale( scale );
			top_invmass_back->Add(top_invmass[a]);
		}
		TCanvas *canvas9a = new TCanvas("name9a","title",2500,0,600,500);
		top_invmass_sig->GetYaxis()->SetRangeUser(0,0.075);
		top_invmass_back->GetYaxis()->SetRangeUser(0,0.075);
		top_invmass_sig->GetXaxis()->SetRangeUser(0,500);
		top_invmass_back->GetXaxis()->SetRangeUser(0,500);
		top_invmass_sig->SetLineColor(2);
		top_invmass_sig->SetLineWidth(3);
		top_invmass_back->SetLineWidth(3);
		top_invmass_sig->Draw();
		top_invmass_sig->SetFillColor(2);
		top_invmass_sig->SetFillStyle(3002);
		top_invmass_back->Draw("same");
		canvas9a->BuildLegend();
		//canvas8a->SetLogy();
		canvas9a->SetLeftMargin(0.127517);
		canvas9a->SetRightMargin(0.02013423);
		top_invmass_back->SetStats(0);
		top_invmass_sig->SetStats(0);
		sprintf(h_name,"8top_invmass_mhc%i.eps", mhc);
		canvas9a->Print(h_name);

		TH1F * met_pz_sig  = new TH1F("met_pz_sig" ,"MET (Pz) Sig;Pz;Relative Occurence",  100, 0, 500);
		TH1F * met_pz_back = new TH1F("met_pz_back","MET (Pz) Back;Pz;Relative Occurence", 100, 0, 500);
		// analysis over the signal
		for (int a=0;a<4;a++)
		{
			if (miss_pz_mhc[a]->GetEntries()==0) continue;
			double scale= (1.0/miss_pz_mhc[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			miss_pz_mhc[a]->Scale( scale );
			met_pz_sig->Add(miss_pz_mhc[a]);
		}
		for (int a=4;a<number_of_datasets;a++)
		{
            if (miss_pz_mhc[a]->GetEntries()==0)continue;
			double scale= (1.0/miss_pz_mhc[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			miss_pz_mhc[a]->Scale( scale );
			met_pz_back->Add(miss_pz_mhc[a]);
		}

		TCanvas *canvas10a = new TCanvas("name10a","title",2500,0,600,500);
		met_pz_sig->GetYaxis()->SetRangeUser(0,0.05);
		met_pz_back->GetYaxis()->SetRangeUser(0,0.05);
		met_pz_sig->GetXaxis()->SetRangeUser(0,500);
		met_pz_back->GetXaxis()->SetRangeUser(0,500);
		met_pz_sig->SetLineColor(2);
		met_pz_sig->SetLineWidth(3);
		met_pz_back->SetLineWidth(3);
		met_pz_sig->Draw();
		met_pz_back->Draw("same");
		canvas10a->BuildLegend();
		//canvas8a->SetLogy();
		canvas10a->SetLeftMargin(0.127517);
		canvas10a->SetRightMargin(0.02013423);
		met_pz_back->SetStats(0);
		met_pz_sig->SetStats(0);
		sprintf(h_name,"8met_pz_mhc%i.eps", mhc);
		canvas10a->Print(h_name);
	}
	
	
	if(comparison_for_jetcut)
	{
		cout << " --> comparison_for_jetcut" << endl;

		TH1F * jet_size_cut8_sig  = new TH1F("jet_size_cut8_sig" ,"Jet Size before cut-8 Sig;Jet Size;Relative Occurence", 10, 0, 10);
		TH1F * jet_size_cut8_back  = new TH1F("jet_size_cut8_back","Jet Size before cut-8  Back;Jet Size;Relative Occurence", 10, 0, 10);
		
		// analysis over the signal
		
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/jet_size_cut8[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			jet_size_cut8[a]->Scale( scale );
			jet_size_cut8_sig->Add(jet_size_cut8[a]);
		}
		
		for (int a=4;a<number_of_datasets;a++)
		{
            if (jet_size_cut8[a]->GetEntries()==0)continue;
			double scale= (1.0/jet_size_cut8[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			jet_size_cut8[a]->Scale( scale );
			jet_size_cut8_back->Add(jet_size_cut8[a]);
		}
		
		TCanvas *canvas8a = new TCanvas("name8b","title",2500,0,600,500);
		jet_size_cut8_sig ->GetYaxis()->SetRangeUser(0.00001,1);
		jet_size_cut8_back->GetYaxis()->SetRangeUser(0.00001,1);
		jet_size_cut8_sig ->GetXaxis()->SetRangeUser(0,250);
		jet_size_cut8_back->GetXaxis()->SetRangeUser(0,250);
		jet_size_cut8_sig->SetLineColor(2);
		jet_size_cut8_sig->SetLineWidth(3);
		jet_size_cut8_back->SetLineWidth(3);
		jet_size_cut8_sig->Draw();
		jet_size_cut8_sig->SetFillColor(2);
		jet_size_cut8_sig->SetFillStyle(3002);
		jet_size_cut8_back->Draw("same");
		canvas8a->BuildLegend();
		//canvas8a->SetLogy();
		canvas8a->SetLeftMargin(0.127517);
		canvas8a->SetRightMargin(0.02013423);
		jet_size_cut8_back->SetStats(0);
		jet_size_cut8_sig->SetStats(0);
		sprintf(h_name,"9jet_size_cut8_mhc%i.eps", mhc);
		canvas8a->Print(h_name);
	
		TH1F * jet_size_cut9_sig  = new TH1F("jet_size_cut9_sig" ,"Jet Size before cut-9 Sig;Jet Size;Relative Occurence", 10, 0, 10);
		TH1F * jet_size_cut9_back  = new TH1F("jet_size_cut9_back","Jet Size before cut-9  Back;Jet Size;Relative Occurence", 10, 0, 10);
		
		// analysis over the signal
		
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/jet_size_cut9[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			jet_size_cut9[a]->Scale( scale );
			jet_size_cut9_sig->Add(jet_size_cut9[a]);
		}
		
		for (int a=4;a<number_of_datasets;a++)
		{
            if (jet_size_cut9[a]->GetEntries()==0)continue;
			double scale= (1.0/jet_size_cut9[a]->GetEntries())*(cross_sections[a])/total_cross_b;
			jet_size_cut9[a]->Scale( scale );
			jet_size_cut9_back->Add(jet_size_cut9[a]);
		}
		
		TCanvas *canvas9a = new TCanvas("name9b","title",2500,0,600,500);
		jet_size_cut9_sig ->GetYaxis()->SetRangeUser(0.00001,1);
		jet_size_cut9_back->GetYaxis()->SetRangeUser(0.00001,1);
		jet_size_cut9_sig ->GetXaxis()->SetRangeUser(0,6);
		jet_size_cut9_back->GetXaxis()->SetRangeUser(0,6);
		jet_size_cut9_sig->SetLineColor(2);
		jet_size_cut9_sig->SetLineWidth(3);
		jet_size_cut9_back->SetLineWidth(3);
		jet_size_cut9_sig->Draw();
		jet_size_cut9_sig->SetFillColor(2);
		jet_size_cut9_sig->SetFillStyle(3002);
		jet_size_cut9_back->Draw("same");
		canvas9a->BuildLegend();
		//canvas9a->SetLogy();
		canvas9a->SetLeftMargin(0.127517);
		canvas9a->SetRightMargin(0.02013423);
		jet_size_cut9_back->SetStats(0);
		jet_size_cut9_sig->SetStats(0);
		sprintf(h_name,"10jet_size_cut9_mhc%i.eps", mhc);
		canvas9a->Print(h_name);

		
		
		
		
		TH1F * jet_pt_sig  = new TH1F("jet_pt_sig" ,"Jet PT Sig;PT;Relative Occurence", 1000, 0.0, 1000.0);
		TH1F * jet_pt_back = new TH1F("jet_pt_back","Jet PT Back;PT;Relative Occurence", 1000, 0.0, 1000.0);

		// analysis over the signal
		
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/jet_pt0[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			jet_pt0[a]->Scale( scale );
			jet_pt_sig->Add(jet_pt0[a]);
		}
		
		for (int a=4;a<number_of_datasets;a++)
		{
            if (jet_pt0[a]->GetEntries()==0|| a==11)continue;
			double scale= (1.0/jet_pt0[a]->GetEntries())*(cross_sections[a])/total_cross_b;
            jet_pt0[a]->Scale( scale );
			jet_pt_back->Add(jet_pt0[a]);
		}
		
		TCanvas *canvas10c = new TCanvas("name10c","title",0,0,600,500);
		jet_pt_sig ->GetYaxis()->SetRangeUser(0.00001,1);
		jet_pt_back->GetYaxis()->SetRangeUser(0.00001,1);
//		jet_pt_sig ->GetXaxis()->SetRangeUser(0,6);
//		jet_pt_back->GetXaxis()->SetRangeUser(0,6);
		jet_pt_sig->Rebin(10);
		jet_pt_back->Rebin(10);
		jet_pt_sig->SetLineColor(2);
		jet_pt_sig->SetLineWidth(3);
		jet_pt_back->SetLineWidth(3);
		jet_pt_sig->Draw();
		jet_pt_sig->SetFillColor(2);
		jet_pt_sig->SetFillStyle(3002);
		jet_pt_back->Draw("same");
		canvas10c->BuildLegend();
		canvas10c->SetLogy();
		canvas10c->SetLeftMargin(0.127517);
		canvas9a->SetRightMargin(0.02013423);
		jet_pt_sig->SetStats(0);
		jet_pt_back->SetStats(0);
		sprintf(h_name,"11jet_pt_mhc%i.eps", mhc);
		canvas10c->Print(h_name);

		
		TH1F * jet_eta_sig  = new TH1F("jet_eta_sig" ,"Jet ETA Sig;ETA;Relative Occurence", 100, -5.0, 5.0);
		TH1F * jet_eta_back = new TH1F("jet_eta_back","Jet ETA Back;ETA;Relative Occurence", 100, -5.0, 5.0);

		// analysis over the signal
		
		for (int a=0;a<4;a++)
		{
			double scale= (1.0/jet_eta0[a]->GetEntries())*(cross_sections[a])/total_cross_s;
			jet_eta0[a]->Scale( scale );
			jet_eta_sig->Add(jet_eta0[a]);
		}
		
		for (int a=4;a<number_of_datasets;a++)
		{
            if (jet_eta0[a]->GetEntries()==0 || a==11)continue;
			double scale= (1.0/jet_eta0[a]->GetEntries())*(cross_sections[a])/total_cross_b;
            jet_eta0[a]->Scale( scale );
			jet_eta_back->Add(jet_eta0[a]);
		}
		
		TCanvas *canvas11c = new TCanvas("name11c","title",0,0,600,500);
		jet_eta_sig ->GetYaxis()->SetRangeUser(0.00001,0.05);
		jet_eta_back->GetYaxis()->SetRangeUser(0.00001,0.05);
		//		jet_pt_sig ->GetXaxis()->SetRangeUser(0,6);
		//		jet_pt_back->GetXaxis()->SetRangeUser(0,6);
		//jet_pt_sig->Rebin(10);
		//jet_pt_back->Rebin(10);
		jet_eta_sig->SetLineColor(2);
		jet_eta_sig->SetLineWidth(3);
		jet_eta_back->SetLineWidth(3);
		jet_eta_sig->Draw();
		jet_eta_sig->SetFillColor(2);
		jet_eta_sig->SetFillStyle(3002);
		jet_eta_back->Draw("same");
		canvas11c->BuildLegend();
		//canvas11c->SetLogy();
		canvas11c->SetLeftMargin(0.127517);
		canvas11c->SetRightMargin(0.02013423);
		jet_eta_sig->SetStats(0);
		jet_eta_back->SetStats(0);
		sprintf(h_name,"jet_eta_mhc%i.eps", mhc);
		canvas11c->Print(h_name);

	}
	
	
    if(comparison_for_jetsize)
    {
		cout << " --> comparison_for_jetcut" << endl;

        TH1F * jet_size_cut8_pt35_sig  = new TH1F("jet_size_cut8_pt35_sig" ,"Jet Size before cut-8 Sig;Jet Size;Relative Occurence", 10, 0, 10);
        TH1F * jet_size_cut8_pt35_back  = new TH1F("jet_size_cut8_pt35_back","Jet Size before cut-8  Back;Jet Size;Relative Occurence", 10, 0, 10);
        
        // analysis over the signal
        
        for (int a=0;a<4;a++)
        {
            double scale= (1.0/jet_size_cut9_50[a]->GetEntries())*(cross_sections[a])/total_cross_s;
            jet_size_cut9_50[a]->Scale( scale );
            jet_size_cut8_pt35_sig->Add(jet_size_cut9_50[a]);
        }
        
        for (int a=4;a<number_of_datasets;a++)
        {
            if (jet_size_cut9_50[a]->GetEntries()==0)continue;
            double scale= (1.0/jet_size_cut9_50[a]->GetEntries())*(cross_sections[a])/total_cross_b;
            jet_size_cut9_50[a]->Scale( scale );
            jet_size_cut8_pt35_back->Add(jet_size_cut9_50[a]);
        }

        TCanvas *canvas_jet = new TCanvas("canvas_jet","title",2500,0,600,500);
        jet_size_cut8_pt35_sig ->GetYaxis()->SetRangeUser(0.00001,1);
        jet_size_cut8_pt35_back->GetYaxis()->SetRangeUser(0.00001,1);
        jet_size_cut8_pt35_sig ->GetXaxis()->SetRangeUser(0,250);
        jet_size_cut8_pt35_back->GetXaxis()->SetRangeUser(0,250);
        jet_size_cut8_pt35_sig->SetLineColor(2);
        jet_size_cut8_pt35_sig->SetLineWidth(3);
        jet_size_cut8_pt35_back->SetLineWidth(3);
        jet_size_cut8_pt35_sig->Draw();
        jet_size_cut8_pt35_sig->SetFillColor(2);
        jet_size_cut8_pt35_sig->SetFillStyle(3002);
        jet_size_cut8_pt35_back->Draw("same");
        canvas_jet->BuildLegend();
        //canvas8a->SetLogy();
        canvas_jet->SetLeftMargin(0.127517);
        canvas_jet->SetRightMargin(0.02013423);
        jet_size_cut8_pt35_back->SetStats(0);
        jet_size_cut8_pt35_sig->SetStats(0);
        sprintf(h_name,"jet_size_cut9_pt50_mhc%i.eps", mhc);
        canvas_jet->Print(h_name);
        

    }
	
	

    
}