{
	gROOT->ProcessLine(".L /Users/nsonmez/root/macros/tdrStyle.C");
	gStyle->SetOptStat("n");
	
	setTDRStyle();
	char h_name[1000];
	/*
	 // number of dataset
	 int number_of_datasets = 14;
	 
	 double cross_sections[] =	{
		0.54920349716010086,	//singletop_s
		15.499481531905113,		//singletop_t2
		13.512667626306223,		//singletop_t3
		0.96083416733012073,	//singletop_tw
	 
		6.7532165336012762,		//back_singletop_s
		176.08374950986266,		//back_singletop_t2
		153.13877227715599,		//back_singletop_t3
		54.445885840029987,		//back_singletop_tw
	 
		6010.56,	//qcd1000-Inf
		8085,		//qcd500-1000
		8085,		//qcd250-500
		8085,		//qcd100-250
		558.58,					//ttjets
		18290,					//wmpjets
	 };
	 
	 */
	
	char* data_name[] =
	{
		"singletop_s",
		"singletop_t2",
		"singletop_t3",
		"singletop_tw",
		
		"back_singletop_s",
		"back_singletop_t2",
		"back_singletop_t3",
		"back_singletop_tw",
		
		"ttjets_lept",
		"ttjets_semilept",
		"ttjets_hadr",
		
		"wmp0jets",
		"wmp1jets",
		"wmp2jets",
		"wmp3jets",
		
		"wc0jets",
		"wc1jets",
		"wc2jets",
		//		"wc3jets",
		
		"wbb0jets",
		"wbb1jets",
		"wbb2jets"
		//		"wmp3jets"
		
	};
	
	// number of dataset
	int number_of_datasets = 21;
	
	//Matched cross section
	double cross_sections[] =
	{
		0.448,					// singletop_s
		13.11,					// singletop_t2
		11.38,					// singletop_t3
		0.812,					// singletop_tw
		
		6.7532165336012762,		//back_singletop_s
		176.08374950986266,		//back_singletop_t2
		153.13877227715599,		//back_singletop_t3
		54.445885840029987,		//back_singletop_tw
		
		107.7,					// ttjets_lept
		364.7,					// ttjets_semilept
		562.7,					// ttjets_hadr
		
		143787.04786199666,		// wmp0jets
		45457.421706146481,		// wmp1jets
		17028.775569156776,    	// wmp2jets
		5800.3387351386527,    	// wmp3jets
		
		9415.0009155022526,		// wc0jets
		7271.9538068456850,		// wc1jets
		3708.1308715138002,    	// wc2jets
								//1622.8438854931915,    // wc3jets
		
		110.84542630992861,		// wbb0jets
		214.67538706479982,		// wbb1jets
		168.06193240935943,    	// wbb2jets
		114.9			    	// wbb3jets
	};
	
	double total_cross_b =234413.1438;
	double total_cross_s =25.75;
	
	//______________________________________________________________
	// read the root files
	
	TFile * rootfile[23];
	for (int a=0;a<number_of_datasets;a++)
	{
		sprintf(h_name,"%s.root", data_name[a]);
		cout << h_name << endl;
		rootfile[a] = new TFile( h_name, "READ" );
	}
	
	//______________________________________________________________
	// get the lepton
	
	TH1F * lepton_pt[number_of_datasets];
	TH1F * met_pt[number_of_datasets];
	
	TH1F * lepton_phi[number_of_datasets];
	TH1F * lepton_eta[number_of_datasets];
	
	TH1F * lepton_numb_before_trig[number_of_datasets];
	TH1F * lepton_numb_pass_trig [number_of_datasets];
	TH1F * lepton_numb_pass_10gev[number_of_datasets];
	
	TH1F * jet_pt[number_of_datasets][10];
	TH1F * jet_eta[number_of_datasets][10];
	TH1F * jet_phi[number_of_datasets][10];
	TH1F * jet_size[number_of_datasets];
	TH1F * jet_btag[number_of_datasets];
	TH1F * jet_btag_pt[number_of_datasets];
	TH1F * lepton_invmass[number_of_datasets];
	TH1F * top_invmass[number_of_datasets];
	
	for (int a=0;a<number_of_datasets;a++)
	{
		lepton_eta[a] =  (TH1F*)rootfile[a]->Get("lepton/lepton_eta");
		lepton_phi[a] =  (TH1F*)rootfile[a]->Get("lepton/lepton_phi");
		lepton_numb_before_trig[a] = (TH1F*)rootfile[a]->Get("lepton/numb_lepton_before_trig");
		lepton_numb_pass_trig[a] = (TH1F*)rootfile[a]->Get("lepton/numb_lepton_pass_trig");
		lepton_numb_pass_10gev[a]  = (TH1F*)rootfile[a]->Get("lepton/numb_lepton_pass_10gev");
		
		jet_size[a] =  (TH1F*)rootfile[a]->Get("jets/jet_size");
		jet_btag[a] =  (TH1F*)rootfile[a]->Get("jets/histJet_btag");
		jet_btag_pt[a] =  (TH1F*)rootfile[a]->Get("jets/histJet_btag_pt");
		
		met_pt[a]    =  (TH1F*)rootfile[a]->Get("met/histMET_et");
		lepton_pt[a] =  (TH1F*)rootfile[a]->Get("lepton/lepton_pt");
		lepton_invmass[a] =  (TH1F*)rootfile[a]->Get("lepton/lepton_invmass");
		top_invmass[a] =  (TH1F*)rootfile[a]->Get("lepton/top_invmass");
		
	}
	
	//______________________________________________________________________
	// LEPTON spectrum for background and the signal
	
	TH1F * lepton_numb_before_trig_sig   = new TH1F("lepton_numb_before_trig_sig" ,"sig;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
	TH1F * lepton_numb_before_trig_back  = new TH1F("lepton_numb_before_trig_back" ,"back;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
	
	for (int a=0;a<4;a++)
	{
		double scale= (1.0/lepton_numb_before_trig[a]->GetEntries())*(cross_sections[a])/total_cross_s;
		lepton_numb_before_trig[a]->Scale( scale );
		lepton_numb_before_trig_sig->Add(lepton_numb_before_trig[a]);
	}
	for (int a=4;a<number_of_datasets;a++)
	{
		double scale= (1.0/lepton_numb_before_trig[a]->GetEntries())*(cross_sections[a])/total_cross_b;
		lepton_numb_before_trig[a]->Scale( scale );
		lepton_numb_before_trig_back->Add(lepton_numb_before_trig[a]);
	}
	TCanvas *canvas1a = new TCanvas("name1abc","title",2500,0,600,500);
	lepton_numb_before_trig_sig->Draw();
	lepton_numb_before_trig_sig->SetLineColor(2);
	lepton_numb_before_trig_back->Draw("same");
	lepton_numb_before_trig_sig->SetLineWidth(2);
	lepton_numb_before_trig_back->SetLineWidth(2);
	canvas1a->BuildLegend();
	//canvas1a->SetLogy();
	lepton_numb_before_trig_sig->SetStats(0);
	lepton_numb_before_trig_back->SetStats(0);
	canvas1a->SetLeftMargin(0.127517);
	canvas1a->SetRightMargin(0.02013423);
	canvas1a->Print("lepton_numb_before_trig.eps");
	
	//______________________________________________________________________
	
	TCanvas *canvas2a = new TCanvas("name2abc","title",2500,0,600,500);
	
	TH1F * signal_added_10gev_sig   = new TH1F("signal_added_10gev_sig" ,"sig;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
	TH1F * signal_added_10gev_back  = new TH1F("signal_added_10gev_back" ,"back;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
	
	for (int a=0;a<4;a++)
	{
		double scale= (1.0/lepton_numb_pass_10gev[a]->GetEntries())*(cross_sections[a])/total_cross_s;
		lepton_numb_pass_10gev[a]->Scale( scale );
		signal_added_10gev_sig->Add(lepton_numb_pass_10gev[a]);
	}
	for (int a=4;a<number_of_datasets;a++)
	{
		double scale= (1.0/lepton_numb_pass_10gev[a]->GetEntries())*(cross_sections[a])/total_cross_b;
		lepton_numb_pass_10gev[a]->Scale( scale );
		signal_added_10gev_back->Add(lepton_numb_pass_10gev[a]);
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
	canvas2a->Print("signal_added_10gev.eps");
	
	//______________________________________________________________________
	
	TCanvas *canvas3a = new TCanvas("name3abc","title",2500,0,600,500);
	
	TH1F * lepton_numb_after_trig_sig   = new TH1F("lepton_numb_after_trig_sig" ,"sig;Number of Leptons;Relative Occurence", 10, 0.0, 10.0);
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
	canvas3a->Print("lepton_numb_after_trig.eps");
	
	/////////////////////////////////////////////////////////////////////////
	// COMPARISON FOR MET
	/////////////////////////////////////////////////////////////////////////
	
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
		//met_pt[a]->Rebin(5);
		double scale= (1.0/met_pt[a]->GetEntries())*(cross_sections[a])/total_cross_b;
		met_pt[a]->Scale( scale );
		back_met->Add(met_pt[a]);
	}
	
	TCanvas *canvas4a = new TCanvas("name4a","title",2500,0,600,500);
	signal_met->Rebin(5);
	back_met->Rebin(5);
	signal_met->SetLineColor(2);
	signal_met->SetLineWidth(2);
	back_met->SetLineWidth(2);
	signal_met->GetYaxis()->SetRangeUser(0.00001,1);
	back_met->GetYaxis()->SetRangeUser(0.00001,1);
	back_met->SetLineWidth(2);
	signal_met->SetLineWidth(2);
	signal_met->Draw();
	back_met->Draw("same");
	canvas4a->BuildLegend();
	canvas4a->SetLogy();
	canvas4a->SetLeftMargin(0.127517);
	canvas4a->SetRightMargin(0.02013423);
	signal_met->SetStats(0);
	back_met->SetStats(0);
	canvas4a->Print("met_signalback.eps");
	
	
	/////////////////////////////////////////////////////////////////////////
	//// COMPARISON FOR LEPTON PT
	/////////////////////////////////////////////////////////////////////////
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
		double scale= (1.0/lepton_pt[a]->GetEntries())*(cross_sections[a])/total_cross_b;
		lepton_pt[a]->Scale( scale );
		lepton_pt_back->Add(lepton_pt[a]);
	}
	TCanvas *canvas5a = new TCanvas("name5a","title",2500,0,600,500);
	lepton_pt_sig->GetYaxis()->SetRangeUser(0,0.30);
	lepton_pt_back->GetYaxis()->SetRangeUser(0,0.30);
	lepton_pt_sig->GetXaxis()->SetRangeUser(0,200);
	lepton_pt_back->GetXaxis()->SetRangeUser(0,200);
	lepton_pt_sig->SetLineColor(2);
	lepton_pt_sig->SetLineWidth(2);
	lepton_pt_back->SetLineWidth(2);
	lepton_pt_sig->Draw();
	lepton_pt_back->Draw("same");
	canvas5a->BuildLegend();
	//canvas5a->SetLogy();
	canvas5a->SetLeftMargin(0.127517);
	canvas5a->SetRightMargin(0.02013423);
	lepton_pt_sig->SetStats(0);
	lepton_pt_back->SetStats(0);
	canvas5a->Print("lepton_signalback.eps");
	
	
	/////////////////////////////////////////////////////////////////////////
	//// COMPARISON FOR BTAG
	/////////////////////////////////////////////////////////////////////////
	
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
		double scale= (1.0/jet_btag_pt[a]->GetEntries())*(cross_sections[a])/total_cross_b;
		jet_btag_pt[a]->Scale( scale );
		btag_pt_back->Add(jet_btag_pt[a]);
	}
	TCanvas *canvas6a = new TCanvas("name6a","title",2500,0,600,500);
	btag_pt_sig->GetYaxis()->SetRangeUser(0.00001,1);
	btag_pt_back->GetYaxis()->SetRangeUser(0.00001,1);
	btag_pt_sig->SetLineColor(2);
	btag_pt_sig->SetLineWidth(2);
	btag_pt_back->SetLineWidth(2);
	btag_pt_sig->Draw();
	btag_pt_back->Draw("same");
	canvas6a->BuildLegend();
	canvas6a->SetLogy();
	canvas6a->SetLeftMargin(0.127517);
	canvas6a->SetRightMargin(0.02013423);
	btag_pt_back->SetStats(0);
	btag_pt_sig->SetStats(0);
	canvas6a->Print("btag_pt.eps");
	
	TH1F * btag_numb_sig  = new TH1F("btag_numb_sig" ,"BTAG Jet Numb Distrib Sig;Number of BTAG Jets;Relative Occurence", 10, 0, 10);
	TH1F * btag_numb_back = new TH1F("btag_numb_back","BTAG Jet Numb Distrib Back;Number of BTAG Jets;Relative Occurence", 10, 0, 10);
	// analysis over the signal
	for (int a=0;a<4;a++)
	{
		double scale= (1.0/jet_btag[a]->GetEntries())*(cross_sections[a])/total_cross_s;
		jet_btag[a]->Scale( scale );
		btag_numb_sig->Add(jet_btag[a]);
	}
	for (int a=4;a<number_of_datasets;a++)
	{
		double scale= (1.0/jet_btag[a]->GetEntries())*(cross_sections[a])/total_cross_b;
		jet_btag[a]->Scale( scale );
		btag_numb_back->Add(jet_btag[a]);
	}
	TCanvas *canvas7a = new TCanvas("name7a","title",2500,0,600,500);
	btag_numb_sig->GetYaxis()->SetRangeUser(0.00001,1);
	btag_numb_back->GetYaxis()->SetRangeUser(0.00001,1);
	btag_numb_sig->SetLineColor(2);
	btag_numb_sig->SetLineWidth(2);
	btag_numb_back->SetLineWidth(2);
	btag_numb_sig->Draw();
	btag_numb_back->Draw("same");
	canvas7a->BuildLegend();
	//canvas7a->SetLogy();
	canvas7a->SetLeftMargin(0.127517);
	canvas7a->SetRightMargin(0.02013423);
	btag_numb_back->SetStats(0);
	btag_numb_sig->SetStats(0);
	canvas7a->Print("btag_numb.eps");
	
	
	
	TH1F * lepton_invmass_sig  = new TH1F("lepton_invmass_sig" ,"Lepton Invariant Mass Sig;M_{T}^{l\nu};Relative Occurence", 100, 0, 500);
	TH1F * lepton_invmass_back  = new TH1F("lepton_invmass_back","Lepton Invariant Mass Back;M_{T}^{l\nu};Relative Occurence", 100, 0, 500);
	// analysis over the signal
	for (int a=0;a<4;a++)
	{
		double scale= (1.0/lepton_invmass[a]->GetEntries())*(cross_sections[a])/total_cross_s;
		lepton_invmass[a]->Scale( scale );
		lepton_invmass_sig->Add(lepton_invmass[a]);
	}
	for (int a=4;a<number_of_datasets;a++)
	{
		double scale= (1.0/lepton_invmass[a]->GetEntries())*(cross_sections[a])/total_cross_b;
		lepton_invmass[a]->Scale( scale );
		lepton_invmass_back->Add(lepton_invmass[a]);
	}
	TCanvas *canvas8a = new TCanvas("name8a","title",2500,0,600,500);
	lepton_invmass_sig->GetYaxis()->SetRangeUser(0.00001,1);
	lepton_invmass_back->GetYaxis()->SetRangeUser(0.00001,1);
	lepton_invmass_sig->SetLineColor(2);
	lepton_invmass_sig->SetLineWidth(2);
	lepton_invmass_back->SetLineWidth(2);
	lepton_invmass_sig->Draw();
	lepton_invmass_back->Draw("same");
	canvas8a->BuildLegend();
	canvas8a->SetLogy();
	canvas8a->SetLeftMargin(0.127517);
	canvas8a->SetRightMargin(0.02013423);
	lepton_invmass_back->SetStats(0);
	lepton_invmass_sig->SetStats(0);
	canvas8a->Print("lepton_invmass.eps");
	
	
	TH1F * top_invmass_sig  = new TH1F("top_invmass_sig" ,"Top Invariant Mass Sig;M_{T}^{bl\nu};Relative Occurence", 100, 0, 500);
	TH1F * top_invmass_back  = new TH1F("top_invmass_back","Top Invariant Mass Back;M_{T}^{bl\nu};Relative Occurence", 100, 0, 500);
	// analysis over the signal
	for (int a=0;a<4;a++)
	{
		double scale= (1.0/top_invmass[a]->GetEntries())*(cross_sections[a])/total_cross_s;
		top_invmass[a]->Scale( scale );
		top_invmass_sig->Add(top_invmass[a]);
	}
	for (int a=4;a<number_of_datasets;a++)
	{
		double scale= (1.0/top_invmass[a]->GetEntries())*(cross_sections[a])/total_cross_b;
		top_invmass[a]->Scale( scale );
		top_invmass_back->Add(top_invmass[a]);
	}
	TCanvas *canvas9a = new TCanvas("name9a","title",2500,0,600,500);
	top_invmass_sig->GetYaxis()->SetRangeUser(0,0.30);
	top_invmass_back->GetYaxis()->SetRangeUser(0,0.30);
	top_invmass_sig->GetXaxis()->SetRangeUser(0,500);
	top_invmass_back->GetXaxis()->SetRangeUser(0,500);
	top_invmass_sig->SetLineColor(2);
	top_invmass_sig->SetLineWidth(2);
	top_invmass_back->SetLineWidth(2);
	top_invmass_sig->Draw();
	top_invmass_back->Draw("same");
	canvas9a->BuildLegend();
	//canvas8a->SetLogy();
	canvas9a->SetLeftMargin(0.127517);
	canvas9a->SetRightMargin(0.02013423);
	top_invmass_back->SetStats(0);
	top_invmass_sig->SetStats(0);
	canvas9a->Print("top_invmass.eps");
	
	
	
	continue;
	
	
	
	
	for (int a=0;a<number_of_datasets;a++)
	{
		jet_btag_pt[a]->Scale((double)total_cross/cross_sections[a]);
		//elec_pt[a]->Rebin(10);
		jet_btag_pt[a]->GetYaxis()->SetRangeUser(0.00001,1);
		if (a==0) jet_btag_pt[a]->Draw("");
		else jet_btag_pt[a]->Draw("same");
		jet_btag_pt[a]->SetLineColor(a+1);
		sprintf(h_name,"%s;P_{T} [GeV];Btag-Jet PT Spectrum #\sigma/N",data_name[a]);
		jet_btag_pt[a]->SetTitle(h_name);
	}
	
	canvas3a->BuildLegend(0.7055058,0.3357271,0.9475032,0.8617594);
	canvas3a->SetLogy();
	canvas3a->Print("btag_pt.eps");
	
	for (int a=0;a<number_of_datasets;a++)
	{
		lepton_pt[a]->Scale((double)total_cross/(cross_sections[a]));
		//elec_pt[a]->Rebin(10);
		lepton_pt[a]->GetYaxis()->SetRangeUser(0.00001,1);
		if (a==0) lepton_pt[a]->Draw("");
		else lepton_pt[a]->Draw("same");
		lepton_pt[a]->SetLineColor(a+1);
		sprintf(h_name,"%s;P_{T} [GeV];Lepton Spectrum #\sigma/N",data_name[a]);
		lepton_pt[a]->SetTitle(h_name);
	}
	canvas2a->BuildLegend(0.7055058,0.3357271,0.9475032,0.8617594);
	canvas2a->SetLogy();
	canvas2a->Print("elec_pt.eps");
	
	TCanvas *canvas4a = new TCanvas("name4abc","title",2500,0,600,500);
	
	for (int a=0;a<number_of_datasets;a++)
	{
		jet_btag[a]->Scale((double)total_cross/cross_sections[a]);
		//elec_pt[a]->Rebin(10);
		jet_btag[a]->GetYaxis()->SetRangeUser(0.00001,1);
		if (a==0) jet_btag[a]->Draw("");
		else jet_btag[a]->Draw("same");
		jet_btag[a]->SetLineColor(a+1);
		sprintf(h_name,"%s;Number of BTAG ;Btag-Jet #\sigma/N",data_name[a]);
		jet_btag[a]->SetTitle(h_name);
	}
	
	canvas4a->BuildLegend(0.7055058,0.3357271,0.9475032,0.8617594);
	canvas4a->SetLogy();
	canvas4a->Print("btag.eps");
	
	
