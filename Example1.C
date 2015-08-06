/*
root -l examples/Example1.C\(\"delphes_output.root\"\)
*/

//------------------------------------------------------------------------------

void Example1(const char *inputFile)
{

  gSystem->Load("libDelphes");

   string OutputFileName="testing_delphes.root";
   TFile *outf = new TFile(OutputFileName.c_str(),"RECREATE");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");
  
  // Book histograms
  TH1 *histJetPT1 = new TH1F("jet_pt1", "1st jet P_{T}", 1000, 0.0, 1000.0);
  TH1 *histJetPT2 = new TH1F("jet_pt2", "2nd jet P_{T}", 500,  0.0, 500.0);
  TH1 *histJetPT3 = new TH1F("jet_pt3", "3rd jet P_{T}", 400, 0.0, 400.0);

  TH1 *histJetPT1_eta = new TH1F("jet_pt1_eta", "1st jet eta", 100, -5.0, 5.0);
  TH1 *histJetPT2_eta = new TH1F("jet_pt2_eta", "2nd jet eta", 100, -5.0, 5.0);

  TH1 *histJetPT1_phi = new TH1F("jet_pt1_phi", "1st jet phi", 100, -5.0, 5.0);
  TH1 *histJetPT2_phi = new TH1F("jet_pt2_phi", "2nd jet phi", 100, -5.0, 5.0);

  TH1 *histMass = new TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0);
  TH1 *jet_size = new TH1F("jet_size", "Number of Jets", 10, 0, 10.0);
  
  TH1 *elec_size = new TH1F("elec_size", "Number of Electrons", 10, 0, 10.0);
  TH1 *histElecPT1 = new TH1F("elec_pt1", "1st elec P_{T}", 100, 0.0, 500.0);
  TH1 *histElecPT2 = new TH1F("elec_pt2", "2nd elec P_{T}", 100,  0.0, 500.0);
  TH1 *histElecPT3 = new TH1F("elec_pt3", "3rd elec P_{T}", 100,  0.0, 500.0);

  TH1 *histElecPT1_phi = new TH1F("elec_pt1_phi", "1st elec Phi", 100, -5.0, 5.0);
  TH1 *histElecPT2_phi = new TH1F("elec_pt2_phi", "2nd elec Phi", 100, -5.0, 5.0);

  TH1 *histElecPT1_eta = new TH1F("elec_pt1_eta", "1st elec Eta", 100, -5.0, 5.0);
  TH1 *histElecPT2_eta = new TH1F("elec_pt2_eta", "2nd elec Eta", 100, -5.0, 5.0);

  TH1 *muon_size = new TH1F("muon_size", "Number of Muons", 10, 0, 10.0);

  TH1 *histMuonPT1 = new TH1F("muon_pt1", "1st mu P_{T}", 100, 0.0, 500.0);
  TH1 *histMuonPT2 = new TH1F("muon_pt2", "2nd mu P_{T}", 100,  0.0, 500.0);
  TH1 *histMuonPT3 = new TH1F("muon_pt3", "3rd mu P_{T}", 100,  0.0, 500.0);
  
  TH1 *histMuonPT1_phi = new TH1F("muon_pt1_phi", "1st mu Phi", 100, -5.0, 5.0);
  TH1 *histMuonPT2_phi = new TH1F("muon_pt2_phi", "2nd mu Phi", 100, -5.0, 5.0);

  TH1 *histMuonPT1_eta = new TH1F("muon_pt1_eta", "1st mu Eta", 100, -5.0, 5.0);
  TH1 *histMuonPT2_eta = new TH1F("muon_pt2_eta", "2nd mu Eta", 100, -5.0, 5.0);


  TH1 *histMET_et = new TH1F("histMET_et", "MET",  500,  0.0, 500.0);
  TH1 *histMET_eta = new TH1F("histMET_eta", "MET Eta", 100, -5.0, 5.0);
  TH1 *histMET_phi = new TH1F("histMET_phi", "histMET_phi Phi", 100, -5.0, 5.0);    
    
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // If event contains at least 1 jet
    if(branchJet->GetEntries() > 0)
    {
    
	jet_size->Fill( branchJet->GetEntries());
      	// Take first jet
      	Jet *jet1 = (Jet*) branchJet->At(0);
      	histJetPT1->Fill(jet1->PT);
      	histJetPT1_eta->Fill(jet1->Eta);
      	histJetPT1_phi->Fill(jet1->Phi);
      
    	if(branchJet->GetEntries() > 1)
    	{
      		Jet *jet2 = (Jet*) branchJet->At(1);
      		histJetPT2->Fill(jet2->PT);
      		histJetPT2_eta->Fill(jet2->Eta);
      		histJetPT2_phi->Fill(jet2->Phi);
    	}
    
     	if(branchJet->GetEntries() > 2)
    	{	
    	  Jet *jet3 = (Jet*) branchJet->At(2);
	jet3->PT

	}
    }



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	JET BRANCH
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
     	Jet *jet1, *jet2, *jet3; 
  // If event contains at least 1 jet
    if(branchJet->GetEntries() > 0)
    {
    
	jet_size->Fill( branchJet->GetEntries());
      	// Take first jet
      	jet1 = (Jet*) branchJet->At(0);
      	histJetPT1->Fill(jet1->PT);
      	histJetPT1_eta->Fill(jet1->Eta);
      	histJetPT1_phi->Fill(jet1->Phi);
      
    	if(branchJet->GetEntries() > 1)
    	{
      		jet2 = (Jet*) branchJet->At(1);
      		histJetPT2->Fill(jet2->PT);
      		histJetPT2_eta->Fill(jet2->Eta);
      		histJetPT2_phi->Fill(jet2->Phi);
    	}
    
     	if(branchJet->GetEntries() > 2)
    	{	
    	  jet3 = (Jet*) branchJet->At(2);
    	  histJetPT3->Fill(jet3->PT);
    	  //histJetPT3_eta->Fill(jet3->Eta);    
    	  //histJetPT3_phi->Fill(jet3->Phi);    
	}
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	ELECTRON BRANCH
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Electron *elec1, *elec2, *elec3;

	if(branchElectron->GetEntries() > 0)
    {
    
      elec_size->Fill( branchElectron->GetEntries());
  
      	elec1 = (Electron *) branchElectron->At(0);
      	histElecPT1->Fill(elec1->PT);
      	histElecPT1_eta->Fill(elec1->Eta);
      	histElecPT1_phi->Fill(elec1->Phi);
      
    	if(branchElectron->GetEntries() > 1)
    	{
      		elec2 = (Electron *) branchElectron->At(1);
      		histElecPT2->Fill(elec2->PT);
      		histElecPT2_eta->Fill(elec2->Eta);
      		histElecPT2_phi->Fill(elec2->Phi);
    	}
    
     	if(branchElectron->GetEntries() > 2)
    	{	
      		elec3 = (Electron *) branchElectron->At(2);
    	  	histElecPT3->Fill(elec3->PT);
    	  //histJetPT3_eta->Fill(jet3->Eta);    
    	  //histJetPT3_phi->Fill(jet3->Phi);    
   	}
    }
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	MUON BRANCH
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
	Muon *mu1, *mu2, *mu3;
	if(branchMuon->GetEntries() > 0)
	{
		muon_leaf_size=branchMuon->GetEntries();

	/*
		for (int i=0; i<muon_leaf_size; i++)
		{
		      	mu1 = (Muon *) branchMuon->At(i);
			histMuonPT[i]->Fill(mu1->PT);
		      	histMuonPT[i]_eta->Fill(mu1->Eta);
		      	histMuonPT[i]_phi->Fill(mu1->Phi);
		}
	*/

	muon_size->Fill( branchMuon->GetEntries());
  
      	mu1 = (Muon *) branchMuon->At(0);
      	histMuonPT1->Fill(mu1->PT);
      	histMuonPT1_eta->Fill(mu1->Eta);
      	histMuonPT1_phi->Fill(mu1->Phi);

    	if(branchMuon->GetEntries() > 1)
    	{
      		mu2 = (Muon *) branchMuon->At(1);
      		histMuonPT2->Fill(mu2->PT);
      		histMuonPT2_eta->Fill(mu2->Eta);
      		histMuonPT2_phi->Fill(mu2->Phi);
    	}

     	if(branchMuon->GetEntries() > 2)
    	{
      		mu3 = (Muon *) branchMuon->At(2);
    	  	histMuonPT3->Fill(mu3->PT);
    	  //histJetPT3_eta->Fill(jet3->Eta);    
    	  //histJetPT3_phi->Fill(jet3->Phi);    
   	}


    }
   
   
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MissingET BRANCH
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
    MissingET *met;

	double this_met=0;
	if(branchMET->GetEntries() > 0)
    {  
      	met = (MissingET *) branchMET->At(0);
      	histMET_et->Fill(met->MET);
      	histMET_eta->Fill(met->Eta);
      	histMET_phi->Fill(met->Phi);

	this_met=met->MET;
    }

	//if(this_met < 50) continue;
      
    // If event contains at least 2 electrons
    if(branchElectron->GetEntries() > 1)
    {
      // Take first two electrons
      elec1 = (Electron *) branchElectron->At(0);
      elec2 = (Electron *) branchElectron->At(1);
      // Plot their invariant mass
      histMass->Fill(((elec1->P4()) + (elec2->P4())).M());
    }
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*

if( met->MET>50 )
{
  if( branchMuon->GetEntries() == 1 || branchElectron->GetEntries() == 1 )
  {
     if( ( mu1->PT > 30 || elec1->PT>20 ) && (abs(mu1->Eta)<2.5 || abs(elec1->Eta)<2.5 ) )
     {
 	if( ( mu1->PT < 55&&elec1->PT<55 )
	{
	  if( jet1->PT <75 && jet1->BTag ==1 )
          {
	     if ( number_of_jets )
             {
		double mtlep=0;
	     if( ( mu1->PT > 30) && (abs(mu1->Eta)<2.5 ) ) mtlep=mtpar(mu1,met);
	     else if( ( elec1->PT>20 ) && (abs(elec1->Eta)<2.5 ) ) mtlep=mtpar(elec1,met);

		if( mtlep > 280 )
	        {


      	histJetPT1_eta->Fill(jet1->Eta);


}
}
}
}

*/

  // Show resulting histograms
	histJetPT1->Write();
	histJetPT2->Write();
	histJetPT3->Write();

	histJetPT1_eta->Write();
	histJetPT2_eta->Write();
	//histJetPT3_eta->Write();

	histJetPT1_phi->Write();
	histJetPT2_phi->Write();
	//histJetPT3_phi->Write();
	histMass->Write();
	jet_size->Write();

	elec_size ->Write();
	histElecPT1 ->Write();
	histElecPT2 ->Write();
	histElecPT3 ->Write();
	histElecPT1_phi ->Write();
	histElecPT2_phi ->Write();
	histElecPT1_eta ->Write();
	histElecPT2_eta ->Write();

	muon_size ->Write();
	histMuonPT1 ->Write();
	histMuonPT2 ->Write();
	histMuonPT3 ->Write();
	histMuonPT1_phi ->Write();
	histMuonPT2_phi ->Write();
	histMuonPT1_eta ->Write();
	histMuonPT2_eta ->Write();

	histMET_et->Write();
	histMET_eta->Write();
	histMET_phi->Write();


	outf->Close(); 

}

double mtpar(*ptlep, *miss)
{
double result=sqrt(2*(ptlep->PT)*(miss->PT)-2*((ptlep->PX)*(miss->PX)+(ptlep->PY)*(miss->PY) ));

return result;
}
