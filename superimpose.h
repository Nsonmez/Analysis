const char* data_name[] =
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
    // "wmp3jets",
    
    "wc0jets",
    "wc1jets",
    "wc2jets",
    // "wc3jets",

    "wbb0jets",
    "wbb1jets",
    "wbb2jets"
    // "wmp3jets"
    
};

// number of dataset
const int number_of_datasets = 20;

//Matched cross section
double cross_sections[] =
{
    
    //mhc100 icin degerler
    0.36*0.549,					// singletop_s
    0.36*15.514771802801450,	// singletop_s
    0.36*13.512667626306223,	// singletop_s
    0.36*0.96083416733012073,	// singletop_s
    
    
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
    //5800.3387351386527,   // wmp3jets
    
    9415.0009155022526,		// wc0jets
    7271.9538068456850,		// wc1jets
    3708.1308715138002,    	// wc2jets
    //1622.8438854931915,   // wc3jets
    
    110.84542630992861,		// wbb0jets
    214.67538706479982,		// wbb1jets
    168.06193240935943,    	// wbb2jets
    //114.9			    	// wbb3jets
};


	double total_cross_b =228613.1851;
	double total_cross_s =30.5372736*0.36; //mhc100 icin degerler

    TH2F * gen_lepton2D[number_of_datasets];
    

    TH1F * gen_lepton[number_of_datasets];
	TH1F * lepton_numb_before_trig[number_of_datasets];
    TH1F * lepton_numb_pass_trig [number_of_datasets];
    TH1F * lepton_before_hardphot [number_of_datasets];
	TH1F * lepton_numb_before_10gev[number_of_datasets];
	
    TH1F * lepton_pt[number_of_datasets];
	TH1F * met_pt[number_of_datasets];
	TH1F * alpha_t[number_of_datasets];
	
    TH1F * lepton_phi[number_of_datasets];
    TH1F * lepton_eta[number_of_datasets];
    TH1F * jet_pt0[number_of_datasets];
    TH1F * jet_pt1[number_of_datasets];
    TH1F * jet_pt2[number_of_datasets];

    TH1F * jet_eta0[number_of_datasets];
    TH1F * jet_eta1[number_of_datasets];
    TH1F * jet_eta2[number_of_datasets];


    TH1F * jet_eta[number_of_datasets][10];
	TH1F * jet_phi[number_of_datasets][10];
	TH1F * jet_size[number_of_datasets];
	TH1F * jet_size_cut8[number_of_datasets];
	TH1F * jet_size_cut9[number_of_datasets];
	TH1F * jet_btag[number_of_datasets];
	TH1F * jet_btag_pt[number_of_datasets];
	TH1F * jet_btag_eta[number_of_datasets];
	TH1F * lepton_invmass[number_of_datasets];
	TH1F * top_invmass[number_of_datasets];
	TH1F * miss_pz_mhc[number_of_datasets];


    TH1F * jet_size_cut8_35[number_of_datasets];
    TH1F * jet_size_cut8_40[number_of_datasets];
    TH1F * jet_size_cut8_45[number_of_datasets];
    TH1F * jet_size_cut8_50[number_of_datasets];

TH1F * jet_size_cut9_35[number_of_datasets];
TH1F * jet_size_cut9_40[number_of_datasets];
TH1F * jet_size_cut9_45[number_of_datasets];
TH1F * jet_size_cut9_50[number_of_datasets];

TH2F * bjet_lepton_dR[number_of_datasets];
TH2F * bjet_lepton_dEta[number_of_datasets];
TH2F * bjet_lepton_dPhi[number_of_datasets];

