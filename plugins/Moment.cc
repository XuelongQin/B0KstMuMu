#include <TTree.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TMath.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistFunc.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TH3.h>
#include <TCanvas.h>

#include <sstream>
#include <time.h>
#include <ctime>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::vector;
using std::stringstream;
using std::string;
using namespace RooFit;
using namespace std;

// #############
// # Variables #
// #############
double q2Min = 0.;
double q2Max = 0.;
int q2Bin = -1;
double mumuMass2 = -1.;

double wight ;
double M_6s ;
double M_6c ;
double f_1s;
double f_3 ;
double f_4 ;
double f_5;
double f_6s ;
double f_6c;
double f_7 ;
double f_8 ;
double f_9 ;

double Mi ;
// ##############
// # Parameters #
// ##############
double value;
double error;
double Fl = 0.;
double AFB = 0.;
double S3 = 0.;
double S4 = 0.;
double S5 = 0.;
double S6 = 0.;
double S7 = 0.;
double S8 = 0.;
double S9 = 0.;

double P1 = 0.;
double P2 = 0.;
double P3 = 0.;
double P4p = 0.;
double P5p = 0.;
double P6p = 0.;
double P8p = 0.;

// ####################
// # angular varables #
// ####################
// #################
// # Reco-Level MC #
// #################
double cos_theta_k = 0;
TBranch *b_cos_theta_k = 0;

double cos_theta_l = 0;
TBranch *b_cos_theta_l = 0;

double phi_kst_mumu = 0;
TBranch *b_phi_kst_mumu = 0;

double mumuMass = 0;
TBranch *b_mumuMass = 0;

double tagB0 = 0;
TBranch *b_tagB0 = 0;

double genSignal = 0;
TBranch *b_genSignal = 0;


double recoB0Mass = 0;
TBranch *b_recoB0Mass = 0;

double PDGB0Mass = 5.27958;
double PDGJpsiMass = 3.096916;
double PDGPsiPrimeMass = 3.686109;


// ################
// # gen-Level MC #
// ################
//float cos_theta_k = 0;
//TBranch *b_cos_theta_k = 0;

//float cos_theta_l = 0;
//TBranch *b_cos_theta_l = 0;

//float phi_kst_mumu = 0;
//TBranch *b_phi_kst_mumu = 0;

double genQ = 0;
TBranch *b_genQ = 0;

void quzhi() {
	char a = '0' + q2Bin;
	switch (a) {
		case '0' :
			q2Min = 1.0;
			q2Max = 2.0;
			q2Bin = 0;
			break;

		case '1' :
			q2Min = 2.0;
			q2Max = 4.3;
			q2Bin = 1;
			break;

		case '2' :
			q2Min = 4.3;
			q2Max = 6.0;
			q2Bin = 2;
			break;

		case '3':
			q2Min = 6.0;
			q2Max = 8.68;
			q2Bin = 3;
			break;

			//        case '4':
			//            q2Min = 8.68;
			//            q2Max = 10.09;
			//            q2Bin = 4;
			//            break;

		case '5':
			q2Min = 10.09;
			q2Max = 12.86;
			q2Bin = 5;
			break;

			//        case '6':
			//            q2Min = 12.86;
			//            q2Max = 14.18;
			//            q2Bin = 6;
			//            break;

		case '7':
			q2Min = 14.18;
			q2Max = 16.0;
			q2Bin = 7;
			break;

			//        case '8':
			//            q2Min = 16.0;
			//            q2Max = 19.0;
			//            q2Bin = 8;
			//            break;

		default:
			break;
	}
}

void GenCalValue(int q2Bin, double q2Min, double q2Max, string Paratype) 
{
	value = 0.0;
	error= 0.0;
	M_6s = 0.;
	M_6c = 0.;
	f_1s = 0.;
	f_3 = 0.;
	f_4 = 0.;
	f_5 = 0.;
	f_6s = 0.;
	f_6c = 0.;
	f_7 = 0.;
	f_8 = 0.;
	f_9 = 0.;

	TH1D *f1s = new TH1D("f1s", "f1s", 100, 0.0, 1.0);
	TH1D *m6s = new TH1D("m6s", "m6s", 200, -5.0, 5.0);
	TH1D *m6c = new TH1D("m6c", "m6c", 200, -5.0, 5.0);
	TH1D *f3 = new TH1D("f3", "f3", 100, -5.0, 5.0);
	TH1D *f4 = new TH1D("f4", "f4", 200, -5.0, 5.0);
	TH1D *f5 = new TH1D("f5", "f5", 200, -5.0, 5.0);
	TH1D *f7 = new TH1D("f7", "f7", 200, -5.0, 5.0);
	TH1D *f8 = new TH1D("f8", "f8", 200, -5.0, 5.0);
	TH1D *f9 = new TH1D("f9", "f9", 200, -5.0, 5.0);
	// ############################
	// # Read tree and set Branch #
	// ############################
	//TFile *f = new TFile("/afs/cern.ch/user/l/llinwei/work2/qinxl/data/2016/skims/GEN/gen_B0_miniaodWithoutGenCuts.root");
	TFile *f = new TFile("/afs/cern.ch/user/x/xuqin/B0KstMuMu/fit/B0KstMuMu2/plugins/data/GEN_BFilter_B0MuMuKstar_p1.root");
	cout << "\n[Moment::GenCalValue]\tTry to open\t" << "GEN_BFilter_B0MuMuKstar_p1.root" << endl;
	//cout << "\n[Moment::GenCalValue]\tTry to open\t" << "gen_B0_miniaodWithoutGenCuts.root" << endl;

	TTree *t = (TTree *) f->Get("ntuple");
	t->SetBranchAddress("cos_theta_k", &cos_theta_k, &b_cos_theta_k);
	t->SetBranchAddress("cos_theta_l", &cos_theta_l, &b_cos_theta_l);
	t->SetBranchAddress("phi_kst_mumu", &phi_kst_mumu, &b_phi_kst_mumu);
	t->SetBranchAddress("genQ", &genQ, &b_genQ);

	Int_t entries = (Int_t) t->GetEntries();
	cout << "\n[Moment::GenCalValue]\tTotal number of events in the tree: " << entries << " @@@" << endl;
	for (Int_t i = 0; i < entries; i++)
	{
		t->GetEntry(i);
		mumuMass2 = genQ * genQ;

		// ###############################
		// # define orthogonal functions #
		// ###############################
		f_1s = 1 - cos_theta_k * cos_theta_k;
		M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
		M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
		f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
		f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
		// ##################################
		// # begin to compute the variables #
		// ##################################
		if (mumuMass2 > q2Min && mumuMass2 < q2Max) 
		{
			f1s->Fill(f_1s);
			m6s->Fill(M_6s);
			m6c->Fill(M_6c);
			f3->Fill(f_3);
			f4->Fill(f_4);
			f5->Fill(f_5);
			f7->Fill(f_7);
			f8->Fill(f_8);
			f9->Fill(f_9);       
		}
	} // end for 
	if  (Paratype == "FlS") 
	{
		value = 2.0 - 2.5 * f1s->GetMean();
		error = 2.5 * f1s->GetRMS() /TMath::Sqrt(f1s->GetEntries() - 1);
	}
	else if  (Paratype == "AFBS")
	{
		value = 3.0/ 4.0 * (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean());
		error = 3.0 / 4.0 * TMath::Sqrt(4.0 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) +
				9.0 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1));
	}
	else if  (Paratype == "S3S")
	{
		value = 25.0 / 8.0 * f3->GetMean();
		error = 25.0 / 8.0 * f3->GetRMS() /TMath::Sqrt(f3->GetEntries() - 1);
	}
	else if  (Paratype == "S4S")
	{
		value = 25.0 / 8.0 * f4->GetMean();
		error = 25.0 / 8.0 * f4->GetRMS() /TMath::Sqrt(f4->GetEntries() - 1);
	}
	else if (Paratype == "S5S")
	{
		value = 2.5 * f5->GetMean();
		error = 2.5 * f5->GetRMS() /TMath::Sqrt(f5->GetEntries() - 1);
	}
	else if (Paratype == "S7S")
	{
		value = 2.5 * f7->GetMean();
		error = 2.5 * f7->GetRMS() /TMath::Sqrt(f7->GetEntries() - 1);
	}
	else if  (Paratype == "S8S")
	{
		value = 25.0 / 8.0 * f8->GetMean();
		error = 25.0 / 8.0 * f8->GetRMS() /TMath::Sqrt(f8->GetEntries() - 1);
	}
	else if  (Paratype == "S9S")
	{
		value = 25.0 / 8.0 * f9->GetMean();
		error = 25.0 / 8.0 * f9->GetRMS() /TMath::Sqrt(f9->GetEntries() - 1);
	}
	else if  (Paratype == "P1S")
	{
		value = 25./4. * f3->GetMean() / ( 1. - (2.0 - 2.5 * f1s->GetMean()) );
		error = TMath::Sqrt(TMath::Power(2/(1.0 - (2.0 - 2.5 * f1s->GetMean())) * (25/8 * f3->GetRMS() /TMath::Sqrt(f3->GetEntries() - 1)), 2 ) + TMath::Power( (25/4 * f3->GetMean() / TMath::Power( (1.0 - (2.0 - 2.5 * f1s->GetMean())), 2) ) * (2.5 * f1s->GetRMS() /TMath::Sqrt(f1s->GetEntries() - 1)), 2 ) );
	}
	else if  (Paratype == "P2S")
	{
		value = 2.5*0.5 * m6s->GetMean() / ( 1. - (2.0 - 2.5 * f1s->GetMean()) );
		error = 2.5*TMath::Sqrt( TMath::Power( 1/(2 * (1- (2.0 - 2.5 * f1s->GetMean()))) ,2) * (4 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) + 9 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1)) + TMath::Power( (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean())/( 2 * TMath::Power( (1-(2.0 - 2.5 * f1s->GetMean())) ,2 )) ,2) * (TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1)) );
	}
	else if (Paratype == "P3S")
	{
		value = - 25.0/8.0 * f9->GetMean() / (1.0 - (2.0 - 2.5 * f1s->GetMean()));
		double a = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double b = TMath::Power( 25/8 * f9->GetRMS() ,2 ) / (f9->GetEntries() - 1);
		error = TMath::Sqrt(TMath::Power( 1.0 / (1.0 - (2.0 - 2.5 * f1s->GetMean())) ,2 ) * b + TMath::Power( 25/8 * f9->GetMean() / TMath::Power((1.0 - (2.0 - 2.5 * f1s->GetMean())) ,2) ,2) * a );
	}
	else if (Paratype == "P4pS")
	{
		value = 2*25/8 * f4->GetMean() /TMath::Sqrt((2.0 - 2.5 * f1s->GetMean()) * ( 1.0 - (2.0 - 2.5 * f1s->GetMean())) );
		double a = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double b = TMath::Power( 25/8 * f4->GetRMS() ,2 ) / (f4->GetEntries() - 1);
		error = 2*TMath::Sqrt( 1/ ((2.0 - 2.5 * f1s->GetMean()) * ( 1- (2.0 - 2.5 * f1s->GetMean()))) * b + TMath::Power((1 - 2 * (2.0 - 2.5 * f1s->GetMean())) * 25/8 * f4->GetMean() ,2) * a / ( 4 * TMath::Power((1 - (2.0 - 2.5 * f1s->GetMean())) ,3)) );
	}
	else if (Paratype == "P5pS")
	{
		double a = 2.0 - 2.5 * f1s->GetMean();
		value = 2.5 * f5->GetMean() /TMath::Sqrt(a * ( 1 - a) );
		double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double c = TMath::Power( 2.5 * f5->GetRMS() ,2 ) / (f5->GetEntries() - 1);
		error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 2.5 * f5->GetMean() ,2) * b /( 4 * TMath::Power((1 - a) ,3)) );
	}
	else if (Paratype == "P6pS")
	{
		double a = 2.0 - 2.5 * f1s->GetMean();
		value = - 2.5 * f7->GetMean() /TMath::Sqrt(a * ( 1 - a) );
		double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double c = TMath::Power( 2.5 * f7->GetRMS() ,2 ) / (f7->GetEntries() - 1);
		error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 2.5 * f7->GetMean() ,2) * b /(4 * TMath::Power((1 - a) ,3)) ); 
	}
	else if (Paratype == "P8pS")
	{
		double a = 2.0 - 2.5 * f1s->GetMean();
		value = 25/8 * f8->GetMean() /TMath::Sqrt(a * ( 1 - a) );
		double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double c = TMath::Power( 25/8 * f8->GetRMS() ,2 ) / (f8->GetEntries() - 1);
		error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 25/8 * f8->GetMean() ,2) * b / ( 4 * TMath::Power((1 - a) ,3)) );
	}
	cout << "\n[Moment::GenCalValue]\t q2bin = " << q2Bin << ":" << Paratype << " = " << value << "+/- " << error  << endl;
}

void ReCalValue(int q2Bin, double q2Min, double q2Max, string TagType, string Paratype) {

	cout << "TagType = " << TagType;
	cout << "q2Bin = " << q2Bin;
	cout << "q2Min = " << q2Min;
	cout << "q2Max = " << q2Max;
	wight = 0.0;
	M_6s = 0.;
	M_6c = 0.;
	f_1s = 0.;
	f_3 = 0.;
	f_4 = 0.;
	f_5 = 0.;
	f_6s = 0.;
	f_6c = 0.;
	f_7 = 0.;
	f_8 = 0.;
	f_9 = 0.;
	// ############################
	// # Read Effcicency Function #
	// ############################
	TH1D *f1s = new TH1D("f1s", "f1s", 100, 0.0, 1.0);
	TH1D *m6s = new TH1D("m6s", "m6s", 200, -5.0, 5.0);
	TH1D *m6c = new TH1D("m6c", "m6c", 200, -5.0, 5.0);
	TH1D *f3 = new TH1D("f3", "f3", 200, -5.0, 5.0);
	TH1D *f4 = new TH1D("f4", "f4", 200, -5.0, 5.0);
	TH1D *f5 = new TH1D("f5", "f5", 200, -5.0, 5.0);
	TH1D *f7 = new TH1D("f7", "f7", 200, -5.0, 5.0);
	TH1D *f8 = new TH1D("f8", "f8", 200, -5.0, 5.0);
	TH1D *f9 = new TH1D("f9", "f9", 200, -5.0, 5.0);

	int test = 4;  // you can choose the Maximum order
	stringstream myString;
	myString.clear();
	myString.str("");
	//myString << "/afs/cern.ch/user/l/llinwei/work2/qinxl/B0KstMuMufull/B0KstMuMu2016/efficiency/effProjection_sh" << test << "o_b" << q2Bin << "ct_25_25_25.root";
	//myString << "/afs/cern.ch/user/a/aboletti/public/Run2-KstarMuMu/KDEeff/latestVersion-integrals/KDEeff_b" << q2Bin << "_od.root";
	myString << "/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/fit/eff-KDE/files_forSimFit/KDEeff_b" << q2Bin << "_od_2017.root";
		cout << "\n[Moment::ReCalValue]\tTry to open" << myString.str().c_str() << endl;

	//    TFile *Eff = new TFile(myString.str().c_str(), "READ");
	//    RooWorkspace *w = (RooWorkspace *) Eff->Get("ws");
	//    RooAbsPdf *EffPdf = (RooAbsPdf *) w->function("projectedFunc");
	//    RooArgSet *param = EffPdf->getVariables();
	int parity = 1;
	TFile *efffile = new TFile(myString.str().c_str(), "READ");
	if (efffile!=0){
		cout<< "file is existed" << endl;
	}

	string effCstring = Form(("effCHist_b%ip%i"),q2Bin,parity); 
	TH3D *effCHist = (TH3D*)efffile->Get(effCstring.c_str()); 
	if (effCHist!=0){
		cout << "correctly tagged efficiency exists" << endl; 
	}

	string effWstring = Form(("effWHist_b%ip%i"),q2Bin,parity);
	TH3D *effWHist = (TH3D*)efffile->Get(effWstring.c_str());
	if (effWHist!=0){
		cout << "wrongly tagged efficiency exists" << endl;
	}

	//    TCanvas c1;
	//    effCHist->Draw();
	//    c1.SaveAs("effc.png") ;  
	RooRealVar ctK("ctK","ctK",-1,1);
	RooRealVar ctL("ctL","ctK",-1,1);
	RooRealVar phi("phi","phi",-TMath::Pi(),TMath::Pi());

	RooArgSet vars (ctK,ctL,phi);

	RooDataHist* effCData = new RooDataHist("effCData","effCData",vars,effCHist);
	RooAbsReal* effC = new RooHistFunc("effC","effC",vars,*effCData,1);

	RooDataHist* effWData = new RooDataHist("effWData","effWData",vars,effWHist);
	RooAbsReal* effW = new RooHistFunc("effW","effW",vars,*effWData,1);

	vector<double> intCVec (0);
	string intCstring = Form(("MCint_b%ip%it1"),q2Bin,parity);
	TH1D* intCHist = (TH1D*)efffile->Get(intCstring.c_str());
	if ( !intCHist || intCHist->IsZombie() ) {

		cout<<"Integral histograms not found in file: "<< myString.str().c_str() <<endl<<"Using rooFit integration"<<endl;
		intCVec.push_back(0);

	} else if ( strcmp( intCHist->GetTitle(), effCHist->GetTitle() )) {
		cout<<"Integral histogram is incoherent with efficiency "<<endl;
		cout<<"Correctly Efficiency conf: "<<effCHist->GetTitle()<<endl;
		cout<<"Wrongly Efficiency conf: "<<effWHist->GetTitle()<<endl;
		cout<<"Integral conf: "<<intCHist->GetTitle()<<endl;
		cout<<"Using rooFit integration"<<endl;
		intCVec.push_back(0);
	}

	else {
		for (int i=1; i<=intCHist->GetNbinsX(); ++i) intCVec.push_back(intCHist->GetBinContent(i));
	}



	// ##################
	// # Read Data & MC #
	// ##################
	cout << "\n[Moment::CalValue]\t @@@ Making datasets @@@ " << endl;
	//    TFile *f = new TFile("/afs/cern.ch/user/l/llinwei/work2/qinxl/data/2016/skims/2016MC_RECO_p1p2p3_newtag_LMNR_addW_add4BDT_addvars_bestBDTv4.root");
	TFile *f = new TFile("/afs/cern.ch/user/x/xuqin/data/2017/skims/2017MC_LMNR.root");
	TTree *t = (TTree *) f->Get("ntuple");
	t->SetBranchAddress("cos_theta_k", &cos_theta_k, &b_cos_theta_k);
	t->SetBranchAddress("cos_theta_l", &cos_theta_l, &b_cos_theta_l);
	t->SetBranchAddress("phi_kst_mumu", &phi_kst_mumu, &b_phi_kst_mumu);
	t->SetBranchAddress("mumuMass", &mumuMass, &b_mumuMass);
	t->SetBranchAddress("tagB0", &tagB0, &b_tagB0);
	t->SetBranchAddress("genSignal", &genSignal, &b_genSignal);
	t->SetBranchAddress( "tagged_mass", &recoB0Mass ,&b_recoB0Mass);
	int entries1 = 0;
	int entries2 = 0;
	int totentries = 0; 
	Int_t entries = (Int_t) t->GetEntries();
	cout << "\n[Moment::CalValue]\tTotal number of events in the tree: " << entries << " @@@" << endl;
	for (Int_t i = 0; i < entries; i++) {

		totentries = totentries+1;
		t->GetEntry(i);
		if ( mumuMass < PDGJpsiMass ) { // below Jpsi
			if ( fabs( recoB0Mass - PDGB0Mass - mumuMass + PDGJpsiMass ) < 0.18 ) continue;
		} else if ( mumuMass > PDGPsiPrimeMass ) { // above PsiPrime
			if ( fabs( recoB0Mass - PDGB0Mass - mumuMass + PDGPsiPrimeMass ) < 0.08 ) continue;
		} else { // between the resonances
			if ( fabs( recoB0Mass - PDGB0Mass - mumuMass + PDGJpsiMass ) < 0.08 ) continue;
			if ( fabs( recoB0Mass - PDGB0Mass - mumuMass + PDGPsiPrimeMass ) < 0.09 ) continue;
		}          
		//std::cout << "okay!" << std::endl;
		mumuMass2 = mumuMass*mumuMass;
		//	if ((tagB0 ==1 && genSignal ==1) || (tagB0 ==0 && genSignal ==2)){
		if (TagType == "good") {

			if (((tagB0 ==1 && genSignal ==1) || (tagB0 ==0 && genSignal ==2))&&(mumuMass2>q2Min) && (mumuMass2<q2Max)){
				ctK.setVal(cos_theta_k);
				ctL.setVal(cos_theta_l);
				phi.setVal(phi_kst_mumu);

				//std::cout << "okayhaha!" << std::endl;
				// ##################
				// # complete wight #
				// ##################
				//         effC->setRealValue("ctK", cos_theta_k);
				//         effC->setRealValue("ctL", cos_theta_l);
				//         effC->setRealValue("phi", phi_kst_mumu);

				// ###############################
				// # define orthogonal functions #
				// ###############################
				f_1s = 1 - cos_theta_k * cos_theta_k;
				M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
				M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
				f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
				f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
				f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
				f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
				f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
				f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
				double eff = effC->getValV();
				//  if (i%100000==0){
				//   cout << "process" << i << "  eff = " << eff << "  costhetak = " << cos_theta_k << "  costhetal = " << cos_theta_l << "  phi = " << phi_kst_mumu << endl;}
				//if (effC->getValV() < 0.0)  
				//cout  << "Eff = " << effC->getValV() << " " << cos_theta_k << "" << cos_theta_l << "" << phi_kst_mumu<< endl;
				wight = std::max(1 / eff,0.);
				//cout << "wight = " << wight << endl;

				f1s->Fill(f_1s, wight);
				m6s->Fill(M_6s, wight);
				m6c->Fill(M_6c, wight);
				f3->Fill(f_3, wight);
				f4->Fill(f_4, wight);
				f5->Fill(f_5, wight);
				f7->Fill(f_7, wight);
				f8->Fill(f_8, wight);
				f9->Fill(f_9, wight);
				entries1 = entries1+1;
			}
			else
				continue;
		}
		if (TagType == "mis") {
			if (((tagB0 ==0 && genSignal ==1) || (tagB0 ==1 && genSignal ==2))&&(mumuMass2>q2Min) && (mumuMass2<q2Max)){
		//		entries2 = entries2+1;
				ctK.setVal(cos_theta_k);
				ctL.setVal(cos_theta_l);
				phi.setVal(phi_kst_mumu);

				f_1s = 1 - cos_theta_k * cos_theta_k;
				M_6s = -(1 - cos_theta_k * cos_theta_k) * cos_theta_l;
				M_6c = -cos_theta_k * cos_theta_k * cos_theta_l;
				f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
				f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
				f_5 = -2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
				f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
				f_8 = -4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
				f_9 = -(1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
				double eff = effW->getValV();
				// if (i%100000==0){
				//   cout << "process" << i << "  eff = " << eff << "  costhetak = " << cos_theta_k << "  costhetal = " << cos_theta_l << "  phi = " << phi_kst_mumu << endl;}
				//           //            if (effC->getValV() < 0.0)  
				//           //            cout  << "Eff = " << effC->getValV() << " " << cos_theta_k << "" << cos_theta_l << "" << phi_kst_mumu<< endl;
				wight = std::max(1 / eff,0.);
				//                                           //cout << "wight = " << wight << endl;
				//
				f1s->Fill(f_1s, wight);
				m6s->Fill(M_6s, wight);
				m6c->Fill(M_6c, wight);
				f3->Fill(f_3, wight); 
				f4->Fill(f_4, wight);
				f5->Fill(f_5, wight);
				f7->Fill(f_7, wight);
				f8->Fill(f_8, wight);
				f9->Fill(f_9, wight);
				//
				entries2 = entries2+1;
			}
			else
				continue;
		}       

	}
	cout << "total entries is "<< totentries << endl;
	cout << "correctly-tagged events is " << entries1 << endl;
	cout << "wrongly-tagged events is " << entries2 << endl;
	if  (Paratype == "FlS")
	{
		cout << f1s->GetMean() << endl;
		value = 2.0 - 2.5 * f1s->GetMean();
		error = 2.5 * f1s->GetRMS() /TMath::Sqrt(f1s->GetEntries() - 1);
	}
	else if  (Paratype == "AFBS")
	{
		value = 3.0/ 4.0 * (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean());
		error = 3.0 / 4.0 * TMath::Sqrt(4.0 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) +
				9.0 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1));
	}
	else if  (Paratype == "S3S")
	{
		value = 25.0 / 8.0 * f3->GetMean();
		error = 25.0 / 8.0 * f3->GetRMS() /TMath::Sqrt(f3->GetEntries() - 1);
	}
	else if  (Paratype == "S4S")
	{
		value = 25.0 / 8.0 * f4->GetMean();
		error = 25.0 / 8.0 * f4->GetRMS() /TMath::Sqrt(f4->GetEntries() - 1);
	}
	else if (Paratype == "S5S")
	{
		value = 2.5 * f5->GetMean();
		error = 2.5 * f5->GetRMS() /TMath::Sqrt(f5->GetEntries() - 1);
	}
	else if (Paratype == "S7S")
	{
		value = 2.5 * f7->GetMean();
		error = 2.5 * f7->GetRMS() /TMath::Sqrt(f7->GetEntries() - 1);
	}
	else if  (Paratype == "S8S")
	{
		value = 25.0 / 8.0 * f8->GetMean();
		error = 25.0 / 8.0 * f8->GetRMS() /TMath::Sqrt(f8->GetEntries() - 1);
	}
	else if  (Paratype == "S9S")
	{
		value = 25.0 / 8.0 * f9->GetMean();
		error = 25.0 / 8.0 * f9->GetRMS() /TMath::Sqrt(f9->GetEntries() - 1);
	}
	else if  (Paratype == "P1S")
	{
		value = 25./4. * f3->GetMean() / ( 1. - (2.0 - 2.5 * f1s->GetMean()) );
		error = TMath::Sqrt(TMath::Power(2/(1.0 - (2.0 - 2.5 * f1s->GetMean())) * (25/8 * f3->GetRMS() /TMath::Sqrt(f3->GetEntries() - 1)), 2 ) + TMath::Power( (25/4 * f3->GetMean() / TMath::Power( (1.0 - (2.0 - 2.5 * f1s->GetMean())), 2) ) * (2.5 * f1s->GetRMS() /TMath::Sqrt(f1s->GetEntries() - 1)), 2 ) );
	}
	else if  (Paratype == "P2S")
	{
		value =2.5* 0.5 * m6s->GetMean() / ( 1. - (2.0 - 2.5 * f1s->GetMean()) );
		error = 2.5*TMath::Sqrt( TMath::Power( 1/(2 * (1- (2.0 - 2.5 * f1s->GetMean()))) ,2) * (4 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) + 9 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1)) + TMath::Power( (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean())/( 2 * TMath::Power( (1-(2.0 - 2.5 * f1s->GetMean())) ,2 )) ,2) * (TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1)) );
	}
	else if (Paratype == "P3S")
	{
		value = - 25.0/8.0 * f9->GetMean() / (1.0 - (2.0 - 2.5 * f1s->GetMean()));
		double a = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double b = TMath::Power( 25/8 * f9->GetRMS() ,2 ) / (f9->GetEntries() - 1);
		error = TMath::Sqrt(TMath::Power( 1.0 / (1.0 - (2.0 - 2.5 * f1s->GetMean())) ,2 ) * b + TMath::Power( 25/8 * f9->GetMean() / TMath::Power((1.0 - (2.0 - 2.5 * f1s->GetMean())) ,2) ,2) * a );
	}
	else if (Paratype == "P4pS")
	{
		value = 2 * 25/8 * f4->GetMean() /TMath::Sqrt((2.0 - 2.5 * f1s->GetMean()) * ( 1.0 - (2.0 - 2.5 * f1s->GetMean())) );
		double a = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double b = TMath::Power( 25/8 * f4->GetRMS() ,2 ) / (f4->GetEntries() - 1);
		error = 2*TMath::Sqrt( 1/ ((2.0 - 2.5 * f1s->GetMean()) * ( 1- (2.0 - 2.5 * f1s->GetMean()))) * b + TMath::Power((1 - 2 * (2.0 - 2.5 * f1s->GetMean())) * 25/8 * f4->GetMean() ,2) * a / ( 4 * TMath::Power((1 - (2.0 - 2.5 * f1s->GetMean())) ,3)) );

		cout << "P4 value is " << value << endl;
		cout << "P4 error is " << value << endl;
	}
	else if (Paratype == "P5pS")
	{
		double a = 2.0 - 2.5 * f1s->GetMean();
		value = 2.5 * f5->GetMean() /TMath::Sqrt(a * ( 1 - a) );
		double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double c = TMath::Power( 2.5 * f5->GetRMS() ,2 ) / (f5->GetEntries() - 1);
		error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 2.5 * f5->GetMean() ,2) * b /( 4 * TMath::Power((1 - a) ,3)) );
	}
	else if (Paratype == "P6pS")
	{
		double a = 2.0 - 2.5 * f1s->GetMean();
		value = - 2.5 * f7->GetMean() /TMath::Sqrt(a * ( 1 - a) );
		double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double c = TMath::Power( 2.5 * f7->GetRMS() ,2 ) / (f7->GetEntries() - 1);
		error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 2.5 * f7->GetMean() ,2) * b /(4 * TMath::Power((1 - a) ,3)) );      
	}
	else if (Paratype == "P8pS")
	{
		double a = 2.0 - 2.5 * f1s->GetMean();
		value = 25/8 * f8->GetMean() /TMath::Sqrt(a * ( 1 - a) );
		double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double c = TMath::Power( 25/8 * f8->GetRMS() ,2 ) / (f8->GetEntries() - 1);
		error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 25/8 * f8->GetMean() ,2) * b / ( 4 * TMath::Power((1 - a) ,3)) );

	}
	cout << "\n[Moment::RecoCalValue]\t q2bin = " << q2Bin << ":" << Paratype << " = " << value << "+/- " << error  << endl;


	}

	int main(int argc, char **argv) {

		string Paratype = "";
		string TagType = "";
		string SampleType = "";

		// ########################
		// # Give format of input #
		// ########################
		if (argc<3){
			cout << "\n[Moment::main]\tParameters missing\t" << endl;
			cout << "\n[Moment::main]\tPlease input correct format like following\t" << endl;
			cout << "\n./Moment\tParatype\tSampletype\tTagtype\t" << endl;
			cout << "\nParatype:FlS,AFBS,P1S,P2S,P3S,P4pS,P5pS,P6pS,P8pS,S3S,S4S,S5S,S6S,S7S,S8S,S9S\t" << endl;
			cout << "\nTagtype:mis,good" << endl;
			cout << "\nSampletype:gen,reco\t" << endl;
			return EXIT_FAILURE;
		}
		// ###################################
		// # Check that Parameter is correct #
		// ###################################
		Paratype = argv[1];
		if ((Paratype != "FlS") && (Paratype != "AFBS") && (Paratype != "P1S") && (Paratype != "P2S") &&
				(Paratype != "P3S") && (Paratype != "P4pS") && (Paratype != "P5pS") && (Paratype != "P6pS") &&
				(Paratype != "P8pS") && (Paratype != "S3S") && (Paratype != "S4S") && (Paratype != "S5S") &&
				(Paratype != "S6S") && (Paratype != "S7S") && (Paratype != "S8S") && (Paratype != "S9S")) {
			cout << "\n[Moment::main]\tIncorrect Parameters\t" << Paratype << endl;
			return EXIT_FAILURE;
		}
		cout << "\n[Moment::main]\tParameter Type = \t" << Paratype << endl;

		// ####################
		// # select Data & MC #
		// ####################
		SampleType = argv[2];
		if ((SampleType != "gen") && (SampleType != "reco")) {
			cout << "\n[Moment::main]\tIncorrect Sample Type\t" << SampleType << endl;
			return EXIT_FAILURE;
		}
		cout << "\n[Moment::main]\tSample Type = \t" << SampleType << endl;
		if (SampleType == "reco")
		{
			// ##################
			// # Check Tag Type #
			// ##################
			TagType = argv[3];
			if ((TagType != "good") && (TagType != "mis")) {
				cout << "\n[Moment::main]\tIncorrect Events Tag Type\t" << TagType << endl;
				return EXIT_FAILURE;
			}
			cout << "\n[Moment::main]\tEvents Tag Type = \t" << TagType << endl;
		}

		// ###########################
		// # calculate final results #
		// ###########################
		for (q2Bin = 0; q2Bin < 8; q2Bin++) {
			if (q2Bin == 4||q2Bin == 6){
				continue;
			}
			quzhi();
			if (SampleType == "gen") 
				GenCalValue(q2Bin, q2Min, q2Max, Paratype);
			else if (SampleType == "reco") 
				ReCalValue(q2Bin, q2Min, q2Max, TagType, Paratype);
		}
		return 0;
	}

