#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "TFile.h"
#include <TCanvas.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include "TH1.h"
#include "TF1.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TStyle.h"
#include <TH2.h>
#include <TH2F.h>
#include <TH3.h>
#include <TH3F.h>
#include "TF2.h"
#include "TGraph2D.h"
#include "TH1F.h"
#include "TMath.h"
#include "TFitter.h"
#include <string>
#include <iomanip>

#include "/user/cezhang/Soft/mue/code/MuEtree.h"//include the header file of the MuEtree class, which is needed to read some mesmer info (number of events, wnorm etc) from the results.root file

TRandom3 r(time(NULL));


//change the numbber of toys
const int  NTOYS = 1000;




//set the angular cuts. the default values correspond to theta_e < 32 mrad, theta_mu > 0.2 mrad
void fitTRLumi_TOYS_Cedric(Double_t cutThetae_low = 0, Double_t cutThetae_high = 32, Double_t cutThetamu_low = 0.2, Double_t cutThetamu_high = 6) {
	
	//set the input parameters for the template generation.
	//They are used to retrieve the value of K for the various templates starting from the template label.
	Double_t Kref = 0.13724;
	Double_t Mref = 0.0525;
	Double_t dKu = 0.045;
	Double_t dMu = 0.;
	Int_t sigmaLim = 5;
	Int_t sigmaStep = 4;
	Int_t ngrid = sigmaLim*sigmaStep*2 + 1;

	//path to the .root files containing the histograms in output from mue.
	TString job = "";//"job_20-10-14_19_41_54";
	TString dir = "results.root";//"~/Scrivania/MUonE/Analysis/mue/test/"+job+"/results.root";
	
	TFile *infile = new TFile(dir);
	if(infile->IsOpen()) cout<<"input file opened successfully"<<endl;
	else { cout<<"where is the input file?"<<endl; return ; }
	
	//info needed to calculate the MonteCarlo equivalent luminosity, and determine the scale factor
	//which is used to set the statistical errors according to the TR Luminosity.
	
	TTree *MuESetup = (TTree*) infile->Get("MuEsetup");
        MuE::Setup *MuEparams = 0;
        MuESetup->SetBranchAddress("MuEparams", &MuEparams);
    	TBranch *MCsums = MuESetup->GetBranch("MCsums");
    	TBranch *MCpargen = MuESetup->GetBranch("MCpargen");
    	MuESetup->GetEntry(0);

	Double_t totalEvents = MuEparams->MCsums.Nevgen;//ntot events from input file
	Double_t sigma0_tot  = MuEparams->MCpargen.Wnorm;//ub (Wnorm from NLO file)
	Double_t TRLumi      = 5e6;//ub-1. Test Run luminosity (5pb-1)
	Double_t MCLumi      = totalEvents/sigma0_tot;//ub-1

	//Test by Cedric 20 Dec 2022
	cout<<"MCLumi: "<<MCLumi<<" totalEvents: "<<totalEvents<<" sigma0_tot: "<<sigma0_tot<<endl;

	Double_t scaleFactor = TRLumi/MCLumi;

	//define the directories
	TString dir_reference;
	TString dir_lep;
	TString dir_templ;
	dir_reference = "det/templ/NLO/full/hsn_thmuVsthe_NLO_ref";
	dir_lep       = "det/templ/NLO/lep/hsn_thmuVsthe_NLO_lep";
	dir_templ     = "det/templ/NLO/had/hsn_thmuVsthe_NLO_templ";

	//get the reference histogram
	TH2D *hn_thmuVsthe_NLO_ref = (TH2D*) infile->Get(dir_reference);
	
	//get the template histograms
	std::vector<TH2D*> hn_thmuVsthe_NLO_templ;
	for(int iK = 0; iK < ngrid; iK++) {
		TString histo_templ = dir_templ + Form("_K%i", iK);
		hn_thmuVsthe_NLO_templ.push_back((TH2D*) infile->Get(histo_templ));
		if(hn_thmuVsthe_NLO_templ.back() == nullptr) { cout<<"PROBLEM: template (K) = ("<<iK<<") does not exist"<<endl; return ;}
	}

	//apply angular cuts
	Int_t binThetae_low    = std::max(hn_thmuVsthe_NLO_ref->ProjectionX()->FindBin(cutThetae_low), 1);
	Int_t binThetae_high   = std::min(hn_thmuVsthe_NLO_ref->ProjectionX()->FindBin(cutThetae_high), hn_thmuVsthe_NLO_ref->GetNbinsX()+1);
	Int_t binThetamu_low   = std::max(hn_thmuVsthe_NLO_ref->ProjectionY()->FindBin(cutThetamu_low), 1);
	Int_t binThetamu_high  = std::min(hn_thmuVsthe_NLO_ref->ProjectionY()->FindBin(cutThetamu_high), hn_thmuVsthe_NLO_ref->GetNbinsY()+1);
	for(int i = 1; i <= hn_thmuVsthe_NLO_ref->GetNbinsX(); i++) {
		for(int j = 1; j <= hn_thmuVsthe_NLO_ref->GetNbinsY(); j++) {
			if((i < binThetae_low || i >= binThetae_high) || (j < binThetamu_low || j >= binThetamu_high)) {
				//data histogram
				hn_thmuVsthe_NLO_ref->SetBinContent(i, j, 0);
				hn_thmuVsthe_NLO_ref->SetBinError(i, j, 0);

				//hn_thmuVsthe_NLO_ref->Draw();
				//c->Modified();
				//c->Update();
				//gSystem->Sleep(1);
				//if( gSystem->ProcessEvents()) break;

				//template histograms
				for(int iK = 0; iK < ngrid; iK++) {
					hn_thmuVsthe_NLO_templ[iK]->SetBinContent(i, j, 0);
					hn_thmuVsthe_NLO_templ[iK]->SetBinError(i, j, 0);
				}
			}
			
		}
	}



	//set the errors according to TRLumi
	TH2D *href = (TH2D*) hn_thmuVsthe_NLO_ref->Clone("href");
	for(int i = 0; i < href->GetNbinsX(); i++) {
		for(int j = 0; j < href->GetNbinsY(); j++) {
			Double_t Nevents = hn_thmuVsthe_NLO_ref->GetBinContent(i+1, j+1);
			//rescale bincontent to the TR statistics
			Nevents = Nevents*scaleFactor;

			Double_t relativeError = 1./TMath::Sqrt(Nevents);
									 
			if(Nevents != 0) {
				//accept only bins with a low statistical error (relativeError < 20%)
				if(relativeError < 0.2) {
					href->SetBinContent(i+1, j+1, Nevents);
					href->SetBinError(i+1, j+1, TMath::Sqrt(Nevents));
				}
				//else, do not take into account these bins in the fit
				else {
					href->SetBinContent(i+1, j+1, 0);
					href->SetBinError(i+1, j+1, 0);
				}
			}
			//if the bin does not have any event, set the error to zero (just to be safe...)
			else href->SetBinError(i+1, j+1, 0);
		}
	}
	
	//also rescale the templates to TRLumi
	std::vector<TH2D*> htempl;
	for(int iK = 0; iK < ngrid; iK++) {
		htempl.push_back((TH2D*) hn_thmuVsthe_NLO_templ[iK]->Clone(Form("htempl_K%i", iK)));

		for(int i = 0; i < htempl.back()->GetNbinsX(); i++) {
			for(int j = 0; j < htempl.back()->GetNbinsY(); j++) {
				Double_t Nevents = hn_thmuVsthe_NLO_templ[iK]->GetBinContent(i+1, j+1);
				//rescale bincontent to the TR statistics
				Nevents = Nevents*scaleFactor;

				Double_t relativeError = 1./TMath::Sqrt(Nevents);
										 
				if(Nevents != 0) {
					//accept only bins with a low statistical error (relativeError < 20%)
					if(relativeError < 0.2) {
						htempl.back()->SetBinContent(i+1, j+1, Nevents);
						htempl.back()->SetBinError(i+1, j+1, TMath::Sqrt(Nevents));
					}
					//else, do not take into account these bins in the fit
					else {
						htempl.back()->SetBinContent(i+1, j+1, 0);
						htempl.back()->SetBinError(i+1, j+1, 0);
					}
				}
				//if the bin does not have any event, set the error to zero (just to be safe...)
				else htempl.back()->SetBinError(i+1, j+1, 0);
			}
		}
	}

	

	
	//TOYS
	TGraphErrors *g[NTOYS];	
	TH1D *hKBest = new TH1D("hKBest", Form("K_{min} distribution (from #chi^{2} minimization of %1.0d toys); K; Entries ;", NTOYS), 
		       50, Kref - 5*dKu, Kref + 5*dKu);

	for(int toy = 0; toy < NTOYS; toy++) {
		
		TH2D *href_toy = (TH2D*) href->Clone("href_toy");
		//smear the reference histogram
		for(int i = 0; i < href_toy->GetNbinsX(); i++) {
			for(int j = 0; j < href_toy->GetNbinsY(); j++) {
				href_toy->SetBinContent(i+1, j+1, r.Gaus(href_toy->GetBinContent(i+1, j+1), href_toy->GetBinError(i+1, j+1)));
			}
		}


		g[toy] = new TGraphErrors();
		g[toy]->SetName(Form("g_toy%d", toy+1));
		g[toy]->SetTitle("#chi^{2} parabola; K; #chi^{2}");

		//build the chi2 parabola as a function of K for the different templates.	
		Double_t chi2 = 0;
		Double_t ndf = 0;
		
		for(int iK = 0; iK < ngrid; iK++) {
			chi2 = 0;
			ndf = 0;
			//calculate the chi2 for a given template.
			for(int i = 0; i < href_toy->GetNbinsX(); i++) {
				for(int j = 0; j < href_toy->GetNbinsY(); j++) {
					Double_t refContent   = href_toy->GetBinContent(i+1, j+1);
					Double_t refError     = href_toy->GetBinError(i+1, j+1);
					
					Double_t templContent = htempl[iK]->GetBinContent(i+1, j+1);
					
					Double_t error        = refError;//assume that template errors are negligible wrt pseudodata errors.
					
					if(error != 0) {
						chi2 += (refContent - templContent)*(refContent - templContent)/error/error;
						ndf++;
					}
				}
			}
			
			ndf = ndf - 1;//1 fit parameter (K)
			
			double idK = double(iK) - double(sigmaLim*sigmaStep); // position in the grid in units of one division 
			double K = Kref + idK*dKu/sigmaStep; // true values of grid points 

			g[toy]->SetPoint(iK, K, chi2);
		}


		//do the template fit: interpolate the chi2 parabola	
		TF1 *func = new TF1("func", "pol2", Kref-10*dKu, Kref+10*dKu);
		func->SetNpx(1e4);
		func->SetParameters(2.5, -1, 1);	
	
		g[toy]->Fit("func", "Q0");
	
		Double_t Kbest = func->GetMinimumX();

		hKBest->Fill(Kbest);

		if(toy%(NTOYS/100) == 0) cout<<"\r"<<toy<<flush;

		href_toy->~TH2D();

	}
	cout<<endl;
	
	TCanvas *c0 = new TCanvas("c0", "", 1080*1.3, 720*1.3);
	c0->SetGrid();
	c0->SetLogz();
	hn_thmuVsthe_NLO_ref->Draw("ZCOL");

	TCanvas *c1 = new TCanvas("c1", "", 1080, 720);
	c1->Draw();
	TFitResultPtr fit = hKBest->Fit("gaus", "MES");

	cout<<"Kbest = "<<fit->Parameter(1)<<" +/- "<<fit->Parameter(2)<<endl;



}










