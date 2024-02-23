R__LOAD_LIBRARY(/afs/cern.ch/user/c/cez/eos/Soft/mue/fairconvt/fair_to_mesmer/libMuEtree.so);
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


#define SAVE 1

const bool SMEARDATA = 1;
const int MINPOINTS = 5;
const int CHI2INTERVAL = 25;//means that I'm doing the parabolic interpolation around +/-5sigma
const int NORMALIZE = 1;


//SIMPLE CHI2 TEMPLATE FIT DELTA ALPHA LEP






void shape_chi2_MultFunc1Par_fitDALep_TOYS_v2(TString angular_cuts, double lumi = 1, int NTOYS = 1000, int iparam = 1211, bool randomSeed = 0) {
	/*lumi: pb-1*/

	if(randomSeed) gRandom->SetSeed(time(NULL));
	
	const double TRLumi    = lumi*1e6;//ub-1

	//some settings related to data selection
	Double_t NentriesCut = 0;
	Double_t RelativeUncertaintyCut = 0.1;
	Int_t    energyCut = 0;


	//determine the selection cuts on the two scattering angles.
	bool PID = 0;
	Double_t cutThL_low  = -10;
	Double_t cutThL_high = -10;
	Double_t cutThR_low  = -10;
	Double_t cutThR_high = -10;


	if(angular_cuts.SubString("thL") == "") {
		PID = 1;
		if(sscanf(angular_cuts, "thmu%lf_the%lf-%lf", &cutThL_low, &cutThR_low, &cutThR_high) == 2) {
			sscanf(angular_cuts, "thmu%lf_the%lf", &cutThL_low, &cutThR_high);
			cutThR_low  = 0;
		}
		cutThL_high = 6;
	}
	else {
		PID = 0;
		sscanf(angular_cuts, "thL%lf-%lf_thR%lf-%lf", &cutThL_low, &cutThL_high, &cutThR_low, &cutThR_high);
	}


//0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0                          MAKE HISTOGRAMS                 0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0	


	//TString mainDir = "/afs/cern.ch/user/r/rpilato/CMSSW_10_2_13/src/MUonE/AlphaLep_Sept23/GenerateSamples_Condor/job0/";
	//TString mainDir = "/eos/user/r/rpilato/Analysis_MuE/AlphaLep_Sept23/OutputHistograms/Test0/";
	TString mainDir = "/afs/cern.ch/user/c/cez/eos/Soft/mue/test/";
	TString fileLocation = "test_MuE_mesmer_events_exampleProductionJob_3cm_1708506156_1_40_results.root";
	//"/MuE/results.root";
	TString dirData = mainDir + fileLocation;//mainDir + Form("/iparam%i/statFit_DALep_Iparam%i_seed1-10/", iparam, iparam) + fileLocation;

	/*
	if(iparam == 1211) dirData = mainDir + Form("/iparam%i_shape/statFit_DALep_Iparam%i_shape_seed1-10/", 121, 121) + fileLocation;
	else if(iparam == 1221) dirData = mainDir + Form("/iparam%i_shape/statFit_DALep_Iparam%i_shape_seed1-10/", 122, 122) + fileLocation;
	*/


//	TString dirData = mainDir + Form("/statFit_DALep_Iparam%i_seed1-10/", iparam, iparam) + fileLocation;
	TFile *infile = new TFile(dirData);
	if(infile->IsOpen()) cout<<"input file opened successfully"<<endl;
	else { cout<<"where is the input file?"<<endl; return ; }	
	
	//get the reference histogram
	//directories in the results.root file where to find the 2D histograms needed for the template fit (DO NOT CHANGE THIS)
	TString dir_reference;
	TString dir_had;
	TString dir_templ;
	if(PID) {
		if(energyCut == 0) {
			dir_reference = "det/templ/NLO/full/hsn_thmuVsthe_NLO_ref";
			dir_had       = "det/templ/NLO/had/hsn_thmuVsthe_NLO_had";
			dir_templ     = "det/templ/NLO/lep/hsn_thmuVsthe_NLO_templ";
		}
		else {
			dir_reference = "det/templ/NLO/full/e20_hsn_thmuVsthe_NLO_ref";
			dir_had       = "det/templ/NLO/had/e20_hsn_thmuVsthe_NLO_had";
			dir_templ     = "det/templ/NLO/lep/e20_hsn_thmuVsthe_NLO_templ";
		}
	}
	else {
		if(energyCut == 0) {
			dir_reference = "det/templ/NLO/full/hsn_thL_vs_thR_NLO_ref";
			dir_had       = "det/templ/NLO/had/hsn_thL_vs_thR_NLO_had";
			dir_templ     = "det/templ/NLO/lep/hsn_thL_vs_thR_NLO_templ";
		}
		else {
			dir_reference = "det/templ/NLO/full/e20_hsn_thL_vs_thR_NLO_ref";
			dir_had       = "det/templ/NLO/had/e20_hsn_thL_vs_thR_NLO_had";
			dir_templ     = "det/templ/NLO/lep/e20_hsn_thL_vs_thR_NLO_templ";
		}
	}

	TTree *MuESetup = (TTree*) infile->Get("MuEsetup");
	MuE::Setup *MuEparams = 0;
	MuESetup->SetBranchAddress("MuEparams", &MuEparams);
	MuESetup->GetEntry(0);
	double Kref   = 2.3223e-3; //MuEparams->Kref;
	double Mref   = 511e-6; //MuEparams->Mref;
	double dKu    = 8e-5; //MuEparams->KerrRef;
	double dMu    = 0; //MuEparams->MerrRef;
	int sigmaLim  = 5; //MuEparams->rangeSigma;
	int sigmaStep = 4; //MuEparams->divSigma;
	int ngrid     = sigmaLim*sigmaStep*2 + 1;

	cout<<"***** Template fit parameters: *****"<<endl;
	cout<<"***** Kref = "<<Kref<<" +/- "<<dKu<<endl;
	cout<<"***** Mref = "<<Mref<<" +/- "<<dMu<<endl;
	cout<<"***** sigmaLim = "<<sigmaLim<<" sigmaStep = "<<sigmaStep<<endl;
	cout<<"***** iparam = "<<iparam<<endl;
	cout<<"************************************"<<endl;

	double FitPar, dFitPar;
	if(iparam == 121 || iparam == 1211) { FitPar = Kref; dFitPar = dKu; }
	else if(iparam == 122 || iparam == 1221) { FitPar = Mref; dFitPar = dMu; }

	//get the pseudodata histogram
	TH2D *hn_thmuVsthe_NLO_ref = (TH2D*) infile->Get(dir_reference);
	//get the hadronic histogram
	TH2D *hn_thmuVsthe_NLO_had = (TH2D*) infile->Get(dir_had);
	//get the templates
	std::vector<TH2D*> hn_thmuVsthe_NLO_templ;
	for (int iK = 0; iK < ngrid; iK++) {
		TString histo_templ = dir_templ + Form("_FitPar%i", iK);
		hn_thmuVsthe_NLO_templ.push_back((TH2D*) infile->Get(histo_templ));
		if(hn_thmuVsthe_NLO_templ.back() == nullptr) { cout<<"PROBLEM: template FitPar = "<<iK<<" does not exist"<<endl; return ;}
	}


	//set the angular cuts
	Int_t binThR_low   = std::max(hn_thmuVsthe_NLO_ref->ProjectionX()->FindBin(cutThR_low), 1);
	Int_t binThR_high  = std::min(hn_thmuVsthe_NLO_ref->ProjectionX()->FindBin(cutThR_high), hn_thmuVsthe_NLO_ref->GetNbinsX()+1);
	Int_t binThL_low   = std::max(hn_thmuVsthe_NLO_ref->ProjectionY()->FindBin(cutThL_low), 1);
	Int_t binThL_high  = std::min(hn_thmuVsthe_NLO_ref->ProjectionY()->FindBin(cutThL_high), hn_thmuVsthe_NLO_ref->GetNbinsY()+1);
	for(int i = 1; i <= hn_thmuVsthe_NLO_ref->GetNbinsX(); i++) {
		for(int j = 1; j <= hn_thmuVsthe_NLO_ref->GetNbinsY(); j++) {
			if((i < binThR_low || i >= binThR_high) || (j < binThL_low || j >= binThL_high)) {
				//data histogram
				hn_thmuVsthe_NLO_ref->SetBinContent(i, j, 0);
				hn_thmuVsthe_NLO_ref->SetBinError(i, j, 0);
				//hadronic histogram
				hn_thmuVsthe_NLO_had->SetBinContent(i, j, 0);
				hn_thmuVsthe_NLO_had->SetBinError(i, j, 0);

				//template histograms
				for(int iK = 0; iK < ngrid; iK++) {
					hn_thmuVsthe_NLO_templ[iK]->SetBinContent(i, j, 0);
					hn_thmuVsthe_NLO_templ[iK]->SetBinError(i, j, 0);
				}
			}
		}
	}

	//calculate integrals
	Double_t Integral_ref = hn_thmuVsthe_NLO_ref->Integral();
	Double_t Integral_had = hn_thmuVsthe_NLO_had->Integral();
	std::vector<Double_t> Integral_templ;
	for(int iK = 0; iK < ngrid; iK++) {
		Integral_templ.push_back(hn_thmuVsthe_NLO_templ[iK]->Integral());
	}

	//get Rhad
	TH2D *href = (TH2D*) hn_thmuVsthe_NLO_ref->Clone("href");
	if(NORMALIZE) href->Divide(hn_thmuVsthe_NLO_ref, hn_thmuVsthe_NLO_had, 1./Integral_ref, 1./Integral_had);
	else href->Divide(hn_thmuVsthe_NLO_ref, hn_thmuVsthe_NLO_had);


	//define the total number of bins for the combine histogram
	Int_t nbins_tot = hn_thmuVsthe_NLO_ref->GetNbinsX()*hn_thmuVsthe_NLO_ref->GetNbinsY();


	//set the errors on data_obs according to FullLumi
	Double_t totalEvents = 0;//ntot events from NLO file
	Double_t sigma0_tot  = 0;//ub (Wnorm from NLO file)
	for(int i = 0; i < MuESetup->GetEntries(); i++) {
		MuESetup->GetEntry(i);
		totalEvents += MuEparams->MCsums.Nevgen;//ntot events from NLO file
		sigma0_tot += MuEparams->MCpargen.Wnorm;//ub (Wnorm from NLO file)
		cout<<i<<") totEvents = "<<MuEparams->MCsums.Nevgen<<"   sigma0tot = "<<MuEparams->MCpargen.Wnorm<<endl;
	}
	sigma0_tot /= MuESetup->GetEntries();
	Double_t MCLumi      = totalEvents/sigma0_tot;//ub-1
	cout<<"MCLumi = "<<MCLumi<<endl;
	Double_t scaleFactor = TRLumi/MCLumi;

	//define the pseudodata 1D histogram which will be used by combine
	TH1D *data_obs_histo = new TH1D("data_obs", hn_thmuVsthe_NLO_ref->GetTitle(), nbins_tot, 0, nbins_tot);//pseudodata histogram
	//fill data_obs and set the errors
	int binCounter = 1;
	int expectedNDF = 0;
	for(int i = 0; i < href->GetNbinsX(); i++) {
		for(int j = 0; j < href->GetNbinsY(); j++) {
			Double_t Nevents = hn_thmuVsthe_NLO_ref->GetBinContent(i+1, j+1);
			Nevents = Nevents*scaleFactor;
			if(Nevents > NentriesCut) {
				double relativeError = 1./TMath::Sqrt(Nevents);
				if(relativeError <= RelativeUncertaintyCut) {
					data_obs_histo->SetBinContent(binCounter, href->GetBinContent(i+1, j+1));
					data_obs_histo->SetBinError(binCounter, relativeError*href->GetBinContent(i+1, j+1));
					expectedNDF++;
				}
				else {
					data_obs_histo->SetBinContent(binCounter, 0);
					data_obs_histo->SetBinError(binCounter, 0);
				}
			}
			else {
				data_obs_histo->SetBinContent(binCounter, 0);
				data_obs_histo->SetBinError(binCounter, 0);
			}
			binCounter++;
		}
	}

	//get Rhad for the templates
	std::vector<TH2D*> htempl_2D;
	for(int iK = 0; iK < ngrid; iK++) {
		htempl_2D.push_back((TH2D*) hn_thmuVsthe_NLO_templ[iK]->Clone(Form("htempl_2D_FitPar%i", iK)));
		if(NORMALIZE) htempl_2D.back()->Divide(hn_thmuVsthe_NLO_templ[iK], hn_thmuVsthe_NLO_had, 1./Integral_templ[iK], 1./Integral_had);
		else htempl_2D.back()->Divide(hn_thmuVsthe_NLO_templ[iK], hn_thmuVsthe_NLO_had);
	}

	std::vector<TH1D*> htempl;//template histograms with the nominal detector modelization

	//fill the templates and set the errors (...although template errors are neglected by combine...)
	for(int iK = 0; iK < ngrid; iK++) {
		htempl.push_back(new TH1D(Form("htempl_FitPar%i", iK), hn_thmuVsthe_NLO_templ[iK]->GetTitle(), nbins_tot, 0, nbins_tot));

		binCounter = 1;
		for(int i = 0; i < href->GetNbinsX(); i++) {
			for(int j = 0; j < href->GetNbinsY(); j++) {
				//template histograms with the nominal detector modelization
				Double_t Nevents = hn_thmuVsthe_NLO_templ[iK]->GetBinContent(i+1, j+1);
				Nevents = Nevents*scaleFactor;
				if(Nevents > NentriesCut) {
					double relativeError = 1./TMath::Sqrt(Nevents);
					if(relativeError <= 1) {
						htempl.back()->SetBinContent(binCounter, htempl_2D[iK]->GetBinContent(i+1, j+1));
						htempl.back()->SetBinError(binCounter, relativeError*htempl_2D[iK]->GetBinContent(i+1, j+1));
					}
					else {
						htempl.back()->SetBinContent(binCounter, 0);
						htempl.back()->SetBinError(binCounter, 0);
					}
				}
				else {
					htempl.back()->SetBinContent(binCounter, 0);
					htempl.back()->SetBinError(binCounter, 0);
				}
				binCounter++;
			}
		}
	}


//0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0     0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0     0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0










//0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0                 CALCULATE CHI2 AND FIT K, M                 0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0

	TStopwatch t;
	t.Start();
	
	TFile *outfile;
	TString outfilename = Form("shape_fitDALep_TRLumi%1.0fpb-1_Iparam%i_NTOYS%i_", lumi, iparam, NTOYS) + angular_cuts + "_v2.root";
	TH1D *hChi2min        = new TH1D("hChi2min", Form("#chi^{2}_{min} distribution; #chi^{2}; Entries"), 1000, 0, expectedNDF*1.5);
	TH1D *hChi2Reducedmin = new TH1D("hChi2Reducedmin", Form("reduced #chi^{2}_{min} distribution; #chi^{2}/ndf; Entries"), 1000, 0, 5);
	TH1D *hBestValues     = new TH1D("hBestValues", Form("FitParBest distribution (#chi^{2} minimzation of %i toys); FitPar", NTOYS),
	                             50, FitPar-4*dFitPar, FitPar+4*dFitPar);

        TH1D *h_chi21D = new TH1D("h_chi21D", "#chi^{2} distribution for the FitPar template grid; #frac{FitPar_{template} - FitPar_{ref}}{#sigma_{FitPar}}; #chi^{2}", ngrid, -sigmaLim, sigmaLim);
	TH1D *h_chi21D_FitPar  = new TH1D("h_chi21D_FitPar", "#chi^{2} distribution for the FitPar template grid; FitPar; #chi^{2}",
                            ngrid, FitPar - sigmaLim*dFitPar, FitPar + sigmaLim*dFitPar);

	TCanvas *cK;

	if(SAVE) outfile = new TFile(outfilename, "UPDATE");

	Double_t chi2 = 0;
	Double_t ndf  = 0;
	Double_t FitParBest = 0;
	TGraphErrors *results_FitPar;
	for(int toy = 0; toy < NTOYS; toy++) {

		results_FitPar = new TGraphErrors();
		results_FitPar->SetName("results_FitPar");
		results_FitPar->SetTitle("FitPar template fit; FitPar; #chi^{2}");
		results_FitPar->SetMarkerStyle(8);
		
		//smear the reference histogram
		TH1D *href_toy = (TH1D*) data_obs_histo->Clone("href_toy");
		if(SMEARDATA) {
			for(int i = 0; i < href_toy->GetNbinsX(); i++) {
				href_toy->SetBinContent(i+1, gRandom->Gaus(href_toy->GetBinContent(i+1), href_toy->GetBinError(i+1)));
			}
		}

		//do the template fit
		for(int iK = 0; iK < ngrid; iK++) {
			//get chi2
			chi2 = 0;
			ndf = 0;
			for(int i = 0; i < htempl[iK]->GetNbinsX(); i++) {
				Double_t refContent   = href_toy->GetBinContent(i+1);//already smeared (if SMEARDATA==1)
				Double_t refError     = href_toy->GetBinError(i+1);//already rescaled
				Double_t templContent = htempl[iK]->GetBinContent(i+1);
				Double_t error        = refError;
				if(error != 0) {
					chi2 += (refContent - templContent)*(refContent - templContent)/error/error;
					ndf++;
				}
			}

			ndf = ndf - 1;//1 free parameter
			double idK = iK - sigmaLim*sigmaStep;
			double FitPar_template = FitPar + idK*dFitPar/sigmaStep;
			results_FitPar->SetPoint(iK, FitPar_template, chi2);
			if(toy == NTOYS-1) {
				h_chi21D->SetBinContent(iK+1, chi2);
				h_chi21D_FitPar->SetBinContent(iK+1, chi2);
			}
		}	

		//find Kbest
		if(toy == NTOYS-1) cK = new TCanvas("cK", "Likelihood as a function of K", 1080, 720);
		TF1 *funK = new TF1("funK", "pol2", FitPar - sigmaLim*dFitPar, FitPar + sigmaLim*dFitPar);
		
		Double_t min_lnL = TMath::MinElement(results_FitPar->GetN(), results_FitPar->GetY());
		Int_t min_index = -1;
		for(int i = 0; i < results_FitPar->GetN(); i++) {
			if(results_FitPar->GetY()[i] == min_lnL) min_index = i;
		}

		Double_t fit_minlim = 0;
		for(int i = min_index; i >= 0; i--) if( (results_FitPar->GetY()[i] - min_lnL) <= CHI2INTERVAL) fit_minlim = results_FitPar->GetX()[i]*(1-0.01);
		Double_t fit_maxlim = 0;
		for(int i = min_index; i < results_FitPar->GetN(); i++) if( (results_FitPar->GetY()[i] - min_lnL) <= CHI2INTERVAL) fit_maxlim = results_FitPar->GetX()[i]*(1+0.01);

		results_FitPar->Fit("funK", "Q0", "", fit_minlim, fit_maxlim);
		if(toy == NTOYS-1) {
			results_FitPar->Draw("APE");
			funK->DrawF1(fit_minlim, fit_maxlim, "SAME");
			//results_FitPar->SetMinimum(-1);
		}
		FitParBest = funK->GetMinimumX();
		double chi2minimum = funK->Eval(FitParBest);

		hBestValues->Fill(FitParBest);
		hChi2min->Fill(chi2minimum);
		hChi2Reducedmin->Fill(chi2minimum/ndf);
		
		if(NTOYS >= 100 && toy%(NTOYS/100) == 0) cout<<toy<<") FitParBest = "<<FitParBest<<endl;

		if(toy != NTOYS-1) {
			href_toy->~TH1D();
			results_FitPar->~TGraphErrors();
			funK->~TF1();
		}
	}
	
	cout<<"Ntot template = "<<htempl.size()<<endl;
	cout<<"Ndf = "<<ndf<<endl;

	TCanvas *c0 = new TCanvas("c0", "", 1080, 720);
	c0->Draw();
	c0->SetLogz();
	href->Draw("ZCOL");
	

	TCanvas *c2 = new TCanvas("c2", "", 1080, 720);
	c2->Draw();
	hBestValues->Draw("ZCOL");
	TF1 *fGaus = new TF1("fGaus", "gaus", FitPar-4*dFitPar, FitPar+4*dFitPar);
	fGaus->SetNpx(1e3);
	fGaus->SetParameters(1, FitPar, dFitPar);
//	fGaus->SetParLimits(0, 0, 100);
	fGaus->SetParLimits(1, FitPar-10*dFitPar, FitPar+10*dFitPar);
	Double_t bestFitPar, errorFitPar;










	if(NTOYS >= 50) {
		TFitResultPtr fitBestValues = hBestValues->Fit("fGaus", "S");//"ME"
		bestFitPar  = fitBestValues->Parameter(1);
		errorFitPar = fitBestValues->Parameter(2);

		cout<<"\nFitParbest = "<<bestFitPar<<" +/- "<<errorFitPar<<endl;
	}


	t.Stop();
	t.Print();
	
	if(SAVE) {
		hn_thmuVsthe_NLO_ref->Write("", TObject::kOverwrite);
		hChi2min->Write("", TObject::kOverwrite);
		hChi2Reducedmin->Write("", TObject::kOverwrite);
		hBestValues->Write("", TObject::kOverwrite);
		h_chi21D->Write("", TObject::kOverwrite);
		h_chi21D_FitPar->Write("", TObject::kOverwrite);
		cK->Write("", TObject::kOverwrite);
		results_FitPar->Write("", TObject::kOverwrite);
		
		outfile->Close();
	}

return ;

}







