/*
 * cutopt.C
 *
 *  Created on: Jan 12, 2018
 *  	Author: mcolak
 */

#include <stdio.h>
#include <iostream>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TString.h>
#include <math.h>
#include <cassert>
#include <cmath>
#include <TGraph2D.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TPaveLabel.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TImage.h>
#include <TVirtualPad.h>
#include <TLatex.h>
#include <fstream>
#include <TLine.h>

#include "MGeomCamMagicTwo.h"
#include "MEnergyEst.h"
#include "MStereoPar.h"
#include "MHadronness.h"
#include "MSrcPosCam.h"
#include "MMath.h"
#include "MPointingPos.h"
using namespace std;

void drawHistogram(TH2D* histogram[], Int_t bin, Int_t cycleno, Double_t energyLevels[], TImage *img, TString title) //
{
  TCanvas *canvasN = new TCanvas("Significances of Nebula","Graph2D example",0,0,1800,1200);
    //canvasN->Divide(1, bin, 0.000001, 0.0000001);
	canvasN->Divide(1, bin);
  for(Int_t i=0; i<bin; i++)
  {
    TString lastTitle= title + Form("for Cycle %d Energy between %f-%fGeV", cycleno, energyLevels[i], energyLevels[i+1]) ;
    cout << lastTitle << endl;
    gStyle->SetOptStat(0);
    canvasN->cd(i+1);
    histogram[i]->SetTitle(Form(lastTitle));
    histogram[i]->SetTitleSize(0.006);
    histogram[i]->GetYaxis()->SetTitle("Theta2");
    histogram[i]->GetYaxis()->CenterTitle();
    histogram[i]->GetYaxis()->SetTitleOffset(00.3);
    histogram[i]->GetYaxis()->SetTitleSize(0.005);
    histogram[i]->GetXaxis()->SetTitle("Hadronness");
    histogram[i]->GetXaxis()->SetTitleOffset(0.07);
    histogram[i]->GetXaxis()->SetTitleSize(0.005);
    histogram[i]->GetXaxis()->CenterTitle();
    histogram[i]->Draw("colz");
  }
  img->FromPad(canvasN);
  img->WriteImage(Form("SignificanceNebulaCycle%d.png",cycleno));
}

void cutopt()
{

	TString file="./20*.root";

	Int_t cycleno = 070535;
	Int_t theta2bin = 50 ;
	Double_t tmax = 0.16 ;
	Double_t tmin = 0.;
	const Int_t hadronnessbin = 50;
	Double_t hmin = 0. ;
	Double_t hmax = 1. ;

	MGeomCamMagicTwo geom;
	Float_t mm2deg = geom.GetConvMm2Deg();

	TChain* chain = new TChain("Events");
	//chain->Add(file56);
	//chain->Add(file71);
	//chain->Add(file72);
	//chain->Add(file81);
	//chain->Add(file82);
	chain->Add(file);


	chain->SetBranchStatus("*", 0); // Improve speed: de-activate reading of all branches

	// now activate the needed branches:
	chain->SetBranchStatus("MHadronness.fHadronness", 1);
	chain->SetBranchStatus("MEnergyEst.fEnergy", 1);
	//chain->SetBranchStatus("MEnergyEst_1.fEnergy", 1);  // M1 estimated energy
	//chain->SetBranchStatus("MEnergyEst_2.fEnergy", 1);  // M2 estimated energy
	chain->SetBranchStatus("MPointingPos_1.fZd", 1);	// M1 telescope pointing
	chain->SetBranchStatus("MStereoParDisp.*", 1);
	chain->SetBranchStatus("MSrcPosCam_1.fX", 1);  // To get the source position
	chain->SetBranchStatus("MSrcPosCam_1.fY", 1);  // To get the source position

	// Pointers for container we want to read:
	MHadronness   *had 	= 0;
	MEnergyEst	*eest	= 0;
	MEnergyEst	*eest1   = 0;
	MEnergyEst	*eest2   = 0;
	MPointingPos  *point1  = 0;
	MStereoPar	*stereo  = 0;
	MSrcPosCam	*srcpos1 = 0;

	// Associate pointers to containers inside the file:
	chain->SetBranchAddress("MHadronness.",	&had);
	chain->SetBranchAddress("MEnergyEst.", 	&eest);
	//chain->SetBranchAddress("MEnergyEst_1.",   &eest1);
	//chain->SetBranchAddress("MEnergyEst_2.",   &eest2);
	chain->SetBranchAddress("MStereoParDisp.", &stereo);
	chain->SetBranchAddress("MPointingPos_1.", &point1);
	chain->SetBranchAddress("MSrcPosCam_1.",   &srcpos1);

	Float_t nentries = chain->GetEntries();

	///////////// This is for test
	//
	//	nentries=1000.;
	//
	////////////////////////////////////

	Int_t P2 = 2;

	Double_t offno = 1.;

	const Int_t maxEnergybin = 30;
	Double_t energyLevels[maxEnergybin + 1]={};
	// Generating levels and mapping values based on the levels//////////
	Double_t minEnergy = 5.54;
	Double_t maxEnergy = 55400;
	Double_t step = (log10(maxEnergy) - log10(minEnergy))/maxEnergybin;
  cout << "step: " << step << endl;
	// Set up the levels
	energyLevels[0] = minEnergy;
	for(Int_t i=1; i<=maxEnergybin; i++)
  {
		energyLevels[i] = energyLevels[i-1] * pow(10,step);
    cout << "EnergyLevels: " << energyLevels[i] << endl;
	}
  	////////////Histograms are in here   set bin numbers wrt hadronnessbin and theta2bin!!!!!
	TH2D* histONe[maxEnergybin];
	TH2D* histOFFe[maxEnergybin];
	for (Int_t e=0 ; e< maxEnergybin ;e++)
	{
		histONe[e] = new TH2D(Form("histONe_%d",e),"", hadronnessbin, hmin, hmax, theta2bin, tmin, tmax);
		histOFFe[e] = new TH2D(Form("histOFF_e%d",e),"", hadronnessbin, hmin, hmax, theta2bin, tmin, tmax);
		//cout << "maxEnergybin loop1: " << e << endl;
	}

	for (Int_t ievent = 0; ievent < nentries; ievent++)
	{
		if (ievent % 1000000 == 0)
			cout << ievent << " / " << nentries << endl;

		chain->GetEvent(ievent); // This will retrieve all the info from event ievent

		// In order to see what "Getter" functions you have in each container, check the MARS documentation:
		// look for the Get* functions of a class, e.g. MEnergyEst::GetEnergy()

		////////////////////////////////////////////////////////////////////////////////////////////////////
		///************************
		/// DON'T FORGET TO SELECT ZENIT ANGLE AND SIZE CUT (IMPACT PARAMETER)

		// Example of impact parameter cuts
		//Float_t impact1 = 0.01*eest1->GetImpact();
		//Float_t impact2 = 0.01*eest2->GetImpact(); // Converted from cm to m
		//if (impact1 > 150 || impact2 > 150)
		//continue;

		// Get the event direction
		Float_t event_dirx = stereo->GetDirectionX(); // degree!
		Float_t event_diry = stereo->GetDirectionY(); // degree!
		//cout << event_dirx << " / " <<  event_diry << endl;

		// Get the nominal source position on the camera:
		Float_t srcx = srcpos1->GetX()*mm2deg;
		Float_t srcy = srcpos1->GetY()*mm2deg; //Converted to degrees
		//cout << srcx << " / " <<  srcy << endl;

		//ON events
		Float_t  ON_square = (event_dirx - srcx)*(event_dirx - srcx)+(event_diry - srcy)*(event_diry - srcy);
		//cout << "ON_square" << ON_square << endl;

		// The "reflected" OFF would be at -srcx, -srcy
		Float_t  OFF1_square = (event_dirx + srcx)*(event_dirx + srcx)+(event_diry + srcy)*(event_diry + srcy);
		//cout << "ON_square vs. OFF_square" << ON_square << " / " <<  OFF1_square << endl;

		//   	 // The OFF2 would be at -srcx, +srcy
		//   	 Float_t  OFF2_square = (event_dirx + srcx)*(event_dirx + srcx)+(event_diry - srcy)*(event_diry - srcy);
		//
		//   	 // The "reflected" OFF would be at +srcx, -srcy
		//   	 Float_t  OFF3_square = (event_dirx - srcx)*(event_dirx - srcx)+(event_diry + srcy)*(event_diry + srcy);

		Float_t hadronness = had->GetHadronness();
		//cout << "hadronness" << hadronness <<endl;

		////////////////////////////////////////////////////////////
		//
		//ENERGY AND ZD SELECTION:

		// Example of a cut in zenith distance. Skip events below 10 or above 40deg:
		Float_t zd = point1->GetZd();
		//Example of energy estimation parameter cuts
		Float_t estimated_energy = eest->GetEnergy();

		//cout << "Energy1: " << estimated_energy << "ZD1: " << zd << endl;

		if (zd<10)
			continue;
		if (zd>30)
			continue;

		// Check which level this estimated_energy falls within
		for(Int_t i=0; i<maxEnergybin; i++){
			if(estimated_energy >= energyLevels[i] && estimated_energy < energyLevels[i+1])
			{
				histONe[i]->Fill(hadronness, ON_square);
				histOFFe[i]->Fill(hadronness, OFF1_square);
			}
		}
		//////////////////////////////////////////////////////////////////////////////////

	}
	///////////////////Summing up Hadronnesses (Hadronness Cut)

	for (Int_t i = 1; i <=histONe[1]->GetNbinsY() ; i++) // loop the Theta2
	{
		for (Int_t j = 2; j <= histONe[1]->GetNbinsX(); j++) // loop the Hadronness
		{
			for(Int_t e=0 ; e<maxEnergybin; e++) // Loop the Energy
			{
				///for Energy E1
				histONe[e]->SetBinContent(j, i, histONe[e]->GetBinContent(j,i)+ histONe[e]->GetBinContent(j-1,i));
				histONe[e]->SetBinError(j, i, sqrt(histONe[e]->GetBinContent(j,i)));

				histOFFe[e]->SetBinContent(j, i, histOFFe[e]->GetBinContent(j,i)+ histOFFe[e]->GetBinContent(j-1,i));
				histOFFe[e]->SetBinError(j, i, sqrt(histOFFe[e]->GetBinContent(j,i)));
			}
		}
	}

	TFile* foute[maxEnergybin];
	TH2D* ONClone[maxEnergybin];
	TH2D* OFFClone[maxEnergybin];
	TH2D* OFFCloneSmooth[maxEnergybin];

	for (Int_t e=0 ; e< maxEnergybin ;e++)
	{
		foute[e] = new TFile(Form("testoute_%d",e),"recreate");
		histONe[e]->Write();
		histOFFe[e]->Write();
		foute[e]->Close();


		ONClone[e] = (TH2D*) histONe[e]->Clone(Form("ON_e%d",e));
		OFFClone[e] = (TH2D*) histOFFe[e]->Clone(Form("OFF_e%d",e));
		OFFCloneSmooth[e] = (TH2D*) histOFFe[e]->Clone(Form("OFFS_e%d",e));
		//cout << "maxEnergybin loop1: " << e << endl;
	}
	//////////////NORMALIZATION FACTOR
	// Do this later &&&&& EXCESS & OFF without any cut -- For Efficiency calculations
	Double_t xminN = 0.15 ;
	Double_t xmaxN = 0.4 ;

	////****************To check
	//cout << "minBinN: " << minBinN << endl;
	//cout << "maxBinN: " << maxBinN << endl;
	////*****************

	//HISTOGRAMS
	//Energy:
	//
	TH2D* PulsarSigP2[maxEnergybin];
	TH2D* NebulaSig[maxEnergybin];
	TH2D* efficiencyGamma[maxEnergybin];
	TH2D* efficiencyCR[maxEnergybin];


	for (Int_t e=0 ; e< maxEnergybin ;e++) //Loop the energy
	{
		NebulaSig[e] = new TH2D(Form("NebulaSig_E%d",e),"", hadronnessbin, hmin, hmax, theta2bin, tmin, tmax);
		efficiencyGamma[e] = new TH2D(Form("GammaEfficiency_E%d",e),"", hadronnessbin, hmin, hmax, theta2bin, tmin, tmax);
		efficiencyCR[e] = new TH2D(Form("CREfficiency_E%d",e),"", hadronnessbin, hmin, hmax, theta2bin, tmin, tmax);
		//cout << "maxEnergybin loop2: " << e << endl;
		PulsarSigP2[e] = new TH2D(Form("PulsarSig_E%d_P2",e),"", hadronnessbin, hmin, hmax, theta2bin, tmin, tmax);
	}


	Double_t maxP2[maxEnergybin]={};
	Double_t maxP2T[maxEnergybin]={};
	Double_t maxP2H[maxEnergybin]={};

	Double_t maxN[maxEnergybin]={};
	Double_t maxNH[maxEnergybin]={};
	Double_t maxNT[maxEnergybin]={};
	Double_t NormFactor[maxEnergybin][100];

	TImage *img = TImage::Create();

	//              	cout << "Normalization factor: " << NormFactor << endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////
	//THETA2 CUT && Significance Calculations
	//
	ofstream SignificancesP2 (Form("SignificancesP2Cycle%d.txt",cycleno)); //Will be saved in a txt file
	SignificancesP2.is_open();

	ofstream SignificancesN (Form("SignificancesNCycle%d.txt",cycleno)); //Will be saved in a txt file
	SignificancesN.is_open();


	for(Int_t e=0 ; e<maxEnergybin; e++) // Loop the energy
	{
		//cout << "maxEnergybin loop3: " << e << endl;
		for (Int_t i = 1; i <= histONe[1]->GetNbinsX(); i++) // Loop the Hadronness
		{
			for (Int_t j = 2; j <= histONe[1]->GetNbinsY(); j++) //Loop the Theta2
			{
				histONe[e]->SetBinContent(i,j,histONe[e]->GetBinContent(i,j)+histONe[e]->GetBinContent(i,j-1));
				histONe[e]->SetBinError(i,j,sqrt(histONe[e]->GetBinContent(i,j)));

				histOFFe[e]->SetBinContent(i,j,histOFFe[e]->GetBinContent(i,j)+histOFFe[e]->GetBinContent(i,j-1));
				histOFFe[e]->SetBinError(i,j,sqrt(histOFFe[e]->GetBinContent(i,j)));
			}

			for (Int_t k = 1; k < histOFFe[e]->GetNbinsY(); k++) // Loop the Theta2
			{
				Int_t tminbin = OFFCloneSmooth[e]->GetYaxis()->FindBin(0.);
				Int_t tmaxbin = OFFCloneSmooth[e]->GetYaxis()->FindBin(0.16);
				Double_t OFFperBin = (OFFCloneSmooth[e]->GetBinContent(i,tmaxbin)-OFFCloneSmooth[e]->GetBinContent(i,tminbin))/ (tmaxbin-tminbin);
				OFFCloneSmooth[e]->SetBinContent(i,k ,OFFperBin*k);
			}

			NormFactor[e][i]=1.;
		}
		///////////////////
		//SIGNIFICANCE CALCULATIONS

		for (Int_t i = 1; i <=  histONe[1]->GetNbinsY() ; i++) // Loop the Theta2
		{
			for (Int_t j = 1 ; j <= histONe[1]->GetNbinsX();  j++) //Loop the Hadronness
			{
				Double_t on = histONe[e]->GetBinContent(j,i);
				Double_t offs = OFFCloneSmooth[e]->GetBinContent(j,i);
				Double_t off = histOFFe[e]->GetBinContent(j,i);;
				Double_t normalizedoff = OFFCloneSmooth[e]->GetBinContent(j,i)*NormFactor[e][j];
				NebulaSig[e]->SetBinContent( j ,i ,MMath::SignificanceLiMa(on,off, NormFactor[e][j])); //MMath::SignificanceLiMa(on,off, NormFactor[e][h])
				//EfficiencyGamma = excess/ Excess no cut
				efficiencyGamma[e]-> SetBinContent (j ,i , (on-off)/(histONe[e]->GetBinContent(histONe[1]->GetNbinsX(),histONe[1]->GetNbinsY()) - histOFFe[e]->GetBinContent(histONe[1]->GetNbinsX(),histONe[1]->GetNbinsY())));
				//EfficiencyCR = OFF/ OFF no cut
				efficiencyCR[e]-> SetBinContent (j ,i , off/histOFFe[e]->GetBinContent(histONe[1]->GetNbinsX(),histONe[1]->GetNbinsY()));

					// To test
					//                    	counterP=counterP+1;
					//                            	cout<< "counterP" << counterP <<endl;
          Double_t F = 0.1 ;
					Double_t s2= (on-off)*F+0.45*(off+(on-off)*(1-F)) ;
					Double_t b2= 0.45*(off+(on-off)*(1-F));

					PulsarSigP2[e]->SetBinContent( j ,i ,MMath::SignificanceLiMa(s2,b2, NormFactor[e][j])); //MMath::SignificanceLiMa(s,b, NormFactor[e][h])

					//To test
					//cout << "e: " << e << " f: " << f <<" PulsarSig: " << endl;
					//cout<< "PulsarSig Bin Content: " << PulsarSigP12[1][1]->GetBinContent(h,t) << " h: " << h << " t: " << t << endl;

					if (PulsarSigP2[e]->GetBinContent(j,i) > maxP2[e])
					{
						maxP2[e] = PulsarSigP2[e]->GetBinContent(j,i);
						maxP2H[e] = j*(hmax-hmin)/hadronnessbin;
						maxP2T[e] = i*(tmax-tmin)/(theta2bin);
					}


					if (NebulaSig[e]->GetBinContent(j,i) > maxN[e])
					{
						maxN[e] = NebulaSig[e]->GetBinContent(j,i);
						maxNH[e] = j*(hmax-hmin)/hadronnessbin;
						maxNT[e] = i*(tmax-tmin)/(theta2bin);
					}
			} //hadronness ends

		} // theta ends
	} //Energy ends

SignificancesP2 << " Energybin " << " MaxSigP2 " << " Theta2 " << " Hadronness "<<'\n' ;
SignificancesN << " Energybin " << " MaxSigN " << " Theta2 " << " Hadronness "<<'\n' ;

	for(Int_t e=0 ; e<maxEnergybin; e++) // Loop the energy
	{
		//cout << "maxEnergybin loop4: " << e << endl;


			////To check!
			//cout<< e << " " << F << " " << "maxP[e][f]: " <<maxP1[e][f]<<endl;

			SignificancesP2 << e << "        "<< maxP2[e] << "           "<< maxP2T[e] << "         "<< maxP2H[e] << " "<< '\n' ;


			SignificancesN << e << "        " << maxN[e] << "            "<< maxNT[e] << "          "<< maxNH[e] << " "<< '\n' ;

	}

	SignificancesP2.close();
	SignificancesN.close();

  //cout << "cycleno: " << cycleno <<endl;
  //TString nebulasigtitle="Nebula Significances ";
	//drawHistogram(NebulaSig, maxEnergybin, cycleno, energyLevels, img, nebulasigtitle); // img

	return;
}
