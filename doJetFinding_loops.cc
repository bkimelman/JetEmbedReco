#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TProfile.h"
//#include "Math/Vector3D.h"
//#include "Math/Vector4D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include <iostream>
using namespace fastjet;
using namespace std;

int main(){
	
	
	TFile *input = new TFile("caloTowerExtract_eventData.root", "READ");
	
	TTree *tree = (TTree*)input->Get("jetRecoTree");

	//zVtx, MBD Energy
	vector<double> *eventData = 0;
	//data, embed, sim
	vector<vector<double>> *sourceTowers[3] = {0,0,0};
	//Px, Py, Pz, E, source
	//source = 0: CEMC, source = 1: IHCAL, source = 2: OHCAL
	vector<vector<double>> *dataTowers = 0;
	vector<vector<double>> *simTowers = 0;
	vector<vector<double>> *embedTowers = 0;

	
	TBranch *evtBr = 0;
	TBranch *dataBr = 0;
	TBranch *simBr = 0;
	TBranch *embBr = 0;

	tree->SetBranchAddress("eventData",&eventData,&evtBr);
	
	tree->SetBranchAddress("dataTowers",&dataTowers,&dataBr);
	tree->SetBranchAddress("simTowers",&simTowers,&simBr);
	tree->SetBranchAddress("embedTowers",&embedTowers,&embBr);

	cout << "set branches" << endl;

	JetDefinition jet_def_kt(kt_algorithm, 0.4);
	JetDefinition jet_def_antikt(antikt_algorithm, 0.4);

	Selector dataSelector = (!SelectorNHardest(1)) * SelectorAbsRapMax(1.0);
	Selector embedSelector = (!SelectorNHardest(3)) * SelectorAbsRapMax(1.0);
	Selector simSelector = (!SelectorNHardest(2)) * SelectorAbsRapMax(1.0);

	Selector embedSelector_drop0 = SelectorAbsRapMax(1.0);
	Selector embedSelector_drop1 = (!SelectorNHardest(1)) * SelectorAbsRapMax(1.0);
	Selector embedSelector_drop2 = (!SelectorNHardest(2)) * SelectorAbsRapMax(1.0);



	AreaDefinition area_def(active_area,GhostedAreaSpec(5.0));
	AreaDefinition area_def_bkg(active_area_explicit_ghosts,GhostedAreaSpec(5.0));


	JetMedianBackgroundEstimator bge_data(dataSelector, jet_def_kt, area_def_bkg);
	JetMedianBackgroundEstimator bge_embed(embedSelector, jet_def_kt, area_def_bkg);
	JetMedianBackgroundEstimator bge_sim(simSelector, jet_def_kt, area_def_bkg);

	JetMedianBackgroundEstimator bge[3] = {bge_data, bge_embed, bge_sim};

	JetMedianBackgroundEstimator bge_embed_drop0(embedSelector_drop0, jet_def_kt, area_def_bkg);
	JetMedianBackgroundEstimator bge_embed_drop1(embedSelector_drop1, jet_def_kt, area_def_bkg);
	JetMedianBackgroundEstimator bge_embed_drop2(embedSelector_drop2, jet_def_kt, area_def_bkg);

	JetMedianBackgroundEstimator bgeEmbed[3] = {bge_embed_drop0, bge_embed_drop1, bge_embed_drop2};

	TH1D *rho[3][3];
	TH1D *sigma[3][3];
	TH1D *pT[3][3][2];
	TProfile *rhoProfile[3][3][2];
	TProfile *sigmaProfile[3][3][2];

	string srcStr[3] = {"Data","Embed","Sim"};
	string centStr[3] = {"central","peripheral","all"};
	string subStr[2] = {"unsubtracted","subtracted"};

	int color[3] = {418,600,632}; //cent
	int marker[3] = {20,21,22}; //src

	//darker, lighter, nominal
	int colorCent[3][3] = {{419,410,417},{602,594,600},{634,626,632}};


	for(int src=0; src<3; src++){

		for(int cent=0; cent<3; cent++){

			rho[src][cent] = new TH1D(Form("rho%s_%sCentrality",srcStr[src].c_str(),centStr[cent].c_str()),Form("%s %s Centrality;#rho",srcStr[src].c_str(),centStr[cent].c_str()),1000,0,100);
			rho[src][cent]->SetMarkerStyle(marker[src]);
			rho[src][cent]->SetMarkerColor(colorCent[src][cent]);
			rho[src][cent]->SetLineColor(colorCent[src][cent]);
			sigma[src][cent] = new TH1D(Form("sigma%s_%sCentrality",srcStr[src].c_str(),centStr[cent].c_str()),Form("%s %s Centrality;#sigma",srcStr[src].c_str(),centStr[cent].c_str()),100,0,20);
			sigma[src][cent]->SetMarkerStyle(marker[src]);
			sigma[src][cent]->SetMarkerColor(colorCent[src][cent]);
			sigma[src][cent]->SetLineColor(colorCent[src][cent]);

			for(int sub=0; sub<2; sub++){

				rhoProfile[src][cent][sub] = new TProfile(Form("rhoProfile%s_%s_%sCentrality",srcStr[src].c_str(),subStr[sub].c_str(),centStr[cent].c_str()),Form("%s %s %s Centrality;%s leading p_{T} [GeV/c];#LT#rho#GT",srcStr[src].c_str(),subStr[sub].c_str(),centStr[cent].c_str(),subStr[sub].c_str()),100,0,50,0,100);
				rhoProfile[src][cent][sub]->SetMarkerStyle(marker[src]);
				rhoProfile[src][cent][sub]->SetMarkerColor(colorCent[src][cent]);
				sigmaProfile[src][cent][sub] = new TProfile(Form("sigmaProfile%s_%s_%sCentrality",srcStr[src].c_str(),subStr[sub].c_str(),centStr[cent].c_str()),Form("%s %s %s Centrality;%s leading p_{T} [GeV/c];#LT#sigma#GT",srcStr[src].c_str(),subStr[sub].c_str(),centStr[cent].c_str(),subStr[sub].c_str()),100,0,50,0,100);
				sigmaProfile[src][cent][sub]->SetMarkerStyle(marker[src]);
				sigmaProfile[src][cent][sub]->SetMarkerColor(colorCent[src][cent]);

				pT[src][cent][sub] = new TH1D(Form("pT%s_%s_%sCentrality",srcStr[src].c_str(),subStr[sub].c_str(),centStr[cent].c_str()),Form("%s %s %s Centrality;%s p_{T} [GeV/c]",srcStr[src].c_str(),subStr[sub].c_str(),centStr[cent].c_str(),subStr[sub].c_str()),100,0,50);
				pT[src][cent][sub]->SetMarkerStyle(marker[src]);
				pT[src][cent][sub]->SetMarkerColor(colorCent[src][cent]);

			}
		}
	}

	//drop 0, 1, 2, then cent
	TH1D *rhoEmbed[4][3];
	TH1D *sigmaEmbed[4][3];
	TProfile *rhoEmbedProfile[4][3][2];
	TProfile *sigmaEmbedProfile[4][3][2];

	int colorEmbed[3] = {591,593,596};

	for(int drop=0; drop<3; drop++){

		for(int cent=0; cent<3; cent++){

			rhoEmbed[drop][cent] = new TH1D(Form("rhoEmbed_drop%d_%sCentrality",drop,centStr[cent].c_str()),Form("Embed %d dropped jets %s Centrality;#rho",drop,centStr[cent].c_str()),1000,0,100);
			rhoEmbed[drop][cent]->SetMarkerStyle(marker[1]);
			rhoEmbed[drop][cent]->SetMarkerColor(colorEmbed[drop]);
			rhoEmbed[drop][cent]->SetLineColor(colorEmbed[drop]);
			sigmaEmbed[drop][cent] = new TH1D(Form("sigmaEmbed_drop%d_%sCentrality",drop,centStr[cent].c_str()),Form("Embed %d dropped jets %s Centrality;#sigma",drop,centStr[cent].c_str()),100,0,20);
			sigmaEmbed[drop][cent]->SetMarkerStyle(marker[1]);
			sigmaEmbed[drop][cent]->SetMarkerColor(colorEmbed[drop]);
			sigmaEmbed[drop][cent]->SetLineColor(colorEmbed[drop]);


			for(int sub=0; sub<2; sub++){

				rhoEmbedProfile[drop][cent][sub] = new TProfile(Form("rhoEmbedProfile_drop%d_%s_%sCentrality",drop,subStr[sub].c_str(),centStr[cent].c_str()),Form("Embed %d dropped jets %s %s Centrality;%s leading p_{T} [GeV/c];#LT#rho#GT",drop,subStr[sub].c_str(),centStr[cent].c_str(),subStr[sub].c_str()),100,0,50,0,100);
				rhoEmbedProfile[drop][cent][sub]->SetMarkerStyle(marker[1]);
				rhoEmbedProfile[drop][cent][sub]->SetMarkerColor(colorEmbed[drop]);
				sigmaEmbedProfile[drop][cent][sub] = new TProfile(Form("sigmaEmbedProfile_drop%d_%s_%sCentrality",drop,subStr[sub].c_str(),centStr[cent].c_str()),Form("Embed %d dropped jets %s %s Centrality;%s leading p_{T} [GeV/c];#LT#sigma#GT",drop,subStr[sub].c_str(),centStr[cent].c_str(),subStr[sub].c_str()),100,0,50,0,100);
				sigmaEmbedProfile[drop][cent][sub]->SetMarkerStyle(marker[1]);
				sigmaEmbedProfile[drop][cent][sub]->SetMarkerColor(colorEmbed[drop]);

			}
		}
	}

	cout << "defined histos" << endl;

	bool skip_CEMC = false;
	bool skip_IHCal = false;
	bool skip_OHCal = false;

	int events = tree->GetEntries();
	//int events = 100;
	cout << "starting event loop" << endl;
	for(int ev=0; ev<events; ev++){

		tree->GetEntry(ev);

		cout << "working on event " << ev << endl;

		
		sourceTowers[0] = dataTowers;
		sourceTowers[1] = embedTowers;
		sourceTowers[2] = simTowers;
	

		double mbdEnergy = eventData->at(1);
		int centBin = -1;
		if(mbdEnergy > 150000) centBin = 0;
		else if(mbdEnergy < 50000) centBin = 1;


		cout << "done with event info" << endl;

		//loop over source DSTs
		for(int src=0; src<3; src++){

			cout << "working on " << srcStr[src] << endl;

			//cout << "sourceTowers[src]: " << sourceTowers[src] << "   size: " << sourceTowers[src]->size() << endl;
			//cout << "sourceTowers[src]: " << "   size: " << sourceTowers[src]->size() << endl;

			vector<PseudoJet> sourcePseudoJets;
			for(int tow=0; tow<sourceTowers[src]->size(); tow++){
				//cout << "working on tower " << tow << endl;
				if((skip_CEMC && sourceTowers[src]->at(tow)[4] == 0.0) || (skip_IHCal && sourceTowers[src]->at(tow)[4] == 1.0) || (skip_OHCal && sourceTowers[src]->at(tow)[4] == 2.0)) continue;
				//cout << "good tower" << endl;
				sourcePseudoJets.push_back( PseudoJet( sourceTowers[src]->at(tow)[0], sourceTowers[src]->at(tow)[1], sourceTowers[src]->at(tow)[2], sourceTowers[src]->at(tow)[3]) );
				//cout << "made PseudoJet" << endl;
			}
			cout << "got inputs into jet reco" << endl;
			if(sourcePseudoJets.size() == 0) continue;
			bge[src].set_particles(sourcePseudoJets);
			Subtractor subtract(&bge[src]);
			BackgroundEstimate bkgd_estimate = bge[src].estimate();
			double rhoValue = bkgd_estimate.rho();
			double sigmaValue = bkgd_estimate.sigma();

			ClusterSequenceArea cs(sourcePseudoJets, jet_def_antikt, area_def);
			double leading_pt = sorted_by_pt(cs.inclusive_jets())[0].pt();
			//double leading_pt_sub = subtract(sorted_by_pt(cs.inclusive_jets())[0]).pt();
			double leading_pt_sub = sorted_by_pt(subtract(cs.inclusive_jets()))[0].pt();

			for(int j=0; j<cs.inclusive_jets().size(); j++){
				double pTVal = sorted_by_pt(cs.inclusive_jets())[j].pt();
				pT[src][2][0]->Fill(pTVal);
				if(centBin != -1) pT[src][centBin][0]->Fill(pTVal);
			}

			for(int j=0; j<subtract(cs.inclusive_jets()).size(); j++){
				double pTVal = sorted_by_pt(subtract(cs.inclusive_jets()))[j].pt();
				if(pTVal > 0){
					pT[src][2][1]->Fill(pTVal);
					if(centBin != -1) pT[src][centBin][1]->Fill(pTVal);
				}
			}
			

			if(centBin != -1){
				rho[src][centBin]->Fill(rhoValue);
				sigma[src][centBin]->Fill(sigmaValue);

				rhoProfile[src][centBin][0]->Fill(leading_pt,rhoValue);
				sigmaProfile[src][centBin][0]->Fill(leading_pt,sigmaValue);

				rhoProfile[src][centBin][1]->Fill(leading_pt_sub,rhoValue);
				sigmaProfile[src][centBin][1]->Fill(leading_pt_sub,sigmaValue);
			}
			rho[src][2]->Fill(rhoValue);
			sigma[src][2]->Fill(sigmaValue);

			rhoProfile[src][2][0]->Fill(leading_pt,rhoValue);
			sigmaProfile[src][2][0]->Fill(leading_pt,sigmaValue);

			rhoProfile[src][2][1]->Fill(leading_pt_sub,rhoValue);
			sigmaProfile[src][2][1]->Fill(leading_pt_sub,sigmaValue);

			//if embed, loop over additional cases with different numbers of jets to drop
			if(src == 1){
				for(int drop=0; drop<3; drop++){
					bgeEmbed[drop].set_particles(sourcePseudoJets);
					Subtractor subtractEmbed(&bgeEmbed[drop]);
					BackgroundEstimate bkgd_estimate_embed = bgeEmbed[drop].estimate();
					double rhoValueEmbed = bkgd_estimate_embed.rho();
					double sigmaValueEmbed = bkgd_estimate_embed.sigma();

					//double leading_pt_sub_embed = subtractEmbed(sorted_by_pt(cs.inclusive_jets())[0]).pt();
					double leading_pt_sub_embed = sorted_by_pt(subtractEmbed(cs.inclusive_jets()))[0].pt();

					if(centBin != -1){
						rhoEmbed[drop][centBin]->Fill(rhoValueEmbed);
						sigmaEmbed[drop][centBin]->Fill(sigmaValueEmbed);

						rhoEmbedProfile[drop][centBin][0]->Fill(leading_pt,rhoValueEmbed);
						sigmaEmbedProfile[drop][centBin][0]->Fill(leading_pt,sigmaValueEmbed);

						rhoEmbedProfile[drop][centBin][1]->Fill(leading_pt_sub_embed,rhoValueEmbed);
						sigmaEmbedProfile[drop][centBin][1]->Fill(leading_pt_sub_embed,sigmaValueEmbed);
					}
					rhoEmbed[drop][2]->Fill(rhoValueEmbed);
					sigmaEmbed[drop][2]->Fill(sigmaValueEmbed);

					rhoEmbedProfile[drop][2][0]->Fill(leading_pt,rhoValueEmbed);
					sigmaProfile[drop][2][0]->Fill(leading_pt,sigmaValueEmbed);

					rhoEmbedProfile[drop][2][1]->Fill(leading_pt_sub_embed,rhoValueEmbed);
					sigmaEmbedProfile[drop][2][1]->Fill(leading_pt_sub_embed,sigmaValueEmbed);

				}
			}

		}
	}

	TCanvas *c1 = new TCanvas();
	c1->SetLogy();
	gStyle->SetOptStat(0);

	string imageDir = "images";

	if(!skip_CEMC && skip_IHCal && skip_OHCal) imageDir += "_onlyCEMC";
	if(skip_CEMC && !skip_IHCal && skip_OHCal) imageDir += "_onlyIHCAL";
	if(skip_CEMC && skip_IHCal && !skip_OHCal) imageDir += "_onlyOHCAL";

	//one cent at a time, all sources
	for(int cent=0; cent<3; cent++){

		c1->SetLogy();

		c1->Clear();
		rho[0][cent]->SetTitle(Form("#rho %s Centrality",centStr[cent].c_str()));
		rho[0][cent]->Draw("PH");
		rho[1][cent]->Draw("PHSAME");
		rho[2][cent]->Draw("PHSAME");

		TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
		leg->AddEntry(rho[0][cent],"Data","p");
		leg->AddEntry(rho[2][cent],"Sim","p");
		leg->AddEntry(rho[1][cent],"Embed","p");
		leg->Draw("same");

		c1->SaveAs(Form("%s/rho_%sCentrality.png",imageDir.c_str(),centStr[cent].c_str()));

		c1->Clear();
		sigma[0][cent]->SetTitle(Form("#sigma %s Centrality",centStr[cent].c_str()));
		sigma[0][cent]->Draw("PH");
		sigma[1][cent]->Draw("PHSAME");
		sigma[2][cent]->Draw("PHSAME");
		leg->Draw("same");

		c1->SaveAs(Form("%s/sigma_%sCentrality.png",imageDir.c_str(),centStr[cent].c_str()));

		for(int sub=0; sub<2; sub++){

			c1->SetLogy();

			c1->Clear();
			pT[0][cent][sub]->SetTitle(Form("p_{T} %s Centrality %s",centStr[cent].c_str(),subStr[sub].c_str()));
			pT[0][cent][sub]->Draw("P");
			pT[1][cent][sub]->Draw("PSAME");
			pT[2][cent][sub]->Draw("PSAME");
			leg->Draw("same");
			c1->SaveAs(Form("%s/pT_%sCentrality_%s.png",imageDir.c_str(),centStr[cent].c_str(),subStr[sub].c_str()));

			c1->SetLogy(0);

			c1->Clear();
			rhoProfile[0][cent][sub]->SetTitle(Form("#rho %s Centrality %s",centStr[cent].c_str(),subStr[sub].c_str()));
			rhoProfile[0][cent][sub]->Draw("P");
			rhoProfile[1][cent][sub]->Draw("PSAME");
			rhoProfile[2][cent][sub]->Draw("PSAME");
			leg->Draw("same");

			c1->SaveAs(Form("%s/rhoProfile_%sCentrality_%s.png",imageDir.c_str(),centStr[cent].c_str(),subStr[sub].c_str()));

			c1->Clear();
			sigmaProfile[0][cent][sub]->SetTitle(Form("#sigma %s Centrality %s",centStr[cent].c_str(),subStr[sub].c_str()));
			sigmaProfile[0][cent][sub]->Draw("P");
			sigmaProfile[1][cent][sub]->Draw("PSAME");
			sigmaProfile[2][cent][sub]->Draw("PSAME");
			leg->Draw("same");

			c1->SaveAs(Form("%s/sigmaProfile_%sCentrality_%s.png",imageDir.c_str(),centStr[cent].c_str(),subStr[sub].c_str()));

		}

	}

	//one source at a time, all cent
	for(int src=0; src<3; src++){

		c1->SetLogy();

		c1->Clear();
		rho[src][2]->SetTitle(Form("#rho %s",srcStr[src].c_str()));
		rho[src][2]->Draw("PH");
		rho[src][0]->Draw("PHSAME");
		rho[src][1]->Draw("PHSAME");

		TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
		leg->AddEntry(rho[src][2],"Minimum Bias","p");
		leg->AddEntry(rho[src][0],"Central","p");
		leg->AddEntry(rho[src][1],"Peripheral","p");
		leg->Draw("same");

		c1->SaveAs(Form("%s/rho%s.png",imageDir.c_str(),srcStr[src].c_str()));

		c1->Clear();
		sigma[src][2]->SetTitle(Form("#sigma %s",srcStr[src].c_str()));
		sigma[src][2]->Draw("PH");
		sigma[src][0]->Draw("PHSAME");
		sigma[src][1]->Draw("PHSAME");
		leg->Draw("same");

		c1->SaveAs(Form("%s/sigma%s.png",imageDir.c_str(),srcStr[src].c_str()));

		for(int sub=0; sub<2; sub++){

			c1->SetLogy();

			c1->Clear();
			pT[src][2][sub]->SetTitle(Form("p_{T} %s %s",srcStr[src].c_str(),subStr[sub].c_str()));
			pT[src][2][sub]->Draw("P");
			pT[src][0][sub]->Draw("PSAME");
			pT[src][1][sub]->Draw("PSAME");
			leg->Draw("same");
			c1->SaveAs(Form("%s/pT%s_%s.png",imageDir.c_str(),srcStr[src].c_str(),subStr[sub].c_str()));

			c1->SetLogy(0);

			c1->Clear();
			rhoProfile[src][2][sub]->SetTitle(Form("#rho %s %s",srcStr[src].c_str(),subStr[sub].c_str()));
			rhoProfile[src][2][sub]->Draw("P");
			rhoProfile[src][0][sub]->Draw("PSAME");
			rhoProfile[src][1][sub]->Draw("PSAME");
			leg->Draw("same");

			c1->SaveAs(Form("%s/rhoProfile%s_%s.png",imageDir.c_str(),srcStr[src].c_str(),subStr[sub].c_str()));

			c1->Clear();
			sigmaProfile[src][2][sub]->SetTitle(Form("#sigma %s %s",srcStr[src].c_str(),subStr[sub].c_str()));
			sigmaProfile[src][2][sub]->Draw("P");
			sigmaProfile[src][0][sub]->Draw("PSAME");
			sigmaProfile[src][1][sub]->Draw("PSAME");
			leg->Draw("same");

			c1->SaveAs(Form("%s/sigmaProfile%s_%s.png",imageDir.c_str(),srcStr[src].c_str(),subStr[sub].c_str()));

		}

	}

	c1->SetLogy();

	//all cent and source
	c1->Clear();
	rho[0][2]->SetTitle("#rho");
	rho[0][2]->Draw("PH");

	TLegend *legAll = new TLegend(0.6,0.6,0.9,0.9);
	for(int src=0; src<3; src++){
		for(int cent=0; cent<3; cent++){
			legAll->AddEntry(rho[src][cent],Form("%s %s Centrality",srcStr[src].c_str(),centStr[cent].c_str()));
			if(src == 0 && cent == 2) continue;
			rho[src][cent]->Draw("PHSAME");
		}
	}
	legAll->Draw("same");
	c1->SaveAs(Form("%s/rho_allSources_allCentrality.png",imageDir.c_str()));

	c1->Clear();
	sigma[0][2]->SetTitle("#sigma");
	sigma[0][2]->Draw("PH");

	for(int src=0; src<3; src++){
		for(int cent=0; cent<3; cent++){
			if(src == 0 && cent == 2) continue;
			sigma[src][cent]->Draw("PHSAME");
		}
	}
	legAll->Draw("same");
	c1->SaveAs(Form("%s/sigma_allSources_allCentrality.png",imageDir.c_str()));

	c1->SetLogy(0);

	c1->Clear();
	rhoProfile[0][2][0]->SetTitle("#rho unsubtracted");
	rhoProfile[0][2][0]->Draw("P");

	for(int src=0; src<3; src++){
		for(int cent=0; cent<3; cent++){
			if(src == 0 && cent == 2) continue;
			rhoProfile[src][cent][0]->Draw("PSAME");
		}
	}
	legAll->Draw("same");
	c1->SaveAs(Form("%s/rhoProfile_allSources_allCentrality_unsubtracted.png",imageDir.c_str()));

	c1->Clear();
	rhoProfile[0][2][1]->SetTitle("#rho subtracted");
	rhoProfile[0][2][1]->Draw("P");

	for(int src=0; src<3; src++){
		for(int cent=0; cent<3; cent++){
			if(src == 0 && cent == 2) continue;
			rhoProfile[src][cent][1]->Draw("PSAME");
		}
	}
	legAll->Draw("same");
	c1->SaveAs(Form("%s/rhoProfile_allSources_allCentrality_subtracted.png",imageDir.c_str()));


	c1->SetLogy();
	c1->Clear();
	sigma[0][2]->SetTitle("#sigma");
	sigma[0][2]->Draw("PH");

	for(int src=0; src<3; src++){
		for(int cent=0; cent<3; cent++){
			if(src == 0 && cent == 2) continue;
			sigma[src][cent]->Draw("PHSAME");
		}
	}
	legAll->Draw("same");
	c1->SaveAs(Form("%s/sigma_allSources_allCentrality.png",imageDir.c_str()));

	c1->SetLogy(0);
	c1->Clear();
	sigmaProfile[0][2][0]->SetTitle("#sigma unsubtracted");
	sigmaProfile[0][2][0]->Draw("P");

	for(int src=0; src<3; src++){
		for(int cent=0; cent<3; cent++){
			if(src == 0 && cent == 2) continue;
			sigmaProfile[src][cent][0]->Draw("PSAME");
		}
	}
	legAll->Draw("same");
	c1->SaveAs(Form("%s/sigmaProfile_allSources_allCentrality_unsubtracted.png",imageDir.c_str()));

	c1->Clear();
	sigmaProfile[0][2][1]->SetTitle("#sigma subtracted");
	sigmaProfile[0][2][1]->Draw("P");

	for(int src=0; src<3; src++){
		for(int cent=0; cent<3; cent++){
			if(src == 0 && cent == 2) continue;
			sigmaProfile[src][cent][1]->Draw("PSAME");
		}
	}
	legAll->Draw("same");
	c1->SaveAs(Form("%s/sigmaProfile_allSources_allCentrality_subtracted.png",imageDir.c_str()));




	//embed drop, one cent at a time, all drop
	for(int cent=0; cent<3; cent++){
		rhoEmbed[3][cent] = rho[1][cent];

		c1->SetLogy();

		c1->Clear();
		rhoEmbed[0][cent]->SetTitle(Form("#rho Embed %s Centrality",centStr[cent].c_str()));
		rhoEmbed[0][cent]->Draw("PH");
		rhoEmbed[1][cent]->Draw("PHSAME");
		rhoEmbed[2][cent]->Draw("PHSAME");
		rhoEmbed[3][cent]->Draw("PHSAME");

		TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
		leg->AddEntry(rhoEmbed[0][cent],"0 dropped jets","p");
		leg->AddEntry(rhoEmbed[2][cent],"1 dropped jet","p");
		leg->AddEntry(rhoEmbed[1][cent],"2 dropped jets","p");
		leg->AddEntry(rhoEmbed[3][cent],"3 dropped jets","p");
		leg->Draw("same");

		c1->SaveAs(Form("%s/rhoEmbed_%sCentrality.png",imageDir.c_str(),centStr[cent].c_str()));

		sigmaEmbed[3][cent] = sigma[1][cent];

		c1->Clear();
		sigmaEmbed[0][cent]->SetTitle(Form("#sigma Embed %s Centrality",centStr[cent].c_str()));
		sigmaEmbed[0][cent]->Draw("PH");
		sigmaEmbed[1][cent]->Draw("PHSAME");
		sigmaEmbed[2][cent]->Draw("PHSAME");
		sigmaEmbed[3][cent]->Draw("PHSAME");
		leg->Draw("same");

		c1->SaveAs(Form("%s/sigmaEmbed_%sCentrality.png",imageDir.c_str(),centStr[cent].c_str()));

		for(int sub=0; sub<2; sub++){

			c1->SetLogy(0);

			rhoEmbedProfile[3][cent][sub] = rhoProfile[1][cent][sub];

			c1->Clear();
			rhoEmbedProfile[0][cent][sub]->SetTitle(Form("#rho Embed %s Centrality %s",centStr[cent].c_str(),subStr[sub].c_str()));
			rhoEmbedProfile[0][cent][sub]->Draw("P");
			rhoEmbedProfile[1][cent][sub]->Draw("PSAME");
			rhoEmbedProfile[2][cent][sub]->Draw("PSAME");
			rhoEmbedProfile[3][cent][sub]->Draw("PSAME");
			leg->Draw("same");

			c1->SaveAs(Form("%s/rhoEmbedProfile_%sCentrality_%s.png",imageDir.c_str(),centStr[cent].c_str(),subStr[sub].c_str()));

			sigmaEmbedProfile[3][cent][sub] = sigmaProfile[1][cent][sub];

			c1->Clear();
			sigmaEmbedProfile[0][cent][sub]->SetTitle(Form("#sigma Embed %s Centrality %s",centStr[cent].c_str(),subStr[sub].c_str()));
			sigmaEmbedProfile[0][cent][sub]->Draw("P");
			sigmaEmbedProfile[1][cent][sub]->Draw("PSAME");
			sigmaEmbedProfile[2][cent][sub]->Draw("PSAME");
			sigmaEmbedProfile[3][cent][sub]->Draw("PSAME");
			leg->Draw("same");

			c1->SaveAs(Form("%s/sigmaEmbedProfile_%sCentrality_%s.png",imageDir.c_str(),centStr[cent].c_str(),subStr[sub].c_str()));

		}
	}

//set rhoEmbed[3] to rho[1] so that it can be called easily, Haven't yet done the loop below
	//embed drop, one drop at a time, all cent
	for(int drop=0; drop<4; drop++){
		c1->Clear();

		c1->SetLogy();
		
		rhoEmbed[drop][2]->SetTitle(Form("#rho Embed %d dropped jets",drop));
		rhoEmbed[drop][2]->Draw("PH");
		rhoEmbed[drop][0]->Draw("PHSAME");
		rhoEmbed[drop][1]->Draw("PHSAME");

		TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
		leg->AddEntry(rhoEmbed[drop][2],"Minimum Bias","p");
		leg->AddEntry(rhoEmbed[drop][0],"Central","p");
		leg->AddEntry(rhoEmbed[drop][1],"Peripheral","p");
		leg->Draw("same");

		c1->SaveAs(Form("%s/rhoEmbed_%dDropped.png",imageDir.c_str(),drop));

		c1->Clear();
		sigmaEmbed[drop][2]->SetTitle(Form("#sigma Embed %d dropped jets",drop));
		sigmaEmbed[drop][2]->Draw("PH");
		sigmaEmbed[drop][0]->Draw("PHSAME");
		sigmaEmbed[drop][1]->Draw("PHSAME");
		leg->Draw("same");

		c1->SaveAs(Form("%s/sigmaEmbed_%dDropped.png",imageDir.c_str(),drop));

		for(int sub=0; sub<2; sub++){


			c1->SetLogy(0);
			c1->Clear();
			rhoEmbedProfile[drop][2][sub]->SetTitle(Form("#rho Embed %d dropped jets %s",drop,subStr[sub].c_str()));
			rhoEmbedProfile[drop][2][sub]->Draw("P");
			rhoEmbedProfile[drop][0][sub]->Draw("PSAME");
			rhoEmbedProfile[drop][1][sub]->Draw("PSAME");
			leg->Draw("same");

			c1->SaveAs(Form("%s/rhoEmbedProfile_%dDropped_%s.png",imageDir.c_str(),drop,subStr[sub].c_str()));

			c1->Clear();
			sigmaEmbedProfile[drop][2][sub]->SetTitle(Form("#sigma Embed %d dropped jets %s",drop,subStr[sub].c_str()));
			sigmaEmbedProfile[drop][2][sub]->Draw("P");
			sigmaEmbedProfile[drop][0][sub]->Draw("PSAME");
			sigmaEmbedProfile[drop][1][sub]->Draw("PSAME");
			leg->Draw("same");

			c1->SaveAs(Form("%s/sigmaEmbedProfile_%dDropped_%s.png",imageDir.c_str(),drop,subStr[sub].c_str()));

		}
	}



//all cent and drop
	c1->SetLogy();
	c1->Clear();
	rhoEmbed[0][2]->SetTitle("#rho Embed");
	rhoEmbed[0][2]->Draw("PH");

	TLegend *legAllEmb = new TLegend(0.6,0.6,0.9,0.9);
	for(int drop=0; drop<4; drop++){
		for(int cent=0; cent<3; cent++){
			legAllEmb->AddEntry(rhoEmbed[drop][cent],Form("%d dropped jets %s Centrality",drop,centStr[cent].c_str()));
			if(drop == 0 && cent == 2) continue;
			rhoEmbed[drop][cent]->Draw("PHSAME");
		}
	}
	legAllEmb->Draw("same");
	c1->SaveAs(Form("%s/rhoEmbed_allDrop_allCentrality.png",imageDir.c_str()));

	c1->Clear();
	sigmaEmbed[0][2]->SetTitle("#sigma Embed");
	sigmaEmbed[0][2]->Draw("PH");

	for(int drop=0; drop<4; drop++){
		for(int cent=0; cent<3; cent++){
			if(drop == 0 && cent == 2) continue;
			sigmaEmbed[drop][cent]->Draw("PHSAME");
		}
	}
	legAllEmb->Draw("same");
	c1->SaveAs(Form("%s/sigmaEmbed_allDrop_allCentrality.png",imageDir.c_str()));


	c1->SetLogy(0);
	c1->Clear();
	rhoEmbedProfile[0][2][0]->SetTitle("#rho Embed unsubtracted");
	rhoEmbedProfile[0][2][0]->Draw("P");

	for(int drop=0; drop<4; drop++){
		for(int cent=0; cent<3; cent++){
			if(drop == 0 && cent == 2) continue;
			rhoEmbedProfile[drop][cent][0]->Draw("PSAME");
		}
	}
	legAllEmb->Draw("same");
	c1->SaveAs(Form("%s/rhoEmbedProfile_allDrop_allCentrality_unsubtracted.png",imageDir.c_str()));

	c1->Clear();
	rhoEmbedProfile[0][2][1]->SetTitle("#rho Embed subtracted");
	rhoEmbedProfile[0][2][1]->Draw("P");

	for(int drop=0; drop<4; drop++){
		for(int cent=0; cent<3; cent++){
			if(drop == 0 && cent == 2) continue;
			rhoEmbedProfile[drop][cent][1]->Draw("PSAME");
		}
	}
	legAllEmb->Draw("same");
	c1->SaveAs(Form("%s/rhoEmbedProfile_allDrop_allCentrality_subtracted.png",imageDir.c_str()));

	c1->Clear();
	sigmaEmbed[0][2]->SetTitle("#sigma Embed");
	sigmaEmbed[0][2]->Draw("PH");

	for(int drop=0; drop<4; drop++){
		for(int cent=0; cent<3; cent++){
			if(drop == 0 && cent == 2) continue;
			sigmaEmbed[drop][cent]->Draw("PHSAME");
		}
	}
	legAllEmb->Draw("same");
	c1->SaveAs(Form("%s/sigmaEmbed_allDrop_allCentrality.png",imageDir.c_str()));

	c1->Clear();
	sigmaEmbedProfile[0][2][0]->SetTitle("#sigma Embed unsubtracted");
	sigmaEmbedProfile[0][2][0]->Draw("P");

	for(int drop=0; drop<4; drop++){
		for(int cent=0; cent<3; cent++){
			if(drop == 0 && cent == 2) continue;
			sigmaEmbedProfile[drop][cent][0]->Draw("PSAME");
		}
	}
	legAllEmb->Draw("same");
	c1->SaveAs(Form("%s/sigmaEmbedProfile_allDrop_allCentrality_unsubtracted.png",imageDir.c_str()));

	c1->Clear();
	sigmaEmbedProfile[0][2][1]->SetTitle("#sigma Embed subtracted");
	sigmaEmbedProfile[0][2][1]->Draw("P");

	for(int drop=0; drop<4; drop++){
		for(int cent=0; cent<3; cent++){
			if(drop == 0 && cent == 2) continue;
			sigmaEmbedProfile[drop][cent][1]->Draw("PSAME");
		}
	}
	legAllEmb->Draw("same");
	c1->SaveAs(Form("%s/sigmaEmbedProfile_allDrop_allCentrality_subtracted.png",imageDir.c_str()));


}