// implementing s_weight to plot comparison plots for the control channel.
//
// \author Rishabh Raturi

#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"


using namespace RooFit;
using namespace RooStats;

void AddModel(RooWorkspace *);
void AddData(RooWorkspace *);
void DoSPlot(RooWorkspace *);

void splot()
{
	//create a workspace
	RooWorkspace *wspace = new RooWorkspace("w_space");

	//create a function to add the real-time data into the workspace
        AddData(wspace);

	//add signal and background model to the workspace
	// create a function AddModel() which has the description of the model
	AddModel(wspace);

	//choice to print the workspac
	// wspace->Print();
	
	//this function creates a new dataset with s_weight applied to every event
	DoSPlot(wspace);

	//this fucntion plots the comparison plots for the control channel
	//MakePlots(wspace);

	//save workspace in a file
	wspace->writeToFile("workspace.root");
	
}

void AddModel(RooWorkspace *ws)
{
	//make model for JPsi - double crystal ball for signal and exponential for the background

	double bm_min(5.2), bm_max(5.6);
	RooRealVar Bmass("Bmass","#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}", bm_min, bm_max);	

	RooRealVar  mean("mean","common means for Crystal Balls", 5.367, 5.2, 5.6);
        RooRealVar  sigma1("sigma1","sigma of CB1",  0.0252789);
        RooRealVar  sigma2("sigma2","sigma of CB2",  0.0474879);
        RooRealVar  sigM_frac("sigM_frac","fraction of CB",  0.557917);
        RooRealVar  n1("n1", "", 1.79897);
        RooRealVar  n2("n2", "", 10.7529);
        RooRealVar  alpha1("alpha1","alpha for CB1", 1.72585);
        RooRealVar  alpha2("alpha2","alpha for CB2", -2.09195);

	RooCBShape CB1("CB1","Crystal Ball-1", Bmass, mean, sigma1, alpha1,n1);
        RooCBShape CB2("CB2","Crystal Ball-2", Bmass, mean, sigma2, alpha2,n2);

	RooAddPdf CB("CB","CB1+CB2", RooArgList(CB1,CB2), RooArgList(sigM_frac));

	RooRealVar lambda("lambda","slope",-10.,0.);
        RooExponential bkgE("bkgE","exponential PDF", Bmass, lambda);

	RooRealVar nsig("nsig","nsig",1E3,100,1E9);
        RooRealVar nbkgE("nbkgE","number of exponential bkg events", 8000, 0, 1E7);

	//create final model for fitting
        RooAddPdf model("model", "CB+e", RooArgList(CB,bkgE), RooArgList(nsig,nbkgE));
        model.Print("t");

	std::cout << "import model to workspace" << std::endl;
	ws->import(model);
}

void AddData(RooWorkspace *ws)
{
	// this function add the real-time data for phimumu
	
	TChain *ch = new TChain("tree");
        ch->Add("./fits/JPsi_Data_ForFit.root");
        TTree *tr = ch;

        int nentries_ = tr->GetEntries();
        cout << "\n=> total entries in data tree = " << nentries_ << endl;

	double bm_min(5.1), bm_max(5.6);

	RooRealVar Bmass("Bmass","#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}", bm_min, bm_max);
        RooRealVar Mumumass("Mumumass","M^{#mu#mu}", 0., 10.);
        RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}", 0., 10.);
        RooRealVar Phimass("Phimass","phi mass", 0.,10.);
        RooRealVar Kmpt("Kmpt","K- pt", 0., 200.);
        RooRealVar Kppt("Kppt","K+ pt", 0., 200.);
        RooRealVar Kmtrkdcasigbs("Kmtrkdcasigbs","K^{-} track DCA/#sigma beam spot", 0., 1000.);
        RooRealVar Kptrkdcasigbs("Kptrkdcasigbs","K^{+} track DCA/#sigma beam spot", 0., 1000.);
        RooRealVar Blxysig("Blxysig","Blxy sig.", 0., 1000.);
        RooRealVar Bcosalphabs2d("Bcosalphabs2d","Bcosalphabs2d", 0., 1.);
        RooRealVar Bvtxcl("Bvtxcl","B vtx cl.", 0., 1.);

	RooArgSet  observables(Bmass,Mumumass,Kmpt,Kppt,Kmtrkdcasigbs,Kptrkdcasigbs,Blxysig,Bcosalphabs2d,Bvtxcl);
        observables.add(RooArgSet(Phimass));

	RooDataSet data("data","dataset with Bmass", ch, observables);

	TCut c1 = Form("Bmass>%f && Bmass<%f",bm_min,bm_max);
        TCut cutTotal = c1  ;

	RooDataSet *redData = (RooDataSet*)data.reduce(cutTotal);   //redData consists data to be fitted

	//importing this dataset to the workspace
	ws->import(*redData, Rename("data"));

	
}

void DoSPlot(RooWorkspace *ws)
{
	TChain *ch = new TChain("tree");
    ch->Add("./fits/JPsi_Data_ForFit.root");
    TTree *tr = ch;
	std::cout << "calculate sWeight" << std::endl;

	// get what you need out of the workspace
	RooAbsPdf *model = ws->pdf("model");
	RooRealVar *nsig =  ws->var("nsig");
	RooRealVar *nbkgE = ws->var("nbkgE");
	RooDataSet *data = (RooDataSet *)ws->data("data");

	//fit model to data
	model->fitTo(*data);

	RooMsgService::instance().setSilentMode(true);

	// now make use of the SPlot class to add s_weight to out dataset
	// based on out model and yields
	RooStats::SPlot *sData = new RooStats::SPlot("sData", "An SPlot", *data, model, RooArgList(*nsig, *nbkgE));

	// create new root file for s_weighted output
    TFile *outFile = new TFile("JPsi_DataSet_with_sWeight.root", "RECREATE");
    outFile->cd(); 

    //cloning original tree
    TTree *tr1 = tr->CloneTree(); 

    double s_wt;
    double b_wt;

    TBranch *_s_wt = tr1->Branch("s_wt", &s_wt, "s_wt/D");
    TBranch *_b_wt = tr1->Branch("b_wt", &b_wt, "b_wt/D");

	// check that out weigt have the desired property
	std::cout << "Check SWeights:" << std::endl;

	std::cout << std::endl
		  << "N_{sig} " << nsig->getVal() << ": From s_weight N_{sig}"
		  << sData->GetYieldFromSWeight("nsig") <<std::endl;

	std::cout << std::endl
                  << "N_{bkg} " << nbkgE->getVal() << ": From s_weight N_{bkg}"
                  << sData->GetYieldFromSWeight("nbkgE") <<std::endl;

	for (Int_t i = 0; i < tr->GetEntries(); i++)
	 {
		if (i%1000 == 0) std::cout << "--- ... Processing event: " << i << std::endl;
      	s_wt = sData->GetSWeight( i, "nsig");
      	b_wt = sData->GetSWeight( i, "nbkgE");

      	_s_wt->Fill();
     	_b_wt->Fill();
    
     }

	 std::cout << "Loop Finished" << std::endl;


	tr1->Write();
    outFile->Write();
    std::cout << "Files written" << std::endl;
    outFile->Close();
}


	






