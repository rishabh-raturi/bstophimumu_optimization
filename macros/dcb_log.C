// This script performs the MC Fit using double sided crystal ball pdf.
// Using the initialization of the paramters as shown in the code, for rare mode MC fitting we are getiing
// the best fit.

// \author- Rishabh Raturi



#include <TLatex.h>
#include <TString.h>
#include <TLegend.h>
#include <RooCBShape.h>

void dcb_log()
{
	gSystem->Load("libRooFit");
    using namespace RooFit;

    TChain *ch = new TChain("tree");
    ch->Add("/home/rishabh/project/root_files/IsoVar/2017/ntuples/sel_BsToPhiMuMu_OfficialMC_signal_2017Mini_Presel_mc.lite_cut0_dcacut_s02.root");
    TTree *tr = ch;

	int nentries_ = tr->GetEntries();
    cout << "\n=> total entries in signal tree = " << nentries_ << endl;

	double bm_min(5.2), bm_max(5.6);

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
    RooRealVar Triggers("Triggers","",0.,10.);

	RooArgSet  observables(Bmass,Mumumass,Kmpt,Kppt,Kmtrkdcasigbs,Kptrkdcasigbs,Blxysig,Bcosalphabs2d,Bvtxcl);
	observables.add(RooArgSet(Phimass,Triggers));

	RooDataSet data("data","dataset with Bmass", ch, observables);

	TCut c1 = Form("Bmass>%f && Bmass<%f",bm_min,bm_max);
    TCut cutTotal = c1  ;

    RooDataSet *redData = (RooDataSet*)data.reduce(cutTotal);   //redData consists data to be fitted
	
	RooRealVar  mean("mean","common means for Crystal Balls", 5.367, bm_min, bm_max);
	RooRealVar  sigma1("sigma1","sigma of CB1",  0.04, 0., 1.);
	RooRealVar  sigma2("sigma2","sigma of CB2",  0.02, 0., 1.);
	RooRealVar  sigM_frac("sigM_frac","fraction of CB",  0., 1.);
	RooRealVar  n1("n1", "", 1., 140.);
	RooRealVar  n2("n2", "", 1., 140.);
	RooRealVar  alpha1("alpha1","alpha for CB1", 0.1, 5.);
	RooRealVar  alpha2("alpha2","alpha for CB2", -0.1, -5.0, 0.);

	RooCBShape CB1("CB1","Crystal Ball-1", Bmass, mean, sigma1, alpha1,n1);
	RooCBShape CB2("CB2","Crystal Ball-2", Bmass, mean, sigma2, alpha2,n2);

	 
	RooAddPdf CB("CB","CB1+CB2", RooArgList(CB1,CB2), RooArgList(sigM_frac));
	RooRealVar nsig("nsig","nsig",1E3,500,1.1*tr->GetEntries());
	RooExtendPdf model("model","model", CB, nsig);


	model.Print("t");

	// fitting the model to the data
	RooFitResult* fitres = model.fitTo(*redData, Extended(true), Minos(true)); 
	

	TCanvas *c = new TCanvas("c","c",1250,1400);
	TPad *p1   = new TPad("p1","p1", 0.1, 0.25, 0.995, 0.97);
  	TPad *p2   = new TPad("p2","p2", 0.1, 0.02, 0.995, 0.24);
  	p1->Draw();
  	p2->Draw();

	p1->cd();
	gPad->SetLogy();

	// creating frame for Bmass invariant mass distribution
	RooPlot *xframe = Bmass.frame(Title(""), Bins(80));
	xframe->SetTitle("");
	xframe->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");	

	// plotting data of the Bmass from file on the frame
	redData->plotOn(xframe,RooFit::Name("data"));

	// plot the model defined which is the addition of two crystal ball for signal and exponential for background on the frame
	model.plotOn(xframe, RooFit::Name("Full PDF"), LineColor(kBlue));

	// create a pull distribution for the fit 
	RooHist* hpull = xframe->pullHist();


	printf("**************************************");
	printf("\n #sigma_{1}: %.7f", sigma1.getVal());
	printf("\n #sigma_{2}: %.7f", sigma2.getVal());
	double eff_sigma = sqrt(sigM_frac.getVal()*pow(sigma1.getVal(),2) + (1 - sigM_frac.getVal()) * pow(sigma2.getVal(),2));
	printf("\n #sigma_{eff}: %.7f", eff_sigma);
	printf("\n************************************\n");

	double chi2dof = xframe->chiSquare();
	std::cout<<"\n"<<std::endl;
	std::cout<<"#chi^{2}/dof= "<< chi2dof << std::endl;

	// plotting individual crystal ball in the frame 
	model.plotOn(xframe, RooFit::Name("CB1"), Components("CB1"), LineStyle(kDashed), LineColor(kRed));
	model.plotOn(xframe, RooFit::Name("CB2"), Components("CB2"), LineStyle(kDashed), LineColor(kGreen));

	xframe->Draw();

	TPaveText* paveText = new TPaveText(0.65,0.60,0.83,0.88,"NDC");
  	paveText->SetBorderSize(0.0);
  	paveText->SetFillColor(kWhite);
  	paveText->SetFillStyle(0);
  	paveText->SetTextSize(0.02);
  	paveText->AddText(Form("Sig. Yield = %.0f #pm %.0f", nsig.getVal(), nsig.getError()));	
  	paveText->AddText(Form("mean  = %.7f #pm %.7f GeV" , mean.getVal() , mean.getError()));
  	paveText->AddText(Form("#sigma_{1} = %.7f #pm %.7f GeV", sigma1.getVal(), sigma1.getError()));
  	paveText->AddText(Form("#sigma_{2} = %.7f #pm %.7f GeV", sigma2.getVal(), sigma2.getError()));
  	paveText->AddText(Form("#sigma_{eff} = %.7f  GeV", eff_sigma));
  	paveText->AddText(Form("frac = %.7f #pm %.7f ", sigM_frac.getVal(), sigM_frac.getError()));
	paveText->AddText(Form("#alpha_{1} = %.7f #pm %.7f ", alpha1.getVal(), alpha1.getError()));
	paveText->AddText(Form("#alpha_{2} = %.7f #pm %.7f ", alpha2.getVal(), alpha2.getError()));	
	paveText->AddText(Form("n_{1} = %.7f #pm %.7f", n1.getVal(), n1.getError()));
        paveText->AddText(Form("n_{2} = %.7f #pm %.7f", n2.getVal(), n2.getError()));
  	paveText->AddText(Form("#chi^{2}/dof  = %.5f ", chi2dof ));

  	paveText->Draw();

	TLatex *mark = new TLatex();
	mark->SetNDC(true);
    double startY = 0.92;
    mark->SetTextFont(42);
	mark->SetTextSize(0.035);
    mark->DrawLatex(0.12,startY,"#bf{CMS} #it{Simulation}");

    mark->DrawLatex(0.75,startY,"#scale[0.8]{66226.56 fb^{-1} (13 TeV)}");
        
	mark->Draw();

  	p2->cd();

    RooPlot* pull = Bmass.frame(Title(""));
    pull->addPlotable(hpull, "BX");
    pull->SetMinimum(-5);
    pull->SetMaximum(5);
	hpull->SetFillColor(kBlue);
    pull->GetYaxis()->SetNdivisions(5);
    pull->SetTitle("#bf{Pull Distribution}");
	pull->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");
	pull->GetYaxis()->SetTitle("");
	pull->GetYaxis()->SetLabelSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
	pull->GetXaxis()->SetLabelSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
    pull->GetXaxis()->SetTitle("");
    pull->GetXaxis()->SetLimits(bm_min, bm_max);
    pull->Draw("P");

	TLine* l1 = new TLine(5.2, 3, 5.6, 3);
	TLine* l2 = new TLine(5.2, -3, 5.6, -3);
	l1->SetLineColor(2);
	l2->SetLineColor(2);
	TLine* l3 = new TLine(5.2, 0, 5.6, 0);
	l3->SetLineColor(2);
	l1->Draw();
	l2->Draw();
	l3->Draw();

	c->SaveAs("rare_2017_MCFit_raw.pdf");


}

