// This script performs the MC Fit using double sided crystal ball pdf.
// Data fit is obtained after gaussian constraining the parameters from the MC fit.
// 

// \author- Rishabh Raturi
// IIT Bhubaneswar



#include <TLatex.h>
#include <TString.h>
#include <TLegend.h>
#include <RooCBShape.h>

void gauss_const()
{
	gSystem->Load("libRooFit");
    	using namespace RooFit;

    	TChain *ch = new TChain("tree");
    	ch->Add("/home/rishabh/project/root_files/DCA/red_ntuples/PhiMuMu_MC_ForFit.root");
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
	RooRealVar  n2("n2", "", 1., 180.);
	RooRealVar  alpha1("alpha1","alpha for CB1", 0.1, 5.);
	RooRealVar  alpha2("alpha2","alpha for CB2", -0.1, -5.0, 0.);

	RooCBShape CB1("CB1","Crystal Ball-1", Bmass, mean, sigma1, alpha1,n1);
	RooCBShape CB2("CB2","Crystal Ball-2", Bmass, mean, sigma2, alpha2,n2);

	 
	RooAddPdf CB("CB","CB1+CB2", RooArgList(CB1,CB2), RooArgList(sigM_frac));
	RooRealVar nsig("nsig","nsig",1E3,500,1.1*tr->GetEntries());
	RooExtendPdf model("model","model", CB, nsig);


	model.Print("t");

	// fitting the model to the data
	RooFitResult* fitres = model.fitTo(*redData); 
	

	TCanvas *c = new TCanvas("c","c",1250, 1500);
	TPad *p1   = new TPad("p1","p1", 0.01, 0.25, 0.995, 0.97);
  	TPad *p2   = new TPad("p2","p2", 0.01, 0.02, 0.995, 0.24);
  	p1->Draw();
  	p2->Draw();

	p1->cd();
	gPad->SetLogy();

	// creating frame for Bmass invariant mass distribution
	RooPlot *xframe = Bmass.frame(Title(""), Bins(400));
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

	TCanvas *c2 = new TCanvas("c2","c2",1250, 1500);
        TPad *P1   = new TPad("P1","P1", 0.01, 0.25, 0.995, 0.97);
        TPad *P2   = new TPad("P2","P2", 0.01, 0.02, 0.995, 0.24);
        P1->Draw();
        P2->Draw();

        P1->cd();
        xframe->Draw();
	paveText->Draw();
	mark->Draw();

	P2->cd();
	pull->Draw("P");
	l1->Draw();
        l2->Draw();
        l3->Draw();

	// using gaussian constrained parameters from MC to fit the data
	TChain *ch1 = new TChain("tree");
	ch1->Add("/home/rishabh/project/root_files/DCA/red_ntuples/PhiMuMu_Data_ForFit.root");
	TTree *tr1 = ch1;
	RooDataSet Data("Data","data sample 2016", ch1, observables);
	RooDataSet *RedData = (RooDataSet*)Data.reduce(cutTotal);

	RooRealVar nsigD("nsigD","nsigD",1E3,100,1.1*tr1->GetEntries());
	RooRealVar lambda("lambda","slope",-10.,0.);
        RooExponential bkgE("bkgE","exponential PDF", Bmass, lambda);
	RooRealVar nbkgE("nbkgE","number of exponential bkg events", 8000, 0, 1E7);

	//create final model for fitting
        RooAddPdf Model("Model", "CB+e", RooArgList(CB,bkgE), RooArgList(nsigD,nbkgE));
        Model.Print("t");

	RooGaussian   gaus_mean("gaus_mean", "gaus_mean", mean, RooConst(mean.getVal()), RooConst(mean.getError()));
	RooGaussian   gaus_sigma1("gaus_sigma1", "gaus_sigma1", sigma1, RooConst(sigma1.getVal()), RooConst(sigma1.getError()));
	RooGaussian   gaus_sigma2("gaus_sigma2", "gaus_sigma2", sigma2, RooConst(sigma2.getVal()), RooConst(sigma2.getError()));
	RooGaussian   gaus_alpha1("gaus_alpha1", "gaus_alpha1", alpha1, RooConst(alpha1.getVal()), RooConst(alpha1.getError()));
	RooGaussian   gaus_alpha2("gaus_alpha2", "gaus_alpha2", alpha2, RooConst(alpha2.getVal()), RooConst(alpha2.getError())); 
	RooGaussian   gaus_n1("gaus_n1", "gaus_n1", n1, RooConst(n1.getVal()), RooConst(n1.getError()));
	RooGaussian   gaus_n2("gaus_n2", "gaus_n2", n2, RooConst(n2.getVal()), RooConst(n2.getError()));
	RooGaussian   gaus_frac("gaus_frac", "gaus_frac", sigM_frac, RooConst(sigM_frac.getVal()), RooConst(sigM_frac.getError())); 

	RooArgSet gausConstraints(gaus_sigma1, gaus_alpha1, gaus_n1, gaus_sigma2, gaus_alpha2, gaus_n2, gaus_frac);
	
	RooFitResult* fitres1 = Model.fitTo(*RedData, Extended(kTRUE), Save(kTRUE), Minimizer("Minuit"), Minos(kTRUE), ExternalConstraints(gausConstraints));

	TCanvas *c3 = new TCanvas("c3", "c3", 1250, 1500);
	TPad *pad1   = new TPad("pad1","pad1", 0.01, 0.25, 0.995, 0.97);
        TPad *pad2   = new TPad("pad2","pad2", 0.01, 0.02, 0.995, 0.24);
        pad1->Draw();
        pad2->Draw();

	pad1->cd();
	RooPlot *xframe1 = Bmass.frame(Title(""), Bins(50));
	xframe1->SetTitle("");
        xframe1->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");
        xframe1->GetYaxis()->SetTitleOffset(1.55);

	// plotting data of the Bmass from file on the frame
        RedData->plotOn(xframe1,RooFit::Name("Data"), MarkerStyle(8));

	//plot the model on the frame
	Model.plotOn(xframe1, RooFit::Name("Full PDF_data"), LineColor(kBlue), LineWidth(3));

	// create a pull distribution for the fit
        RooHist* hpull1 = xframe1->pullHist();

	printf("**************************************");
        printf("\n sigma1: %.5f", sigma1.getVal());
        printf("\n sigma2: %.5f", sigma2.getVal());
	double gaus_eff_sigma      = sqrt(sigM_frac.getVal()*pow(sigma1.getVal(),2) + (1 - sigM_frac.getVal()) * pow(sigma2.getVal(),2));
        printf("\n eff_sigma: %.5f", gaus_eff_sigma);
        printf("\n************************************\n");

	double chi2dof1 = xframe1->chiSquare();
        std::cout<<"\n"<<std::endl;
        std::cout<<"Data Fit-> #chi^{2}/dof= "<< chi2dof1 << std::endl;

	// plotting total crystal ball in the frame
        Model.plotOn(xframe1, RooFit::Name("CBall_data"), Components("CB"), LineColor(kGreen), LineWidth(3));
        Model.plotOn(xframe1, RooFit::Name("bkgE_data"), Components("bkgE"), LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3003));

	Bmass.setRange("signal", 5.293135, 5.440865);
	RooAbsReal* fsig_bkg = bkgE.createIntegral(Bmass, NormSet(Bmass), Range("signal"));
	
	double bkg = fsig_bkg->getVal()*nbkgE.getVal();
        double err_bkg = nbkgE.getError()*fsig_bkg->getVal();

	cout<< "expected number of background events in signal region-> "<< bkg<<endl;
        cout<< "error-> "<< err_bkg<<endl;

	TLine* line1 =  new TLine(5.293135, 0, 5.293135, 30);
        line1->SetLineColor(2);
        line1->SetLineWidth(2);

        TLine* line2 =  new TLine(5.440865, 0, 5.440865, 30);
        line2->SetLineColor(2);
        line2->SetLineWidth(2);

        xframe1->addObject(line1);
        xframe1->addObject(line2);

	xframe1->Draw();

	TLegend *leg = new TLegend(0.15, 0.65, 0.28, 0.88);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(0.0);
        leg->AddEntry(xframe->FindObject("Data"), "Data", "ep");
        leg->AddEntry("Full PDF_data", "Fit", "L");
        leg->AddEntry("CBall_data", "Signal", "L");
        leg->AddEntry("bkgE_data", "Background", "F");
        xframe1->addObject(leg);
        xframe1->Draw();

        TPaveText* pT = new TPaveText(0.65,0.60,0.83,0.88,"NDC");
	pT->SetBorderSize(0.0);
        pT->SetFillColor(kWhite);
        pT->SetFillStyle(0);
        pT->SetTextSize(0.02);

	pT->AddText(Form("N_{sig}  = %.0f #pm %.0f", nsigD.getVal(), nsigD.getError()));
	pT->AddText(Form("N_{bkg}  = %.0f #pm %.0f", nbkgE.getVal(), nbkgE.getError()));
	pT->AddText(Form("#mu  = %.6f #pm %.6f GeV" , mean.getVal() , mean.getError()));
	pT->AddText(Form("#sigma_{1} = %.6f #pm %.6f GeV", sigma1.getVal(), sigma1.getError()));
	pT->AddText(Form("#sigma_{2} = %.6f #pm %.6f GeV", sigma2.getVal(), sigma2.getError()));
	pT->AddText(Form("#sigma_{eff} = %.6f  GeV", gaus_eff_sigma));
	pT->AddText(Form("frac = %.6f #pm %.6f ", sigM_frac.getVal(), sigM_frac.getError()));
	pT->AddText(Form("slope = %.6f #pm %.6f", lambda.getVal(), lambda.getError()));
	pT->AddText(Form("#alpha_{1} = %.7f #pm %.7f ", alpha1.getVal(), alpha1.getError()));
        pT->AddText(Form("#alpha_{2} = %.7f #pm %.7f ", alpha2.getVal(), alpha2.getError()));
        pT->AddText(Form("n_{1} = %.7f #pm %.7f", n1.getVal(), n1.getError()));
        pT->AddText(Form("n_{2} = %.7f #pm %.7f", n2.getVal(), n2.getError()));
	pT->AddText(Form("#chi^{2}/dof  = %.5f ", chi2dof1 ));

	pT->Draw();

	TLatex *mark1 = new TLatex();
        mark1->SetNDC(true);

        
        mark1->SetTextFont(42);
        mark1->SetTextSize(0.035);
        mark1->DrawLatex(0.12,startY,"#bf{CMS} #it{Preliminary}");

        mark1->DrawLatex(0.65,startY,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
        mark1->Draw();

	pad2->cd();

	RooPlot* pull1 = Bmass.frame(Title(""));
        pull1->addPlotable(hpull1, "P");
        pull1->SetMinimum(-5);
        pull1->SetMaximum(5);
        pull1->GetYaxis()->SetNdivisions(5);
        pull1->SetTitle("#bf{Pull Distribution}");
        pull1->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");
        pull1->GetXaxis()->SetTitleOffset(1.2);
        pull1->GetYaxis()->SetTitle("");
        pull1->GetXaxis()->SetTitleSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
        pull1->GetYaxis()->SetLabelSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
        pull1->GetXaxis()->SetLabelSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
        pull1->GetXaxis()->SetLimits(bm_min, bm_max);
        pull1->Draw("P");

	TLine* L1 = new TLine(5.2, 3, 5.6, 3);
        TLine* L2 = new TLine(5.2, -3, 5.6, -3);
        L1->SetLineColor(2);
        L2->SetLineColor(2);
        TLine* L3 = new TLine(5.2, 0, 5.6, 0);
        L3->SetLineColor(2);
        L1->Draw();
        L2->Draw();
        L3->Draw();

}

