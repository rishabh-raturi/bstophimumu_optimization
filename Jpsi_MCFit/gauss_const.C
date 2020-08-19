#include <TLatex.h>
#include <TString.h>
#include <TLegend.h>
#include <RooCBShape.h>
#include "StdFitter.cc"

void constt()
{
	gROOT->SetBatch(1);
	gSystem->Load("libRooFit");
	using namespace RooFit;

	TChain *tr = new TChain("tree");
	tr->Add("/home/rishabh/project/root_files/val/2016/red_ntuples/JPsi_MC_ForFit.root");

	int nentries_ = tr->GetEntries();
	cout << "\n=> total entries in signal tree = " << nentries_ << endl;

	double bm_min(5.2), bm_max(5.6);

	RooRealVar Bmass("Bmass","#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}", bm_min, bm_max);
	RooArgSet  observables(Bmass);

	//create dataset from the input file
	RooDataSet data("data","dataset with Bmass", tr, observables);

	TCut c1 = Form("Bmass>%f && Bmass<%f",bm_min,bm_max);
        TCut cutTotal = c1  ;

	//reduced data set- apply the Bmass cut(fitting range)
	RooDataSet *redData = (RooDataSet*)data.reduce(cutTotal);
	std::cout<<"After final cut: "<<redData->sumEntries()<<std::endl;

	//define double crystal ball pdf for fitting
	RooRealVar  mean("mean","common means for Crystal Balls", 5.367, bm_min, bm_max);
	RooRealVar  sigma1("sigma1","sigma of CB1",  0.02, 0., 0.28);
	RooRealVar  sigma2("sigma2","sigma of CB2",  0.04, 0.00039, 0.49);
	RooRealVar  sigM_frac("sigM_frac","fraction of CB", 0.4, 0.01, 1.);
	RooRealVar  n1("n1", "",50, 1., 140.);
	RooRealVar  n2("n2", "",40, 1., 140.);
	RooRealVar  alpha1("alpha1","alpha for CB1", 2., 0.1, 4.);
	RooRealVar  alpha2("alpha2","alpha for CB2", -1., -5.0, 0.);

	RooCBShape CB1("CB1","Crystal Ball-1", Bmass, mean, sigma1, alpha1,n1);
        RooCBShape CB2("CB2","Crystal Ball-2", Bmass, mean, sigma2, alpha2,n2);

	RooAddPdf CB("CB","CB1+CB2", RooArgList(CB1,CB2), RooArgList(sigM_frac));
        RooRealVar nsig("nsig","nsig",1E3,500,1.1*tr->GetEntries());

	 //final model used for fitting
        RooExtendPdf model("model","model", CB, nsig);

	model.Print();

	cout<< "pdf evaluation for MC fitting" << endl;
        cout<< model.getLogVal() <<endl;

	RooFitResult* fitres = model.fitTo(*redData);

	TCanvas *c = new TCanvas("c","c",1250, 1500);
        TPad *p1   = new TPad("p1","p1", 0.01, 0.25, 0.995, 0.97);
        TPad *p2   = new TPad("p2","p2", 0.01, 0.02, 0.995, 0.24);
        p1->Draw();
        p2->Draw();

        p1->cd();
        gPad->SetLogy();

	// creating frame for Bmass invariant mass distribution
        RooPlot *xframe = Bmass.frame(Title(""), Bins(500));
        xframe->SetTitle("");
        xframe->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");

	// plotting data of the Bmass from file on the frame
        redData->plotOn(xframe,RooFit::Name("data"));

	// plotting final pdf on the frame
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

	TPaveText* paveText = new TPaveText(0.65,0.40,0.83,0.88,"NDC");
	paveText->SetBorderSize(0.0);
        paveText->SetFillColor(kWhite);
        paveText->SetFillStyle(0);
        paveText->SetTextSize(0.02);
	paveText->AddText(Form("N_{sig} = %.0f #pm %.0f", nsig.getVal(), nsig.getError()));
	paveText->AddText(Form("#mu  = %.7f #pm %.7f GeV" , mean.getVal() , mean.getError()));
	paveText->AddText(Form("#sigma_{1} = %.7f #pm %.7f GeV", sigma1.getVal(), sigma1.getError()));
	paveText->AddText(Form("#sigma_{2} = %.7f #pm %.7f GeV", sigma2.getVal(), sigma2.getError()));
	paveText->AddText(Form("#sigma_{eff} = %.7f  GeV", eff_sigma));
	paveText->AddText(Form("f = %.7f #pm %.7f ", sigM_frac.getVal(), sigM_frac.getError()));
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

	printf("\n #################################");
	printf("\n #################################");
	printf("\n #################################");

	cout << "JPsiKstar MC fitting started" << endl;

	printf("\n #################################");
	printf("\n #################################");
	printf("\n #################################");

	TChain *tr4 =new TChain("tree");
	tr4->Add("/home/rishabh/project/root_files/val/2016/red_ntuples/JKstar_MC_ForFit.root");
	RooDataSet MC("MC","dataset with Bmass", tr4, observables);
	RooDataSet *redMC = (RooDataSet*)MC.reduce(cutTotal);
	std::cout<<"After final cut: "<<redMC->sumEntries()<<std::endl;

	//define double crystal ball pd for JPsiKstar
	RooRealVar  md("md","common means for Crystal Balls", 5.41, 5.35, 5.45);
	RooRealVar  sig1("sig1","sigma of CB1",  0.04, 0., 1.);
	RooRealVar  sig2("sig2","sigma of CB2",  0.02, 0., 1.);
	RooRealVar  frac("frac","fraction of CB",  0.5, 0., 1.);
	RooRealVar  N1("N1", "",  1., 140.);
	RooRealVar  N2("N2", "",  1., 140.);
	RooRealVar  alp1("alp1","alpha for CB1", 0.1 , 5.);
        RooRealVar  alp2("alp2","alpha for CB2", -0.1, -5.0, 0.);

	RooCBShape cb1("cb1","Crystal Ball-1", Bmass, md, sig1, alp1, N1);
	RooCBShape cb2("cb2","Crystal Ball-1", Bmass, md, sig2, alp2, N2);

	RooAddPdf cb("cb","cb1+cb2", RooArgList(cb1,cb2), RooArgList(frac));
	RooRealVar nsigPB("nsigPB","nsigPB",1E3,500,1.1*tr4->GetEntries());

	//final model for peaking background 
	RooExtendPdf mod("mod","mod", cb, nsigPB);

	mod.Print();

	// fitting the model to the data
        RooFitResult* fitresults = mod.fitTo(*redMC);

	TCanvas *c2 = new TCanvas("c1","c1",1500,1500);
        TPad *pad1   = new TPad("p1","p1", 0.1, 0.25, 0.995, 0.97);
        TPad *pad2   = new TPad("p2","p2", 0.1, 0.02, 0.995, 0.24);
        pad1->Draw();
        pad2->Draw();

	pad1->cd();

	// creating frame for Bmass invariant mass distribution
        RooPlot *xframe1 = Bmass.frame(Title(""), Bins(80));
        xframe1->SetTitle("");
        xframe1->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");

	// plotting data of the Bmass from file on the frame
        redMC->plotOn(xframe1,RooFit::Name("MC_PB"));

	mod.plotOn(xframe1, RooFit::Name("MC_PDF"), LineColor(kBlue));

	// create a pull distribution for the fit
	RooHist* hpull1 = xframe1->pullHist();

	printf("**************************************");
        printf("\n #sigma_{1}: %.7f", sig1.getVal());
        printf("\n #sigma_{2}: %.7f", sig2.getVal());
        double sigma = sqrt(frac.getVal()*pow(sig1.getVal(),2) + (1 - frac.getVal()) * pow(sig2.getVal(),2));
        printf("\n #sigma_{eff}: %.7f", sigma);
        printf("\n************************************\n");

	double chi2dof1 = xframe1->chiSquare();
        std::cout<<"\n"<<std::endl;
        std::cout<<"#chi^{2}/dof= "<< chi2dof1 << std::endl;

	// plotting individual crystal ball in the frame
        mod.plotOn(xframe1, RooFit::Name("cb1"), Components("cb1"), LineStyle(kDashed), LineColor(kRed));
        mod.plotOn(xframe1, RooFit::Name("cb2"), Components("cb2"), LineStyle(kDashed), LineColor(kGreen));

	xframe1->Draw();

	TPaveText* paveT = new TPaveText(0.65,0.45,0.83,0.88,"NDC");
        paveT->SetBorderSize(0.0);
        paveT->SetFillColor(kWhite);
        paveT->SetFillStyle(0);
        paveT->SetTextSize(0.02);
        paveT->AddText(Form("N_{sig} = %.0f #pm %.0f", nsigPB.getVal(), nsigPB.getError()));
        paveT->AddText(Form("#mu  = %.7f #pm %.7f GeV" , md.getVal() , md.getError()));
        paveT->AddText(Form("#sigma_{1} = %.7f #pm %.7f GeV", sig1.getVal(), sig1.getError()));
        paveT->AddText(Form("#sigma_{2} = %.7f #pm %.7f GeV", sig2.getVal(), sig2.getError()));
        paveT->AddText(Form("#sigma_{eff} = %.7f  GeV", sigma));
        paveT->AddText(Form("f = %.7f #pm %.7f ", frac.getVal(), frac.getError()));
        paveT->AddText(Form("#alpha_{1} = %.7f #pm %.7f ", alp1.getVal(), alp1.getError()));
        paveT->AddText(Form("#alpha_{2} = %.7f #pm %.7f ", alp2.getVal(), alp2.getError()));
        paveT->AddText(Form("n_{1} = %.7f #pm %.7f", N1.getVal(), N1.getError()));
        paveT->AddText(Form("n_{2} = %.7f #pm %.7f", N2.getVal(), N2.getError()));
        paveT->AddText(Form("#chi^{2}/dof  = %.5f ", chi2dof1 ));

        paveT->Draw();

	TLatex *mk = new TLatex();
        mk->SetNDC(true);
    	mk->SetTextFont(42);
        mk->SetTextSize(0.035);
   	mk->DrawLatex(0.12,startY,"#bf{CMS} #it{Simulation}");
	mk->DrawLatex(0.75,startY,"#scale[0.8]{45.6 fb^{-1} (13 TeV)}");

        mk->Draw();

	pad2->cd();

    	RooPlot* pull1 = Bmass.frame(Title(""));
    	pull1->addPlotable(hpull1, "BX");
   	pull1->SetMinimum(-5);
   	pull1->SetMaximum(5);
        hpull1->SetFillColor(kBlue);
    	pull1->GetYaxis()->SetNdivisions(5);
    	pull1->SetTitle("#bf{Pull Distribution}");
        pull1->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");
        pull1->GetYaxis()->SetTitle("");
        pull1->GetYaxis()->SetLabelSize( (0.7/0.3)*xframe1->GetXaxis()->GetLabelSize() );
        pull1->GetXaxis()->SetLabelSize( (0.7/0.3)*xframe1->GetXaxis()->GetLabelSize() );
    	pull1->GetXaxis()->SetTitle("");
    	pull1->GetXaxis()->SetLimits(bm_min, bm_max);
    	pull1->Draw("P");

	l1->Draw();
        l2->Draw();
        l3->Draw();

	c->SaveAs("./plots/JPsi_MCFit_2016_Log.pdf");
	c->SaveAs("./plots/JPsi_MCFit_2016_Log.png");

	c2->SaveAs("./plots/Kstar_MCFit_2016_Log.pdf");
	c2->SaveAs("./plots/Kstar_MCFit_2016_Log.png");

	// Data fitting 
	TChain *trd = new TChain("tree");
	trd->Add("/home/rishabh/project/root_files/val/2016/red_ntuples/JPsi_Data_ForFit.root");
	RooDataSet Data("Data","data sample 2016", trd, observables);
	RooDataSet *RedData = (RooDataSet*)Data.reduce(cutTotal);

	//defining model for signal in data

	RooRealVar nsigD("nsigD","nsigD",1E3,500,1.1*trd->GetEntries());

	RooRealVar nsigpb("nsigpb","nsigpb",1E3,500,1.1*tr4->GetEntries());
	
	//defining exponential pdf for combinatorial bakcground
	RooRealVar lambda("lambda","slope",-10.,0.);
        RooExponential bkgE("bkgE","exponential PDF", Bmass, lambda);

	RooRealVar nbkgE("nbkgE","number of exponential bkg events", 8000, 0, 1E7);

	RooAddPdf Model("Model", "CB+cb+e", RooArgList(CB, cb, bkgE), RooArgList(nsigD, nsigpb, nbkgE));

	Model.Print("t");

	//gaussian constraining parameters for signal shape
	RooGaussian   gaus_mean("gaus_mean", "gaus_mean", mean, RooConst(mean.getVal()), RooConst(mean.getError()));
        RooGaussian   gaus_sigma1("gaus_sigma1", "gaus_sigma1", sigma1, RooConst(sigma1.getVal()), RooConst(sigma1.getError()));
        RooGaussian   gaus_sigma2("gaus_sigma2", "gaus_sigma2", sigma2, RooConst(sigma2.getVal()), RooConst(sigma2.getError()));
        RooGaussian   gaus_alpha1("gaus_alpha1", "gaus_alpha1", alpha1, RooConst(alpha1.getVal()), RooConst(alpha1.getError()));
        RooGaussian   gaus_alpha2("gaus_alpha2", "gaus_alpha2", alpha2, RooConst(alpha2.getVal()), RooConst(alpha2.getError()));
        RooGaussian   gaus_n1("gaus_n1", "gaus_n1", n1, RooConst(n1.getVal()), RooConst(n1.getError()));
        RooGaussian   gaus_n2("gaus_n2", "gaus_n2", n2, RooConst(n2.getVal()), RooConst(n2.getError()));
        RooGaussian   gaus_frac("gaus_frac", "gaus_frac", sigM_frac, RooConst(sigM_frac.getVal()), RooConst(sigM_frac.getError()));

	//gaussian constraining the parameters for peaking background shape
	RooGaussian   g_npb("g_npb", "gaus_npb", nsigPB, RooConst(nsigPB.getVal()*0.78728), RooConst(nsigPB.getError()*0.78728));
        RooGaussian   g_sigma1("g_sigma1", "gaus_sigma1", sig1, RooConst(sig1.getVal()), RooConst(sig1.getError()));
        RooGaussian   g_sigma2("g_sigma2", "gaus_sigma2", sig2, RooConst(sig2.getVal()), RooConst(sig2.getError()));
        RooGaussian   g_alpha1("g_alpha1", "gaus_alpha1", alp1, RooConst(alp1.getVal()), RooConst(alp1.getError()));
        RooGaussian   g_alpha2("g_alpha2", "gaus_alpha2", alp2, RooConst(alp2.getVal()), RooConst(alp2.getError()));
        RooGaussian   g_n1("g_n1", "gaus_n1", N1, RooConst(N1.getVal()), RooConst(N1.getError()));
        RooGaussian   g_n2("g_n2", "gaus_n2", N2, RooConst(N2.getVal()), RooConst(N2.getError()));
        RooGaussian   g_frac("g_frac", "gaus_frac", frac, RooConst(frac.getVal()), RooConst(frac.getError()));

	RooArgSet gausConstraints(gaus_sigma1, gaus_alpha1, gaus_n1, gaus_sigma2, gaus_alpha2, gaus_n2, gaus_frac);
	gausConstraints.add(RooArgSet(g_sigma1, g_sigma2, g_alpha1, g_alpha2, g_n1, g_n2, g_frac, g_npb));

	RooFitResult* fitres1 = Model.fitTo(*RedData, Extended(kTRUE), Save(kTRUE), ExternalConstraints(gausConstraints));

	cout << "fit result with constraint" << endl;
        fitres1->Print("v");

	TCanvas *c3 = new TCanvas("c3", "c3", 1250, 1500);
        TPad *P1   = new TPad("P1","P1", 0.01, 0.25, 0.995, 0.97);
        TPad *P2   = new TPad("P2","P2", 0.01, 0.02, 0.995, 0.24);
        P1->Draw();
        P2->Draw();

	P1->cd();
        RooPlot *xframe2 = Bmass.frame(Title(""), Bins(100));
        xframe2->SetTitle("");
        xframe2->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");
        xframe2->GetYaxis()->SetTitleOffset(1.55);

	// plotting data of the Bmass from file on the frame
        RedData->plotOn(xframe2,RooFit::Name("Data"), MarkerStyle(8));

	//plot the model on the frame
        Model.plotOn(xframe2, RooFit::Name("Full PDF_data"), LineColor(kBlue), LineWidth(3));

	// create a pull distribution for the fit
        RooHist* hpull2 = xframe2->pullHist();

	printf("\n **************************************");
        printf("\n sigma1: %.5f", sigma1.getVal());
        printf("\n sigma2: %.5f", sigma2.getVal());
        double gaus_eff_sigma2      = sqrt(sigM_frac.getVal()*pow(sigma1.getVal(),2) + (1 - sigM_frac.getVal()) * pow(sigma2.getVal(),2));
        printf("\n eff_sigma: %.5f", gaus_eff_sigma2);
        printf("\n************************************\n");

	double chi2dof2 = xframe2->chiSquare();
        std::cout<<"\n"<<std::endl;
        std::cout<<"Data Fit-> #chi^{2}/dof= "<< chi2dof2 << std::endl;

	// plotting total crystal ball in the frame
        Model.plotOn(xframe2, RooFit::Name("CBall_data"), Components("CB"), LineColor(kGreen), LineWidth(3));
        Model.plotOn(xframe2, RooFit::Name("bkgE_data"), Components("bkgE"), LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3003));
        Model.plotOn(xframe2, RooFit::Name("cb"), Components("cb"), LineColor(kMagenta), DrawOption("F"), FillColor(kMagenta), FillStyle(3008));

	xframe2->Draw();

	TLegend *leg = new TLegend(0.15, 0.65, 0.28, 0.88);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(0.0);
        leg->AddEntry(xframe2->FindObject("Data"), "Data", "ep");
        leg->AddEntry("Full PDF_data", "Fit", "L");
        leg->AddEntry("CBall_data", "Signal", "L");
        leg->AddEntry("bkgE_data", "Background", "F");
	leg->AddEntry("cb", "Peaking backgound", "F");
        xframe2->addObject(leg);
        xframe2->Draw();
	
	TPaveText* pT = new TPaveText(0.65,0.30,0.83,0.88,"NDC");
        pT->SetBorderSize(0.0);
        pT->SetFillColor(kWhite);
        pT->SetFillStyle(0);
        pT->SetTextSize(0.02);

        pT->AddText(Form("N_{sig}  = %.0f #pm %.0f", nsigD.getVal(), nsigD.getError()));
        pT->AddText(Form("N_{comb}  = %.0f #pm %.0f", nbkgE.getVal(), nbkgE.getError()));
        pT->AddText(Form("N_{peak. bkg}  = %.0f #pm %.0f", nsigpb.getVal(), nsigpb.getError()));
        pT->AddText(Form("#mu  = %.6f #pm %.6f GeV" , mean.getVal() , mean.getError()));
        pT->AddText(Form("#sigma_{1} = %.6f #pm %.6f GeV", sigma1.getVal(), sigma1.getError()));
        pT->AddText(Form("#sigma_{2} = %.6f #pm %.6f GeV", sigma2.getVal(), sigma2.getError()));
        pT->AddText(Form("#sigma_{eff} = %.6f  GeV", gaus_eff_sigma2));
        pT->AddText(Form("frac = %.6f #pm %.6f ", sigM_frac.getVal(), sigM_frac.getError()));
        pT->AddText(Form("slope = %.6f #pm %.6f", lambda.getVal(), lambda.getError()));
        pT->AddText(Form("#alpha_{1} = %.7f #pm %.7f ", alpha1.getVal(), alpha1.getError()));
        pT->AddText(Form("#alpha_{2} = %.7f #pm %.7f ", alpha2.getVal(), alpha2.getError()));
        pT->AddText(Form("n_{1} = %.7f #pm %.7f", n1.getVal(), n1.getError()));
        pT->AddText(Form("n_{2} = %.7f #pm %.7f", n2.getVal(), n2.getError()));
        pT->AddText(Form("#chi^{2}/dof  = %.5f  ", chi2dof2 ));

	pT->Draw();

	TLatex *mark1 = new TLatex();
        mark1->SetNDC(true);


        mark1->SetTextFont(42);
        mark1->SetTextSize(0.035);
        mark1->DrawLatex(0.12,startY,"#bf{CMS} #it{Preliminary}");

        mark1->DrawLatex(0.65,startY,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
        mark1->Draw();

	P2->cd();

        RooPlot* pull2 = Bmass.frame(Title(""));
        pull2->addPlotable(hpull2, "P");
        pull2->SetMinimum(-5);
        pull2->SetMaximum(5);
        pull2->GetYaxis()->SetNdivisions(5);
        pull2->SetTitle("#bf{Pull Distribution}");
        pull2->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");
        pull2->GetXaxis()->SetTitleOffset(1.2);
        pull2->GetYaxis()->SetTitle("");
        pull2->GetXaxis()->SetTitleSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
        pull2->GetYaxis()->SetLabelSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
        pull2->GetXaxis()->SetLabelSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
        pull2->GetXaxis()->SetLimits(bm_min, bm_max);
        pull2->Draw("P");

	TLine* L1 = new TLine(5.2, 3, 5.6, 3);
        TLine* L2 = new TLine(5.2, -3, 5.6, -3);
        L1->SetLineColor(2);
        L2->SetLineColor(2);
        TLine* L3 = new TLine(5.2, 0, 5.6, 0);
        L3->SetLineColor(2);
        L1->Draw();
        L2->Draw();
        L3->Draw();
	
	c3->SaveAs("./plots/JPsi_DataFit_2016.pdf");	
	c3->SaveAs("./plots/JPsi_DataFit_2016.png");


}







