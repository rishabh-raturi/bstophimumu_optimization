#include <TLatex.h>
#include <TString.h>
#include <TLegend.h>
#include <RooCBShape.h>

void Jpsi_MCFit_Log()
{
	gSystem->Load("libRooFit");
        using namespace RooFit;

        TChain *ch = new TChain("tree");
	ch->Add("/home/rishabh/project/root_files/val/2016/red_ntuples/JPsi_MC_ForFit.root");
	TTree* tr=ch;

	int nentries_ = tr->GetEntries();
        cout << "\n=> total entries in resonance MC = " << nentries_ << endl;

        double bm_min(5.2), bm_max(5.6);

	RooRealVar Bmass("Bmass","#bf{m(KK#mu#mu) [GeV]}", bm_min, bm_max);
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

	RooDataSet *redData = (RooDataSet*)data.reduce(cutTotal);  // redData consists data to be fitted

	RooRealVar  mean("mean","common means for Crystal Balls", 5.367, bm_min, bm_max);
        RooRealVar  sigma1("sigma1","sigma of CB1", 0.02, 0., 0.029);
        RooRealVar  sigma2("sigma2","sigma of CB2", 0.04, 0., 0.049);
        RooRealVar  sigM_frac("sigM_frac","fraction of CB", 0.5, 0.4, 0.7);
        RooRealVar  n1("n1", "", 1., 140.);
	RooRealVar  n2("n2", "", 1., 100.);
        RooRealVar  alpha1("alpha1","alpha for CB1", 0.1, 4.);
        RooRealVar  alpha2("alpha2","alpha for CB2", -0.1, -5.0, 0.);

	RooCBShape CB1("CB1","Crystal Ball-1", Bmass, mean, sigma1, alpha1,n1);
        RooCBShape CB2("CB2","Crystal Ball-2", Bmass, mean, sigma2, alpha2,n2);

	RooAddPdf CB("CB","CB1+CB2", RooArgList(CB1,CB2), RooArgList(sigM_frac));
        RooRealVar nsig("nsig","nsig",1E3,500,1.1*tr->GetEntries());
        RooExtendPdf model("model","model", CB, nsig);

	model.Print("t");

	// fitting model to the data
	RooFitResult* fitres = model.fitTo(*redData, Extended(true));

	TCanvas *c = new TCanvas("c","c",1400,1400);
        TPad *p1   = new TPad("p1","p1", 0.1, 0.25, 0.995, 0.97);
        TPad *p2   = new TPad("p2","p2", 0.1, 0.02, 0.995, 0.24);
        p1->Draw();
        p2->Draw();

        p1->cd();
	gPad->SetLogy();

	// creating frame for Bmass invariant mass distribution
        RooPlot *xframe = Bmass.frame(Title(""), Bins(500) );
        xframe->SetTitle("");
	xframe->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");
        // plotting data of the Bmass from file on the frame
        redData->plotOn(xframe,RooFit::Name("data"));

        // plot the model defined which is the addition of two crystal ball for signal and exponential for background on the frame
        model.plotOn(xframe, RooFit::Name("Full PDF"), LineColor(kBlue), MarkerStyle(20));

        // create a pull distribution for the fit 
        RooHist* hpull = xframe->pullHist();

	printf("**************************************");
        printf("\n sigma1: %.5f", sigma1.getVal());
        printf("\n sigma2: %.5f", sigma2.getVal());
        double eff_sigma = sqrt(sigM_frac.getVal()*pow(sigma1.getVal(),2) + (1 - sigM_frac.getVal()) * pow(sigma2.getVal(),2));
        printf("\n eff_sigma: %.5f", eff_sigma);
        printf("\n************************************\n");

        double chi2dof = xframe->chiSquare();
        std::cout<<"\n"<<std::endl;
        std::cout<<"#chi^{2}/dof= "<< chi2dof << std::endl;

	// plotting individual crystal ball in the frame
        model.plotOn(xframe, RooFit::Name("CB1"), Components("CB1"), LineStyle(kDashed), LineColor(kRed));
	model.plotOn(xframe, RooFit::Name("CB2"), Components("CB2"), LineStyle(kDashed), LineColor(kGreen));

	xframe->Draw();

	TPaveText* paveText = new TPaveText(0.65,0.50,0.83,0.88,"NDC");
        paveText->SetBorderSize(0.0);
        paveText->SetFillColor(kWhite);
        paveText->SetFillStyle(0);
        paveText->SetTextSize(0.02);
        paveText->AddText(Form("N_{sig} = %.0f #pm %.0f", nsig.getVal(), nsig.getError()));
	paveText->AddText(Form("mean  = %.5f #pm %.5f GeV" , mean.getVal() , mean.getError()));
        paveText->AddText(Form("#sigma_{1} = %.5f #pm %.5f GeV", sigma1.getVal(), sigma1.getError()));
        paveText->AddText(Form("#sigma_{2} = %.5f #pm %.5f GeV", sigma2.getVal(), sigma2.getError()));
        paveText->AddText(Form("#sigma_{eff} = %.5f  GeV", eff_sigma));
        paveText->AddText(Form("frac = %.5f #pm %.5f ", sigM_frac.getVal(), sigM_frac.getError()));
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
        mark->DrawLatex(0.15,startY,"#bf{CMS} #it{Simulation}");

        mark->DrawLatex(0.65,startY,"#scale[0.70]{539.3205 fb^{-1} (13 TeV)}");
        mark->Draw();

	p2->cd();

        RooPlot* pull = Bmass.frame(Title(""));
        pull->addPlotable(hpull, "BX");
        pull->SetMinimum(-5);
        pull->SetMaximum(5);
        hpull->SetFillColor(kBlue);
        pull->GetYaxis()->SetNdivisions(5);
        pull->SetTitle("Pull Distribution");
	pull->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");
        pull->GetYaxis()->SetTitle("");
	pull->GetYaxis()->SetLabelSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
        pull->GetXaxis()->SetLabelSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
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

}
