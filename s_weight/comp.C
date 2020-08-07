void comp()
{
    TChain *ch = new TChain("tree");
    ch->Add("/home/rishabh/project/root_files/2016/ntuples/reduces_ntuples/JPsi_MC_ForFit.root");
    TTree *sig_tr =ch;

    TChain *ch1 = new TChain("tree");
    ch1->Add("/home/rishabh/project/root_files/2016/ntuples/s_wt/JPsi_DataSet_with_sWeight.root");
    TTree *bkg_tr = ch1;

    

    TH1D *h1 = new TH1D("h1", "h1", 50, -1., 1.);
    TH1D *h2 = new TH1D("h2", "h2", 50, -1., 1.);
    TH1D *t1 = new TH1D("t1", "t1", 50, -1., 1.);
    TH1D *t2 = new TH1D("t2", "t2", 50, -1., 1.);
    

    double k1, k2;
    double d1, d2;

    sig_tr->SetBranchAddress("CosThetaK", &k1);
    bkg_tr->SetBranchAddress("CosThetaK", &k2);
    sig_tr->SetBranchAddress("CosThetaL", &d1);
    bkg_tr->SetBranchAddress("CosThetaL", &d2);



    sig_tr->Project("h1", "CosThetaK");
    bkg_tr->Project("h2", "CosThetaK", "s_wt");
    sig_tr->Project("t1", "CosThetaL");
    bkg_tr->Project("t2", "CosThetaL", "s_wt");
    

    
    double norm = 1;
    double sc1 = norm/h1->Integral();
    h1->Scale(sc1);

    double sc2 = norm/h2->Integral();
    h2->Scale(sc2);

    double sc3 = norm/t1->Integral();
    t1->Scale(sc3);

    double sc4 = norm/t2->Integral();
    t2->Scale(sc4);


    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 900);
    c1->Divide(2,1);

    c1->cd(1);
    TPad *p1   = new TPad("p1","p1", 0.05, 0.25, 0.995, 0.97);
    TPad *p2   = new TPad("p2","p2", 0.05, 0.02, 0.995, 0.24);
    p1->Draw();
    p2->Draw();


    p1->cd();
    h1->SetStats(0);
    h2->SetStats(0);
    h1->SetLineWidth(2);
    h1->SetLineColor(4);
    h1->SetFillColor(kBlue+2);
    h1->SetFillStyle(3001);
    h2->SetTitle("#bf{cos#theta_{K}}");
    h2->GetYaxis()->SetTitle("#bf{A.U.}");
    h2->GetXaxis()->SetTitle("#bf{cos#theta_{K}}");
    h2->GetYaxis()->SetTitleOffset(1.6);
    h2->SetMarkerStyle(8);
    h2->SetMarkerColor(2);
    h2->Draw("e1");
    h1->Draw("histsame");

    auto legend=new TLegend(0.70,0.70,0.85,0.85);
    legend->AddEntry(h1,"MC", "F");
    legend->AddEntry(h2,"Data", "P");
    legend->Draw();

    p2->cd();

    p2->SetBottomMargin(0.2);
    TH1D *h3 = (TH1D*)h1->Clone("h3");
    h3->Sumw2();
    h3->GetYaxis()->SetRangeUser(0.5,1.5);
    h3->SetTitle("");
    h3->GetYaxis()->SetTitle("");
    h3->GetXaxis()->SetTitle("#bf{cos#theta_{K}}");
    h3->GetYaxis()->SetTitleSize( (0.7/0.3)*h1->GetXaxis()->GetLabelSize() );
    h3->GetXaxis()->SetTitleSize( (0.7/0.3)*h1->GetXaxis()->GetLabelSize() );
    h3->GetYaxis()->SetTitleOffset(0.5);
    h3->GetYaxis()->SetLabelSize( (0.7/0.3)*h1->GetXaxis()->GetLabelSize() );
    h3->GetXaxis()->SetLabelSize( (0.7/0.3)*h1->GetXaxis()->GetLabelSize() );
    h3->GetYaxis()->SetNdivisions(505);
    h3->SetStats(0);
    h3->Divide(h2);
    h3->SetMarkerStyle(8);
    h3->SetMarkerColor(1);
    h3->Draw("e1");

    
    

    int n=21;
    double x[n];
    double flag = -1.;
    for(int i=0;i<n;i++)
    {
        x[i] = flag;
        flag = flag + 0.1;
    }
    TGraph *grshade =  new TGraph(2*n);
    for (int i =0;i<n;i++)
	{
		grshade->SetPoint(i, x[i], 0.8);
		grshade->SetPoint(n+i, x[n-i-1], 1.2);
	}

    grshade->SetFillColor(kGreen);
	grshade->SetFillStyle(3008);
	grshade->Draw("f");

    

    c1->cd(2);
    TPad *p_1   = new TPad("p_1","p_1", 0.05, 0.25, 0.995, 0.97);
    TPad *p_2   = new TPad("p_2","p_2", 0.05, 0.02, 0.995, 0.24);
    p_1->Draw();
    p_2->Draw();

    p_1->cd();

    t1->SetStats(0);
    t2->SetStats(0);
    t1->SetLineWidth(2);
    t1->SetLineColor(4);
    t1->SetFillColor(kBlue+2);
    t1->SetFillStyle(3001);
    t2->SetTitle("#bf{cos#theta_{L}}");
    t2->GetYaxis()->SetTitle("#bf{A.U.}");
    t2->GetXaxis()->SetTitle("#bf{cos#theta_{L}}");
    t2->GetYaxis()->SetTitleOffset(1.6);
    t2->SetMarkerStyle(8);
    t2->SetMarkerColor(2);
    t2->Draw("e1");
    t1->Draw("histsame");

    
    auto leg=new TLegend(0.70,0.70,0.85,0.85);
    leg->AddEntry(t1,"MC", "F");
    leg->AddEntry(t2,"Data", "P");
    leg->Draw();

    p_2->cd();

    p_2->SetBottomMargin(0.2);
    TH1D *h33 = (TH1D*)t1->Clone("h33");
    h33->Sumw2();
    h33->GetYaxis()->SetRangeUser(0.5,1.5);
    h33->SetTitle("");
    h33->GetYaxis()->SetTitle("");
    h33->GetXaxis()->SetTitle("#bf{cos#theta_{L}}");
    h33->GetYaxis()->SetTitleSize( (0.7/0.3)*t1->GetXaxis()->GetLabelSize() );
    h33->GetXaxis()->SetTitleSize( (0.7/0.3)*t1->GetXaxis()->GetLabelSize() );
    h33->GetYaxis()->SetTitleOffset(0.5);
    h33->GetYaxis()->SetLabelSize( (0.7/0.3)*t1->GetXaxis()->GetLabelSize() );
    h33->GetXaxis()->SetLabelSize( (0.7/0.3)*t1->GetXaxis()->GetLabelSize() );
    h33->GetYaxis()->SetNdivisions(505);
    h33->SetStats(0);
    h33->Divide(t2);
    h33->SetMarkerStyle(8);
    h33->SetMarkerColor(1);
    h33->Draw("e1");

    int n1 = 21;
    double x1[n];
    double flag1 = -1.;
    for(int i=0;i<n1;i++)
    {
        x[i] = flag1;
        flag1 = flag1 + 0.1;
    }
    TGraph *grshade1 =  new TGraph(2*n1);
    for (int i =0;i<n1;i++)
	{
		grshade1->SetPoint(i, x[i], 0.8);
		grshade1->SetPoint(n1+i, x[n1-i-1], 1.2);
	}

    grshade1->SetFillColor(kGreen);
	grshade1->SetFillStyle(3008);
	grshade1->Draw("f");





}
