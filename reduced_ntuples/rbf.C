double frac(double x, double x1, double y, double y1)
{
    double a = x;
    double del_a = x1;
    double b = y;
    double del_b = y1;

    double t1 = b*del_a;
    double t2 = a*del_b;

    double n = TMath::Sqrt(t1*t1+t2*t2);
    double d = b*b;

    return n/d; 

}

double prod(double x, double x1, double y, double y1)
{
    double a = x;
    double del_a = x1;
    double b = y;
    double del_b = y1;

    double t1 = b*del_a;
    double t2 = a*del_b;

    double delta = TMath::Sqrt(t1*t1+t2*t2);
    return delta;

}


void rbf()
{
    double r_ae = 0.05645;
    double j_ae = 0.04590;
    double p_ae = 0.05830;

    double r_mc, r_dt, j_mc, j_dt, p_mc, p_dt;
    cout<<"Enter rare mode MC yield-> "<<endl;
    cin>> r_mc;
    cout<<"Enter rare mode Data yield-> "<<endl;
    cin>> r_dt;
    cout<<"Enter JPsi mode MC yield-> "<<endl;
    cin>> j_mc;
    cout<<"Enter JPsi mode Data yield-> "<<endl;
    cin>> j_dt;
    cout<<"Enter Psi mode MC yield-> "<<endl;
    cin>> p_mc;
    cout<<"Enter Psi mode Data yield-> "<<endl;
    cin>> p_dt;

    double r_re = (double)  r_mc/(87912702+88098785);
    double j_re = (double)  j_mc/94945299;
    double p_re = (double)  p_mc/35560367;

    double r_eff = (double) r_re*r_ae;
    double j_eff = (double) j_re*j_ae;
    double p_eff = (double) p_re*p_ae;

    double rbf_rj = (double) (r_dt/j_dt)*(j_eff/r_eff)*0.05961;
    double rbf_rp = (double) (r_dt/p_dt)*(p_eff/r_eff)*(8e-03);
    double rbf_pj = (double) (p_dt/j_dt)*(j_eff/p_eff)*7.45125;

    //caclulation for error in rare to jpsi mode

    double r_dte;
    double j_dte;
    double p_dte;
    cout<<"Enter error in rare mode data yield-> "<<endl;
    cin>> r_dte;
    cout<<"Enter error in Jpsi mode data yield-> "<<endl;
    cin>> j_dte;
    cout<< "Enter error in psi mode data yield-> "<<endl;
    cin>> p_dte;
    double t1 = frac(r_dt, r_dte, j_dt, j_dte);
    //final error in rare mode brancing fraction
    double del_rj = prod(r_dt/j_dt, t1, 0.05961, 0.033e-2);

    double pj1 = frac(r_dt, r_dte, p_dt,p_dte);
    //final error in rare to psi mode brancing fraction
    double del_rp = prod(r_dt/p_dt, pj1, 8e-03, 0.6e-03);

    double t2 = frac(p_dt, p_dte, j_dt, j_dte);

    //final error in reso to reso
    double del_pj = prod(p_dt/j_dt, t2, 7.45125, 0.56065);

    cout<<"********************************************"<<endl;
        

    cout<<"Rare Mode Data Yield-> "<< r_dt <<endl;
    cout<<"Jpsi Mode Data Yield-> "<< j_dt <<endl;
    cout<< "Pso Mode Data Yield-> "<< p_dt <<endl;

    cout<<"********************************************"<<endl;

    cout<<"Rare Mode MC yield-> "<< r_mc <<endl;
    cout<<"Jpsi Mode MC yield-> "<< j_mc <<endl;
    cout<<" Psi Mode MC yield-> "<< p_mc <<endl;

    cout<<"**************************************************"<<endl;

    cout<<"Recontruction efficiencies->"<<endl;
    cout<<"Rare mode-> "<< r_re <<endl;
    cout<<"JPsi mode-> "<< j_re <<endl;
    cout<<"Psi mode-> " << p_re <<endl;

    cout<<"**************************************************"<<endl;

    cout<<"**************************************************"<<endl;

    cout<<"Total efficiencies->"<<endl;
    cout<<"Rare mode-> "<< r_eff <<endl;
    cout<<"JPsi mode-> "<< j_eff <<endl;
    cout<<"Psi mode-> " << p_eff <<endl;

    cout<<"**************************************************"<<endl;
    cout<<"**************************************************"<<endl;
    cout<<"**************************************************"<<endl;
    cout<<"**************************************************"<<endl;

    cout<<"Relative Branching Fractions->"<<endl;
    
    cout<<"Rare to Jpsi mode-> " << rbf_rj<<" +/- "<<del_rj <<endl;
    cout<<"Rare to Psi mode-> " <<  rbf_rp<<" +/- "<<del_rp <<endl;
    cout<<"Rare to Jpsi mode-> " << rbf_pj<<" +/- "<<del_pj <<endl;



}
