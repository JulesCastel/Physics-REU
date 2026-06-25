double gdipole(double qsq) {
    //Gdipole parameterization
    double dipole = TMath::Power((1.+qsq/0.71),-2.);
    return dipole;
}
double gepkelly(double tau) {
    //G_E^p from Kelly fit
    double fracnum = 1.0-0.24*tau;
    double fracden = 1+tau*(10.98+tau*(12.82+tau*21.97));
    return fracnum/fracden;
}
double gmpkelly(double tau) {
    //G_M^p from Kelly fit
    double fracnum = 2.793 * (1.0+0.12*tau);
    double fracden = 1+tau*(10.97+tau*(18.86+tau*6.55));
    return fracnum/fracden;
}
double geparrington2004r(double qsq) {
    //G_E^p from Arrington Rosenbluth fit of PHYSICAL REVIEW C 69, 022201(R) (2004)
    double gea = 1. + qsq * (3.226 + qsq * (1.508 + qsq * (-0.3773 + qsq * (0.611 + qsq * (-0.1853 + qsq * 0.01596)))));
    return 1./gea;
}
double gmparrington2004r(double qsq) {
    //G_M^p from Arrington Rosenbluth fit of PHYSICAL REVIEW C 69, 022201(R) (2004)
    double gma = 1. + qsq * (3.19 + qsq * (1.355 + qsq * (0.151 + qsq * (-0.0114 + qsq * (0.000533 + qsq * -0.000009)))));
    gma = 2.793 / gma;
    return gma;
}
double geparrington2004p(double qsq) {
    //G_E^p from Arrington polarization fit of PHYSICAL REVIEW C 69, 022201(R) (2004)
    double gea = 1. + qsq * (2.94 + qsq * (3.04 + qsq * (-2.255 + qsq * (2.002 + qsq * (-0.5338 + qsq * 0.04875)))));
    return 1./gea;
}
double gmparrington2004p(double qsq) {
    //G_M^p from Arrington polarization fit of PHYSICAL REVIEW C 69, 022201(R) (2004)
    double gma = 1. + qsq * (3.00 + qsq * (1.39 + qsq * (0.122 + qsq * (-0.00834 + qsq * (0.000425 + qsq * -0.00000779)))));
    gma = 2.793 / gma;
    return gma;
}
double gepamt(double tau) {
    //G_E^p from AMT fit PHYSICAL REVIEW C 76, 035205 (2007)
    double fracnum = 1.0 + tau * (3.439 + tau * (-1.602 + tau * 0.068));
    double fracden = 1.0 + tau * (15.055 + tau * (48.061 + tau * (99.304 + tau * (0.012 + tau * 8.650))));
    return fracnum/fracden;
}
double gmpamt(double tau) {
    //G_M^p from AMT fit PHYSICAL REVIEW C 76, 035205 (2007)
    double fracnum = 1.0 + tau * (-1.465 + tau * (1.260 + tau * 0.262));
    fracnum = 2.793 * fracnum;
    double fracden = 1.0 + tau * ( 9.627 + tau * (0.000 + tau * (0.000 + tau * (11.179 + tau * 13.245))));
    return fracnum/fracden;
}
double gepbosted(double qsq) {
    //G_E^p from Bosted fit of PHYSICAL REVIEW C 51, 409 (1995)
    double q = TMath::Power(qsq,0.5);
    double geb = 1. / (1. + q * (0.62 + q * (0.68 + q * (2.80 + q * 0.83))));
    return geb;
}
double gmpbosted(double qsq) {
    //G_M^p from Bosted fit of PHYSICAL REVIEW C 51, 409 (1995)
    double q = TMath::Power(qsq,0.5);
    double gmb = 2.793 / (1. + q * (0.35 + q * (2.44 + q * (0.50 + q * (1.04 + q * 0.34)))));
    return gmb;
}

void gratios()
{
    // look at mu p / ep cross section ratio
    gROOT->Reset();
    gSystem->Load("libPhysics"); 
 
    //constants 
    double degr=0.01745329252, pi=3.1415926536, eeee=2.71828, xhbarc=0.197327, alfa=1./137.036, xmup  = 5.586/2.;
    const int n = 101;
    double xplot[n];
  
    //masses
    double xmp=0.938272, xmn=0.939566, xmgam=0.0, xmpi0=0.13495, xme=0.000511;
    double xamu=0.931494, xmmu=0.105658366;
    double xmc=12.*xamu, zee=6.;
    double xm;
    //beam momentum 
    int ipin = 210;
    double pin = 0.01*ipin;
    
    //choice of formfactors
    double gd[n], gk[n], gap[n], gar[n], gamt[n], gb[n];

    //loop over qsq, evaluate each f.f., calculate ratios to dipole
    
    //loop over qsq, evaluate all f.f.
    for (int iq=0; iq<n; iq++) {
        double qsq = 0.01 * iq;
        double tau    = qsq/(4*TMath::Power(xmp,2));
        xplot[iq] = qsq;
        gd[iq]    = gdipole(qsq);
        gk[iq]    = gepkelly(tau) / gd[iq]; //gmpk   = gmpkelly(tau); etc.
        gap[iq]   = geparrington2004p(qsq) / gd[iq]; 
        gar[iq]   = geparrington2004r(qsq) / gd[iq]; 
        gamt[iq]  = gepamt(tau) / gd[iq]; 
        gb[iq]    = gepbosted(qsq) / gd[iq]; 
    } // qsq loop
    
    auto c1 = new TCanvas("GE_wrt_dipole","GE_wrt_dipole"); //,200,10,wwide,wtall);
    //c1->SetCanvasSize(1200, 400);
    //c1->SetWindowSize(1210, 450);
    //TLegend *legend = new TLegend(0.70, 0.65, 0.9, 0.90);

    TMultiGraph *multigraph = new TMultiGraph();
    
    TLegend *legend = new TLegend(0.67,0.65,0.9,0.9); //Create the TLegend object and define it's position 
    legend->SetHeader("Form Factor Param.", "C"); //"C" Center alignment for the header ("L" Left and "R" Right)
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(1);
    legend->SetTextSize(0.04);
    //legend->AddEntry(gr1, "Data points", "lp");  // "p" for point marker, "l" for line, "e" for error bars if TGraphError is used.

    auto gr = new TGraph(n,xplot,gk);
    gr->SetLineColor(2);
    gr->SetLineWidth(1);
    legend->AddEntry(gr, "Kelly" , "l");
    multigraph->Add(gr,"");

    auto gr1 = new TGraph(n,xplot,gap);
    gr1->SetLineColor(3);
    gr1->SetLineWidth(2);
    legend->AddEntry(gr1, "Arr_pol" , "l");
    multigraph->Add(gr1,"");

    auto gr2 = new TGraph(n,xplot,gar);
    gr2->SetLineColor(4);
    gr2->SetLineWidth(3);
    legend->AddEntry(gr2,  "Arr_Ros" , "l");
    multigraph->Add(gr2,"");

    auto gr3 = new TGraph(n,xplot,gamt);
    gr3->SetLineColor(6);
    gr3->SetLineWidth(4);
    legend->AddEntry(gr3,  "AMT" , "l");
    multigraph->Add(gr3,"");
    
    auto gr4 = new TGraph(n,xplot,gb);
    gr4->SetLineColor(8);
    gr4->SetLineWidth(5);
    legend->AddEntry(gr4,  "Bosted" , "l");
    multigraph->Add(gr4,"");
    
    multigraph->GetXaxis()->SetTitle("Q^2 (GeV^2)");
    multigraph->GetYaxis()->SetTitle("G_E / G_{dipole}");
    multigraph->GetXaxis()->CenterTitle();
    multigraph->GetYaxis()->CenterTitle();
    //multigraph->GetXaxis()->SetRangeUser(0., 120.);
    //multigraph->GetYaxis()->SetRangeUser(1., 3.);
    multigraph->GetYaxis()->SetTitleOffset(1.2);
    multigraph->Draw("AL");
    legend->Draw();
    //c1->Modified();
}
