double gepkelly(double tau) {
    //G_E^p from Kelly fit
    //double gepk = (1.0-0.24*tau)/(1+(10.98*tau)+(12.82*tau**2)+(21.97*tau**3));
    double fracnum = 1.0-0.24*tau;
    double fracden = 1+tau*(10.98+tau*(12.82+tau*21.97));
    return fracnum/fracden;
}
double gmpkelly(double tau) {
    //G_M^p from Kelly fit
    //double gmpk = xmup * (1. + 0.12*tau)/(1.+(10.97*tau)+(18.86*tau**2)+(6.55*tau**3));
    double fracnum = (5.586/2.)*(1.0+0.12*tau);
    double fracden = 1+tau*(10.97+tau*(18.86+tau*6.55));
    return fracnum/fracden;
}
double geparrington2004r(double tau) {
    //G_E^p from Arrington Rosenbluth fit of PHYSICAL REVIEW C 69, 022201(R) (2004)
    double xmp = 0.938272;
    double qsq = tau * (4*TMath::Power(xmp,2));
    double gea = 1. + qsq * (3.226 + qsq * (1.508 + qsq * (-0.3773 + qsq * (0.611 + qsq * (-0.1853 + qsq * 0.01596)))));
    return gea;
}
double gmparrington2004r(double tau) {
    //G_M^p from Arrington Rosenbluth fit of PHYSICAL REVIEW C 69, 022201(R) (2004)
    double xmp = 0.938272;
    double qsq = tau * (4*TMath::Power(xmp,2));
    double gma = 1. + qsq * (3.19 + qsq * (1.355 + qsq * (0.151 + qsq * (-0.0114 + qsq * (0.000533 + qsq * -0.000009)))));
    gma = (5.586/2.) * gma;
    return gma;
}
double geparrington2004p(double tau) {
    //G_E^p from Arrington polarization fit of PHYSICAL REVIEW C 69, 022201(R) (2004)
    double xmp = 0.938272;
    double qsq = tau * (4*TMath::Power(xmp,2));
    double gea = 1. + qsq * (2.94 + qsq * (3.04 + qsq * (-2.255 + qsq * (2.002 + qsq * (-0.5338 + qsq * 0.04875)))));
    return gea;
}
double gmparrington2004p(double tau) {
    //G_M^p from Arrington polarization fit of PHYSICAL REVIEW C 69, 022201(R) (2004)
    double xmp = 0.938272;
    double qsq = tau * (4*TMath::Power(xmp,2));
    double gma = 1. + qsq * (3.00 + qsq * (1.39 + qsq * (0.122 + qsq * (-0.00834 + qsq * (0.000425 + qsq * -0.00000779)))));
    gma = (5.586/2.) * gma;
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
    fracnum = (5.586/2.) * fracnum;
    double fracden = 1.0 + tau * ( 9.627 + tau * (0.000 + tau * (0.000 + tau * (11.179 + tau * 13.245))));
    return fracnum/fracden;
}
double gepbosted(double tau) {
    //G_E^p from Bosted fit of PHYSICAL REVIEW C 51, 409 (1995)
    double xmp = 0.938272;
    double qsq = tau * (4*TMath::Power(xmp,2));
    double q = TMath::Power(qsq,0.5);
    double geb = 1. / (1. + q * (0.62 + q * (0.68 + q * (2.80 + q * 0.83))));
    return geb;
}
double gmpbosted(double tau) {
    //G_M^p from Bosted fit of PHYSICAL REVIEW C 51, 409 (1995)
    double xmp = 0.938272;
    double qsq = tau * (4*TMath::Power(xmp,2));
    double q = TMath::Power(qsq,0.5);
    double gmb = (5.586/2.) / (1. + q * (0.35 + q * (2.44 + q * (0.50 + q * (1.04 + q * 0.34)))));
    return gmb;
}


double poutofeth(double pin, double xm, double xmt, double theta) {
  //calculate outgoing lepton momentum for lepton + proton elastic scattering
  //pin: input momentum;     xm: e or mu mass;     theta: angle in radians
  //double xamu=0.931494;
  //double xmt     = 12.*xamu;
  //double xmp     = 0.938272;
  double ein     = TMath::Power(TMath::Power(pin,2)+TMath::Power(xm,2),0.5);
  double xmssq   = TMath::Power(ein+xmt,2) - TMath::Power(pin,2) + TMath::Power(xm,2) - TMath::Power(xmt,2);
  double a       = 4. * (TMath::Power(ein+xmt,2) - TMath::Power(pin*cos(theta),2));
  double b       = -4. * xmssq * pin * cos(theta);
  double c       = 4. * TMath::Power(ein+xmt,2) * TMath::Power(xm,2) - TMath::Power(xmssq,2);
  double bsqm4ac = TMath::Power(b,2) - 4.*a*c;
  double pout    = (-1.*b + TMath::Power(bsqm4ac,0.5)) / (2.*a);
  return pout;
}

void muerat2026()
{
    // look at mu p / ep cross section ratio
    gROOT->Reset();
    gSystem->Load("libPhysics"); 
 
    //constants 
    double degr=0.01745329252, pi=3.1415926536, eeee=2.71828, xhbarc=0.197327, alfa=1./137.036, xmup  = 5.586/2.;
    const int n = 88;
    double dang = 1., angstep = 0.5*degr;
    double sig0[n], srat0[n], srat1[n], srat2[n], thplot[n];
  
    //masses
    double xmp=0.938272, xmn=0.939566, xmgam=0.0, xmpi0=0.13495, xme=0.000511;
    double xamu=0.931494, xmmu=0.105658366;
    double xmc=12.*xamu, zee=6.;
    double xm;
    //beam momentum, set up for 115 
    double pin = 0.115;
    
    //choice of formfactors
    bool kelly = false;
    bool arrington2004r = false;
    bool arrington2004p = false;
    bool amt = false;
    bool bosted = true;

    // define some histograms, in 2-d fill, it is Fill(horiz,vert)


    //loop over momenta
    for (int imom=0; imom<3; imom++) {
        if (imom==1) {pin = 0.160;}
        if (imom==2) {pin = 0.210;}
        //loop over particles
        for (int ipart=0; ipart<2; ipart++) {
            if (ipart==0) {xm = xmmu; printf("muons\n");}
            if (ipart==1) {xm = xme; printf("electrons\n");}
            double ein     = TMath::Power(TMath::Power(pin,2)+TMath::Power(xm,2),0.5);
            double betain   = pin/ein;
            //TLorentzVector lvftmuin;  lvftmuin.SetPxPyPzE(0.,0.,pin,ein);    
            // loop over angles (!) from 15 - 102 degrees in 1 degree bins
            for (int i=0; i<n; i++) {
                //determine angle and solid angle
                double theta = (15 + dang*i) * degr;
                //some H kinematics - our qsq = Q^2, a positive quantity, not q^2 = -Q^2, a negative quantity
                double pmuo   = poutofeth(pin, xm, xmp, theta);
                double emuo   = TMath::Power(TMath::Power(pmuo,2)+TMath::Power(xm,2),0.5);
                double qsq    = 2.*TMath::Power(xm,2) - 2.*ein*emuo + 2.*pin*pmuo*cos(theta); qsq = -1.*qsq;
                double sinsq  = TMath::Power(TMath::Sin(0.5*theta),2);
                double cossqc = 1. - TMath::Power(betain,2)*sinsq;
                // form factors
                double tau    = qsq/(4*TMath::Power(xmp,2));
                double gepk, gmpk;
                if (kelly) {
                    gepk   = gepkelly(tau);
                    gmpk   = gmpkelly(tau);
                }
                else if (arrington2004r) {
                    gepk   = geparrington2004r(tau);
                    gmpk   = gmparrington2004r(tau);
                }
                else if (arrington2004p) {
                    gepk   = geparrington2004p(tau);
                    gmpk   = gmparrington2004p(tau);
                }
                else if (amt) {
                    gepk   = gepamt(tau);
                    gmpk   = gmpamt(tau);
                }
                else if (bosted) {
                    gepk   = gepbosted(tau);
                    gmpk   = gmpbosted(tau);
                }
                //cross section
                //sm (sm_mu) is the lepton m=0 (m!=0) pointlike cross section
                //use massive lepton formula from Preedom and Tegen, PRC36
                //eta from that article is tau here, eta here is the P&T -q^2/4EE'
                //add in hbarc to convert cross section to fm^2/sr
                double eta    = qsq/(4.*ein*emuo);
                double d      = TMath::Power((1.-TMath::Power(xm/ein,2))/(1.-TMath::Power(xm/emuo,2)),0.5);
                //double sm     = TMath::Power(xhbarc*alfa/(2*ein*sinsq),2) * (emuo*cossqc/ein);
                double sm_mu_den = 1. + 2.*ein*d*sinsq/xmp + ein*(1.-d)/xmp;
                double sm_mu  = TMath::Power(xhbarc*alfa/(2*ein),2) * ((1.-eta)/TMath::Power(eta,2)) / (d*sm_mu_den);
                double a      = TMath::Power(gepk,2)/ (1.+tau) + TMath::Power(gmpk,2)*tau / (1.+tau);
                //double b      = 2.*tau * TMath::Power(gmpk,2) * TMath::Power(TMath::Tan(0.5*theta),2);
                //double sig    = sm * (a + b);
                double b_mu    = (2.*tau-TMath::Power(xm/xmp,2)) * TMath::Power(gmpk,2) * (eta/(1.-eta));
                double sig_mu  = sm_mu * (a + b_mu);
                if (ipart==0) {sig0[i] = sig_mu;}
                if (ipart==1) {
                    if (imom==0) {srat0[i] = sig0[i]/sig_mu;  thplot[i] = theta/degr; 
                    //printf("th = %f, rat = %f\n",thplot[i],srat[i]);
                    }
                    if (imom==1) {srat1[i] = sig0[i]/sig_mu;}
                    if (imom==2) {srat2[i] = sig0[i]/sig_mu;}
                }
            } // angle
        } // particle
    } // momentum

    auto c1 = new TCanvas("muesigratio","muesigratio"); //,200,10,wwide,wtall);
    //c1->SetCanvasSize(1200, 400);
    //c1->SetWindowSize(1210, 450);
    //TLegend *legend = new TLegend(0.70, 0.65, 0.9, 0.90);

    TMultiGraph *multigraph = new TMultiGraph();

    auto gr = new TGraph(n,thplot,srat0);
    gr->SetLineColor(2);
    gr->SetLineWidth(6);
    //legend->AddEntry(multigraph, 115.0 , "l");
    multigraph->Add(gr,"");


    auto gr1 = new TGraph(n,thplot,srat1);
    gr1->SetLineColor(3);
    gr1->SetLineWidth(4);
    //legend->AddEntry(multigraph, 160.0 , "l");
    multigraph->Add(gr1,"");

    auto gr2 = new TGraph(n,thplot,srat2);
    gr2->SetLineColor(4);
    gr2->SetLineWidth(2);
    //legend->AddEntry(multigraph, 210.0 , "l");
    multigraph->Add(gr2,"");
    
    if (kelly) {
        multigraph->SetTitle("#sigma_{#mu p} / #sigma_{e p} at 115, 160, 210 MeV/c, Kelly form factors");
    }
    else if (arrington2004r) {
        multigraph->SetTitle("#sigma_{#mu p} / #sigma_{e p} at 115, 160, 210 MeV/c, Arrington 2004r form factors");
    }
    else if (arrington2004p) {
        multigraph->SetTitle("#sigma_{#mu p} / #sigma_{e p} at 115, 160, 210 MeV/c, Arrington 2004p form factors");
    }
    else if (amt) {
        multigraph->SetTitle("#sigma_{#mu p} / #sigma_{e p} at 115, 160, 210 MeV/c, AMT form factors");
    }
    else if (bosted) {
        multigraph->SetTitle("#sigma_{#mu p} / #sigma_{e p} at 115, 160, 210 MeV/c, Bosted form factors");
    }
    multigraph->GetXaxis()->SetTitle("#theta (deg)");
    multigraph->GetYaxis()->SetTitle("#sigma_{#mu p} / #sigma_{e p}");
    multigraph->GetXaxis()->CenterTitle();
    multigraph->GetYaxis()->CenterTitle();
    //multigraph->GetXaxis()->SetRangeUser(0., 120.);
    //multigraph->GetYaxis()->SetRangeUser(1., 3.);
    multigraph->GetYaxis()->SetTitleOffset(1.2);
    multigraph->Draw("AL");
    c1->Modified();
    //legend->Draw();
}
