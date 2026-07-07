
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
//#include "TArrow.h"
//#include "TLatex.h"

void BHStabForReport2024deltas()
{
    double runnum[2500], zeros[2500], te[2500], teunc[2500], tmu[2500], tmuunc[2500], tpi[2500], tpiunc[2500];
    double prunnum[2500],      pte[2500], pteunc[2500], ptmu[2500], ptmuunc[2500], ptpi[2500], ptpiunc[2500];
    double fe[2500], feunc[2500], fmu[2500], fmuunc[2500], fpi[2500], fpiunc[2500];
	int j = 0;
    //int p=160;
    int p=115;
    //int bar = 16;
    //std::string plane = "BHC";
    //int bar = 6;
    std::string plane = "BHD";	///// load files
    //int i1 = 15681, i2 = 15715;
    //int i1 = 15412, i2 = 15655;
    //int i1 = 15324, i2 = 15500;
    //int i1 = 15650, i2 = 15717;
    //int i1 = 15324, i2 = 15717;
    //+210 for report
    //int i1 = 15688, i2 = 15738;
    //all +210
    //int i1 = 15324, i2 = 15821;
    //all -210
    //int i1 = 15490, i2 = 15600;
    //all 2021
    //int i1 = 8940, i2 = 11334;
    //int i1 = 11200, i2 = 11250;
    int i1 = 35566, i2 = 35572;
    for (int irun=i1; irun<i2; irun++) {
        if(
            //remove calo calibration
            irun != 26502 && irun != 26503 && irun != 26504 && irun != 26505 && irun != 26506 &&
            irun != 26507 && irun != 26508 && irun != 26509 && irun != 26510 && irun != 26511 &&
            irun != 26512 && irun != 26513 && irun != 26514
            //remove Collimator scan
            &&
            irun != 26555 && irun != 26556 && irun != 26557 && irun != 26558 && irun != 26559 && irun != 26560
            //remove TOF scan
            &&
            irun != 26609 && irun != 26610 && irun != 26611 && irun != 26612 && irun != 26613 &&
            irun != 26614 && irun != 26615 && irun != 26616 && irun != 26617 && irun != 26618
            //remove some positron trigger runs that mess up muon-electron RF time
            &&
            irun != 26462 && irun != 26463 && irun != 26464 && irun != 26465 && irun != 26466 &&
            irun != 26467 && irun != 26468 && irun != 26469 && irun != 26470 && irun != 26471 &&
            irun != 26472

            //
            //(irun<26938 || irun > 26873) &&

              /*&& irun != 14955 && irun != 14956 && irun != 14957
               && irun != 15154 && irun != 15156 && irun != 15157 && irun != 15158 */
            )
        {
        //printf("j = %d  irun = %d \n",j,irun);
		TFile *file = new TFile(Form("run%d.root",irun));

		if ( file != NULL )  {

            ///// load histograms
            TH1F* h_rf = (TH1F*)file->Get(Form("%s/RF/RF Avg of Plane",plane.c_str()));
            //TH1F* h_rf = (TH1F*)file->Get(Form("%s/RF/Out of Time RF Avg of Plane",plane.c_str()));
            if ( h_rf != NULL ){
                runnum[j] = irun;
                zeros[j] = 0.;

                ///// fit histogram to determine rf peak positions
                TF1* f_gaus = new TF1(Form("f_%s",plane.c_str()),"gaus",6.7,8.7);
                h_rf->Fit(f_gaus,"R");
                //210
                //TF1* f_gaus_mu = new TF1(Form("f_%s_mu",plane.c_str()),"gaus",15.8,17.6);
                //TF1* f_gaus_pi = new TF1(Form("f_%s_pi",plane.c_str()),"gaus",2,4.5);
                //160
                //TF1* f_gaus_mu = new TF1(Form("f_%s_mu",plane.c_str()),"gaus",2.0,4.5);
                //TF1* f_gaus_pi = new TF1(Form("f_%s_pi",plane.c_str()),"gaus",12.0,14.5);
                //115
                TF1* f_gaus_mu = new TF1(Form("f_%s_mu",plane.c_str()),"gaus",13.8,15.8);
                TF1* f_gaus_pi = new TF1(Form("f_%s_pi",plane.c_str()),"gaus",10.5,12.5);
                h_rf->Fit(f_gaus_mu,"R+");
                h_rf->Fit(f_gaus_pi,"R+");

                ///// print fited info
                printf("Electron RF mean for run %d is %f\n", irun, f_gaus->GetParameter(1));
                //histmean[j] = h_all->GetMean();
                //meanerr[j] = h_all->GetMeanError();

                fe[j] = f_gaus->GetParameter(0);
                fmu[j] = f_gaus_mu->GetParameter(0);
                fpi[j] = f_gaus_pi->GetParameter(0);
                feunc[j] = f_gaus->GetParError(0);
                fmuunc[j] = f_gaus_mu->GetParError(0);
                fpiunc[j] = f_gaus_pi->GetParError(0);
                double tot = fe[j] + fmu[j] + fpi[j];
                double dtot = TMath::Power(TMath::Power(feunc[j],2)+TMath::Power(fmuunc[j],2)+TMath::Power(fpiunc[j],2),0.5);
                dtot = dtot/tot;
                feunc[j] = feunc[j]/fe[j];
                fmuunc[j] = fmuunc[j]/fmu[j];
                fpiunc[j] = fpiunc[j]/fpi[j];
                fe[j] = fe[j] / tot;
                fmu[j] = fmu[j] / tot;
                fpi[j] = fpi[j] / tot;
                feunc[j] = fe[j]*TMath::Power(TMath::Power(feunc[j],2)+TMath::Power(dtot,2),0.5);
                fmuunc[j] = fmu[j]*TMath::Power(TMath::Power(fmuunc[j],2)+TMath::Power(dtot,2),0.5);
                fpiunc[j] = fpi[j]*TMath::Power(TMath::Power(fpiunc[j],2)+TMath::Power(dtot,2),0.5);

                te[j] = f_gaus->GetParameter(1);
                teunc[j] = f_gaus->GetParError(1);
                tmu[j] = f_gaus_mu->GetParameter(1) - te[j];
                tmuunc[j] = f_gaus_mu->GetParError(1);
                tmuunc[j] = TMath::Power( TMath::Power(teunc[j],2) + TMath::Power(tmuunc[j],2), 0.5);
                tpi[j] = f_gaus_pi->GetParameter(1) - te[j];
                tpiunc[j] = f_gaus_pi->GetParError(1);
                tpiunc[j] = TMath::Power( TMath::Power(teunc[j],2) + TMath::Power(tpiunc[j],2), 0.5);

                //	cout << "irun: " << runnum[j]
                //		<< " Mean = " << histmean[j] << " +- " << meanerr[j]
                //		<< "\n";
                //remove collimator scans...
                //if (tmuunc[j] > 0.05) {j--;}
                //if ( (irun > 15434 && irun < 15446) || (irun > 15555 && irun < 15566) || (tmuunc[j] > 1.0)) {
                j++;
			}
            else {
                printf("Hist not found for run %d \n", irun);
            }
		}}
	}

    const int n = j;
    int nprune = 0;
    printf(" total number of points = %d \n",n);
    for (int il=0; il<n; il++){
        te[il] = 1000*te[il];
        tmu[il] = 1000*tmu[il];
        tpi[il] = 1000*tpi[il];
        teunc[il] = 1000*teunc[il];
        tmuunc[il] = 1000*tmuunc[il];
        tpiunc[il] = 1000*tpiunc[il];
        if (teunc[il]<30. && tmuunc[il]<30. && tpiunc[il]<30.){
            prunnum[nprune] = runnum[il];
            pte[nprune] = te[il];
            ptmu[nprune] = tmu[il];
            ptpi[nprune] = tpi[il];
            pteunc[nprune] = teunc[il];
            ptmuunc[nprune] = tmuunc[il];
            ptpiunc[nprune] = tpiunc[il];
            nprune++;
        }
    }

	///// create canvas and divide the canvas into two parts
    //TCanvas *c = new TCanvas(Form("c_%s%02d_RF",plane.c_str(),bar),Form("%s%02d RF",plane.c_str(),bar), 1200, 600);
    //TCanvas *c = new TCanvas("RFTimeDiffs","RFTimeDiffs", 1200, 300);
    //TCanvas *c = new TCanvas(Form("Plane_%s_InTimeStability_-210",plane.c_str()),Form("Plane_%s_InTimeStability_-210",plane.c_str()), 1200, 800);
    TCanvas *c0 = new TCanvas(Form("Plane_%s_OOT_e_RF_Stability_%04d",plane.c_str(),p),Form("Plane_%s_OOT_e_RF_Stability_%04d",plane.c_str(),p), 550, 400);
    //c->Divide(3,2);
    //c->Divide(3,1);

	//c->cd(1);
    TGraphErrors *gr0 = new TGraphErrors(nprune,prunnum,pte,zeros,pteunc);
    gr0->SetMarkerColor(kBlue);
    //gr1->SetMarkerStyle(21);
    gr0->SetMarkerStyle(kOpenCircle);
    gr0->SetTitle(Form("t_{RF e} %s;Run Number;Time (ps)",plane.c_str()));
    gr0->DrawClone("APE");
    // Define a linear function
    //TF1 f("Linear law","[0]+x*[1]",15681,15715);
    //TF1 f("Linear fit","[0]+x*[1]",i1,i2);
    TF1 f("Linear fit","[0]",i1,i2);
    // Let's make the function line nicer
    f.SetLineColor(kRed); f.SetLineStyle(2);
    // Fit the graph and draw it
    gr0->Fit(&f);
    f.DrawClone("Same");
    double p0 = f.GetParameter(0);
    //double p1 = f.GetParameter(1);
    double e0 = f.GetParError(0);
    //double e1 = f.GetParError(1);
    double chi2 = f.GetChisquare();
    //printf("mu const = %f +/- %f   linear term = %f +/- %f --> %f ps/100 runs\n",p0,e0,p1,e1,100*p1);
    printf("e RF const = %f +/- %f \n",p0,e0);
    double avetimee = p0;
    double stddev = 0.;
    for (int ii = 0; ii<n; ii++) {
        //stddev += TMath::Power(p0+p1*runnum[ii]-tmu[ii],2.);
        stddev += TMath::Power(p0-te[ii],2.);
    }
    stddev = stddev/(n-1);
    stddev = TMath::Power(stddev,0.5);
    printf("e RF chi^2 = %f for %d points, sigma = %f \n",chi2,n,stddev);

    TCanvas *c1 = new TCanvas(Form("Plane_%s_OOT_mu-e_RF_Stability_%04d",plane.c_str(),p),Form("Plane_%s_OOT_mu-e_RF_Stability_%04d",plane.c_str(),p), 550, 400);
    //c->cd(2);
	TGraphErrors *gr1 = new TGraphErrors(nprune,prunnum,ptmu,zeros,ptmuunc);
 	gr1->SetMarkerColor(kBlue);
 	//gr1->SetMarkerStyle(21);
 	gr1->SetMarkerStyle(kOpenCircle);
 	gr1->SetTitle(Form("t_{RF #mu} - t_{RF e} %s;Run Number;Time (ps)",plane.c_str()));
	gr1->DrawClone("APE");
    // Fit the graph and draw it
    gr1->Fit(&f);
    f.DrawClone("Same");
    p0 = f.GetParameter(0);
    //double p1 = f.GetParameter(1);
    e0 = f.GetParError(0);
    //double e1 = f.GetParError(1);
    chi2 = f.GetChisquare();
    //printf("mu const = %f +/- %f   linear term = %f +/- %f --> %f ps/100 runs\n",p0,e0,p1,e1,100*p1);
    printf("mu RF const = %f +/- %f \n",p0,e0);
    double avetimemu = p0;
    stddev = 0.;
    for (int ii = 0; ii<n; ii++) {
        //stddev += TMath::Power(p0+p1*runnum[ii]-tmu[ii],2.);
        stddev += TMath::Power(p0-tmu[ii],2.);
    }
    stddev = stddev/(n-1);
    stddev = TMath::Power(stddev,0.5);
    printf("mu RF chi^2 = %f for %d points, sigma = %f \n",chi2,n,stddev);

    TCanvas *c2 = new TCanvas(Form("Plane_%s_OOT_pi-e_RF_Stability_%04d",plane.c_str(),p),Form("Plane_%s_OOT_pi-e_RF_Stability_%04d",plane.c_str(),p), 550, 400);
    //c->cd(3);
    TGraphErrors *gr2 = new TGraphErrors(nprune,prunnum,ptpi,zeros,ptpiunc);
    gr2->SetMarkerColor(kBlue);
    //gr1->SetMarkerStyle(21);
    gr2->SetMarkerStyle(kOpenCircle);
    gr2->SetTitle(Form("t_{RF #pi} - t_{RF e} %s;Run Number;Time (ps)",plane.c_str()));
    gr2->DrawClone("APE");
    //Fit the graph and draw it
    gr2->Fit(&f);
    f.DrawClone("Same");
    p0 = f.GetParameter(0);
    //p1 = f.GetParameter(1);
    e0 = f.GetParError(0);
    //e1 = f.GetParError(1);
    chi2 = f.GetChisquare();
    //printf("pi const = %f +/- %f   linear term = %f +/- %f --> %f ps/100 runs\n",p0,e0,p1,e1,100*p1);
    printf("pi RF const = %f +/- %f \n",p0,e0);
    double avetimepi = p0;
    stddev = 0.;
    for (int ii = 0; ii<n; ii++) {
        //stddev += TMath::Power(p0+p1*runnum[ii]-tpi[ii],2.);
        stddev += TMath::Power(p0-tpi[ii],2.);
    }
    stddev = stddev/(n-1);
    stddev = TMath::Power(stddev,0.5);
    printf("pi RF chi^2 = %f for %d points, sigma = %f \n",chi2,n,stddev);

    for (int il=0; il<nprune; il++){
        pte[il] = pte[il]-avetimee;
        ptmu[il] = ptmu[il]-avetimemu;
        ptpi[il] = ptpi[il]-avetimepi;
    }

    TCanvas *etimeave = new TCanvas(Form("Plane_%s_OOT_e_RF_Stability_delta_%04d",plane.c_str(),p),Form("Plane_%s_OOT_e_RF_Stability_delta_%04d",plane.c_str(),p), 550, 400);
    TGraphErrors *gr0ave = new TGraphErrors(nprune,prunnum,pte,zeros,pteunc);
    gr0ave->SetMarkerColor(kBlue);
    gr0ave->SetMarkerStyle(kOpenCircle);
    gr0ave->SetTitle(Form("t_{RF e} %s;Run Number;Time (ps)",plane.c_str()));
    gr0ave->DrawClone("APE");

    TCanvas *mutimeave = new TCanvas(Form("Plane_%s_OOT_mu-e_RF_Stability_delta_%04d",plane.c_str(),p),Form("Plane_%s_OOT_mu-e_RF_Stability_delta_%04d",plane.c_str(),p), 550, 400);
    TGraphErrors *gr1ave = new TGraphErrors(nprune,prunnum,ptmu,zeros,ptmuunc);
    gr1ave->SetMarkerColor(kBlue);
    gr1ave->SetMarkerStyle(kOpenCircle);
    gr1ave->SetTitle(Form("t_{RF #mu} - t_{RF e} %s;Run Number;Time (ps)",plane.c_str()));
    gr1ave->DrawClone("APE");

    TCanvas *pitimaeave = new TCanvas(Form("Plane_%s_OOT_pi-e_RF_Stability_delta_%04d",plane.c_str(),p),Form("Plane_%s_OOT_pi-e_RF_Stability_delta_%04d",plane.c_str(),p), 550, 400);
    TGraphErrors *gr2ave = new TGraphErrors(nprune,prunnum,ptpi,zeros,ptpiunc);
    gr2ave->SetMarkerColor(kBlue);
    gr2ave->SetMarkerStyle(kOpenCircle);
    gr2ave->SetTitle(Form("t_{RF #pi} - t_{RF e} %s;Run Number;Time (ps)",plane.c_str()));
    gr2ave->DrawClone("APE");

    /*
    double fracttot = 0.;
    c->cd(4);
    TGraphErrors *gf1 = new TGraphErrors(n,runnum,fe,zeros,feunc);
    gf1->SetMarkerColor(kBlue);
    //gr1->SetMarkerStyle(21);
    gf1->SetMarkerStyle(kOpenCircle);
    gf1->SetTitle(Form("e fraction %s;Run Number;Fraction",plane.c_str()));
    gf1->DrawClone("APE");
    //f.SetLineColor(kRed); f.SetLineStyle(2);
    // Fit it to the graph and draw it
    gf1->Fit(&f);
    f.DrawClone("Same");
    p0 = f.GetParameter(0);
    //p1 = f.GetParameter(1);
    e0 = f.GetParError(0);
    //e1 = f.GetParError(1);
    chi2 = f.GetChisquare();
    //printf("mu const = %f +/- %f   linear term = %f +/- %f --> %f ps/100 runs\n",p0,e0,p1,e1,100*p1);
    printf("e fraction = %f +/- %f \n",p0,e0);
    stddev = 0.;
    for (int ii = 0; ii<n; ii++) {
        //stddev += TMath::Power(p0+p1*runnum[ii]-tmu[ii],2.);
        stddev += TMath::Power(p0-fe[ii],2.);
    }
    stddev = stddev/(n-1);
    stddev = TMath::Power(stddev,0.5);
    printf("e fraction chi^2 = %f for %d points, sigma = %f \n",chi2,n,stddev);
    fracttot += p0;

    c->cd(5);
    TGraphErrors *gf2 = new TGraphErrors(n,runnum,fmu,zeros,fmuunc);
    gf2->SetMarkerColor(kBlue);
    //gr1->SetMarkerStyle(21);
    gf2->SetMarkerStyle(kOpenCircle);
    gf2->SetTitle(Form("#mu fraction %s;Run Number;Fraction",plane.c_str()));
    gf2->DrawClone("APE");
    //f.SetLineColor(kRed); f.SetLineStyle(2);
    // Fit it to the graph and draw it
    gf2->Fit(&f);
    f.DrawClone("Same");
    p0 = f.GetParameter(0);
    //p1 = f.GetParameter(1);
    e0 = f.GetParError(0);
    //e1 = f.GetParError(1);
    chi2 = f.GetChisquare();
    //printf("mu const = %f +/- %f   linear term = %f +/- %f --> %f ps/100 runs\n",p0,e0,p1,e1,100*p1);
    printf("mu fraction = %f +/- %f \n",p0,e0);
    stddev = 0.;
    for (int ii = 0; ii<n; ii++) {
        //stddev += TMath::Power(p0+p1*runnum[ii]-tmu[ii],2.);
        stddev += TMath::Power(p0-fmu[ii],2.);
    }
    stddev = stddev/(n-1);
    stddev = TMath::Power(stddev,0.5);
    printf("mu fraction chi^2 = %f for %d points, sigma = %f \n",chi2,n,stddev);
    fracttot += p0;

    c->cd(6);
    TGraphErrors *gf3 = new TGraphErrors(n,runnum,fpi,zeros,fpiunc);
    gf3->SetMarkerColor(kBlue);
    //gr1->SetMarkerStyle(21);
    gf3->SetMarkerStyle(kOpenCircle);
    gf3->SetTitle(Form("#pi fraction %s;Run Number;Fraction",plane.c_str()));
    gf3->DrawClone("APE");
    //f.SetLineColor(kRed); f.SetLineStyle(2);
    // Fit it to the graph and draw it
    gf3->Fit(&f);
    f.DrawClone("Same");
    p0 = f.GetParameter(0);
    //p1 = f.GetParameter(1);
    e0 = f.GetParError(0);
    //e1 = f.GetParError(1);
    chi2 = f.GetChisquare();
    //printf("mu const = %f +/- %f   linear term = %f +/- %f --> %f ps/100 runs\n",p0,e0,p1,e1,100*p1);
    printf("pi fraction = %f +/- %f \n",p0,e0);
    stddev = 0.;
    for (int ii = 0; ii<n; ii++) {
        //stddev += TMath::Power(p0+p1*runnum[ii]-tmu[ii],2.);
        stddev += TMath::Power(p0-fpi[ii],2.);
    }
    stddev = stddev/(n-1);
    stddev = TMath::Power(stddev,0.5);
    printf("pi fraction chi^2 = %f for %d points, sigma = %f \n",chi2,n,stddev);
    fracttot += p0;
    printf(" total particle fraction = %f \n",fracttot);
    */

     // Build and Draw a legend
    //TLegend leg(.1,.7,.3,.9,"Mean Paddle Position");
    //leg.SetFillColor(0);
    //gr1->SetFillColor(0);
    //leg.AddEntry(&gr1,"Exp. Points");
    //leg.AddEntry(&f,"Average");
    //leg.DrawClone("Same");

	//gPad->Update();
	//gPad->Modified();

	return;
}
	///// reference
	///// TCanvas: https://root.cern.ch/doc/master/classTCanvas.html
	///// Histograms and Plotting: https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
	///// TF1: https://root.cern.ch/doc/master/classTF1.html
	///// TGraph: https://root.cern.ch/doc/master/classTGraph.html
	///// TAxis: https://root.cern.ch/doc/master/classTAxis.html

    // calculate peak positions for qdc_hit_BHC[paddle#][right] and BHD and tdc_coin_left_BHC[paddle#]-RF
    // calculate uncertainity for each