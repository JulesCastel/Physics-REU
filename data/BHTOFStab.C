
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
//#include "TArrow.h"
//#include "TLatex.h"

void BHTOFStab()
{
    Double_t runnum[2500], zeros[2500];
    Double_t fe[2500], feunc[2500], ce[2500], ceunc[2500], we[2500], weunc[2500]; 
	Int_t j = 0;
    //Int_t p=160;
    //Int_t bar = 16;
    //std::string plane = "BHC";
    //Int_t bar = 6;
    std::string plane = "BHD";	///// load files
    //Int_t i1 = 15681, i2 = 15715;
    //Int_t i1 = 15412, i2 = 15655;
    //Int_t i1 = 15324, i2 = 15500;
    //Int_t i1 = 15650, i2 = 15717;
    //Int_t i1 = 15324, i2 = 15717;
    //+210 for report
    //Int_t i1 = 15688, i2 = 15738;
    //all +210
    //Int_t i1 = 15324, i2 = 15821;
    //all -210
    //Int_t i1 = 15490, i2 = 15600;
    //all 2021
    //Int_t i1 = 8940, i2 = 11334;
    //Int_t i1 = 11200, i2 = 11250;
    //Int_t i1 = 14714, i2 = 15190;
    Int_t i1 = 14000, i2 = 15600;
    for (Int_t irun=i1; irun<i2; irun++) {
        if( 
            irun!=15074 && irun != 15075 && irun != 15076 && irun != 15077 && irun != 15078
             && irun != 15079 && irun != 15080 && irun != 15081 && irun != 15082
              && irun != 14955 && irun != 14956 && irun != 14957
               && irun != 15154 && irun != 15156 && irun != 15157 && irun != 15158  
            )
        {

		TFile *file = new TFile(Form("/data2/processed_fast/BH_detail/run%d.root",irun));

		if ( file != NULL )  {

            ///// load histograms
            TH1F* h_rf = (TH1F*)file->Get("TOF/BHD06 - /BHC07/Electron TOF");
            //TH1F* h_rf = (TH1F*)file->Get("TOF/BHD06 - /BHC08/Electron TOF");
            //TH1F* h_rf = (TH1F*)file->Get(Form("%s/RF/Out of Time RF Avg of Plane",plane.c_str()));
            if ( h_rf != NULL ){
                runnum[j] = irun;
                zeros[j] = 0.;

                ///// fit histogram to determine rf peak positions
                TF1* f_gaus = new TF1("ElectronTOF","gaus",-0.25,0.25);
                h_rf->Fit(f_gaus,"R");
                //save fit values for plotting
                fe[j] = f_gaus->GetParameter(0);
                ce[j] = f_gaus->GetParameter(1);
                we[j] = f_gaus->GetParameter(2);
                feunc[j] = f_gaus->GetParError(0);
                ceunc[j] = f_gaus->GetParError(1);
                weunc[j] = f_gaus->GetParError(2);
                //remove large uncertainty results to allow seeing what is happening
                if ((ceunc[j]>0.02) || (we[j]<0.02)) {j--;}
                j++;
			}
            else {
                printf("Hist not found for run %d \n", irun);
            }
		}
        file->Close();
    }
	}

    const Int_t n = j;
    printf(" total number of points = %d \n",n);
    //convert from ns to ps
    for (Int_t il=0; il<n; il++){
        ce[il] = 1000*ce[il];
        we[il] = 1000*we[il];
        ceunc[il] = 1000*ceunc[il];
        weunc[il] = 1000*weunc[il];
    }

	///// create canvas and divide the canvas into two parts
    TCanvas *c0 = new TCanvas("BHD6BHC7eTOF","BHD6BHC7eTOF", 550, 400);
    //TCanvas *c0 = new TCanvas("BHD6BHC8eTOF","BHD6BHC8eTOF", 550, 400);
    //c->Divide(3,2); 

	//c->cd(1);
    TGraphErrors *gr0 = new TGraphErrors(n,runnum,ce,zeros,ceunc);
    gr0->SetMarkerColor(kBlue);
    //gr1->SetMarkerStyle(21);
    gr0->SetMarkerStyle(kOpenCircle);
    gr0->SetTitle("TOF_{e};Run Number;Time (ps)");
    gr0->DrawClone("APE");
    // Define a linear function
    //TF1 f("Linear law","[0]+x*[1]",15681,15715);
    TF1 f("Linear fit","[0]+x*[1]",i1,i2);
    //TF1 f("Linear fit","[0]",i1,i2);
    // Let's make the function line nicer
    f.SetLineColor(kRed); f.SetLineStyle(2);
    // Fit the graph and draw it
    gr0->Fit(&f);
    f.DrawClone("Same");
    Double_t p0 = f.GetParameter(0);
    Double_t p1 = f.GetParameter(1);
    Double_t e0 = f.GetParError(0);
    Double_t e1 = f.GetParError(1);
    Double_t chi2 = f.GetChisquare();
    Double_t stddev = 0.;
    for (Int_t ii = 0; ii<n; ii++) {
        stddev += TMath::Power(p0+p1*runnum[ii]-ce[ii],2.);
        //stddev += TMath::Power(p0-ce[ii],2.);
    }
    stddev = stddev/(n-1);
    stddev = TMath::Power(stddev,0.5);
    printf("TOF const = %f +/- %f slope %f +/- %f  chi^2 = %f for %d points, sigma = %f \n",p0,e0,p1,e1,chi2,n,stddev);

    TCanvas *c1 = new TCanvas("BHD6BHC7eTOFwid","BHD6BHC7eTOFwid", 550, 400);
    //c->cd(2);
	TGraphErrors *gr1 = new TGraphErrors(n,runnum,we,zeros,weunc);
 	gr1->SetMarkerColor(kBlue);
 	//gr1->SetMarkerStyle(21);
 	gr1->SetMarkerStyle(kOpenCircle);
    gr1->SetTitle("TOF_{e} rms Width;Run Number;Width (ps)");
	gr1->DrawClone("APE");
    // Fit the graph and draw it
    gr1->Fit(&f);
    f.DrawClone("Same");
    p0 = f.GetParameter(0);
    p1 = f.GetParameter(1);
    e0 = f.GetParError(0);
    e1 = f.GetParError(1);
    chi2 = f.GetChisquare();
    stddev = 0.;
    for (Int_t ii = 0; ii<n; ii++) {
        stddev += TMath::Power(p0+p1*runnum[ii]-we[ii],2.);
        //stddev += TMath::Power(p0-we[ii],2.);
    }
    stddev = stddev/(n-1);
    stddev = TMath::Power(stddev,0.5);
    printf("width const = %f +/- %f slope %f +/- %f  chi^2 = %f for %d points, sigma = %f \n",p0,e0,p1,e1,chi2,n,stddev);

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
