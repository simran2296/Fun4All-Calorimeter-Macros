#include <TFile.h>
#include <TH2.h>
#include <TPad.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <utility> // for std::pair

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::pair;


// forward declaration, implemented at the bottom
TProfile* ProfileAndPlot ( TH2* hin);
vector< pair<double,double> > Get_ge_te( TH2* hin );
  
int eres_forward()
{
  //read uncalibrated delta e/e plots from root file //i dont need this step
  
  TFile * fin = new TFile("tree2.root", "READ");
  TH2D* femc_raw = (TH2D*) fin->Get("hist_energy_FEMC_1"); 
  femc_raw->SetName("femc_raw");
  TH2D* fhcal_raw = (TH2D*) fin->Get("hist_energy_FHCAL_1");
  fhcal_raw->SetName("fhcal_raw");
  new TCanvas;
  // Profile and plot // make a profile plot on the main plot (not separately)
  TProfile* femc_raw_p = ProfileAndPlot( femc_raw ); 
  TProfile* fhcal_raw_p = ProfileAndPlot( fhcal_raw ); 

  // Get calibration fits
  new TCanvas;
  TF1* femc_fit = new TF1("femc_fit","pol4",0.1,30);
  femc_fit->SetLineColor(kGreen+1);   femc_fit->SetLineWidth(4); 
  new TCanvas;
  femc_raw_p->Fit(femc_fit); //fitting the profile+raw plot here with pol4
  
  TF1* fhcal_fit = new TF1("fhcal_fit","pol4",0.1,30);
  // hcalin_fit->SetParameters( -1,0.01,0.1);
  fhcal_fit->SetLineColor(kGreen+1);   fhcal_fit->SetLineWidth(4); 
  new TCanvas;
  fhcal_raw_p->Fit(fhcal_fit);

 
  //i have all this upto this point

  // These fit functions approximate the expectation value of measured e ("te") for a given ge
  // f (ge) = mu( (te-ge) / ge ) = mu(te) / ge - 1
  // <-> mu( te ) = ( 1+f(ge) ) * ge
  // IF we only had one detector, we could now introduce ce (for calibrated or corrected)
  // ce == te / (f(ge)+1) --> mu(ce) = mu( te ) / (f(ge)+1) = ( 1+f(ge) )*ge/ ( 1+f(ge) ) = ge
  // and thus have calibrated ce centered on truth energy
  // However, we want to calibrate three detectors and add those up, so we'd end up with 3*ge in the final version
  // One way would be to just divide by three, but that gives uneven weight.
  // In our case, the HCalOut catches about 45% of the ge, the Hcalin less than 10% and the EMC between 20 and 70(!) %
  // I don't love the strong E-dependence of the EMC, but nevertheless, the simple approach we take
  // is to scale the ce from above by the _average_ value of the given system,
  // simply by using (one plus) the y-average from the plot
  // so:
  // ce_tot= 0.44* ce_hcalout + 0.09*ce_hcalin + 0.45*ce_cemc
  // doesn't quite add up to 1, partly because I rounded and partly because it's inexact,
  // but we can run a small afterburner calibration to get rid of the rest
  // Let's get those numbers from the histo though
  double p_femc = 1 + femc_raw->GetMean(2);
  double p_fhcal = 1 + fhcal_raw->GetMean(2);
  cout << "Sum of weights is " << p_femc << " + "
       << p_fhcal << " = "
       << p_femc + p_fhcal << endl;
  // Missing 5%, that's okay. We could reshuffle here, but better to do it at the very end with another fit

  // NOTE: My original idea was to integrate and rescale the fit functions, but just multiplying
  // them with weights like this should be equivalent and conceptually cleaner.

  // Now I need ge and te for the different calos
  TTree * t = (TTree*) fin->Get("t");
  Float_t truth_e, te_FEMC, te_FHCAL;

  t->SetBranchAddress("truth_e",&truth_e);
  t->SetBranchAddress("te_FEMC",&te_FEMC);
  t->SetBranchAddress("te_FHCAL",&te_FHCAL);

  TH2D* allcalos_precalibrated = (TH2D*) femc_raw->Clone("allcalos_precalibrated");
  allcalos_precalibrated->Reset();
  allcalos_precalibrated->SetTitle("#frac{te-ge}{ge} vs ge (weighted sum)");  
  for(int i=0; i<t->GetEntries(); i++)  {
    t->GetEntry(i);

    // This is what we get without precalibration. Useless, just for bookkeeping
    float etot_raw = 0;
    if ( te_FEMC >-1 ){ etot_raw+=te_FEMC; }
    if ( te_FHCAL >-1 ){ etot_raw+=te_FHCAL; }
    
    // This is what we want! Same as before but individually precalibrated and then scaled 
    float etot_precalib = 0;
    if ( te_FEMC >-1 ){
      auto ce = te_FEMC / (femc_fit->Eval(truth_e)+1);
      etot_precalib+=ce*p_femc;
    }
    if ( te_FHCAL >-1 ){
      auto ce = te_FHCAL / (fhcal_fit->Eval(truth_e)+1);
      etot_precalib+=ce*p_fhcal;
    }
    //   cout<<"etot_precalib="<<etot_precalib<<endl;
    if(etot_precalib>0.0)
      allcalos_precalibrated->Fill (truth_e, (etot_precalib-truth_e)/truth_e);
  }
  TProfile* allcalos_precalibrated_p = ProfileAndPlot( allcalos_precalibrated ); 
  
  // Final pass. Because the weights don't quite add up, do another simpler calibration
  TF1* final_fit = new TF1("final_fit","pol3",2.0,30); // pol0 might be enough, but there is still some shape
  final_fit->SetRange(2.0,30);
  final_fit->SetLineColor(kGreen+1);   final_fit->SetLineWidth(4); 
  new TCanvas;
  allcalos_precalibrated_p->Fit(final_fit);

  
  TH2D* allcalos_fullcalibrated = (TH2D*) femc_raw->Clone("allcalos_fullcalibrated");
  allcalos_fullcalibrated->Reset();
  allcalos_fullcalibrated->SetTitle("#frac{te-ge}{ge} vs ge (final calibration)");  
  //for resolution
  const Int_t NBINS = 15;
  double edges[NBINS + 1] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0};
  TH2D *hist_energy_FEMC_FHCAL_1 = new TH2D("hist_energy_FEMC_FHCAL_1","#frac{te-ge}{ge} vs ge (FEMC+FHCAL); ge; #frac{te-ge}{ge}",200,0,30,200,-2.6,2.6);
  TH2D *binned_energy_FEMC_FHCAL_1 = new TH2D("binned_energy_FEMC_FHCAL_1","#frac{te-ge}{ge} vs ge (FEMC+FHCAL); ge; #frac{te-ge}{ge}",NBINS,edges,200,-2.6,2.6);


  for(int i=0; i<t->GetEntries(); i++)  {
    t->GetEntry(i);
    
    // Same as above but also divided by the last shape
    float etot_fullcalib = 0;
    if ( te_FEMC >-1 ){
      auto ce = te_FEMC / (femc_fit->Eval(truth_e)+1);
      etot_fullcalib+=ce*p_femc;
    }
    if ( te_FHCAL >-1 ){
      auto ce = te_FHCAL / (fhcal_fit->Eval(truth_e)+1);
      etot_fullcalib+=ce*p_fhcal;
    }

    
    etot_fullcalib /= ( final_fit->Eval(truth_e) +1 ); // <-- This is the change
    
    allcalos_fullcalibrated->Fill (truth_e, (etot_fullcalib-truth_e)/truth_e);
    hist_energy_FEMC_FHCAL_1->Fill (truth_e, (etot_fullcalib-truth_e)/truth_e);
    binned_energy_FEMC_FHCAL_1->Fill (truth_e, (etot_fullcalib-truth_e)/truth_e);
  }
  TProfile* allcalos_fullcalibrated_p = ProfileAndPlot( allcalos_fullcalibrated ); 
  
  //resolution
  TH1D *e_hist2 = new TH1D("e_hist2"," ", NBINS, edges); //FEMC +FHCAL
  TCanvas *c1 = new TCanvas("c1");

  TString arr[NBINS]; // Used for naming fitted slices used in sigma_e vs ge 
  for(int i = 0; i < NBINS; i++){
    arr[i] = TString::Itoa(i + 1, 10);
  }
 

  ofstream fout2("ERes_sigma_FEMC_FHCAL.dat");

  //FEMC+FHCAL
  for(int i=0; i<NBINS; i++)
    {
      int j = i + 1; 
      TString name = "slice"+arr[i];
      TH1D *h2 =  binned_energy_FEMC_FHCAL_1->ProjectionY(name, j, j);
      //rebin here if needed
      h2->Rebin(2);
      float maxy2 = h2->GetMaximum();
      float max2 = maxy2 + 10;
      //    h2->GetYaxis()->SetRange(0.,max2);
      // h2->GetXaxis()->SetRange(-2,2);
      //gaussian fit for slices 
      
      TF1 *Fgauss2=new TF1(name,"gaus",-2,2);
      // if(j==2)
      // 	Fgauss2->SetRange(-0.1, 1); // for ymean fit
      
      Fgauss2->SetLineStyle(2);
      Fgauss2->SetLineColor(4); 
      h2->SetLineColor(2);
     
      c1->cd(1);
      gStyle->SetOptStat(11);
      gStyle->SetOptFit(102);

      h2->Draw("ELP");
      h2->Fit(Fgauss2,"RM");
      
      e_hist2->SetBinContent(j,Fgauss2->GetParameter(2));    
      e_hist2->SetBinError(j,Fgauss2->GetParError(2));

      fout2<<"slice"<<j<<"    "<<Fgauss2->GetParameter(2)<<"    "<<Fgauss2->GetParError(2)<<endl;
      c1->Print("slices_FEMC_FHCAL/slice_"+arr[i]+"_FEMC_FHCAL.gif");
      


      }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  e_hist2->SetTitle("#sigma_{E}/E vs ge (FEMC+FHCAL)");
  e_hist2->GetXaxis()->SetTitle(" Generated energy (GeV)");
  e_hist2->GetXaxis()->SetLabelFont(42);
  e_hist2->GetXaxis()->SetLabelSize(0.035);
  e_hist2->GetXaxis()->SetTitleSize(0.05);
  e_hist2->GetXaxis()->SetTitleOffset(0.91);
  e_hist2->GetXaxis()->SetTitleFont(42);
  e_hist2->GetYaxis()->SetTitle("#sigma_{E}/E");
  e_hist2->GetYaxis()->SetNdivisions(507);
  e_hist2->GetYaxis()->SetLabelFont(42);
  e_hist2->GetYaxis()->SetLabelSize(0.035);
  e_hist2->GetYaxis()->SetTitleSize(0.05);
  e_hist2->GetYaxis()->SetTitleOffset(0.93);
  e_hist2->GetYaxis()->SetTitleFont(42);
  e_hist2->SetMinimum(0);
  e_hist2->SetMaximum(1);
  e_hist2->SetMarkerColor(2); 
  e_hist2->SetLineColor(2);
  e_hist2->SetMarkerStyle(20);

  TF1 *Efit2 = new TF1("Efit2","[0] + [1] / sqrt(x)", 0.1, 30);
  Efit2->SetParameter(0, 0.02); // constant term
  Efit2->SetParameter(1, 0.1); // 1/sqrt(E)
  Efit2->SetLineStyle(2);
  Efit2->SetLineColor(2);
  e_hist2->Fit("Efit2","","",0.1,30);
  TF1 *TrueFit2 = new TF1("TrueFit2","0.10 + 0.50 / sqrt(x)", 0.1, 30);
  TrueFit2->SetLineStyle(2);
  TrueFit2->SetLineColor(4);
  TrueFit2->Draw("same");

  TLegend *leg2 = new TLegend(0.4300885,0.6604775,0.8513274,0.8355438,NULL,"brNDC");
   leg2->SetBorderSize(0);
   leg2->SetTextSize(0.033);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(1001);
   leg2->AddEntry(e_hist2, "FEMC+FHCAL");
   leg2->AddEntry(TrueFit2, "Required: #sigma_{E}/E = 10% + 50%/#sqrt{E}");
   leg2->AddEntry(Efit2, "Obtained: #sigma_{E}/E = p0 + p1/#sqrt{E}");
   leg2->Draw();
   c1->Print("Resolution_FEMC_FHCAL.gif");

   c1->Close();


  cout<<"FEMC+FHCAL: p0="<<Efit2->GetParameter(0)<<"; p1="<<Efit2->GetParameter(1)<<endl;
  TFile *f = new TFile("Pion_forward_plots.root","RECREATE");
  f->GetList()->Add(hist_energy_FEMC_FHCAL_1);
  f->GetList()->Add(binned_energy_FEMC_FHCAL_1);
  f->Write();

  
  
  return 0;
}

// =============================================================================

TProfile* ProfileAndPlot ( TH2* hin){
  // new TCanvas;
  // hin->DrawCopy("surf1");

  // Profiles - Could use option "s" to get standard deviation, but for this "" works
  TString profname=hin->GetName();  profname += "_p";
  TProfile* ret = hin->ProfileX(profname, 1, -1, "");
  
  new TCanvas;
  hin->DrawCopy("colz");
  ret->SetLineColor(kRed);
  ret->SetLineWidth(1);
  ret->DrawCopy("same");

  
  return ret;
}
// =============================================================================

vector< pair<double,double> > Get_ge_te( TH2* hin ){
  vector< pair<double,double> > ge_te;
  for (int i=1; i<=hin->GetNbinsX() ; ++i){
    auto ge = hin->GetXaxis()->GetBinCenter(i);
    for ( int j=1; j<=hin->GetNbinsY() ; ++j){
      if ( hin->GetBinContent(i,j) == 0 ) continue;
      auto h = hin->GetYaxis()->GetBinCenter(j); // = (te-ge) / ge <--> te = h*ge+ge
      auto te = h*ge + ge;
      // how many of these?
      for ( int k=0; k<hin->GetBinContent(i,j) ; ++k){
	ge_te.push_back( pair<double,double>(ge, te));
      }
    }    
  }

  // // test
  // TH2D* h = (TH2D*) hin->Clone("h"); h->Reset();
  // for (auto p : ge_te) {
  //   auto ge=p.first; auto te=p.second;
  //   h->Fill(ge, (te-ge)/ge);
  // }
  // new TCanvas;  hin->Draw("colz");
  // new TCanvas;  h->Draw("colz");
  
  return ge_te;
}
