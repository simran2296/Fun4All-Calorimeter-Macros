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
  
int eres_barrel()
{
  //read uncalibrated delta e/e plots from root file //i dont need this step
  
  TFile * fin = new TFile("tree1.root", "READ");
  TH2D* cemc_raw = (TH2D*) fin->Get("hist_energy_CEMC_1"); 
  cemc_raw->SetName("cemc_raw");
  TH2D* hcalin_raw = (TH2D*) fin->Get("hist_energy_HCALIN_1");
  hcalin_raw->SetName("hcalin_raw");
  TH2D* hcalout_raw = (TH2D*) fin->Get("hist_energy_HCALOUT_1");
  hcalout_raw->SetName("hcalout_raw");
  new TCanvas;
  // Profile and plot // make a profile plot on the main plot (not separately)
  TProfile* cemc_raw_p = ProfileAndPlot( cemc_raw ); 
  TProfile* hcalin_raw_p = ProfileAndPlot( hcalin_raw ); 
  TProfile* hcalout_raw_p = ProfileAndPlot( hcalout_raw );

  // Get calibration fits
  TCanvas *c3 =  new TCanvas("c3");;
  TF1* cemc_fit = new TF1("cemc_fit","pol4",0.1,30);
  cemc_fit->SetLineColor(kGreen+1);   cemc_fit->SetLineWidth(4); 
  new TCanvas;
  cemc_raw_p->Fit(cemc_fit); //fitting the profile+raw plot here with pol4
  
  TF1* hcalin_fit = new TF1("hcalin_fit","[0]+[1]/sqrt(x) + [2] /x",0.1,30);
  hcalin_fit->SetParameters( -1,0.01,0.1);
  hcalin_fit->SetLineColor(kGreen+1);   hcalin_fit->SetLineWidth(4); 
  new TCanvas;
  hcalin_raw_p->Fit(hcalin_fit);

  TF1* hcalout_fit = new TF1("hcalout_fit","pol3",0.1,30);
  hcalout_fit->SetLineColor(kGreen+1);   hcalout_fit->SetLineWidth(4); 
  new TCanvas;
  hcalout_raw_p->Fit(hcalout_fit);
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
  double p_cemc = 1 + cemc_raw->GetMean(2);
  double p_hcalin = 1 + hcalin_raw->GetMean(2);
  double p_hcalout = 1 + hcalout_raw->GetMean(2);
  cout << "Sum of weights is " << p_cemc << " + "
       << p_hcalin << " + "
       << p_hcalout << " = "
       << p_cemc + p_hcalin + p_hcalout << endl;
  // Missing 5%, that's okay. We could reshuffle here, but better to do it at the very end with another fit

  // NOTE: My original idea was to integrate and rescale the fit functions, but just multiplying
  // them with weights like this should be equivalent and conceptually cleaner.

  // Now I need ge and te for the different calos
  TTree * t = (TTree*) fin->Get("t");
  Float_t truth_e, te_CEMC, te_HCALIN, te_HCALOUT;

  t->SetBranchAddress("truth_e",&truth_e);
  t->SetBranchAddress("te_CEMC",&te_CEMC);
  t->SetBranchAddress("te_HCALIN",&te_HCALIN);
  t->SetBranchAddress("te_HCALOUT",&te_HCALOUT);

  TH2D* allcalos_precalibrated = (TH2D*) cemc_raw->Clone("allcalos_precalibrated");
  allcalos_precalibrated->Reset();
  allcalos_precalibrated->SetTitle("#frac{te-ge}{ge} vs ge (weighted sum)");  
  for(int i=0; i<t->GetEntries(); i++)  {
    t->GetEntry(i);

    // This is what we get without precalibration. Useless, just for bookkeeping
    float etot_raw = 0;
    if ( te_CEMC >-1 ){ etot_raw+=te_CEMC; }
    if ( te_HCALIN >-1 ){ etot_raw+=te_HCALIN; }
    if ( te_HCALOUT >-1 ){ etot_raw+=te_HCALOUT; }
    
    // This is what we want! Same as before but individually precalibrated and then scaled 
    float etot_precalib = 0;
    if ( te_CEMC >-1 ){
      auto ce = te_CEMC / (cemc_fit->Eval(truth_e)+1);
      etot_precalib+=ce*p_cemc;
    }
    if ( te_HCALIN >-1 ){
      auto ce = te_HCALIN / (hcalin_fit->Eval(truth_e)+1);
      etot_precalib+=ce*p_hcalin;
    }
    if ( te_HCALOUT >-1 ){
      auto ce = te_HCALOUT / (hcalout_fit->Eval(truth_e)+1);
      etot_precalib+=ce*p_hcalout;
    }
    if(etot_precalib>0.0)
      allcalos_precalibrated->Fill (truth_e, (etot_precalib-truth_e)/truth_e);
  }
  TProfile* allcalos_precalibrated_p = ProfileAndPlot( allcalos_precalibrated ); 
  
  // Final pass. Because the weights don't quite add up, do another simpler calibration
  TF1* final_fit = new TF1("final_fit","pol3",0.1,30); // pol0 might be enough, but there is still some shape
  final_fit->SetLineColor(kGreen+1);   final_fit->SetLineWidth(4); 
  new TCanvas;
  allcalos_precalibrated_p->Fit(final_fit);

  
  TH2D* allcalos_fullcalibrated = (TH2D*) cemc_raw->Clone("allcalos_fullcalibrated");
  allcalos_fullcalibrated->Reset();
  allcalos_fullcalibrated->SetTitle("#frac{te-ge}{ge} vs ge (final calibration)");  
  //for resolution
  const Int_t NBINS = 15;
  double edges[NBINS + 1] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0};
  TH2D *hist_energy_CEMC_HCALIN_HCALOUT_1 = new TH2D("hist_energy_CEMC_HCALIN_HCALOUT_1","#frac{te-ge}{ge} vs ge (CEMC+HCALIN+HCALOUT); ge; #frac{te-ge}{ge}",200,0,30,200,-2.6,2.6);
  TH2D *binned_energy_CEMC_HCALIN_HCALOUT_1 = new TH2D("binned_energy_CEMC_HCALIN_HCALOUT_1","#frac{te-ge}{ge} vs ge (CEMC+HCALIN+HCALOUT); ge; #frac{te-ge}{ge}",NBINS,edges,200,-2.6,2.6);

  for(int i=0; i<t->GetEntries(); i++)  {
    t->GetEntry(i);
    
    // Same as above but also divided by the last shape
    float etot_fullcalib = 0;
    if ( te_CEMC >-1 ){
      auto ce = te_CEMC / (cemc_fit->Eval(truth_e)+1);
      etot_fullcalib+=ce*p_cemc;
    }
    if ( te_HCALIN >-1 ){
      auto ce = te_HCALIN / (hcalin_fit->Eval(truth_e)+1);
      etot_fullcalib+=ce*p_hcalin;
    }
    if ( te_HCALOUT >-1 ){
      auto ce = te_HCALOUT / (hcalout_fit->Eval(truth_e)+1);
      etot_fullcalib+=ce*p_hcalout;
    }
    
    etot_fullcalib /= ( final_fit->Eval(truth_e) +1 ); // <-- This is the change
    
    allcalos_fullcalibrated->Fill (truth_e, (etot_fullcalib-truth_e)/truth_e);
    hist_energy_CEMC_HCALIN_HCALOUT_1->Fill (truth_e, (etot_fullcalib-truth_e)/truth_e);
    binned_energy_CEMC_HCALIN_HCALOUT_1->Fill (truth_e, (etot_fullcalib-truth_e)/truth_e);
  }
  TProfile* allcalos_fullcalibrated_p = ProfileAndPlot( allcalos_fullcalibrated ); 
  //resolution
  
  TH1D *e_hist1 = new TH1D("e_hist1"," ", NBINS, edges); //CEMC+HCALIN+HCALOUT
  TCanvas *c1 = new TCanvas("c1");

  TString arr[NBINS]; // Used for naming fitted slices used in sigma_e vs ge 
  for(int i = 0; i < NBINS; i++){
    arr[i] = TString::Itoa(i + 1, 10);
  }
 

  ofstream fout1("ERes_sigma_CEMC_HCALIN_HCALOUT.dat");

  //CEMC+HCALIN+HCALOUT
  for(int i=0; i<NBINS; i++)
    {
      int j = i + 1; 
      TString name = "slice"+arr[i];
      TH1D *h1 =  binned_energy_CEMC_HCALIN_HCALOUT_1->ProjectionY(name, j, j);
      //rebin here if needed
      h1->Rebin(2);
      float maxy1 = h1->GetMaximum();
      float max1 = maxy1 + 10;
      //   h1->GetYaxis()->SetRange(0.,max1);
      // h1->GetXaxis()->SetRange(-2,2);
      //gaussian fit for slices
      TF1 *Fgauss1=new TF1(name,"gaus",-2,2);
      if(i==1)
      	Fgauss1->SetRange(-1,0.9); //for y mean fit
      if(i==2)
	Fgauss1->SetRange(-1,0.9);
      Fgauss1->SetLineStyle(2);
      Fgauss1->SetLineColor(4); 
      h1->SetLineColor(2);

      c1->cd(1);
      gStyle->SetOptStat(11);
      gStyle->SetOptFit(102);
  
      h1->Draw("ELP");
      h1->Fit(Fgauss1,"RM");

      e_hist1->SetBinContent(j,Fgauss1->GetParameter(2));    
      e_hist1->SetBinError(j,Fgauss1->GetParError(2));

      fout1<<"slice"<<j<<"    "<<Fgauss1->GetParameter(2)<<"    "<<Fgauss1->GetParError(2)<<endl;
      c1->Print("slices_CEMC_HCALIN_HCALOUT/slice_"+arr[i]+"_CEMC_HCALIN_HCALOUT.gif");


      }
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  e_hist1->SetTitle("#sigma_{E}/E vs ge (CEMC+HCALIN+HCALOUT)");
  e_hist1->GetXaxis()->SetTitle(" Generated energy (GeV)");
  e_hist1->GetXaxis()->SetLabelFont(42);
  e_hist1->GetXaxis()->SetLabelSize(0.035);
  e_hist1->GetXaxis()->SetTitleSize(0.05);
  e_hist1->GetXaxis()->SetTitleOffset(0.91);
  e_hist1->GetXaxis()->SetTitleFont(42);
  e_hist1->GetYaxis()->SetTitle("#sigma_{E}/E");
  e_hist1->GetYaxis()->SetNdivisions(507);
  e_hist1->GetYaxis()->SetLabelFont(42);
  e_hist1->GetYaxis()->SetLabelSize(0.035);
  e_hist1->GetYaxis()->SetTitleSize(0.05);
  e_hist1->GetYaxis()->SetTitleOffset(0.93);
  e_hist1->GetYaxis()->SetTitleFont(42);
  e_hist1->SetMinimum(0);
  e_hist1->SetMaximum(1);
  e_hist1->SetMarkerColor(2); 
  e_hist1->SetLineColor(2);
  e_hist1->SetMarkerStyle(20);

  TF1 *Efit1 = new TF1("Efit1","[0] + [1] / sqrt(x)", 0.1, 30);
  Efit1->SetParameter(0, 0.02); // constant term
  Efit1->SetParameter(1, 0.1); // 1/sqrt(E)
  Efit1->SetLineStyle(2);
  Efit1->SetLineColor(2);
  e_hist1->Fit("Efit1","","",0.1,30);
  TF1 *TrueFit1 = new TF1("TrueFit1","0.135 + 0.649 / sqrt(x)", 0.1, 30);
  TrueFit1->SetLineStyle(2);
  TrueFit1->SetLineColor(4);
  TrueFit1->Draw("same");

  TLegend *leg1 = new TLegend(0.4300885,0.6604775,0.8513274,0.8355438,NULL,"brNDC");
  //TLegend *leg1 = new TLegend(1.75,1.75);
   leg1->SetBorderSize(0);
   leg1->SetTextSize(0.033);
   leg1->SetLineColor(1);
   leg1->SetLineStyle(1);
   leg1->SetLineWidth(1);
   leg1->SetFillColor(0);
   leg1->SetFillStyle(1001);
   leg1->AddEntry(e_hist1, "CEMC+HCALIN+HCALOUT");
   leg1->AddEntry(TrueFit1, "Required: #sigma_{E}/E = 13.5% + 64.9%/#sqrt{E}");
   leg1->AddEntry(Efit1, "Obtained:  #sigma_{E}/E = p0 + p1/#sqrt{E}");
   leg1->Draw();
   c1->Print("Resolution_CEMC_HCALIN_HCALOUT.gif");

   c1->Close();
   cout<<"CEMC+HCALIN+HCALOUT: p0="<<Efit1->GetParameter(0)<<"; p1="<<Efit1->GetParameter(1)<<endl;

  TFile *f = new TFile("Pion_barrel_plots.root","RECREATE");
  f->GetList()->Add(hist_energy_CEMC_HCALIN_HCALOUT_1);
  f->GetList()->Add(binned_energy_CEMC_HCALIN_HCALOUT_1);
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
