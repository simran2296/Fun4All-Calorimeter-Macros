/*
Macro: Energy resolution for electrons
Calorimeters: CEMC, FEMC, EEMC
To run: root -l electron_aggcut.C
Outputs: 
1. Electron_calorimeter_plots.root
2. Sigma vs generated energy plot: resolution_CEMC.gif, resolution_FEMC.gif, resolution_EEMC.gif
3. Fitted slices for sigma: create directories to save - slices_CEMC, slices_FEMC, slices_EEMC

*/

#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"


R__LOAD_LIBRARY(libeicqa_modules.so)

void electron_aggcut()
{
  TFile *f1 = new TFile("Eval_CEMC.root","READ");
  TFile *f2 = new TFile("Eval_FEMC.root","READ");
  TFile *f3 = new TFile("Eval_EEMC.root","READ");

  TTree* T1 = (TTree*)f1->Get("T");
  TTree* T2 = (TTree*)f2->Get("T");
  TTree* T3 = (TTree*)f3->Get("T");
 
  EvalRootTTree *evaltree1 = nullptr;
  EvalRootTTree *evaltree2 = nullptr;
  EvalRootTTree *evaltree3 = nullptr;

  const Int_t NBINS = 15;
  double edges[NBINS + 1] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0};
  float ecut = 0.1;

  auto prof1 = new TProfile("prof1","CEMC ",200,0,30,-2.6,2.6);
  auto prof2 = new TProfile("prof2","FEMC ",200,0,30,-2.6,2.6);
  auto prof3 = new TProfile("prof3","EEMC ",200,0,30,-2.6,2.6);

  //Energy resolution plots - towers
 
  TH2D *hist_energy_CEMC_1 = new TH2D("hist_energy_CEMC_1","#frac{te-ge}{ge} vs ge (CEMC); ge; #frac{te-ge}{ge}",200,0,30,200,-2.6,2.6);
  TH2D *hist_energy_FEMC_1 = new TH2D("hist_energy_FEMC_1","#frac{te-ge}{ge} vs ge (FEMC); ge; #frac{te-ge}{ge}",200,0,30,200,-2.6,2.6);
  TH2D *hist_energy_EEMC_1 = new TH2D("hist_energy_EEMC_1","#frac{te-ge}{ge} vs ge (EEMC); ge; #frac{te-ge}{ge}",200,0,30,200,-2.6,2.6);

  // binned plots
  TH2D *binned_energy_CEMC_1 = new TH2D("binned_energy_CEMC_1","#frac{te-ge}{ge} vs ge (CEMC); ge; #frac{te-ge}{ge}",NBINS,edges,200,-2.6,2.6);
  TH2D *binned_energy_FEMC_1 = new TH2D("binned_energy_FEMC_1","#frac{te-ge}{ge} vs ge (FEMC); ge; #frac{te-ge}{ge}",NBINS,edges,200,-2.6,2.6);
  // TH2D *binned_energy_EEMC_1 = new TH2D("binned_energy_EEMC_1","#frac{te-ge}{ge} vs ge (EEMC); ge; #frac{te-ge}{ge}",NBINS,edges,50,-1.5,0.5);//50 for orig
  TH2D *binned_energy_EEMC_1 = new TH2D("binned_energy_EEMC_1","#frac{te-ge}{ge} vs ge (EEMC); ge; #frac{te-ge}{ge}",NBINS,edges,400,-0.5,0.5);//for cali

  //Energy conservation plots
  TH2D *te_ge_CEMC = new TH2D("te_ge_CEMC","#frac{te}{ge} vs ge (CEMC); ge; #frac{te}{ge}",200,0,30,200,-1,5);
  TH2D *te_ge_FEMC = new TH2D("te_ge_FEMC","#frac{te}{ge} vs ge (FEMC); ge; #frac{te}{ge}",200,0,30,200,-1,5);
  TH2D *te_ge_EEMC = new TH2D("te_ge_EEMC","#frac{te}{ge} vs ge (EEMC); ge; #frac{te}{ge}",200,0,30,200,-1,5);

  // phi theta plots using tower
  TH2D *dphi_dtheta_CEMC = new TH2D("dphi_dtheta_CEMC","dphi vs dtheta (CEMC); ttheta - gtheta; tphi - gphi", 200, -1.5, 1.5, 200, -1.5, 1.5);
  TH2D *dphi_dtheta_FEMC = new TH2D("dphi_dtheta_FEMC","dphi vs dtheta (FEMC); ttheta - gtheta; tphi - gphi", 200, -1.5, 1.5, 200, -1.5, 1.5);
  TH2D *dphi_dtheta_EEMC = new TH2D("dphi_dtheta_EEMC","dphi vs dtheta (EEMC); ttheta - gtheta; tphi - gphi", 200, -1.5, 1.5, 200, -1.5, 1.5);
 
  T1->SetBranchAddress("DST#EvalTTree_CEMC",&evaltree1);
  T2->SetBranchAddress("DST#EvalTTree_FEMC",&evaltree2);
  T3->SetBranchAddress("DST#EvalTTree_EEMC",&evaltree3);
  for(int i=0; i<T1->GetEntries(); i++)
    {
      T1->GetEntry(i);
      T2->GetEntry(i);
      T3->GetEntry(i);
      
      float te1 = 0, te2 = 0, te3 = 0, te4 = 0, te5 = 0, te6 = 0, ge, ge1, ge2, ge3, ge4, ge5, ge6, geta1, geta2, geta3, geta4 ,geta5, geta6, gphi1, gphi2, gphi3, gphi4, gphi5, gphi6, tower1, tower2, tower3, tower4, tower5, tower6, ttheta1, ttheta2, ttheta3, ttheta4, ttheta5, ttheta6, gtheta1, gtheta2, gtheta3, gtheta4, gtheta5, gtheta6, tphi1, tphi2, tphi3, tphi4, tphi5, tphi6;
      
      ge1 = evaltree1->get_ge();      
      ge2 = evaltree2->get_ge();      
      ge3 = evaltree3->get_ge();
      
      geta1 = evaltree1->get_geta();
      geta2 = evaltree2->get_geta();
      geta3 = evaltree3->get_geta();

      gphi1 = evaltree1->get_gphi();
      gphi2 = evaltree2->get_gphi();
      gphi3 = evaltree3->get_gphi();

      gtheta1 = evaltree1->get_gtheta();
      gtheta2 = evaltree2->get_gtheta();
      gtheta3 = evaltree3->get_gtheta();
      
      //CEMC
      if( geta1 >= -1.5 && geta1 <= 1.2 ) 
	{
	  for (int i=0; i<evaltree1->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree1->get_tower(i);
	      if (twr)
		{
		  tower1 = twr->get_te();
		  ttheta1 = twr->get_ttheta();
		  tphi1 = twr->get_tphi();
	        
		  if(tphi1 >= - sqrt(0.25*0.25 - (ttheta1-gtheta1)*(ttheta1-gtheta1)) + gphi1 && tphi1 <= sqrt(0.25*0.25 - (ttheta1-gtheta1)*(ttheta1-gtheta1)) + gphi1 )
		    {
		      te1 = tower1 + te1;
		  	
		      dphi_dtheta_CEMC->Fill(ttheta1 - gtheta1, tphi1 - gphi1);
		    }
		}
	    }
	  
	  if(te1 > ecut)
	    {
	      te1 = te1/0.991468; //pol0
	      float rt_CEMC = te1/ge1;
	      te_ge_CEMC->Fill(ge1, rt_CEMC);
	      prof1->Fill(ge1, rt_CEMC);

	      float frac_CEMC_1 = (te1-ge1)/ge1;
	      hist_energy_CEMC_1->Fill(ge1, frac_CEMC_1);
	      binned_energy_CEMC_1->Fill(ge1, frac_CEMC_1);
	    }
	}
      

      
      //FEMC
      
      if(geta2 >= 1.3 && geta2 <= 3.3)
	{
	  for (int i=0; i<evaltree2->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree2->get_tower(i);
	      if (twr)
		{
		  tower2 = twr->get_te();
		  ttheta2 = twr->get_ttheta();
		  tphi2 = twr->get_tphi();

		  if(tphi2 >= - sqrt(0.40*0.40 - (ttheta2-gtheta2)*(ttheta2-gtheta2)) + gphi2 && tphi2 <= sqrt(0.40*0.40 - (ttheta2-gtheta2)*(ttheta2-gtheta2)) + gphi2 ) //circle
		    {
		      te2 = tower2 + te2;
		      
		      dphi_dtheta_FEMC->Fill(ttheta2 - gtheta2, tphi2 - gphi2);
		      
		    }
		}
	    }
	  if(te2 > ecut)
	    {
	      te2 = te2/0.830733; //ol0
	      float rt_FEMC = te2/ge2;
	      te_ge_FEMC->Fill(ge2, rt_FEMC);
	      prof2->Fill(ge2, rt_FEMC);

	      float frac_FEMC_1 = (te2-ge2)/ge2;
	      hist_energy_FEMC_1->Fill(ge2, frac_FEMC_1);
	      binned_energy_FEMC_1->Fill(ge2, frac_FEMC_1);
	    }
	  
	}      
      
      
      
      //EEMC

      if(geta3 >= -3.5 && geta3 <= -1.7)
	{
	  for (int i=0; i<evaltree3->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree3->get_tower(i);
	      if (twr)
		{
		  tower3 = twr->get_te();
		  ttheta3 = twr->get_ttheta();
		  tphi3 = twr->get_tphi();

		  if(tphi3 >= - sqrt(0.50*0.50 - (ttheta3-gtheta3)*(ttheta3-gtheta3)) + gphi3 && tphi3 <= sqrt(0.50*0.50 - (ttheta3-gtheta3)*(ttheta3-gtheta3)) + gphi3 )
		    {
		      te3 = tower3 + te3;
		      
		      dphi_dtheta_EEMC->Fill(ttheta3 - gtheta3, tphi3 - gphi3);

		    }
		}
	    }
	  if(te3 > ecut)
	    {
	      te3 = te3/0.939041; //pol0
	      float rt_EEMC = te3/ge3;
	      te_ge_EEMC->Fill(ge3, rt_EEMC);
	      prof3->Fill(ge3, rt_EEMC);
	      
	      float frac_EEMC_1 = (te3-ge3)/ge3;
	      hist_energy_EEMC_1->Fill(ge3, frac_EEMC_1);
	      binned_energy_EEMC_1->Fill(ge3, frac_EEMC_1);
	    }
	}    
    }
  



  TH1D *e_hist1 = new TH1D("e_hist1"," ", NBINS, edges); //CEMC
  TH1D *e_hist2 = new TH1D("e_hist2"," ", NBINS, edges); //FEMC
  TH1D *e_hist3 = new TH1D("e_hist3"," ", NBINS, edges); //EEMC
  TCanvas *c1 = new TCanvas("c1");

  TString arr[NBINS]; // Used for naming fitted slices used in sigma_e vs ge 
  for(int i = 0; i < NBINS; i++){
    arr[i] = TString::Itoa(i + 1, 10);
  }
 

  ofstream fout1("ERes_sigma_CEMC.dat");
  ofstream fout2("ERes_sigma_FEMC.dat");
  ofstream fout3("ERes_sigma_EEMC.dat");

  //CEMC
  for(int i=0; i<NBINS; i++)
    {
      int j = i + 1; 
      TString name = "slice"+arr[i];
      TH1D *h1 =  binned_energy_CEMC_1->ProjectionY(name, j, j);
      //rebin here if needed
      float maxy1 = h1->GetMaximum();
      float max1 = maxy1 + 10;
   
      h1->SetAxisRange(-1,1,"X");
      //gaussian fit for slices
      TF1 *Fgauss1=new TF1(name,"gaus",-1,1);
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
      c1->Print("slices_CEMC/slice_"+arr[i]+"_CEMC.gif");


      }
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  e_hist1->SetTitle("#sigma_{E}/E vs ge (CEMC)");
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
  e_hist1->SetMaximum(0.4);
  e_hist1->SetMarkerColor(2); 
  e_hist1->SetLineColor(2);
  e_hist1->SetMarkerStyle(20);

  TF1 *Efit1 = new TF1("Efit1","[0] + [1] / sqrt(x) + [2] / x", 0.1, 30);
  Efit1->SetParameter(0, 0.02); // constant term
  Efit1->SetParameter(1, 0.1); // 1/sqrt(E)
  Efit1->SetParameter(2, 0.02); // 1/sqrt(E)
  Efit1->SetLineStyle(2);
  Efit1->SetLineColor(2);
  e_hist1->Fit("Efit1","","",0.1,30); //minimum req
  TF1 *TrueFit1 = new TF1("TrueFit1","0.02 + 0.13 / sqrt(x) + 0.02 / x", 0.1, 30);
  TrueFit1->SetLineStyle(2);
  TrueFit1->SetLineColor(4);
  TrueFit1->Draw("same");

  TLegend *leg1 = new TLegend(0.3300885,0.6604775,0.8513274,0.8355438,NULL,"brNDC");
  //TLegend *leg1 = new TLegend(1.75,1.75);
   leg1->SetBorderSize(0);
   leg1->SetTextSize(0.033);
   leg1->SetLineColor(1);
   leg1->SetLineStyle(1);
   leg1->SetLineWidth(1);
   leg1->SetFillColor(0);
   leg1->SetFillStyle(1001);
   leg1->AddEntry(e_hist1, "CEMC");
   leg1->AddEntry(TrueFit1, "Required: #sigma_{E}/E = 2% + 13%/#sqrt{E} + 2%/E");
   leg1->AddEntry(Efit1, "Obtained:  #sigma_{E}/E = p0 + p1/#sqrt{E} + p2/E");
   leg1->Draw();
   c1->Print("Resolution_CEMC.gif");



  //FEMC
  for(int i=0; i<NBINS; i++)
    {
      int j = i + 1; 
      TString name = "slice"+arr[i];
      TH1D *h2 =  binned_energy_FEMC_1->ProjectionY(name, j, j);
      //rebin here if needed
      float maxy2 = h2->GetMaximum();
      float max2 = maxy2 + 10;
      //    h2->GetYaxis()->SetRange(0.,max2);
      //   h2->GetXaxis()->SetRange(-1,1);
      h2->SetAxisRange(-0.5,0.4,"X");
      //gaussian fit for slices 
      
      TF1 *Fgauss2=new TF1(name,"gaus",-1,1);
      if(j==1)
	Fgauss2->SetRange(-0.1,0.1);
      /*   if(j==1)
	   Fgauss2->SetRange(-0.5, 0.0);*/ //for orignal
      
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
      c1->Print("slices_FEMC/slice_"+arr[i]+"_FEMC.gif");
      


      }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  e_hist2->SetTitle("#sigma_{E}/E vs ge (FEMC)");
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
  e_hist2->SetMaximum(0.2);
  e_hist2->SetMarkerColor(2); 
  e_hist2->SetLineColor(2);
  e_hist2->SetMarkerStyle(20);

  TF1 *Efit2 = new TF1("Efit2","[0] + [1] / sqrt(x) + [2] / x", 0.1, 30);
  Efit2->SetParameter(0, 0.02); // constant term
  Efit2->SetParameter(1, 0.1); // 1/sqrt(E)
  Efit1->SetParameter(2, 0.02); // 1/E
  Efit2->SetLineStyle(2);
  Efit2->SetLineColor(2);
  e_hist2->Fit("Efit2","","",0.1,30);
  TF1 *TrueFit2 = new TF1("TrueFit2","0.02 + 0.06 / sqrt(x) + 0.02 / x", 0.1, 30);
  TrueFit2->SetLineStyle(2);
  TrueFit2->SetLineColor(4);
  TrueFit2->Draw("same");

  TLegend *leg2 = new TLegend(0.3300885,0.6604775,0.8513274,0.8355438,NULL,"brNDC");
  //  TLegend *leg2 = new TLegend(1.75,1.75);
   leg2->SetBorderSize(0);
   leg2->SetTextSize(0.033);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(1001);
   leg2->AddEntry(e_hist2, "FEMC");
   leg2->AddEntry(TrueFit2, "Required: #sigma_{E}/E = 2% + (4-12%)/#sqrt{E} + 2%/E");
   leg2->AddEntry(Efit2, "Obtained: #sigma_{E}/E = p0 + p1/#sqrt{E} + p2/E");
   leg2->Draw();
   c1->Print("Resolution_FEMC.gif");

   /* delete c1;
 c1 = new TCanvas("c1"); 
   */

  //EEMC
  for(int i=0; i<NBINS; i++)
    {
      int j = i + 1; 
      TString name = "slice"+arr[i];
      TH1D *h3 =  binned_energy_EEMC_1->ProjectionY(name, j, j);
      //rebin here if needed
      //  h3->Rebin(2);
      float maxy3 = h3->GetMaximum();
      float max3 = maxy3 + 10;
      //    h3->SetAxisRange(-1,0.1,"X");//for original
      h3->SetAxisRange(-0.2,0.2,"X");//for cali

      //    h2->GetYaxis()->SetRange(0.,max2);
      //  h3->GetXaxis()->SetRange(-1,1);
      //gaussian fit for slices 
      
      //  TF1 *Fgauss3=new TF1(name,"gaus",0.02,0.1);
      TF1 *Fgauss3=new TF1(name,"crystalball",-0.2,0.2); //crystal ball fuction
      Fgauss3->SetParameter(0,65.36);
      Fgauss3->SetParameter(1,0.0436);
      Fgauss3->SetParameter(2,0.009424);
      Fgauss3->SetParameter(3,0.414);
      Fgauss3->SetParameter(4,4.993e+08);
      //  if(j==10)
      //	   Fgauss3->SetRange(-0.04, 0.1); //for original
      
      Fgauss3->SetLineStyle(2);
      Fgauss3->SetLineColor(4); 
      h3->SetLineColor(2);
     
      c1->cd(1);
      gStyle->SetOptStat(11);
      gStyle->SetOptFit(102);

      h3->Draw("ELP");
      h3->Fit(Fgauss3,"RM");
      
      e_hist3->SetBinContent(j,Fgauss3->GetParameter(2));    
      e_hist3->SetBinError(j,Fgauss3->GetParError(2));

      fout3<<"slice"<<j<<"    "<<Fgauss3->GetParameter(2)<<"    "<<Fgauss3->GetParError(2)<<endl;

      /*    e_hist3->SetBinContent(j,Fgauss3->GetParameter(2));    
      e_hist3->SetBinError(j,Fgauss3->GetParError(3));

      fout3<<"slice"<<j<<"    "<<Fgauss3->GetParameter(3)<<"    "<<Fgauss3->GetParError(3)<<endl;*/
      c1->Print("slices_EEMC/slice_"+arr[i]+"_EEMC.gif");
      


      }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  e_hist3->SetTitle("#sigma_{E}/E vs ge (EEMC)");
  e_hist3->GetXaxis()->SetTitle(" Generated energy (GeV)");
  e_hist3->GetXaxis()->SetLabelFont(42);
  e_hist3->GetXaxis()->SetLabelSize(0.035);
  e_hist3->GetXaxis()->SetTitleSize(0.05);
  e_hist3->GetXaxis()->SetTitleOffset(0.91);
  e_hist3->GetXaxis()->SetTitleFont(42);
  e_hist3->GetYaxis()->SetTitle("#sigma_{E}/E");
  e_hist3->GetYaxis()->SetNdivisions(507);
  e_hist3->GetYaxis()->SetLabelFont(42);
  e_hist3->GetYaxis()->SetLabelSize(0.035);
  e_hist3->GetYaxis()->SetTitleSize(0.05);
  e_hist3->GetYaxis()->SetTitleOffset(0.93);
  e_hist3->GetYaxis()->SetTitleFont(42);
  e_hist3->SetMinimum(0);
  e_hist3->SetMaximum(0.1);
  e_hist3->SetMarkerColor(2); 
  e_hist3->SetLineColor(2);
  e_hist3->SetMarkerStyle(20);

  TF1 *Efit3 = new TF1("Efit3","[0] + [1] / sqrt(x) + [2] / x", 0.1, 30);
  Efit3->SetParameter(0, 0.02); // constant term
  Efit3->SetParameter(1, 0.1); // 1/sqrt(E)
  Efit3->SetParameter(2, 0.02); // 1/E
  Efit3->SetLineStyle(2);
  Efit3->SetLineColor(2);
  e_hist3->Fit("Efit3","","",0.1,30);
  TF1 *TrueFit3 = new TF1("TrueFit3","0.01 + 0.025 / sqrt(x) + 0.01 / x", 0.1, 30);
  TrueFit3->SetLineStyle(2);
  TrueFit3->SetLineColor(4);
  TrueFit3->Draw("same");

  TLegend *leg3 = new TLegend(0.3300885,0.6604775,0.8513274,0.8355438,NULL,"brNDC");
  //  TLegend *leg2 = new TLegend(1.75,1.75);
   leg3->SetBorderSize(0);
   leg3->SetTextSize(0.033);
   leg3->SetLineColor(1);
   leg3->SetLineStyle(1);
   leg3->SetLineWidth(1);
   leg3->SetFillColor(0);
   leg3->SetFillStyle(1001);
   leg3->AddEntry(e_hist3, "EEMC");
   leg3->AddEntry(TrueFit3, "Required: #sigma_{E}/E = 1% + 2.5%/#sqrt{E} + 1%/E");
   leg3->AddEntry(Efit3, "Obtained: #sigma_{E}/E = p0 + p1/#sqrt{E} + p2/E");
   leg3->Draw();
   c1->Print("Resolution_EEMC.gif");

   c1->Close();


  cout<<"CEMC: p0="<<Efit1->GetParameter(0)<<"; p1="<<Efit1->GetParameter(1)<<"; p2="<<Efit1->GetParameter(2)<<endl;
  cout<<"FEMC: p0="<<Efit2->GetParameter(0)<<"; p1="<<Efit2->GetParameter(1)<<"; p2="<<Efit2->GetParameter(2)<<endl;
  cout<<"EEMC: p0="<<Efit3->GetParameter(0)<<"; p1="<<Efit3->GetParameter(1)<<"; p2="<<Efit3->GetParameter(2)<<endl;




  TF1 *fit = new TF1("fit","[0]*exp(-[1]/x)", 0.0, 30.0);
  //  fit->SetParameter(0, 0.0); // constant term
  //  fit->SetParameter(1, 0.0); // 1/sqrt(E)

  prof1->Fit("pol0");
  prof1->GetXaxis()->SetTitle("ge");
  prof1->GetYaxis()->SetTitle("Mean of te/ge");
  /* TF1 *func1 = prof1->GetFunction("pol3"); 
  double int1 = func1->Integral(0.0, 30.0);
  cout<<"integral="<<int1<<endl;*/

  prof2->Fit("pol0");
  prof2->GetXaxis()->SetTitle("ge");
  prof2->GetYaxis()->SetTitle("Mean of te/ge");
  /*  TF1 *func2 = prof5->GetFunction("pol3"); 
  double int2 = func2->Integral(0.0, 30.0);
  cout<<"integral="<<int2<<endl;*/

  prof3->Fit("pol0");
  prof3->GetXaxis()->SetTitle("ge");
  prof3->GetYaxis()->SetTitle("Mean of te/ge");
  /*  TF1 *func2 = prof5->GetFunction("pol3"); 
  double int2 = func2->Integral(0.0, 30.0);
  cout<<"integral="<<int2<<endl;*/


  TFile *f = new TFile("Electron_calorimeter_plots.root","RECREATE");
  f->GetList()->Add(te_ge_CEMC);
  f->GetList()->Add(hist_energy_CEMC_1);
  f->GetList()->Add(binned_energy_CEMC_1);
  f->GetList()->Add(te_ge_FEMC);
  f->GetList()->Add(hist_energy_FEMC_1);
  f->GetList()->Add(binned_energy_FEMC_1);
  f->GetList()->Add(te_ge_EEMC);
  f->GetList()->Add(hist_energy_EEMC_1);                        
  f->GetList()->Add(binned_energy_EEMC_1);                        
  f->GetList()->Add(dphi_dtheta_CEMC);
  f->GetList()->Add(dphi_dtheta_FEMC);
  f->GetList()->Add(dphi_dtheta_EEMC);
  f->GetList()->Add(prof1);
  f->GetList()->Add(prof2);
  f->GetList()->Add(prof3);
  f->Write();
}

