//note: eemc tower cluster not switched, eemc not updated

#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"


R__LOAD_LIBRARY(libeicqa_modules.so)

void eres_muon_para()
{
  TFile *f1 = new TFile("/gpfs02/eic/iiti/July21/muon_characterization/merged_Eval_CEMC.root","READ");
  TFile *f2 = new TFile("/gpfs02/eic/iiti/July21/muon_characterization/merged_Eval_HCALIN.root","READ");
  TFile *f3 = new TFile("/gpfs02/eic/iiti/July21/muon_characterization/merged_Eval_HCALOUT.root","READ");
  TFile *f4 = new TFile("/gpfs02/eic/iiti/July21/muon_characterization/merged_Eval_FEMC.root","READ");
  TFile *f5 = new TFile("/gpfs02/eic/iiti/July21/muon_characterization/merged_Eval_FHCAL.root","READ");
  TFile *f6 = new TFile("/gpfs02/eic/iiti/July21/muon_characterization/merged_Eval_EEMC.root","READ");

  TTree* T1 = (TTree*)f1->Get("T");
  TTree* T2 = (TTree*)f2->Get("T");
  TTree* T3 = (TTree*)f3->Get("T");
  TTree* T4 = (TTree*)f4->Get("T");
  TTree* T5 = (TTree*)f5->Get("T");
  TTree* T6 = (TTree*)f6->Get("T");
 
  EvalRootTTree *evaltree1 = nullptr;
  EvalRootTTree *evaltree2 = nullptr;
  EvalRootTTree *evaltree3 = nullptr;
  EvalRootTTree *evaltree4 = nullptr;
  EvalRootTTree *evaltree5 = nullptr;
  EvalRootTTree *evaltree6 = nullptr;

 


  TH1D *edep_FEMC = new TH1D("edep_FEMC", "Energy Deposition per event in FEMC",200,0,0.5);
  TH1D *edep_CEMC = new TH1D("edep_CEMC", "Energy Deposition per event in CEMC",200,0,1);
  TH1D *edep_EEMC = new TH1D("edep_EEMC", "Energy Deposition per event in EEMC",200,0,1);
  TH1D *edep_FHCAL = new TH1D("edep_FHCAL", "Energy Deposition per event in FHCAL",200,0,10);
  TH1D *edep_HCALIN = new TH1D("edep_HCALIN", "Energy Deposition per event in HCALIN",200,0,1);
  TH1D *edep_HCALOUT = new TH1D("edep_HCALOUT", "Energy Deposition per event in HCALOUT",200,0,5);


  TH2D *te_geta_CEMC = new TH2D("te_geta_CEMC", "Agg tower energy vs generated #eta (CEMC)",200,-4,+4,200,0,10);
  TH2D *te_geta_EEMC = new TH2D("te_geta_EEMC", "Agg tower energy vs generated #eta (EEMC)",200,-4,+4,200,0,10);
  TH2D *te_geta_FEMC = new TH2D("te_geta_FEMC", "Agg tower energy vs generated #eta (FEMC)",200,-4,+4,200,0,10);
  TH2D *te_geta_HCALIN = new TH2D("te_geta_HCALIN", "Agg tower energy vs generated #eta (HCALIN)",200,-4,+4,200,0,10);
  TH2D *te_geta_HCALOUT = new TH2D("te_geta_HCALOUT", "Agg tower energy vs generated #eta (HCALOUT)",200,-4,+4,200,0,10);
  TH2D *te_geta_FHCAL = new TH2D("te_geta_FHCAL", "Agg tower energy vs generated #eta (FHCAL)",200,-4,+4,200,0,10);

  gStyle->SetOptStat(0);
  auto prof1 = new TProfile("prof1","Mean of te_{agg} vs geta (CEMC)",200,-1.5,1.5,"");
  auto prof2 = new TProfile("prof2","Mean of te_{agg} vs geta (HCALIN)",200,-1.5,1.5,"");
  auto prof3 = new TProfile("prof3","Mean of te_{agg} vs geta (HCALOUT)",200,-1.5,1.5,"");
  auto prof4 = new TProfile("prof4","Mean of te_{agg} vs geta (FEMC)",200,1.2,3.4,"");
  auto prof5 = new TProfile("prof5","Mean of te_{agg} vs geta (FHCAL)",200,1.2,3.4,"");
  auto prof6 = new TProfile("prof6","Mean of te_{agg} vs geta (EEMC)",200,-3.6,-1.6,"");
 
  auto prof11 = new TProfile("prof11","Mean of te_{agg} vs gtheta (CEMC)",200,0.5,2.6,"");
  auto prof22 = new TProfile("prof22","Mean of te_{agg} vs gtheta (HCALIN)",200,0.5,2.6,"");
  auto prof33 = new TProfile("prof33","Mean of te_{agg} vs gtheta (HCALOUT)",200,0.5,2.6,"");
  auto prof44 = new TProfile("prof44","Mean of te_{agg} vs gtheta (FEMC)",200,0.,0.6,"");
  auto prof55 = new TProfile("prof55","Mean of te_{agg} vs gtheta (FHCAL)",200,0.,0.6,"");
  auto prof66 = new TProfile("prof66","Mean of te_{agg} vs gtheta (EEMC)",200,2.7,3.5,"");

  T1->SetBranchAddress("DST#EvalTTree_CEMC",&evaltree1);
  T2->SetBranchAddress("DST#EvalTTree_HCALIN",&evaltree2);
  T3->SetBranchAddress("DST#EvalTTree_HCALOUT",&evaltree3);
  T4->SetBranchAddress("DST#EvalTTree_FEMC",&evaltree4);
  T5->SetBranchAddress("DST#EvalTTree_FHCAL",&evaltree5);
  T6->SetBranchAddress("DST#EvalTTree_EEMC",&evaltree6);
  for(int i=0; i<T1->GetEntries(); i++)
    {
      T1->GetEntry(i);
      T2->GetEntry(i);
      T3->GetEntry(i);
      T4->GetEntry(i);
      T5->GetEntry(i);
      T6->GetEntry(i);
      
      float te1 = 0, te2 = 0, te3 = 0, te4 = 0, te5 = 0, te6 = 0, ce1 = 0, ce2 = 0, ce3 = 0, ce4 = 0, ce5 = 0, ce6 = 0, ge, ge1, ge2, ge3, ge4, ge5, ge6, geta1, geta2, geta3, geta4 ,geta5, geta6, gphi1, gphi2, gphi3, gphi4, gphi5, gphi6, cluster1, cluster2, cluster3, cluster4, cluster5, cluster6, ceta1, ceta2, ceta3, ceta4, ceta5, ceta6, cphi1, cphi2, cphi3, cphi4, cphi5, cphi6,  tower1, tower2, tower3, tower4, tower5, tower6, ttheta1, ttheta2, ttheta3, ttheta4, ttheta5, ttheta6, gtheta1, gtheta2, gtheta3, gtheta4, gtheta5, gtheta6, r1, r2, r3, r4,r5,r6, tphi1, tphi2, tphi3, tphi4, tphi5, tphi6;
      
      ge1 = evaltree1->get_ge();
      ge2 = evaltree2->get_ge();
      ge3 = evaltree3->get_ge();
      ge4 = evaltree4->get_ge();
      ge5 = evaltree5->get_ge();     
      ge6 = evaltree6->get_ge();
     
      
      geta1 = evaltree1->get_geta();
      geta2 = evaltree2->get_geta();
      geta3 = evaltree3->get_geta();
      geta4 = evaltree4->get_geta();
      geta5 = evaltree5->get_geta();
      geta6 = evaltree6->get_geta();

      gphi1 = evaltree1->get_gphi();
      gphi2 = evaltree2->get_gphi();
      gphi3 = evaltree3->get_gphi();
      gphi4 = evaltree4->get_gphi();
      gphi5 = evaltree5->get_gphi();
      gphi6 = evaltree6->get_gphi();

      gtheta1 = evaltree1->get_gtheta();
      gtheta2 = evaltree2->get_gtheta();
      gtheta3 = evaltree3->get_gtheta();
      gtheta4 = evaltree4->get_gtheta();
      gtheta5 = evaltree5->get_gtheta();
      gtheta6 = evaltree6->get_gtheta();
      
      //CEMC
     
      if(geta1 >= -0.98 && geta1 <= 0.99) // after checking muon parametrization plots
	{
	  
	  
	  
	
	  

	  for (int i=0; i<evaltree1->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree1->get_tower(i);
	      if (twr)
		{
		  tower1 = twr->get_te();
		  ttheta1 = twr->get_ttheta();
		  tphi1 = twr->get_tphi();


		  //   if(tower1 >=0.2 && geta1 >= -1.5 && geta1 <= 1.2)
		  //  if(tphi1 >= - sqrt(0.4*0.4 - (ttheta1-gtheta1)*(ttheta1-gtheta1)) + gphi1 && tphi1 <= sqrt(0.4*0.4 - (ttheta1-gtheta1)*(ttheta1-gtheta1)) + gphi1 )
		    {
		      te1 = tower1 + te1;
		      
		     
		    }
		}
	    }
	  
	  //  if(te1 > ecut)
	    {
	      edep_CEMC->Fill(te1);
	      te_geta_CEMC->Fill(geta1,te1);
	      prof1->Fill(geta1,te1);
	      prof11->Fill(gtheta1,te1);
	    }

	}
      
      
      //HCALIN
      
      if(geta2 >= -0.98 && geta2 <= 0.99)
	{
	  

	

	  for (int i=0; i<evaltree2->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree2->get_tower(i);
	      if (twr)
		{
		  tower2 = twr->get_te();
		  ttheta2 = twr->get_ttheta();
		  tphi2 = twr->get_tphi();


		  //    if(tower2 >=0.2 && geta2 >= -1.1 && geta2 <= 1.1)
		  //	  if(tphi2 >= - sqrt(0.6*0.6 - (ttheta2-gtheta2)*(ttheta2-gtheta2)) + gphi2 && tphi2 <= sqrt(0.6*0.6 - (ttheta2-gtheta2)*(ttheta2-gtheta2)) + gphi2 )
		    {
		      te2 = tower2 + te2;
		      
		    
		    }
		}
	    }
	  //	  if(te2>ecut)
	    {
	      edep_HCALIN->Fill(te2);
	      te_geta_HCALIN->Fill(geta2,te2);
	      prof2->Fill(geta2,te2);
	      prof22->Fill(gtheta2,te2);
	    }
	  
	}
      
      
      //HCALOUT
      
      if(geta3 >= -0.98 && geta3 <= 0.99)
	{
	  

	  
	


	  for (int i=0; i<evaltree3->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree3->get_tower(i);
	      if (twr)
		{
		  tower3 = twr->get_te();
		  ttheta3 = twr->get_ttheta();
		  tphi3 = twr->get_tphi();


		  //     if(tower3 >=0.2 && geta3 >= -1.1 && geta3 <= 1.1)
		  //	  if(tphi3 >= - sqrt(0.6*0.6 - (ttheta3-gtheta3)*(ttheta3-gtheta3)) + gphi3 && tphi3 <= sqrt(0.6*0.6 - (ttheta3-gtheta3)*(ttheta3-gtheta3)) + gphi3 )
		    {
		      te3 = tower3 + te3;
		      
		    }
		}
	    }
	  //	  if(te3 > ecut)
	    {
	      edep_HCALOUT->Fill(te3);
	      te_geta_HCALOUT->Fill(geta3,te3);
	      prof3->Fill(geta3,te3);
	      prof33->Fill(gtheta3,te3);

	    }
	} 
      
      //FEMC
      
      if(geta4 >= 1.32 && geta4 <= 3.14)
	{
	  

	


	  for (int i=0; i<evaltree4->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree4->get_tower(i);
	      if (twr)
		{
		  tower4 = twr->get_te();
		  ttheta4 = twr->get_ttheta();
		  tphi4 = twr->get_tphi();


		  //    if(tower4 >=0.2 && geta4 >= 1.3 && geta4 <= 3.3)
		  //  if(tphi4 >= - sqrt(0.3*0.3 - (ttheta4-gtheta4)*(ttheta4-gtheta4)) + gphi4 && tphi4 <= sqrt(0.3*0.3 - (ttheta4-gtheta4)*(ttheta4-gtheta4)) + gphi4 )
		    { 
		      te4 = tower4 + te4;
		  
		     

		    }
		}
	    }
	  //  if(te4 > ecut)
	    {
	      //   te4 = te4/(0.250805 - 0.00396617*ge4 + 0.0000331404*ge4*ge4);//energy calibration pol2
	      //    te4 = te4/(0.239815 + 0.00000220243*ge4 - 0.000282063*ge4*ge4 + 0.00000681101*ge4*ge4*ge4);//energy calibration pol3
	      edep_FEMC->Fill(te4);
	      te_geta_FEMC->Fill(geta4,te4);
	      prof4->Fill(geta4,te4);
	      prof44->Fill(gtheta4,te4);

	    }

	}      
      
      
      
      //FHCAL
      if(geta5 >= 1.32 && geta5 <= 3.14)//matched with femc
	{
	  





	  for (int i=0; i<evaltree5->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree5->get_tower(i);
	      if (twr)
		{
		  tower5 = twr->get_te();
		  ttheta5 = twr->get_ttheta();
		  tphi5 = twr->get_tphi();


		  //   if(tower5 >=0.2 && geta5 >= 1.3 && geta5 <= 3.3)//matched with femc
		  //  if(tphi5 >= - sqrt(0.4*0.4 - (ttheta5-gtheta5)*(ttheta5-gtheta5)) + gphi5 && tphi5 <= sqrt(0.4*0.4 - (ttheta5-gtheta5)*(ttheta5-gtheta5)) + gphi5 )
		    {
		      te5 = tower5 + te5;
		      
		    

		    }
		}
	    }
	  //	  if(te5>ecut)
	    {
	      //    te5 = te5/(0.230773 + 0.00100258*ge5 + 0.0000124815*ge5*ge5);//energy calibration pol2
	      //   te5 = te5/(0.237499 - 0.00119102*ge5 + 0.000180067*ge5*ge5 - 0.00000354526*ge5*ge5*ge5);//energy calibration pol3
	      
	      edep_FHCAL->Fill(te5);
	      te_geta_FHCAL->Fill(geta5,te5);
	      prof5->Fill(geta5,te5);
	      prof55->Fill(gtheta5,te5);

	    }
	}
      
      
      
      //EEMC
      
      if(geta5 >= -3.4 && geta5 <= -1.77)
	{
	  
	  for (int i=0; i<evaltree6->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree6->get_tower(i);
	      if (twr)
		{
		  tower6 = twr->get_te();
		  ttheta6 = twr->get_ttheta();
		  tphi6 = twr->get_tphi();


		  //	      if(tower6 >=0.2 && geta5 >= -3.5 && geta5 <= -1.7)
		  te6 = tower6 + te6;


		  
		}
	    }
	  //    if(te6 >= 0.2)
	  //{
	  edep_EEMC->Fill(te6);
	  te_geta_EEMC->Fill(geta6,te6);
	  prof6->Fill(geta6,te6);
	  prof66->Fill(gtheta6,te6);

	  //	}
	  
	  
	
	
     
       
       }
      
 
    

  
    }

 

  cout<<"CEMC eta profile fit parameters:"<<endl;
  TF1 *fit1 = new TF1("fit1","[0]*cosh([1]*x)",-1.1,1.1);
  //  fit1->SetRange(-1.1,1.1);
  fit1->SetParameters(5,10);
  prof1->Fit("fit1","R");  //CEMC
  prof1->GetXaxis()->SetTitle("geta");
  prof1->GetYaxis()->SetTitle("Mean of te_{agg}");
  cout<<"CEMC theta profile fit parameters:"<<endl;
  TF1 *fit11 = new TF1("fit11","pol4",0.65,2.5);
  //  TF1 *fit11 = new TF1("fit11","[0]*cosh([1]*x)",0.7,2.5);
  //  fit11->SetParameters(0.1,1);
  prof11->Fit("fit11","R");  
  prof11->GetXaxis()->SetTitle("gtheta");
  prof11->GetYaxis()->SetTitle("Mean of te_{agg}");
 


  cout<<"HCALIN eta profile fit parameters:"<<endl;
  TF1 *fit2 = new TF1("fit2","[0]*cosh([1]*x)",-1.1,1.1);
  fit2->SetParameters(5,10);
  prof2->Fit("fit2","R");  //HCALIN
  prof2->GetXaxis()->SetTitle("geta");
  prof2->GetYaxis()->SetTitle("Mean of te_{agg}");
  cout<<"HCALIN theta profile fit parameters:"<<endl;
  TF1 *fit22 = new TF1("fit22","pol3",0.65,2.5);
  //  TF1 *fit22 = new TF1("fit22","[0]*cosh([1]*x)",0.7,2.4);
  // fit22->SetParameters(0.1,1);
  prof22->Fit("fit22","R"); 
  prof22->GetXaxis()->SetTitle("gtheta");
  prof22->GetYaxis()->SetTitle("Mean of te_{agg}");

  cout<<"HCALOUT eta profile fit parameters:"<<endl;
  TF1 *fit3a = new TF1("fit3a","pol4",-0.7,0.7);
  prof3->Fit("fit3a","R");  //HCALOUT
  TF1 *fit3b = new TF1("fit3b","pol4",-1.,-0.72);
  prof3->Fit("fit3b","R+");  //HCALOUT
  TF1 *fit3c = new TF1("fit3c","pol4",0.72,1.);
  prof3->Fit("fit3c","R+");  //HCALOUT
  prof3->GetXaxis()->SetTitle("geta");
  prof3->GetYaxis()->SetTitle("Mean of te_{agg}");
  TF1 *total = new TF1("total","pol4(0)+pol4(5)+pol4(10)",-1.,1.); 
  Double_t par[15];
  fit3a->GetParameters(&par[0]);
  fit3b->GetParameters(&par[5]);
  fit3c->GetParameters(&par[10]);
  total->SetParameters(par);
  prof3->Fit("total","R+");
  cout<<"HCALOUT theta profile fit parameters:"<<endl;
  TF1 *fit33a = new TF1("fit33a","pol4",0.9,2.25);
  prof33->Fit("fit33a","R");  
  TF1 *fit33b = new TF1("fit33b","pol4",0.69,0.89);
  prof33->Fit("fit33b","R+");  
  TF1 *fit33c = new TF1("fit33c","pol4",2.24,2.45);
  prof33->Fit("fit33c","R+");  
  prof33->GetXaxis()->SetTitle("gtheta");
  prof33->GetYaxis()->SetTitle("Mean of te_{agg}");
  TF1 *total1 = new TF1("total1","pol4(0)+pol4(5)+pol4(10)",0.69,2.45); 
  Double_t par1[15];
  fit33a->GetParameters(&par1[0]);
  fit33b->GetParameters(&par1[5]);
  fit33c->GetParameters(&par1[10]);
  total1->SetParameters(par1);
  prof33->Fit("total1","R+");
 
  cout<<"FEMC eta profile fit parameters:"<<endl;
  TF1 *fit4 = new TF1("fit4","pol4",1.35,3.15);
  prof4->Fit("fit4","R");  //FEMC
  prof4->GetXaxis()->SetTitle("geta");
  prof4->GetYaxis()->SetTitle("Mean of te_{agg}");
  cout<<"FEMC theta profile fit parameters:"<<endl;
  TF1 *fit44 = new TF1("fit44","pol4",0.08,0.5);
  prof44->Fit("fit44","R");  
  prof44->GetXaxis()->SetTitle("gtheta");
  prof44->GetYaxis()->SetTitle("Mean of te_{agg}");

  cout<<"FHCAL eta profile fit parameters:"<<endl;
  TF1 *fit5 = new TF1("fit5","pol4",1.3,3.3);
  prof5->Fit("fit5","R");  //FHCAL
  prof5->GetXaxis()->SetTitle("geta");
  prof5->GetYaxis()->SetTitle("Mean of te_{agg}");
  cout<<"FHCAL theta profile fit parameters:"<<endl;
  TF1 *fit55 = new TF1("fit55","pol4",0.1,0.51);
  prof55->Fit("fit55","R");  
  prof55->GetXaxis()->SetTitle("gtheta");
  prof55->GetYaxis()->SetTitle("Mean of te_{agg}");

  cout<<"EEMC eta profile fit parameters:"<<endl;
  TF1 *fit6 = new TF1("fit6","pol1",-3.4,-1.8);
  prof6->Fit("fit6","R");  //EEMC
  prof6->GetXaxis()->SetTitle("geta");
  prof6->GetYaxis()->SetTitle("Mean of te_{agg}");
  cout<<"EEMC theta profile fit parameters:"<<endl;
  TF1 *fit66 = new TF1("fit66","pol1",2.8,3.07);
  prof66->Fit("fit66","R");  
  prof66->GetXaxis()->SetTitle("gtheta");
  prof66->GetYaxis()->SetTitle("Mean of te_{agg}");
  


 


  TFile *f = new TFile("Muon_calorimeter_plots.root","RECREATE");
  /*  f->GetList()->Add(edep_FEMC);
  f->GetList()->Add(edep_CEMC);
  f->GetList()->Add(edep_EEMC);
  f->GetList()->Add(edep_FHCAL);
  f->GetList()->Add(edep_HCALIN);
  f->GetList()->Add(edep_HCALOUT);
  f->GetList()->Add(te_geta_FEMC);
  f->GetList()->Add(te_geta_CEMC);
  f->GetList()->Add(te_geta_EEMC);
  f->GetList()->Add(te_geta_FHCAL);
  f->GetList()->Add(te_geta_HCALIN);
  f->GetList()->Add(te_geta_HCALOUT);*/
  f->GetList()->Add(prof1);
  f->GetList()->Add(prof11);
  f->GetList()->Add(prof2);
  f->GetList()->Add(prof22);
  f->GetList()->Add(prof3);
  f->GetList()->Add(prof33);
  f->GetList()->Add(prof4);
  f->GetList()->Add(prof44);
  f->GetList()->Add(prof5);
  f->GetList()->Add(prof55);
  f->GetList()->Add(prof6);
  f->GetList()->Add(prof66);

  f->Write();

}



