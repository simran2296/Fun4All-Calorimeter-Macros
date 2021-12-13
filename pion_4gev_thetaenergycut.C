//note: only ecal eta cuts 

#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"


R__LOAD_LIBRARY(libeicqa_modules.so)

void eres_pion_newcut()
{
  TFile *f1 = new TFile("/gpfs02/eic/iiti/July21/pion_4gev_uni/merged_Eval_CEMC.root","READ");
  TFile *f2 = new TFile("/gpfs02/eic/iiti/July21/pion_4gev_uni/merged_Eval_HCALIN.root","READ");
  TFile *f3 = new TFile("/gpfs02/eic/iiti/July21/pion_4gev_uni/merged_Eval_HCALOUT.root","READ");
  TFile *f4 = new TFile("/gpfs02/eic/iiti/July21/pion_4gev_uni/merged_Eval_FEMC.root","READ");
  TFile *f5 = new TFile("/gpfs02/eic/iiti/July21/pion_4gev_uni/merged_Eval_FHCAL.root","READ");
  TFile *f6 = new TFile("/gpfs02/eic/iiti/July21/pion_4gev_uni/merged_Eval_EEMC.root","READ");

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
      //   cout<<"ge1: "<<ge1<<endl;
      //   gen_energy_CEMC->Fill(ge1);
      
      ge2 = evaltree2->get_ge();
      //   cout<<"ge2: "<<ge1<<endl;
      //    gen_energy_HCALIN->Fill(ge2);
      
      ge3 = evaltree3->get_ge();
      //    cout<<"ge3: "<<ge1<<endl;
      //    gen_energy_HCALOUT->Fill(ge3);

      ge4 = evaltree4->get_ge();
      //    cout<<"ge4: "<<ge1<<endl;
      //    gen_energy_FEMC->Fill(ge4);
      
      ge5 = evaltree5->get_ge();
      //   cout<<"ge5: "<<ge1<<endl;
      //    gen_energy_FHCAL->Fill(ge5);
      
      ge6 = evaltree6->get_ge();
      //   cout<<"ge6: "<<ge1<<endl;
      //   gen_energy_EEMC->Fill(ge6);
      
      
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
     
      if(geta1 >= -0.98 && geta1 <= 0.99) // set to hcal limit
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
	  if(te1 > 0.850493 - 1.3569*gtheta1 + 1.10232*pow(gtheta1,2) - 0.430019*pow(gtheta1,3) + 0.0692608*pow(gtheta1,4) ) //mip cut
	    {
	      edep_CEMC->Fill(te1);
	      te_geta_CEMC->Fill(geta1,te1);
	    }

	}
      
      
      //HCALIN
      
      if(geta2 >= -1.1 && geta2 <= 1.1)
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

	    }
	  
	}
      
      
      //HCALOUT
      
      if(geta3 >= -1.1 && geta3 <= 1.1)
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
	  if(te4 > 0.0433039 - 0.0591906*gtheta4 + 0.33481*pow(gtheta4,2) - 0.33481*pow(gtheta4,3) + 0.444401*pow(gtheta4,4) ) //mip cut

	    {
	      //   te4 = te4/(0.250805 - 0.00396617*ge4 + 0.0000331404*ge4*ge4);//energy calibration pol2
	      //    te4 = te4/(0.239815 + 0.00000220243*ge4 - 0.000282063*ge4*ge4 + 0.00000681101*ge4*ge4*ge4);//energy calibration pol3
	      edep_FEMC->Fill(te4);
	      te_geta_FEMC->Fill(geta4,te4);

	    }

	}      
      
      
      
      //FHCAL
      if(geta5 >= 1.3 && geta5 <= 3.3)//matched with femc
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
	  if(te6 > 0.36477 - 0.046759*gtheta6)
	    {
	      edep_EEMC->Fill(te6);
	      te_geta_EEMC->Fill(geta6,te6);

	    }
	  
	  
	
	
     
       
       }
      
 
    

  
    }





 


  TFile *f = new TFile("Pion_calorimeter_plots.root","RECREATE");
  f->GetList()->Add(edep_FEMC);
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
  f->GetList()->Add(te_geta_HCALOUT);

  f->Write();

}



