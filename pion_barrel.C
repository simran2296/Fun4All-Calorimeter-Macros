#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"


R__LOAD_LIBRARY(libeicqa_modules.so)

void eres_barrel()
{
  TFile *f1 = new TFile("Eval_CEMC.root","READ");
  TFile *f2 = new TFile("Eval_HCALIN.root","READ");
  TFile *f3 = new TFile("Eval_HCALOUT.root","READ");
 
  TTree* T1 = (TTree*)f1->Get("T");
  TTree* T2 = (TTree*)f2->Get("T");
  TTree* T3 = (TTree*)f3->Get("T"); 
 
  EvalRootTTree *evaltree1 = nullptr;
  EvalRootTTree *evaltree2 = nullptr;
  EvalRootTTree *evaltree3 = nullptr; 

  float ecut = 0.1;
  float ecut1 = 0.2;

  TH2D *hist_energy_CEMC_1 = new TH2D("hist_energy_CEMC_1","#frac{te-ge}{ge} vs ge (CEMC); ge; #frac{te-ge}{ge}",200,0,30,200,-2.6,2.6);
  TH2D *hist_energy_HCALIN_1 = new TH2D("hist_energy_HCALIN_1","#frac{te-ge}{ge} vs ge (HCALIN); ge; #frac{te-ge}{ge}",200,0,30,200,-2.6,2.6);
  TH2D *hist_energy_HCALOUT_1 = new TH2D("hist_energy_HCALOUT_1","#frac{te-ge}{ge} vs ge (HCALOUT); ge; #frac{te-ge}{ge}",200,0,30,200,-2.6,2.6);

  T1->SetBranchAddress("DST#EvalTTree_CEMC",&evaltree1);
  T2->SetBranchAddress("DST#EvalTTree_HCALIN",&evaltree2);
  T3->SetBranchAddress("DST#EvalTTree_HCALOUT",&evaltree3);

  TFile *f = new TFile("tree1.root","RECREATE");

  TTree *t = new TTree("t","resolution tree");
  Float_t truth_e, te_CEMC, te_HCALIN, te_HCALOUT;

  t->Branch("truth_e",&truth_e,"truth_e/F");
  t->Branch("te_CEMC",&te_CEMC,"te_CEMC/F");
  t->Branch("te_HCALIN",&te_HCALIN,"te_HCALIN/F");
  t->Branch("te_HCALOUT",&te_HCALOUT,"te_HCALOUT/F");
 
  for(int i=0; i<T1->GetEntries(); i++)
    {
      T1->GetEntry(i);
      T2->GetEntry(i);
      T3->GetEntry(i);
            
      float te1 = 0, te2 = 0, te3 = 0, te4 = 0, te5 = 0, te6 = 0, ge, ge1, ge2, ge3, ge4, ge5, ge6, geta1, geta2, geta3, geta4 ,geta5, geta6, gphi1, gphi2, gphi3, gphi4, gphi5, gphi6, tower1, tower2, tower3, tower4, tower5, tower6, ttheta1, ttheta2, ttheta3, ttheta4, ttheta5, ttheta6, gtheta1, gtheta2, gtheta3, gtheta4, gtheta5, gtheta6, tphi1, tphi2, tphi3, tphi4, tphi5, tphi6;
      
      // initialize with unphysical values
      te_CEMC=-1;
      te_HCALIN=-1;
      te_HCALOUT=-1;

      ge1 = evaltree1->get_ge();
      ge2 = evaltree2->get_ge();
      ge3 = evaltree3->get_ge();

      truth_e = ge1;

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
     
      if(geta1 >= -1.1 && geta1 <= 1.1) // set to hcal limit
	{
	  for (int i=0; i<evaltree1->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree1->get_tower(i);
	      if (twr)
		{
		  tower1 = twr->get_te();
		  ttheta1 = twr->get_ttheta();
		  tphi1 = twr->get_tphi();


		  if(tphi1 >= - sqrt(0.4*0.4 - (ttheta1-gtheta1)*(ttheta1-gtheta1)) + gphi1 && tphi1 <= sqrt(0.4*0.4 - (ttheta1-gtheta1)*(ttheta1-gtheta1)) + gphi1 )//circular cut
		    {
		      if(tower1 > ecut1)
			{
			  te1 = tower1 + te1;
			}
		    }
		}
	    }
	  
	  if(te1 > 0.0)
	    {
	        te_CEMC=te1;
	      float frac_CEMC_1 = (te1-ge1)/ge1;
	      hist_energy_CEMC_1->Fill(ge1, frac_CEMC_1);
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


		  if(tphi2 >= - sqrt(0.6*0.6 - (ttheta2-gtheta2)*(ttheta2-gtheta2)) + gphi2 && tphi2 <= sqrt(0.6*0.6 - (ttheta2-gtheta2)*(ttheta2-gtheta2)) + gphi2 )
		    {
		      te2 = tower2 + te2;
		      
		    }
		}
	    }
	  if(te2>ecut)
	    {
	
	      te_HCALIN=te2;

	      float frac_HCALIN_1 = (te2-ge2)/ge2;
	      hist_energy_HCALIN_1->Fill(ge2, frac_HCALIN_1);

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

		  if(tphi3 >= - sqrt(0.6*0.6 - (ttheta3-gtheta3)*(ttheta3-gtheta3)) + gphi3 && tphi3 <= sqrt(0.6*0.6 - (ttheta3-gtheta3)*(ttheta3-gtheta3)) + gphi3 )
		    {
		      te3 = tower3 + te3;
		      
		    }
		}
	    }
	  if(te3 > ecut)
	    {
	        te_HCALOUT=te3;
	      float frac_HCALOUT_1 = (te3-ge3)/ge3;
	      hist_energy_HCALOUT_1->Fill(ge3, frac_HCALOUT_1);

	    }
	} 
      


  


      // only fill the tree if we recovered something in at least one detector after cuts      
      if ( te_CEMC>-1 ||  te_HCALIN>-1 || te_HCALOUT>-1)
	{
	  t->Fill();
	}
      
    }

  t->Write();

  f->GetList()->Add(hist_energy_CEMC_1);
  f->GetList()->Add(hist_energy_HCALIN_1);
  f->GetList()->Add(hist_energy_HCALOUT_1);

  f->Write();



}



