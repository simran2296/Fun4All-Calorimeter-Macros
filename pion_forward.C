#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"


R__LOAD_LIBRARY(libeicqa_modules.so)

void eres_forward()
{
  TFile *f4 = new TFile("Eval_FEMC.root","READ");
  TFile *f5 = new TFile("Eval_FHCAL.root","READ");
 
  TTree* T4 = (TTree*)f4->Get("T");
  TTree* T5 = (TTree*)f5->Get("T");
 
  EvalRootTTree *evaltree4 = nullptr;
  EvalRootTTree *evaltree5 = nullptr;
 

  float ecut = 0.1;
  float ecut1 = 0.2;

  TH1D *edep_FEMC = new TH1D("edep_FEMC", "Energy Deposition per event in FEMC",200,0,0.5);
  TH1D *edep_FHCAL = new TH1D("edep_FHCAL", "Energy Deposition per event in FHCAL",200,0,0.5);

  TH2D *hist_energy_FEMC_1 = new TH2D("hist_energy_FEMC_1","#frac{te-ge}{ge} vs ge (FEMC); ge; #frac{te-ge}{ge}",200,0,30,200,-1.5,1.5);
  TH2D *hist_energy_FHCAL_1 = new TH2D("hist_energy_FHCAL_1","#frac{te-ge}{ge} vs ge (FHCAL); ge; #frac{te-ge}{ge}",200,0,30,200,-1.5,1.5);

  T4->SetBranchAddress("DST#EvalTTree_FEMC",&evaltree4);
  T5->SetBranchAddress("DST#EvalTTree_FHCAL",&evaltree5);


  TFile *f = new TFile("tree2.root","RECREATE");

  TTree *t = new TTree("t","resolution tree");
  Float_t truth_e, te_FEMC, te_FHCAL;

  t->Branch("truth_e",&truth_e,"truth_e/F");

  t->Branch("te_FEMC",&te_FEMC,"te_FEMC/F");
  t->Branch("te_FHCAL",&te_FHCAL,"te_FHCAL/F");
 
  for(int i=0; i<T4->GetEntries(); i++)
    {
      T4->GetEntry(i);
      T5->GetEntry(i);
      
      float te1 = 0, te2 = 0, te3 = 0, te4 = 0, te5 = 0, te6 = 0, ge, ge1, ge2, ge3, ge4, ge5, ge6, geta1, geta2, geta3, geta4 ,geta5, geta6, gphi1, gphi2, gphi3, gphi4, gphi5, gphi6, tower1, tower2, tower3, tower4, tower5, tower6, ttheta1, ttheta2, ttheta3, ttheta4, ttheta5, ttheta6, gtheta1, gtheta2, gtheta3, gtheta4, gtheta5, gtheta6, tphi1, tphi2, tphi3, tphi4, tphi5, tphi6;
      
      // initialize with unphysical values
      te_FEMC=-1;
      te_FHCAL=-1;

      ge4 = evaltree4->get_ge();
      ge5 = evaltree5->get_ge();

      truth_e = ge4;

      geta4 = evaltree4->get_geta();
      geta5 = evaltree5->get_geta();
 
      gphi4 = evaltree4->get_gphi();
      gphi5 = evaltree5->get_gphi();

      gtheta4 = evaltree4->get_gtheta();
      gtheta5 = evaltree5->get_gtheta();

      


      //FEMC
      
      if(geta4 >= 1.3 && geta4 <= 3.3)
	{

	  for (int i=0; i<evaltree4->get_ntowers(); i++)
	    {
	      EvalTower *twr = evaltree4->get_tower(i);
	      if (twr)
		{
		  tower4 = twr->get_te();
		  ttheta4 = twr->get_ttheta();
		  tphi4 = twr->get_tphi();

		  //  if(tphi4 >= - sqrt(0.25*0.25 - (ttheta4-gtheta4)*(ttheta4-gtheta4)) + gphi4 && tphi4 <= sqrt(0.25*0.25 - (ttheta4-gtheta4)*(ttheta4-gtheta4)) + gphi4 )
		    { 
		      //  if(tower4 > 0.4)
			{
			  te4 = tower4 + te4;
		  	}	  
		    }
		}
	    }
	  //  if(te4 > 0.1)
	    {
	      edep_FEMC->Fill(te4);
	      te_FEMC = te4;

	      float frac_FEMC_1 = (te4-ge4)/ge4;
	      hist_energy_FEMC_1->Fill(ge4, frac_FEMC_1);
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

		  if(tphi5 >= - sqrt(0.4*0.4 - (ttheta5-gtheta5)*(ttheta5-gtheta5)) + gphi5 && tphi5 <= sqrt(0.4*0.4 - (ttheta5-gtheta5)*(ttheta5-gtheta5)) + gphi5 )
		    {
		      te5 = tower5 + te5;
		     		      
		    }
		}
	    }
	  if(te5>ecut)
	    {
	      
	      edep_FHCAL->Fill(te5);

	      te_FHCAL = te5;

	      float frac_FHCAL_1 = (te5-ge5)/ge5;
	      hist_energy_FHCAL_1->Fill(ge5, frac_FHCAL_1);
	    }
	}
      
      
      


      // only fill the tree if we recovered something in at least one detector after cuts      
      if (  te_FEMC>-1 || te_FHCAL>-1 )
	{
	  t->Fill();
	}
    
    }

  t->Write();

  f->GetList()->Add(hist_energy_FEMC_1);
  f->GetList()->Add(edep_FEMC);
  f->GetList()->Add(hist_energy_FHCAL_1);
  f->GetList()->Add(edep_FHCAL);

  f->Write();

}



