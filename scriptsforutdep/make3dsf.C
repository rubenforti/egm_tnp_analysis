void make3dsf() {
	std::string step("isonotrig");
	std::string Step, Title;
	if (step==std::string("triggerplus")) {
		Step=std::string("mu_trigger_plus");
		Title=std::string("TriggerPlus");
	}
	else if (step==std::string("triggerminus")) {
		Step=std::string("mu_trigger_minus");
		Title=std::string("TriggerMinus");
	}
	else if (step==std::string("isolation")) {
		Step=std::string("mu_iso_both");
		Title=std::string("Isolation");
	}
	else if (step==std::string("isonotrig")) {
		Step=std::string("mu_isonotrig_both");
		Title=std::string("IsolationNoTrigger");
	}
	std::string Filename("3dsf");
	Filename+=step+std::string(".root");
	TFile *output=new TFile(Filename.c_str(),"RECREATE");
	double ptarray[16] = {24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.};
	double utarray[18] = {-3000000000,-30,-15,-10,-5,0,5,10,15,30,40,50,60,70,80,90,100,30000000000};
	double etaarray[48];
	for (unsigned int i=0; i!=49; i++) {
		etaarray[i] = -2.4+i*0.1;
	}
	TH3D *histo=new TH3D(Title.c_str(),"",48,etaarray,15,ptarray,17,utarray);
	std::vector<TH2D*> sfs, effdata, effmc;
	for (unsigned int i=1; i!=18; i++) {
		std::string filename("/scratch/bruschin/cmsasymow/2018isotrig/egm_tnp_analysis/bin");
		filename+=std::to_string(i)+std::string("/efficiencies_GtoH/")+Step+std::string("/allEfficiencies_2D.root");
		TFile *file=new TFile(filename.c_str());
		TH2D* Histo=(TH2D*)file->Get("SF2D_nominal")->Clone((std::string("SF2D_nominal")+std::to_string(i)).c_str());
		TH2D* Histo2=(TH2D*)file->Get("EffData2D")->Clone((std::string("EffData2D")+std::to_string(i)).c_str());
		sfs.push_back(Histo);
		effdata.push_back(Histo2);
		//file->Close();
	}
	for (unsigned int i=1; i!=18; i++) {
		std::string filename("/scratch/bruschin/cmsasymow/2018isotrig/Steve_Marc_Raj/");
		filename+=step+std::string("mc_")+std::to_string(i)+std::string(".root");
		TFile *file=new TFile(filename.c_str());
		TH3D* Histo1=(TH3D*)file->Get("pass_mu_DY_postVFP")->Clone((std::string("pass_mu_DY_postVFP")+std::to_string(i)).c_str());
		TH3D* Histo2=(TH3D*)file->Get("fail_mu_DY_postVFP")->Clone((std::string("fail_mu_DY_postVFP")+std::to_string(i)).c_str());
		TH2D* Histo=new TH2D((std::string("EffMC2D")+std::to_string(i)).c_str(),"",48,etaarray,15,ptarray);
		for (unsigned int I=0; I!=48; I++) {
			for (unsigned int J=0; J!=15; J++) {
				double pass=0, fail=0;
				for (unsigned int H=0; H!=Histo1->GetXaxis()->GetNbins(); H++) {
					pass+=Histo1->GetBinContent(H+1,J+1,I+1);
					fail+=Histo2->GetBinContent(H+1,J+1,I+1);
				}
				pass=pass/(pass+fail);
				Histo->SetBinContent(I+1,J+1,pass);
			}
		}
		effmc.push_back(Histo);
		//delete Histo1;
		//delete Histo2;
		//file->Close();
	}
	for (unsigned int h=0; h!=17; h++) {
		for (unsigned int i=0; i!=48; i++) {
			for (unsigned int j=0; j!=15; j++) {
				histo->SetBinContent(i+1,j+1,h+1,effdata[h]->GetBinContent(i+1,j+1)/effmc[h]->GetBinContent(i+1,j+1));
				histo->SetBinError(i+1,j+1,h+1,sfs[h]->GetBinError(i+1,j+1));
			}
		}
	}
	output->cd();
	histo->Write();
	output->Close();
}
