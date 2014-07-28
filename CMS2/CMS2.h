// -*- C++ -*-
#ifndef CMS2_H
#define CMS2_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector> 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std; 
class CMS2 {
private: 
protected: 
	unsigned int index;
	float	pfmet_;
	TBranch *pfmet_branch;
	bool pfmet_isLoaded;
	float	metphi_;
	TBranch *metphi_branch;
	bool metphi_isLoaded;
	float	pfmet_type1cor_;
	TBranch *pfmet_type1cor_branch;
	bool pfmet_type1cor_isLoaded;
	float	gen_met_;
	TBranch *gen_met_branch;
	bool gen_met_isLoaded;
	float	evt_scale1fb_;
	TBranch *evt_scale1fb_branch;
	bool evt_scale1fb_isLoaded;
	int	evt_isRealData_;
	TBranch *evt_isRealData_branch;
	bool evt_isRealData_isLoaded;
	unsigned int	evt_event_;
	TBranch *evt_event_branch;
	bool evt_event_isLoaded;
	unsigned int	evt_lumiBlock_;
	TBranch *evt_lumiBlock_branch;
	bool evt_lumiBlock_isLoaded;
	unsigned int	evt_run_;
	TBranch *evt_run_branch;
	bool evt_run_isLoaded;
	unsigned int	evt_nvtxs_;
	TBranch *evt_nvtxs_branch;
	bool evt_nvtxs_isLoaded;
	bool	dielectronTrigger_;
	TBranch *dielectronTrigger_branch;
	bool dielectronTrigger_isLoaded;
	bool	dimuonTrigger_;
	TBranch *dimuonTrigger_branch;
	bool dimuonTrigger_isLoaded;
	bool	electronmuonTrigger_;
	TBranch *electronmuonTrigger_branch;
	bool electronmuonTrigger_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *els_p4_;
	TBranch *els_p4_branch;
	bool els_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *mus_p4_;
	TBranch *mus_p4_branch;
	bool mus_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *pfjets_p4_;
	TBranch *pfjets_p4_branch;
	bool pfjets_p4_isLoaded;
	vector<int> *genps_id_;
	TBranch *genps_id_branch;
	bool genps_id_isLoaded;
	vector<int> *genps_id_mother_;
	TBranch *genps_id_mother_branch;
	bool genps_id_mother_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *genps_p4_;
	TBranch *genps_p4_branch;
	bool genps_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *genjets_p4_;
	TBranch *genjets_p4_branch;
	bool genjets_p4_isLoaded;
	vector<bool> *looseEl_;
	TBranch *looseEl_branch;
	bool looseEl_isLoaded;
	vector<bool> *looseMu_;
	TBranch *looseMu_branch;
	bool looseMu_isLoaded;
	vector<bool> *tightEl_;
	TBranch *tightEl_branch;
	bool tightEl_isLoaded;
	vector<bool> *tightMu_;
	TBranch *tightMu_branch;
	bool tightMu_isLoaded;
	vector<bool> *passesLoosePFJetID_;
	TBranch *passesLoosePFJetID_branch;
	bool passesLoosePFJetID_isLoaded;
	vector<float> *pfjets_corL1FastL2L3_;
	TBranch *pfjets_corL1FastL2L3_branch;
	bool pfjets_corL1FastL2L3_isLoaded;
	vector<float> *pfjets_combinedSecondaryVertexBJetTag_;
	TBranch *pfjets_combinedSecondaryVertexBJetTag_branch;
	bool pfjets_combinedSecondaryVertexBJetTag_isLoaded;
public: 
void Init(TTree *tree) {
	els_p4_branch = 0;
	if (tree->GetBranch("els_p4") != 0) {
		els_p4_branch = tree->GetBranch("els_p4");
		els_p4_branch->SetAddress(&els_p4_);
	}
	mus_p4_branch = 0;
	if (tree->GetBranch("mus_p4") != 0) {
		mus_p4_branch = tree->GetBranch("mus_p4");
		mus_p4_branch->SetAddress(&mus_p4_);
	}
	pfjets_p4_branch = 0;
	if (tree->GetBranch("pfjets_p4") != 0) {
		pfjets_p4_branch = tree->GetBranch("pfjets_p4");
		pfjets_p4_branch->SetAddress(&pfjets_p4_);
	}
	genps_p4_branch = 0;
	if (tree->GetBranch("genps_p4") != 0) {
		genps_p4_branch = tree->GetBranch("genps_p4");
		genps_p4_branch->SetAddress(&genps_p4_);
	}
	genjets_p4_branch = 0;
	if (tree->GetBranch("genjets_p4") != 0) {
		genjets_p4_branch = tree->GetBranch("genjets_p4");
		genjets_p4_branch->SetAddress(&genjets_p4_);
	}
  tree->SetMakeClass(1);
	pfmet_branch = 0;
	if (tree->GetBranch("pfmet") != 0) {
		pfmet_branch = tree->GetBranch("pfmet");
		pfmet_branch->SetAddress(&pfmet_);
	}
	metphi_branch = 0;
	if (tree->GetBranch("metphi") != 0) {
		metphi_branch = tree->GetBranch("metphi");
		metphi_branch->SetAddress(&metphi_);
	}
	pfmet_type1cor_branch = 0;
	if (tree->GetBranch("pfmet_type1cor") != 0) {
		pfmet_type1cor_branch = tree->GetBranch("pfmet_type1cor");
		pfmet_type1cor_branch->SetAddress(&pfmet_type1cor_);
	}
	gen_met_branch = 0;
	if (tree->GetBranch("gen_met") != 0) {
		gen_met_branch = tree->GetBranch("gen_met");
		gen_met_branch->SetAddress(&gen_met_);
	}
	evt_scale1fb_branch = 0;
	if (tree->GetBranch("evt_scale1fb") != 0) {
		evt_scale1fb_branch = tree->GetBranch("evt_scale1fb");
		evt_scale1fb_branch->SetAddress(&evt_scale1fb_);
	}
	evt_isRealData_branch = 0;
	if (tree->GetBranch("evt_isRealData") != 0) {
		evt_isRealData_branch = tree->GetBranch("evt_isRealData");
		evt_isRealData_branch->SetAddress(&evt_isRealData_);
	}
	evt_event_branch = 0;
	if (tree->GetBranch("evt_event") != 0) {
		evt_event_branch = tree->GetBranch("evt_event");
		evt_event_branch->SetAddress(&evt_event_);
	}
	evt_lumiBlock_branch = 0;
	if (tree->GetBranch("evt_lumiBlock") != 0) {
		evt_lumiBlock_branch = tree->GetBranch("evt_lumiBlock");
		evt_lumiBlock_branch->SetAddress(&evt_lumiBlock_);
	}
	evt_run_branch = 0;
	if (tree->GetBranch("evt_run") != 0) {
		evt_run_branch = tree->GetBranch("evt_run");
		evt_run_branch->SetAddress(&evt_run_);
	}
	evt_nvtxs_branch = 0;
	if (tree->GetBranch("evt_nvtxs") != 0) {
		evt_nvtxs_branch = tree->GetBranch("evt_nvtxs");
		evt_nvtxs_branch->SetAddress(&evt_nvtxs_);
	}
	dielectronTrigger_branch = 0;
	if (tree->GetBranch("dielectronTrigger") != 0) {
		dielectronTrigger_branch = tree->GetBranch("dielectronTrigger");
		dielectronTrigger_branch->SetAddress(&dielectronTrigger_);
	}
	dimuonTrigger_branch = 0;
	if (tree->GetBranch("dimuonTrigger") != 0) {
		dimuonTrigger_branch = tree->GetBranch("dimuonTrigger");
		dimuonTrigger_branch->SetAddress(&dimuonTrigger_);
	}
	electronmuonTrigger_branch = 0;
	if (tree->GetBranch("electronmuonTrigger") != 0) {
		electronmuonTrigger_branch = tree->GetBranch("electronmuonTrigger");
		electronmuonTrigger_branch->SetAddress(&electronmuonTrigger_);
	}
	genps_id_branch = 0;
	if (tree->GetBranch("genps_id") != 0) {
		genps_id_branch = tree->GetBranch("genps_id");
		genps_id_branch->SetAddress(&genps_id_);
	}
	genps_id_mother_branch = 0;
	if (tree->GetBranch("genps_id_mother") != 0) {
		genps_id_mother_branch = tree->GetBranch("genps_id_mother");
		genps_id_mother_branch->SetAddress(&genps_id_mother_);
	}
	looseEl_branch = 0;
	if (tree->GetBranch("looseEl") != 0) {
		looseEl_branch = tree->GetBranch("looseEl");
		looseEl_branch->SetAddress(&looseEl_);
	}
	looseMu_branch = 0;
	if (tree->GetBranch("looseMu") != 0) {
		looseMu_branch = tree->GetBranch("looseMu");
		looseMu_branch->SetAddress(&looseMu_);
	}
	tightEl_branch = 0;
	if (tree->GetBranch("tightEl") != 0) {
		tightEl_branch = tree->GetBranch("tightEl");
		tightEl_branch->SetAddress(&tightEl_);
	}
	tightMu_branch = 0;
	if (tree->GetBranch("tightMu") != 0) {
		tightMu_branch = tree->GetBranch("tightMu");
		tightMu_branch->SetAddress(&tightMu_);
	}
	passesLoosePFJetID_branch = 0;
	if (tree->GetBranch("passesLoosePFJetID") != 0) {
		passesLoosePFJetID_branch = tree->GetBranch("passesLoosePFJetID");
		passesLoosePFJetID_branch->SetAddress(&passesLoosePFJetID_);
	}
	pfjets_corL1FastL2L3_branch = 0;
	if (tree->GetBranch("pfjets_corL1FastL2L3") != 0) {
		pfjets_corL1FastL2L3_branch = tree->GetBranch("pfjets_corL1FastL2L3");
		pfjets_corL1FastL2L3_branch->SetAddress(&pfjets_corL1FastL2L3_);
	}
	pfjets_combinedSecondaryVertexBJetTag_branch = 0;
	if (tree->GetBranch("pfjets_combinedSecondaryVertexBJetTag") != 0) {
		pfjets_combinedSecondaryVertexBJetTag_branch = tree->GetBranch("pfjets_combinedSecondaryVertexBJetTag");
		pfjets_combinedSecondaryVertexBJetTag_branch->SetAddress(&pfjets_combinedSecondaryVertexBJetTag_);
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		pfmet_isLoaded = false;
		metphi_isLoaded = false;
		pfmet_type1cor_isLoaded = false;
		gen_met_isLoaded = false;
		evt_scale1fb_isLoaded = false;
		evt_isRealData_isLoaded = false;
		evt_event_isLoaded = false;
		evt_lumiBlock_isLoaded = false;
		evt_run_isLoaded = false;
		evt_nvtxs_isLoaded = false;
		dielectronTrigger_isLoaded = false;
		dimuonTrigger_isLoaded = false;
		electronmuonTrigger_isLoaded = false;
		els_p4_isLoaded = false;
		mus_p4_isLoaded = false;
		pfjets_p4_isLoaded = false;
		genps_id_isLoaded = false;
		genps_id_mother_isLoaded = false;
		genps_p4_isLoaded = false;
		genjets_p4_isLoaded = false;
		looseEl_isLoaded = false;
		looseMu_isLoaded = false;
		tightEl_isLoaded = false;
		tightMu_isLoaded = false;
		passesLoosePFJetID_isLoaded = false;
		pfjets_corL1FastL2L3_isLoaded = false;
		pfjets_combinedSecondaryVertexBJetTag_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (pfmet_branch != 0) pfmet();
	if (metphi_branch != 0) metphi();
	if (pfmet_type1cor_branch != 0) pfmet_type1cor();
	if (gen_met_branch != 0) gen_met();
	if (evt_scale1fb_branch != 0) evt_scale1fb();
	if (evt_isRealData_branch != 0) evt_isRealData();
	if (evt_event_branch != 0) evt_event();
	if (evt_lumiBlock_branch != 0) evt_lumiBlock();
	if (evt_run_branch != 0) evt_run();
	if (evt_nvtxs_branch != 0) evt_nvtxs();
	if (dielectronTrigger_branch != 0) dielectronTrigger();
	if (dimuonTrigger_branch != 0) dimuonTrigger();
	if (electronmuonTrigger_branch != 0) electronmuonTrigger();
	if (els_p4_branch != 0) els_p4();
	if (mus_p4_branch != 0) mus_p4();
	if (pfjets_p4_branch != 0) pfjets_p4();
	if (genps_id_branch != 0) genps_id();
	if (genps_id_mother_branch != 0) genps_id_mother();
	if (genps_p4_branch != 0) genps_p4();
	if (genjets_p4_branch != 0) genjets_p4();
	if (looseEl_branch != 0) looseEl();
	if (looseMu_branch != 0) looseMu();
	if (tightEl_branch != 0) tightEl();
	if (tightMu_branch != 0) tightMu();
	if (passesLoosePFJetID_branch != 0) passesLoosePFJetID();
	if (pfjets_corL1FastL2L3_branch != 0) pfjets_corL1FastL2L3();
	if (pfjets_combinedSecondaryVertexBJetTag_branch != 0) pfjets_combinedSecondaryVertexBJetTag();
}

	float &pfmet()
	{
		if (not pfmet_isLoaded) {
			if (pfmet_branch != 0) {
				pfmet_branch->GetEntry(index);
			} else { 
				printf("branch pfmet_branch does not exist!\n");
				exit(1);
			}
			pfmet_isLoaded = true;
		}
		return pfmet_;
	}
	float &metphi()
	{
		if (not metphi_isLoaded) {
			if (metphi_branch != 0) {
				metphi_branch->GetEntry(index);
			} else { 
				printf("branch metphi_branch does not exist!\n");
				exit(1);
			}
			metphi_isLoaded = true;
		}
		return metphi_;
	}
	float &pfmet_type1cor()
	{
		if (not pfmet_type1cor_isLoaded) {
			if (pfmet_type1cor_branch != 0) {
				pfmet_type1cor_branch->GetEntry(index);
			} else { 
				printf("branch pfmet_type1cor_branch does not exist!\n");
				exit(1);
			}
			pfmet_type1cor_isLoaded = true;
		}
		return pfmet_type1cor_;
	}
	float &gen_met()
	{
		if (not gen_met_isLoaded) {
			if (gen_met_branch != 0) {
				gen_met_branch->GetEntry(index);
			} else { 
				printf("branch gen_met_branch does not exist!\n");
				exit(1);
			}
			gen_met_isLoaded = true;
		}
		return gen_met_;
	}
	float &evt_scale1fb()
	{
		if (not evt_scale1fb_isLoaded) {
			if (evt_scale1fb_branch != 0) {
				evt_scale1fb_branch->GetEntry(index);
			} else { 
				printf("branch evt_scale1fb_branch does not exist!\n");
				exit(1);
			}
			evt_scale1fb_isLoaded = true;
		}
		return evt_scale1fb_;
	}
	int &evt_isRealData()
	{
		if (not evt_isRealData_isLoaded) {
			if (evt_isRealData_branch != 0) {
				evt_isRealData_branch->GetEntry(index);
			} else { 
				printf("branch evt_isRealData_branch does not exist!\n");
				exit(1);
			}
			evt_isRealData_isLoaded = true;
		}
		return evt_isRealData_;
	}
	unsigned int &evt_event()
	{
		if (not evt_event_isLoaded) {
			if (evt_event_branch != 0) {
				evt_event_branch->GetEntry(index);
			} else { 
				printf("branch evt_event_branch does not exist!\n");
				exit(1);
			}
			evt_event_isLoaded = true;
		}
		return evt_event_;
	}
	unsigned int &evt_lumiBlock()
	{
		if (not evt_lumiBlock_isLoaded) {
			if (evt_lumiBlock_branch != 0) {
				evt_lumiBlock_branch->GetEntry(index);
			} else { 
				printf("branch evt_lumiBlock_branch does not exist!\n");
				exit(1);
			}
			evt_lumiBlock_isLoaded = true;
		}
		return evt_lumiBlock_;
	}
	unsigned int &evt_run()
	{
		if (not evt_run_isLoaded) {
			if (evt_run_branch != 0) {
				evt_run_branch->GetEntry(index);
			} else { 
				printf("branch evt_run_branch does not exist!\n");
				exit(1);
			}
			evt_run_isLoaded = true;
		}
		return evt_run_;
	}
	unsigned int &evt_nvtxs()
	{
		if (not evt_nvtxs_isLoaded) {
			if (evt_nvtxs_branch != 0) {
				evt_nvtxs_branch->GetEntry(index);
			} else { 
				printf("branch evt_nvtxs_branch does not exist!\n");
				exit(1);
			}
			evt_nvtxs_isLoaded = true;
		}
		return evt_nvtxs_;
	}
	bool &	dielectronTrigger()
	{
		if (not dielectronTrigger_isLoaded) {
			if (dielectronTrigger_branch != 0) {
				dielectronTrigger_branch->GetEntry(index);
			} else { 
				printf("branch dielectronTrigger_branch does not exist!\n");
				exit(1);
			}
			dielectronTrigger_isLoaded = true;
		}
		return dielectronTrigger_;
	}
	bool &	dimuonTrigger()
	{
		if (not dimuonTrigger_isLoaded) {
			if (dimuonTrigger_branch != 0) {
				dimuonTrigger_branch->GetEntry(index);
			} else { 
				printf("branch dimuonTrigger_branch does not exist!\n");
				exit(1);
			}
			dimuonTrigger_isLoaded = true;
		}
		return dimuonTrigger_;
	}
	bool &	electronmuonTrigger()
	{
		if (not electronmuonTrigger_isLoaded) {
			if (electronmuonTrigger_branch != 0) {
				electronmuonTrigger_branch->GetEntry(index);
			} else { 
				printf("branch electronmuonTrigger_branch does not exist!\n");
				exit(1);
			}
			electronmuonTrigger_isLoaded = true;
		}
		return electronmuonTrigger_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4()
	{
		if (not els_p4_isLoaded) {
			if (els_p4_branch != 0) {
				els_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_p4_branch does not exist!\n");
				exit(1);
			}
			els_p4_isLoaded = true;
		}
		return *els_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_p4()
	{
		if (not mus_p4_isLoaded) {
			if (mus_p4_branch != 0) {
				mus_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_p4_branch does not exist!\n");
				exit(1);
			}
			mus_p4_isLoaded = true;
		}
		return *mus_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_p4()
	{
		if (not pfjets_p4_isLoaded) {
			if (pfjets_p4_branch != 0) {
				pfjets_p4_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_p4_branch does not exist!\n");
				exit(1);
			}
			pfjets_p4_isLoaded = true;
		}
		return *pfjets_p4_;
	}
	vector<int> &genps_id()
	{
		if (not genps_id_isLoaded) {
			if (genps_id_branch != 0) {
				genps_id_branch->GetEntry(index);
			} else { 
				printf("branch genps_id_branch does not exist!\n");
				exit(1);
			}
			genps_id_isLoaded = true;
		}
		return *genps_id_;
	}
	vector<int> &genps_id_mother()
	{
		if (not genps_id_mother_isLoaded) {
			if (genps_id_mother_branch != 0) {
				genps_id_mother_branch->GetEntry(index);
			} else { 
				printf("branch genps_id_mother_branch does not exist!\n");
				exit(1);
			}
			genps_id_mother_isLoaded = true;
		}
		return *genps_id_mother_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_p4()
	{
		if (not genps_p4_isLoaded) {
			if (genps_p4_branch != 0) {
				genps_p4_branch->GetEntry(index);
			} else { 
				printf("branch genps_p4_branch does not exist!\n");
				exit(1);
			}
			genps_p4_isLoaded = true;
		}
		return *genps_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genjets_p4()
	{
		if (not genjets_p4_isLoaded) {
			if (genjets_p4_branch != 0) {
				genjets_p4_branch->GetEntry(index);
			} else { 
				printf("branch genjets_p4_branch does not exist!\n");
				exit(1);
			}
			genjets_p4_isLoaded = true;
		}
		return *genjets_p4_;
	}
	vector<bool> &looseEl()
	{
		if (not looseEl_isLoaded) {
			if (looseEl_branch != 0) {
				looseEl_branch->GetEntry(index);
			} else { 
				printf("branch looseEl_branch does not exist!\n");
				exit(1);
			}
			looseEl_isLoaded = true;
		}
		return *looseEl_;
	}
	vector<bool> &looseMu()
	{
		if (not looseMu_isLoaded) {
			if (looseMu_branch != 0) {
				looseMu_branch->GetEntry(index);
			} else { 
				printf("branch looseMu_branch does not exist!\n");
				exit(1);
			}
			looseMu_isLoaded = true;
		}
		return *looseMu_;
	}
	vector<bool> &tightEl()
	{
		if (not tightEl_isLoaded) {
			if (tightEl_branch != 0) {
				tightEl_branch->GetEntry(index);
			} else { 
				printf("branch tightEl_branch does not exist!\n");
				exit(1);
			}
			tightEl_isLoaded = true;
		}
		return *tightEl_;
	}
	vector<bool> &tightMu()
	{
		if (not tightMu_isLoaded) {
			if (tightMu_branch != 0) {
				tightMu_branch->GetEntry(index);
			} else { 
				printf("branch tightMu_branch does not exist!\n");
				exit(1);
			}
			tightMu_isLoaded = true;
		}
		return *tightMu_;
	}
	vector<bool> &passesLoosePFJetID()
	{
		if (not passesLoosePFJetID_isLoaded) {
			if (passesLoosePFJetID_branch != 0) {
				passesLoosePFJetID_branch->GetEntry(index);
			} else { 
				printf("branch passesLoosePFJetID_branch does not exist!\n");
				exit(1);
			}
			passesLoosePFJetID_isLoaded = true;
		}
		return *passesLoosePFJetID_;
	}
	vector<float> &pfjets_corL1FastL2L3()
	{
		if (not pfjets_corL1FastL2L3_isLoaded) {
			if (pfjets_corL1FastL2L3_branch != 0) {
				pfjets_corL1FastL2L3_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_corL1FastL2L3_branch does not exist!\n");
				exit(1);
			}
			pfjets_corL1FastL2L3_isLoaded = true;
		}
		return *pfjets_corL1FastL2L3_;
	}
	vector<float> &pfjets_combinedSecondaryVertexBJetTag()
	{
		if (not pfjets_combinedSecondaryVertexBJetTag_isLoaded) {
			if (pfjets_combinedSecondaryVertexBJetTag_branch != 0) {
				pfjets_combinedSecondaryVertexBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_combinedSecondaryVertexBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_combinedSecondaryVertexBJetTag_isLoaded = true;
		}
		return *pfjets_combinedSecondaryVertexBJetTag_;
	}

  static void progress( int nEventsTotal, int nEventsChain ){
    int period = 1000;
    if(nEventsTotal%1000 == 0) {
      // xterm magic from L. Vacavant and A. Cerri
      if (isatty(1)) {
        if( ( nEventsChain - nEventsTotal ) > period ){
          float frac = (float)nEventsTotal/(nEventsChain*0.01);
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
          fflush(stdout);
        }
        else {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", 100.);
          cout << endl;
        }
      }
    }
  }
  
};

#ifndef __CINT__
extern CMS2 cms2;
#endif

namespace tas {
	float &pfmet();
	float &metphi();
	float &pfmet_type1cor();
	float &gen_met();
	float &evt_scale1fb();
	int &evt_isRealData();
	unsigned int &evt_event();
	unsigned int &evt_lumiBlock();
	unsigned int &evt_run();
	unsigned int &evt_nvtxs();
	bool &dielectronTrigger();
	bool &dimuonTrigger();
	bool &electronmuonTrigger();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_p4();
	vector<int> &genps_id();
	vector<int> &genps_id_mother();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genjets_p4();
	vector<bool> &looseEl();
	vector<bool> &looseMu();
	vector<bool> &tightEl();
	vector<bool> &tightMu();
	vector<bool> &passesLoosePFJetID();
	vector<float> &pfjets_corL1FastL2L3();
	vector<float> &pfjets_combinedSecondaryVertexBJetTag();
}
#endif
