#include "CMS2.h"
CMS2 cms2;
namespace tas {
	float &pfmet() { return cms2.pfmet(); }
	float &metphi() { return cms2.metphi(); }
	float &pfmet_type1cor() { return cms2.pfmet_type1cor(); }
	float &gen_met() { return cms2.gen_met(); }
	float &evt_scale1fb() { return cms2.evt_scale1fb(); }
	int &evt_isRealData() { return cms2.evt_isRealData(); }
	unsigned int &evt_event() { return cms2.evt_event(); }
	unsigned int &evt_lumiBlock() { return cms2.evt_lumiBlock(); }
	unsigned int &evt_run() { return cms2.evt_run(); }
	bool &dielectronTrigger() { return cms2.dielectronTrigger(); }
	bool &dimuonTrigger() { return cms2.dimuonTrigger(); }
	bool &electronmuonTrigger() { return cms2.electronmuonTrigger(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4() { return cms2.els_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_p4() { return cms2.mus_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_p4() { return cms2.pfjets_p4(); }
	vector<int> &genps_id() { return cms2.genps_id(); }
	vector<int> &genps_id_mother() { return cms2.genps_id_mother(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_p4() { return cms2.genps_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genjets_p4() { return cms2.genjets_p4(); }
	vector<bool> &looseEl() { return cms2.looseEl(); }
	vector<bool> &looseMu() { return cms2.looseMu(); }
	vector<bool> &tightEl() { return cms2.tightEl(); }
	vector<bool> &tightMu() { return cms2.tightMu(); }
	vector<bool> &passesLoosePFJetID() { return cms2.passesLoosePFJetID(); }
	vector<float> &pfjets_corL1FastL2L3() { return cms2.pfjets_corL1FastL2L3(); }
}
