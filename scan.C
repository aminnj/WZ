#include "TTreeCache.h"

#include "TChain.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TRegexp.h"
#include "TStyle.h"
#include "Math/VectorUtil.h"
#include "TDatabasePDG.h"

#include "./CMS2.h"

#include "/home/users/namin/macros/utils.C"

#include <iostream>

#include <math.h>

#include "Include.C"
#include "goodrun.cc"

using namespace std;
using namespace tas;

vector<int> findZPair(std::vector<LorentzVector> goodEls, std::vector<LorentzVector> goodMus) {
    int nEl = goodEls.size();
    int nMu = goodMus.size();
    float ZMass = 91.1876;
    // 1st elem is el (0) or mu (1), next 2 are the Z lep idxs
    // 4th is el/mu of W lep, last is W lep idx
    vector<int> v(5, -1);
    // Z decays to pair of same flavour leptons, so we can check
    // electrons and muons separately
    if(nEl >= 2) {
        v[0] = 0; // pair must be electrons
        if(nEl == 2) {
            v[1] = 0; v[2] = 1; // found the only possible pair
            v[3] = 1; v[4] = 0; // W lep must be mu 
            return v;
        } else if(nEl == 3) { // 3 possible pairs to check now: 0,1 and 0,2 and 1,2
            int iElZ = -1, jElZ = -1;
            float minDiffFromZMass = 9999.0;
            for(int i = 0; i < 3; i++) {       // loop through idx pairs 0,1 0,2 1,2
                for(int j = i+1; j < 3; j++) { // retain indices with best masses
                    float mass = (goodEls[i] + goodEls[j]).mass();
                    if(fabs(mass - ZMass) < minDiffFromZMass) {
                        iElZ = i; jElZ = j;
                        minDiffFromZMass = fabs(mass - ZMass);
                    }
                }
            }
            // now we've found the pair with mass closest to ZMass
            v[1] = iElZ; v[2] = jElZ; v[3] = 0; // W lep must be el
            v[4] = (0+1+2) - iElZ - jElZ; // index of W lepton
            if(v[1] != -1 && v[2] != -1) return v;
        }
    } else if(nMu >= 2) {
        v[0] = 1; // pair must be muons
        if(nMu == 2) {
            v[1] = 0; v[2] = 1; // found the only possible pair
            v[3] = 0; v[4] = 0; // W lep must be el 
            return v;
        } else if(nMu == 3) { // 3 possible pairs to check now: 0,1 and 0,2 and 1,2
            int iMuZ = -1, jMuZ = -1;
            float minDiffFromZMass = 9999.0;
            for(int i = 0; i < 3; i++) {       // loop through idx pairs 0,1 0,2 1,2
                for(int j = i+1; j < 3; j++) { // retain indices with best masses
                    float mass = (goodMus[i] + goodMus[j]).mass();
                    if(fabs(mass - ZMass) < minDiffFromZMass) {
                        iMuZ = i; jMuZ = j;
                        minDiffFromZMass = fabs(mass - ZMass);
                    }
                }
            }
            // now we've found the pair with mass closest to ZMass
            v[1] = iMuZ; v[2] = jMuZ; v[3] = 1; // W lep must be mu
            v[4] = (0+1+2) - iMuZ - jMuZ; // index of W lepton
            if(v[1] != -1 && v[2] != -1) return v;
        }
    }
    cout << "shouldn't end up at end of zpair function" << endl;
    return v;
}

int scan(){

    TH1F* h1D_njets_data = new TH1F("njets", "Njets", 15, 0, 15); 
    TH1F* h1D_ht_data = new TH1F("ht", "H_{T};GeV", 100, 0, 600); 
    TH1F* h1D_met_data = new TH1F("met", "#slash{E}_{T};GeV", 100, 0, 300); 
    TH1F* h1D_mt_data = new TH1F("mt", "M_{T};GeV", 100, 0, 400); 
    TH1F* h1D_nleps_data = new TH1F("nleps", "# of leptons (no cuts)", 10, 0, 10); 
    TH1F* h1D_nleps_good_data = new TH1F("nleps_good", "# of leptons (jet, quality cuts)", 10, 0, 10); 
    TH1F* h1D_zmass_data = new TH1F("zmass", "Z Mass;GeV", 100, 80, 110); 
    TH1F* h1D_leppt_all_data = new TH1F("leppt_all", "lepton p_{T} (no cuts);p_{T} [GeV]", 100, 0, 150); 
    TH1F* h1D_lepeta_all_data = new TH1F("lepeta_all", "lepton #eta (no cuts)", 100, -3.0, 3.0); 
    TH1F* h1D_leppt_data = new TH1F("leppt", "lepton p_{T} after jet, lep quality cuts;p_{T} [GeV]", 100, 0, 150); 
    TH1F* h1D_lepeta_data = new TH1F("lepeta", "lepton #eta after jet, lep quality cuts", 100, -3.0, 3.0); 

    float luminosity = 19.49;
    // DATA
    {
        TChain *ch = new TChain("tree");

        ch->Add("/hadoop/cms/store/user/namin/wz/baby_2012A*.root");
        ch->Add("/hadoop/cms/store/user/namin/wz/baby_2012B*.root");
        ch->Add("/hadoop/cms/store/user/namin/wz/baby_2012C*.root");
        ch->Add("/hadoop/cms/store/user/namin/wz/baby_2012D*.root");

        int nEventsTotal = ch->GetEntries();
        int nEventsSoFar = 0;
        int nGoodEvents = 0;

        TFile *currentFile = 0;
        TObjArray *listOfFiles = ch->GetListOfFiles();
        TIter fileIter(listOfFiles);

        // File Loop
        while ( (currentFile = (TFile*)fileIter.Next()) ) { 
            TFile *file = new TFile( currentFile->GetTitle() );
            TTree *tree = (TTree*)file->Get("tree");
            cms2.Init(tree);
            // Set Good Run List
            if(evt_isRealData()) set_goodrun_file("final_19p49fb_cms2.txt");

            TString filename(currentFile->GetTitle());

            // Loop over Events in current file
            unsigned int nEventsTree = tree->GetEntriesFast();
            for( unsigned int event = 0; event < nEventsTree; ++event) {

                // XXX
                //if(event > 300000) break;

                // Get Event Content
                cms2.GetEntry(event);
                nEventsSoFar++;

                // Progress
                CMS2::progress( nEventsSoFar, nEventsTotal );

                // Select Good Runs
                if( evt_isRealData() && !goodrun( evt_run(), evt_lumiBlock() ) ) continue;

                if(evt_isRealData()){
                    DorkyEventIdentifier id = { evt_run(), evt_event(), evt_lumiBlock() };
                    if ( is_duplicate(id) ){
                        continue;
                    }
                }

                std::vector<LorentzVector> goodJets;
                std::vector<LorentzVector> goodEls;
                std::vector<LorentzVector> goodMus;
                std::map<int, int> goodToP4MapEl; // map indices in good{Els,Mus} to {els,mus}_p4 indices
                std::map<int, int> goodToP4MapMu;
                float ht = 0, met = 0;

                h1D_nleps_data->Fill(els_p4().size() + mus_p4().size());

                // make electron quality cuts
                for(unsigned int iEl = 0; iEl < els_p4().size(); iEl++) {

                    h1D_leppt_all_data->Fill(els_p4().at(iEl).pt());
                    h1D_lepeta_all_data->Fill(els_p4().at(iEl).eta());

                    if(!looseEl().at(iEl)) continue;
                    if(els_p4().at(iEl).pt() < 25) continue;
                    if(fabs(els_p4().at(iEl).eta()) > 2.4) continue;

                    h1D_leppt_data->Fill(els_p4().at(iEl).pt());
                    h1D_lepeta_data->Fill(els_p4().at(iEl).eta());

                    goodEls.push_back(els_p4().at(iEl));
                    goodToP4MapEl[goodEls.size()] = iEl;
                }
                // mirror for muons
                for(unsigned int iMu = 0; iMu < mus_p4().size(); iMu++) {
                    h1D_leppt_all_data->Fill(mus_p4().at(iMu).pt());
                    h1D_lepeta_all_data->Fill(mus_p4().at(iMu).eta());

                    if(!looseMu().at(iMu)) continue;
                    if(mus_p4().at(iMu).pt() < 25) continue;
                    if(fabs(mus_p4().at(iMu).eta()) > 2.4) continue;

                    h1D_leppt_data->Fill(mus_p4().at(iMu).pt());
                    h1D_lepeta_data->Fill(mus_p4().at(iMu).eta());

                    goodMus.push_back(mus_p4().at(iMu));
                    goodToP4MapMu[goodMus.size()] = iMu;
                }

                // select good jets
                for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++){
                    if (pfjets_p4().at(iJet).pt()*pfjets_corL1FastL2L3().at(iJet) < 40) continue;
                    if (fabs(pfjets_p4().at(iJet).eta()) > 2.4) continue;
                    if (!passesLoosePFJetID().at(iJet)) continue;

                    // deltaR match jets and leptons; if match, get rid of jet
                    bool isIsolatedFromLepton = true;
                    for(unsigned int iMu = 0; iMu < goodMus.size(); iMu++) {
                        if( deltaR(goodMus.at(iMu),pfjets_p4().at(iJet)) < 0.3 ) {
                            isIsolatedFromLepton = false; break;
                        }
                    }
                    for(unsigned int iEl = 0; iEl < goodEls.size(); iEl++) {
                        if( deltaR(goodEls.at(iEl),pfjets_p4().at(iJet)) < 0.3 ) {
                            isIsolatedFromLepton = false; break;
                        }
                    }
                    if(!isIsolatedFromLepton) continue;

                    ht += pfjets_p4().at(iJet).pt()*pfjets_corL1FastL2L3().at(iJet);
                    goodJets.push_back(pfjets_p4().at(iJet));
                }

                met = pfmet_type1cor();
                if(met < 30) continue;

                h1D_nleps_good_data->Fill(goodMus.size() + goodEls.size());
                // require that we have 3 good leptons
                if(goodMus.size() + goodEls.size() != 3) continue;

                vector<int> pair = findZPair(goodEls, goodMus);
                vector<LorentzVector> leps;
                float mass = 0.0;
                if( pair[0] != -1 && pair[1] != -1 && pair[2] != -1 && 
                        pair[3] != -1 && pair[4] != -1 ) {
                    if(pair[0] == 0) { // el
                        mass += (els_p4().at(pair[1])+els_p4().at(pair[2])).mass();
                        leps.push_back(els_p4().at(pair[1]));
                        leps.push_back(els_p4().at(pair[2]));
                    } else { // mu
                        mass += (mus_p4().at(pair[1])+mus_p4().at(pair[2])).mass();
                        leps.push_back(mus_p4().at(pair[1]));
                        leps.push_back(mus_p4().at(pair[2]));
                    }

                    if(pair[3] == 0) { // W lep is el
                        if( !tightEl().at(goodToP4MapEl[pair[4]]) ) continue;
                        leps.push_back(els_p4().at(pair[4]));

                    } else { // W lep is mu
                        if( !tightMu().at(goodToP4MapMu[pair[4]]) ) continue;
                        leps.push_back(mus_p4().at(pair[4]));
                    }
                } else {
                    cout << "shouldn't end up with pair[i] == -1" << endl;
                }

                nGoodEvents++;

                double mt = MT(leps[0]+leps[1]+leps[2], pfmet_type1cor(), metphi());
                h1D_mt_data->Fill(mt);
                h1D_zmass_data->Fill(mass); 
                h1D_njets_data->Fill(goodJets.size());
                h1D_ht_data->Fill(ht);
                h1D_met_data->Fill(met);


            }//event loop

        }//file loop


        std::cout << " nGoodEvents: " << nGoodEvents << " nEventsTotal: " << nEventsTotal << std::endl;

        std::cout << "This dataset (A+B+C+D) has 19.49 fb^-1 of data" << std::endl;
        std::cout << " nGoodEvents scaled to 1/fb: " << nGoodEvents*(1.0/luminosity) << std::endl;
        std::cout << " nGoodEvents scaled to 19.49/fb: " << nGoodEvents*(19.49/luminosity) << std::endl;

    } // DATA

    // MC
    TChain *ch = new TChain("tree");

    ch->Add("/hadoop/cms/store/user/namin/wz/baby_DYJetsToLL.root");
    ch->Add("/hadoop/cms/store/user/namin/wz/baby_TBZToLL.root");
    ch->Add("/hadoop/cms/store/user/namin/wz/baby_TTJets.root");
    ch->Add("/hadoop/cms/store/user/namin/wz/baby_TTW.root");
    ch->Add("/hadoop/cms/store/user/namin/wz/baby_TTZJets.root");
    ch->Add("/hadoop/cms/store/user/namin/wz/baby_VVVJets.root");
    ch->Add("/hadoop/cms/store/user/namin/wz/baby_WJetsToLNu.root");
    ch->Add("/hadoop/cms/store/user/namin/wz/baby_WZ.root");
    ch->Add("/hadoop/cms/store/user/namin/wz/baby_ZZJetsTo2L2Nu.root");
    ch->Add("/hadoop/cms/store/user/namin/wz/baby_ZZJetsTo2L2Q.root");
    ch->Add("/hadoop/cms/store/user/namin/wz/baby_ZZJetsTo4L.root");

    int nEventsTotal = ch->GetEntries();
    int nEventsSoFar = 0;
    int nGoodEvents = 0;
    float nGoodEventsWeighted = 0;

    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);


    vector<TH1F*> h1D_njets_vec;
    vector<TH1F*> h1D_ht_vec;
    vector<TH1F*> h1D_met_vec;
    vector<TH1F*> h1D_mt_vec;
    vector<TH1F*> h1D_nleps_vec;
    vector<TH1F*> h1D_nleps_good_vec;
    vector<TH1F*> h1D_zmass_vec;
    vector<TH1F*> h1D_lepeta_all_vec;
    vector<TH1F*> h1D_leppt_vec;
    vector<TH1F*> h1D_leppt_all_vec;
    vector<TH1F*> h1D_lepeta_vec;

    // File Loop
    int iFile = 0;
    while ( (currentFile = (TFile*)fileIter.Next()) ) { 

        // Get File Content
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("tree");
        cms2.Init(tree);

        iFile++; 
        TString filename(currentFile->GetTitle());


        TH1F* h1D_njets_file = new TH1F("njets"+filename, "Njets", 15, 0, 15); 
        TH1F* h1D_ht_file = new TH1F("ht"+filename, "H_{T};GeV", 100, 0, 600); 
        TH1F* h1D_met_file = new TH1F("met"+filename, "#slash{E}_{T};GeV", 100, 0, 300); 
        TH1F* h1D_mt_file = new TH1F("mt"+filename, "M_{T};GeV", 100, 0, 400); 
        TH1F* h1D_nleps_file = new TH1F("nleps"+filename, "# of leptons (no cuts)", 10, 0, 10); 
        TH1F* h1D_nleps_good_file = new TH1F("nleps_good"+filename, "# of leptons (jet, quality cuts)", 10, 0, 10); 
        TH1F* h1D_zmass_file = new TH1F("zmass"+filename, "Z Mass;GeV", 100, 80, 110); 
        TH1F* h1D_leppt_all_file = new TH1F("leppt_all"+filename, "lepton p_{T} (no cuts);p_{T} [GeV]", 100, 0, 150); 
        TH1F* h1D_lepeta_all_file = new TH1F("lepeta_all"+filename, "lepton #eta (no cuts)", 100, -3.0, 3.0); 
        TH1F* h1D_leppt_file = new TH1F("leppt"+filename, "lepton p_{T} after jet, lep quality cuts;p_{T} [GeV]", 100, 0, 150); 
        TH1F* h1D_lepeta_file = new TH1F("lepeta"+filename, "lepton #eta after jet, lep quality cuts", 100, -3.0, 3.0); 

        h1D_njets_vec.push_back(h1D_njets_file); 
        h1D_ht_vec.push_back(h1D_ht_file); 
        h1D_met_vec.push_back(h1D_met_file); 
        h1D_mt_vec.push_back(h1D_mt_file); 
        h1D_nleps_vec.push_back(h1D_nleps_file);
        h1D_nleps_good_vec.push_back(h1D_nleps_good_file);
        h1D_zmass_vec.push_back(h1D_zmass_file); 
        h1D_leppt_all_vec.push_back(h1D_leppt_all_file);
        h1D_lepeta_all_vec.push_back(h1D_lepeta_all_file);
        h1D_leppt_vec.push_back(h1D_leppt_file); 
        h1D_lepeta_vec.push_back(h1D_lepeta_file); 

        // Loop over Events in current file
        unsigned int nEventsTree = tree->GetEntriesFast();
        for( unsigned int event = 0; event < nEventsTree; ++event) {

            // XXX
            //if(event > 300000) break;

            // Get Event Content
            cms2.GetEntry(event);
            nEventsSoFar++;

            if(event == 0) {
                std::cout << " evt_scale1fb(): " << evt_scale1fb() << " filename: " << filename << std::endl;
            }
            // Progress
            CMS2::progress( nEventsSoFar, nEventsTotal );

            std::vector<LorentzVector> goodJets;
            std::vector<LorentzVector> goodEls;
            std::vector<LorentzVector> goodMus;
            std::map<int, int> goodToP4MapEl; // map indices in good{Els,Mus} to {els,mus}_p4 indices
            std::map<int, int> goodToP4MapMu;
            float ht = 0, met = 0;

            float scale = evt_scale1fb()  * luminosity;

            h1D_nleps_file->Fill(els_p4().size() + mus_p4().size(), scale);

            // make electron quality cuts
            for(unsigned int iEl = 0; iEl < els_p4().size(); iEl++) {
                float elPt = els_p4().at(iEl).pt();
                h1D_leppt_all_file->Fill(elPt, scale);

                h1D_lepeta_all_file->Fill(els_p4().at(iEl).eta(), scale);

                if(!looseEl().at(iEl)) continue;
                if(els_p4().at(iEl).pt() < 25) continue;
                if(fabs(els_p4().at(iEl).eta()) > 2.4) continue;

                h1D_leppt_file->Fill(els_p4().at(iEl).pt(), scale);
                h1D_lepeta_file->Fill(els_p4().at(iEl).eta(), scale);

                goodEls.push_back(els_p4().at(iEl));
                goodToP4MapEl[goodEls.size()] = iEl;
            }
            // mirror for muons
            for(unsigned int iMu = 0; iMu < mus_p4().size(); iMu++) {
                float muPt = mus_p4().at(iMu).pt();
                h1D_leppt_all_file->Fill(muPt, scale);

                h1D_lepeta_all_file->Fill(mus_p4().at(iMu).eta(), scale);

                if(!looseMu().at(iMu)) continue;
                if(mus_p4().at(iMu).pt() < 25) continue;
                if(fabs(mus_p4().at(iMu).eta()) > 2.4) continue;

                h1D_leppt_file->Fill(mus_p4().at(iMu).pt(), scale);
                h1D_lepeta_file->Fill(mus_p4().at(iMu).eta(), scale);

                goodMus.push_back(mus_p4().at(iMu));
                goodToP4MapMu[goodMus.size()] = iMu;
            }

            // select good jets
            for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++){
                if (pfjets_p4().at(iJet).pt()*pfjets_corL1FastL2L3().at(iJet) < 40) continue;
                if (fabs(pfjets_p4().at(iJet).eta()) > 2.4) continue;
                if (!passesLoosePFJetID().at(iJet)) continue;

                // deltaR match jets and leptons; if match, get rid of jet
                bool isIsolatedFromLepton = true;
                for(unsigned int iMu = 0; iMu < goodMus.size(); iMu++) {
                    if( deltaR(goodMus.at(iMu),pfjets_p4().at(iJet)) < 0.3 ) {
                        isIsolatedFromLepton = false; break;
                    }
                }
                for(unsigned int iEl = 0; iEl < goodEls.size(); iEl++) {
                    if( deltaR(goodEls.at(iEl),pfjets_p4().at(iJet)) < 0.3 ) {
                        isIsolatedFromLepton = false; break;
                    }
                }
                if(!isIsolatedFromLepton) continue;


                ht += pfjets_p4().at(iJet).pt()*pfjets_corL1FastL2L3().at(iJet);
                goodJets.push_back(pfjets_p4().at(iJet));
            }

            met = pfmet_type1cor();
            if(met < 30) continue;


            h1D_nleps_good_file->Fill(goodMus.size() + goodEls.size(), scale);
            // require that we have 3 good leptons
            if(goodMus.size() + goodEls.size() != 3) continue;

            vector<int> pair = findZPair(goodEls, goodMus);
            vector<LorentzVector> leps;
            float mass = 0.0;
            if( pair[0] != -1 && pair[1] != -1 && pair[2] != -1 && 
                    pair[3] != -1 && pair[4] != -1 ) {
                if(pair[0] == 0) { // el
                    mass += (els_p4().at(pair[1])+els_p4().at(pair[2])).mass();
                    leps.push_back(els_p4().at(pair[1]));
                    leps.push_back(els_p4().at(pair[2]));
                } else { // mu
                    mass += (mus_p4().at(pair[1])+mus_p4().at(pair[2])).mass();
                    leps.push_back(mus_p4().at(pair[1]));
                    leps.push_back(mus_p4().at(pair[2]));
                }

                if(pair[3] == 0) { // W lep is el
                    if( !tightEl().at(goodToP4MapEl[pair[4]]) ) continue;
                    leps.push_back(els_p4().at(pair[4]));

                } else { // W lep is mu
                    if( !tightMu().at(goodToP4MapMu[pair[4]]) ) continue;
                    leps.push_back(mus_p4().at(pair[4]));

                }
            } else {
                cout << "shouldn't end up with pair[i] == -1" << endl;
            }

            nGoodEvents++;
            nGoodEventsWeighted+=scale*1.0;

            double mt = MT(leps[0]+leps[1]+leps[2], pfmet_type1cor(), metphi());
            h1D_mt_file->Fill(mt, scale);
            h1D_zmass_file->Fill(mass, scale); 
            h1D_njets_file->Fill(goodJets.size(), scale);
            h1D_ht_file->Fill(ht, scale);
            h1D_met_file->Fill(met, scale);

        }//event loop

    }//file loop

    std::cout << " nGoodEvents: " << nGoodEvents << " nEventsTotal: " << nEventsTotal << std::endl;
    std::cout << " nGoodEventsWeighted to 19.49 1/fb: " << nGoodEventsWeighted << std::endl;

    TString prefix("plots/");

    std::string common = " --luminosity 19.49 --percentages ";

    drawStacked(h1D_leppt_all_data, h1D_leppt_all_vec, prefix+"h1D_leppt_all.pdf", "--logscale"+common);
    drawStacked(h1D_njets_data, h1D_njets_vec, prefix+"h1D_njets.pdf", ""+common);
    drawStacked(h1D_ht_data, h1D_ht_vec,prefix+"h1D_ht.pdf","--logscale"+common);
    drawStacked(h1D_leppt_data, h1D_leppt_vec,prefix+"h1D_leppt.pdf","--logscale"+common);
    drawStacked(h1D_nleps_good_data, h1D_nleps_good_vec,prefix+"h1D_nleps_good.pdf","--logscale"+common);
    drawStacked(h1D_nleps_data, h1D_nleps_vec,prefix+"h1D_nleps.pdf","--logscale"+common);
    drawStacked(h1D_met_data, h1D_met_vec,prefix+"h1D_met.pdf",""+common);
    drawStacked(h1D_mt_data, h1D_mt_vec,prefix+"h1D_mt.pdf",""+common);
    drawStacked(h1D_lepeta_all_data, h1D_lepeta_all_vec,prefix+"h1D_lepeta_all.pdf",""+common);
    drawStacked(h1D_lepeta_data, h1D_lepeta_vec,prefix+"h1D_lepeta.pdf",""+common);
    drawStacked(h1D_zmass_data, h1D_zmass_vec,prefix+"h1D_zmass.pdf",""+common);


    // std::string common = " --luminosity 19.49 --percentages";
    // TH1F* dog = new TH1F("dog","dog",1,0,1);
    // vector<TH1F*> h1D_njets_vec;
    // TH1F* test = new TH1F("test","test",15,0,15);
    // h1D_njets_vec.push_back(test);
    // test->Fill(1); test->Fill(2);
    // test->Fill(2); test->Fill(3);
    // drawStacked(h1D_njets_data, h1D_njets_vec, prefix+"h1D_njets.pdf", ""+common);

    return 0;
}
