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
#include <sstream>

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

int scan(unsigned int njetsCut=0, TString tag=""){

    TH1F* h1D_njets_data = new TH1F("njets", "", 15, 0, 15); 
    TH1F* h1D_ht_data = new TH1F("ht", "", 20, 0, 600); 
    TH1F* h1D_met_data = new TH1F("met", "", 20, 0, 300); 
    TH1F* h1D_mt_data = new TH1F("mt", "", 20, 0, 400); 
    TH1F* h1D_nleps_good_data = new TH1F("nleps_good", "", 10, 0, 10); 
    TH1F* h1D_zmass_data = new TH1F("zmass", "", 30, 70, 120); 
    TH1F* h1D_lepeta_data = new TH1F("lepeta", "", 50, -3.0, 3.0); 
    TH1F* h1D_leppt_data = new TH1F("leppt", "", 50, 0, 150); 
    TH1F* h1D_Wleppt_data = new TH1F("Wleppt", "", 30, 0, 150); 
    TH1F* h1D_Wmetet_data = new TH1F("Wmetet", "", 40, 0, 400); 
    TH1F* h1D_Wmetphi_data = new TH1F("Wmetphi", "", 25, 0, 3.141592); 

    float luminosity = 19.407;
    float ptCut = 20;
    map<int, vector<int> > runLumi;
    cout << "njetsCut will exclude evts with njets<" << njetsCut << endl;
    // DATADATA
    {
        TChain *ch = new TChain("tree");

        ch->Add("/home/users/namin/sandbox/condorTest/output/baby_2012ABCD.root");

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
                // if(event > 30000) break;

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

                // int run = evt_run();
                // int lumi = evt_lumiBlock();
                // if(runLumi.find(run) != runLumi.end()) {
                //     if(find(runLumi[run].begin(),runLumi[run].end(),lumi) != runLumi[run].end()) {
                //         // duplicate - don't do anything
                //     } else {
                //         runLumi[run].push_back(lumi);
                //     }
                // } else {
                //     runLumi[run] = vector<int>();
                //     runLumi[run].push_back(lumi);
                // }

                std::vector<LorentzVector> goodJets;
                std::vector<LorentzVector> goodEls;
                std::vector<LorentzVector> goodMus;
                std::map<int, int> goodToP4MapEl; // map indices in good{Els,Mus} to {els,mus}_p4 indices
                std::map<int, int> goodToP4MapMu;
                float ht = 0, met = 0;

                met = pfmet_type1cor();

                if(met < 30) continue;

                // make electron quality cuts
                for(unsigned int iEl = 0; iEl < els_p4().size(); iEl++) {

                    if(!looseEl().at(iEl)) continue;
                    if(els_p4().at(iEl).pt() < ptCut) continue;
                    if(fabs(els_p4().at(iEl).eta()) > 2.4) continue;

                    goodToP4MapEl[goodEls.size()] = iEl;
                    goodEls.push_back(els_p4().at(iEl));
                }
                // mirror for muons
                for(unsigned int iMu = 0; iMu < mus_p4().size(); iMu++) {

                    if(!looseMu().at(iMu)) continue;
                    if(mus_p4().at(iMu).pt() < ptCut) continue;
                    if(fabs(mus_p4().at(iMu).eta()) > 2.4) continue;

                    goodToP4MapMu[goodMus.size()] = iMu;
                    goodMus.push_back(mus_p4().at(iMu));
                }

                h1D_nleps_good_data->Fill(goodMus.size() + goodEls.size());
                // require that we have 3 good leptons
                if(goodMus.size() + goodEls.size() != 3) continue;

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

                // XXX
                if(goodJets.size() < njetsCut) continue;


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
                        h1D_Wleppt_data->Fill(els_p4().at(pair[4]).pt());

                        // h1D_Wmetpt_data->Fill(els_p4().at(pair[4]).pt()+pfmet_type1cor());
                        h1D_Wmetet_data->Fill(els_p4().at(pair[4]).Et()+pfmet_type1cor());
                        float dphi = fabs(els_p4().at(pair[4]).phi() - metphi());
                        if(dphi > M_PI) dphi -= M_PI;
                        h1D_Wmetphi_data->Fill(dphi);

                    } else { // W lep is mu
                        if( !tightMu().at(goodToP4MapMu[pair[4]]) ) continue;
                        leps.push_back(mus_p4().at(pair[4]));
                        h1D_Wleppt_data->Fill(mus_p4().at(pair[4]).pt());

                        // h1D_Wmetpt_data->Fill(mus_p4().at(pair[4]).pt()+pfmet_type1cor());
                        h1D_Wmetet_data->Fill(mus_p4().at(pair[4]).Et()+pfmet_type1cor());
                        float dphi = fabs(mus_p4().at(pair[4]).phi() - metphi());
                        if(dphi > M_PI) dphi -= M_PI;
                        h1D_Wmetphi_data->Fill(dphi);
                    }
                } else {
                    cout << "shouldn't end up with pair[i] == -1" << endl;
                }


                // We are now in the region of interest

                for(unsigned int iMu = 0; iMu < goodMus.size(); iMu++) {
                    h1D_leppt_data->Fill(goodMus.at(iMu).pt());
                    h1D_lepeta_data->Fill(goodMus.at(iMu).eta());
                }
                for(unsigned int iEl = 0; iEl < goodEls.size(); iEl++) {
                    h1D_leppt_data->Fill(goodEls.at(iEl).pt());
                    h1D_lepeta_data->Fill(goodEls.at(iEl).eta());
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

        std::cout << "This dataset (A+B+C+D) has 19.4 fb^-1 of data" << std::endl;
        std::cout << " nGoodEvents scaled to 1/fb: " << nGoodEvents*(1.0/luminosity) << std::endl;
        std::cout << " nGoodEvents scaled to 19.4/fb: " << nGoodEvents*(19.4/luminosity) << std::endl;

    } // DATADATA

    // MCMC
    TChain *ch = new TChain("tree");

    ch->Add("/home/users/namin/sandbox/condorTest/output/baby_DYJetsToLL.root");
    ch->Add("/home/users/namin/sandbox/condorTest/output/baby_TBZToLL.root");
    ch->Add("/home/users/namin/sandbox/condorTest/output/baby_TTJets.root");
    ch->Add("/home/users/namin/sandbox/condorTest/output/baby_TTW.root");
    ch->Add("/home/users/namin/sandbox/condorTest/output/baby_TTZJets.root");
    ch->Add("/home/users/namin/sandbox/condorTest/output/baby_VVV.root");
    ch->Add("/home/users/namin/sandbox/condorTest/output/baby_WJetsToLNu.root");
    ch->Add("/home/users/namin/sandbox/condorTest/output/baby_WZ.root");
    ch->Add("/home/users/namin/sandbox/condorTest/output/baby_ZZJets.root");

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
    vector<TH1F*> h1D_nleps_good_vec;
    vector<TH1F*> h1D_zmass_vec;
    vector<TH1F*> h1D_lepeta_vec;
    vector<TH1F*> h1D_leppt_vec;
    vector<TH1F*> h1D_Wleppt_vec;

    vector<TH1F*> h1D_Wmetet_vec;
    vector<TH1F*> h1D_Wmetphi_vec;
    
    TH1F* error = new TH1F("error","",1,0,1);
    error->Sumw2();

    // File Loop
    int iFile = 0;
    while ( (currentFile = (TFile*)fileIter.Next()) ) { 

        // Get File Content
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("tree");
        cms2.Init(tree);

        iFile++; 
        TString filename(currentFile->GetTitle());


        TH1F* h1D_njets_file = new TH1F("njets"+filename, "Njets;;Entries", 15, 0, 15); 
        TH1F* h1D_ht_file = new TH1F("ht"+filename, "H_{T};GeV;Entries", 20, 0, 600); 
        TH1F* h1D_met_file = new TH1F("met"+filename, "#slash{E}_{T};GeV;Entries", 20, 0, 300); 
        TH1F* h1D_mt_file = new TH1F("mt"+filename, "M_{T};GeV;Entries", 20, 0, 400); 
        TH1F* h1D_nleps_good_file = new TH1F("nleps_good"+filename, "# of leptons (jet, quality cuts);;Entries", 10, 0, 10); 
        TH1F* h1D_zmass_file = new TH1F("zmass"+filename, "Z Mass;GeV;Entries", 30, 70, 120); 
        TH1F* h1D_lepeta_file = new TH1F("lepeta"+filename, "lepton #eta;#eta;Entries", 50, -3.0, 3.0); 
        TH1F* h1D_leppt_file = new TH1F("leppt"+filename, "lepton p_{T};p_{T} [GeV];Entries", 50, 0, 150); 
        TH1F* h1D_Wleppt_file = new TH1F("Wleppt"+filename, "W lepton p_{T};p_{T} [GeV];Entries", 30, 0, 150); 
        TH1F* h1D_Wmetet_file = new TH1F("Wmetet"+filename, "#slash{E}_{T} + E_{T} of W lepton;[GeV];Entries", 40, 0, 400); 
        TH1F* h1D_Wmetphi_file = new TH1F("Wmetphi"+filename, "|#phi_{MET}-#phi(W lep)|", 25, 0, 3.141592); 

        h1D_njets_vec.push_back(h1D_njets_file); 
        h1D_ht_vec.push_back(h1D_ht_file); 
        h1D_met_vec.push_back(h1D_met_file); 
        h1D_mt_vec.push_back(h1D_mt_file); 
        h1D_nleps_good_vec.push_back(h1D_nleps_good_file);
        h1D_zmass_vec.push_back(h1D_zmass_file); 
        h1D_lepeta_vec.push_back(h1D_lepeta_file);
        h1D_leppt_vec.push_back(h1D_leppt_file); 
        h1D_Wleppt_vec.push_back(h1D_Wleppt_file); 
        h1D_Wmetet_vec.push_back(h1D_Wmetet_file); 
        h1D_Wmetphi_vec.push_back(h1D_Wmetphi_file); 

        // Loop over Events in current file
        unsigned int nEventsTree = tree->GetEntriesFast();
        for( unsigned int event = 0; event < nEventsTree; ++event) {

            // XXX
            // if(event > 30000) break;

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

            float scale = evt_scale1fb()  * luminosity * 1674.0/1415.4; // ptCut 20
            // float scale = evt_scale1fb()  * luminosity * 1674.0/1415.4 * 259/204.5; // ptCut 20 njets >= 2
            // float scale = evt_scale1fb()  * luminosity * 1202.0/1072.2; // ptCut 25
            // float scale = evt_scale1fb()  * luminosity;

            met = pfmet_type1cor();

            if(met < 30) continue;

            // make electron quality cuts
            for(unsigned int iEl = 0; iEl < els_p4().size(); iEl++) {

                if(!looseEl().at(iEl)) continue;
                if(els_p4().at(iEl).pt() < ptCut) continue;
                if(fabs(els_p4().at(iEl).eta()) > 2.4) continue;

                goodToP4MapEl[goodEls.size()] = iEl;
                goodEls.push_back(els_p4().at(iEl));
            }
            // mirror for muons
            for(unsigned int iMu = 0; iMu < mus_p4().size(); iMu++) {

                if(!looseMu().at(iMu)) continue;
                if(mus_p4().at(iMu).pt() < ptCut) continue;
                if(fabs(mus_p4().at(iMu).eta()) > 2.4) continue;

                goodToP4MapMu[goodMus.size()] = iMu;
                goodMus.push_back(mus_p4().at(iMu));
            }


            h1D_nleps_good_file->Fill(goodMus.size() + goodEls.size(), scale);
            // require that we have 3 good leptons
            if(goodMus.size() + goodEls.size() != 3) continue;

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


            // XXX
            if(goodJets.size() < njetsCut) continue;

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
                    h1D_Wleppt_file->Fill(els_p4().at(pair[4]).pt(), scale);

                    h1D_Wmetet_file->Fill(els_p4().at(pair[4]).Et()+pfmet_type1cor(), scale);
                    float dphi = fabs(els_p4().at(pair[4]).phi() - metphi());
                    if(dphi > M_PI) dphi -= M_PI;
                    h1D_Wmetphi_file->Fill(dphi,scale);

                } else { // W lep is mu
                    if( !tightMu().at(goodToP4MapMu[pair[4]]) ) continue;
                    leps.push_back(mus_p4().at(pair[4]));
                    h1D_Wleppt_file->Fill(mus_p4().at(pair[4]).pt(), scale);

                    h1D_Wmetet_file->Fill(mus_p4().at(pair[4]).Et()+pfmet_type1cor(), scale);
                    float dphi = fabs(mus_p4().at(pair[4]).phi() - metphi());
                    if(dphi > M_PI) dphi -= M_PI;
                    h1D_Wmetphi_file->Fill(dphi,scale);

                }
            } else {
                cout << "shouldn't end up with pair[i] == -1" << endl;
            }


            // We are now in the region of interest

            for(unsigned int iMu = 0; iMu < goodMus.size(); iMu++) {
                h1D_leppt_file->Fill(goodMus.at(iMu).pt(), scale);
                h1D_lepeta_file->Fill(goodMus.at(iMu).eta(), scale);
            }
            for(unsigned int iEl = 0; iEl < goodEls.size(); iEl++) {
                h1D_leppt_file->Fill(goodEls.at(iEl).pt(), scale);
                h1D_lepeta_file->Fill(goodEls.at(iEl).eta(), scale);
            }



            nGoodEvents++;
            nGoodEventsWeighted+=scale*1.0;

            double mt = MT(leps[0]+leps[1]+leps[2], pfmet_type1cor(), metphi());
            h1D_mt_file->Fill(mt, scale);
            h1D_zmass_file->Fill(mass, scale); 
            h1D_njets_file->Fill(goodJets.size(), scale);
            h1D_ht_file->Fill(ht, scale);
            h1D_met_file->Fill(met, scale);
            error->Fill(0.0,scale);

        }//event loop

    }//file loop MCMC

    std::cout << " nGoodEvents: " << nGoodEvents << " nEventsTotal: " << nEventsTotal << std::endl;
    std::cout << " nGoodEventsWeighted to 19.4 1/fb: " << nGoodEventsWeighted << std::endl;

    std::cout << " error->GetBinContent(1): " << error->GetBinContent(1) << " error->GetBinError(1): " << error->GetBinError(1) << std::endl;

    TString prefix("plots");
    prefix += tag;
    prefix += "/";

    stringstream ss; ss << njetsCut;
    std::string common = " --luminosity 19.4 --percentages --label njets>="+ss.str();

    drawStacked(h1D_njets_data, h1D_njets_vec, prefix+"h1D_njets.pdf", "--centerlabel"+common);
    drawStacked(h1D_ht_data, h1D_ht_vec,prefix+"h1D_ht.pdf","--logscale --binsize --centerlabel"+common);
    drawStacked(h1D_leppt_data, h1D_leppt_vec,prefix+"h1D_leppt.pdf","--logscale --centerlabel"+common);
    drawStacked(h1D_Wleppt_data, h1D_Wleppt_vec,prefix+"h1D_Wleppt.pdf","--logscale --centerlabel"+common);
    drawStacked(h1D_nleps_good_data, h1D_nleps_good_vec,prefix+"h1D_nleps_good.pdf","--logscale --centerlabel"+common);
    drawStacked(h1D_met_data, h1D_met_vec,prefix+"h1D_met.pdf","--binsize --centerlabel"+common);
    drawStacked(h1D_mt_data, h1D_mt_vec,prefix+"h1D_mt.pdf","--binsize --centerlabel"+common);
    drawStacked(h1D_lepeta_data, h1D_lepeta_vec,prefix+"h1D_lepeta.pdf",""+common);
    drawStacked(h1D_zmass_data, h1D_zmass_vec,prefix+"h1D_zmass.pdf","--binsize"+common);
    drawStacked(h1D_Wmetet_data, h1D_Wmetet_vec,prefix+"h1D_Wmetet.pdf","--binsize --centerlabel"+common);
    drawStacked(h1D_Wmetphi_data, h1D_Wmetphi_vec,prefix+"h1D_Wmetphi.pdf","-binsize --centerlabel"+common);

    // for(map<int, vector<int> >::iterator it = runLumi.begin(); it != runLumi.end(); it++) {
    //     for(int ilumi = 0; ilumi < it->second.size(); ilumi++) {
    //         cout << it->first << " " << it->second.at(ilumi) << endl;
    //     }
    // }

    return 0;
}
