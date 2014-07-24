#include "TTreeCache.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "Math/VectorUtil.h"
#include "TDatabasePDG.h"

#include <iostream>
#include <math.h>
#include <sstream>

#include "Include.C"
#include "goodrun.cc"
#include "./CMS2.h"
#include "/home/users/namin/macros/utils.C"

using namespace std;
using namespace tas;

struct JetStruct {
    LorentzVector jet;
    float pt = 0.0;
    unsigned int idx = 0;
};

bool jetCompare(JetStruct j1, JetStruct j2) {
    return j1.pt > j2.pt;
}

float scan(unsigned int njetsCut=0, int btagCut=999, TString tag="", float manualScale=-1.0){

    TH1F* h1D_dummy_data = new TH1F("dummy", "dummyhisto", 10, 0, 10); 

    TH1F* h1D_njets_data = new TH1F("njets", "", 15, 0, 15); 
    TH1F* h1D_ht_data = new TH1F("ht", "", 20, 0, 600); 
    TH1F* h1D_met_data = new TH1F("met", "", 20, 0, 300); 
    TH1F* h1D_mt_data = new TH1F("mt", "", 20, 0, 400); 
    TH1F* h1D_mtW_data = new TH1F("mtW", "", 40, 0, 200); 
    TH1F* h1D_zmass_data = new TH1F("zmass", "", 30, 70, 120); 
    TH1F* h1D_lepeta_data = new TH1F("lepeta", "", 50, -3.0, 3.0); 
    TH1F* h1D_leppt_data = new TH1F("leppt", "", 50, 0, 150); 
    TH1F* h1D_Wleppt_data = new TH1F("Wleppt", "", 30, 0, 150); 

    TH1F* h1D_nbtags_data = new TH1F("nbtags", "", 15, 0, 15); 
    TH1F* h1D_btagval_data = new TH1F("btagval", "", 20, 0, 1); 
    TH1F* h1D_ptZ_data = new TH1F("ptZ", "", 25, 0, 400);
    TH1F* h1D_st_data = new TH1F("st", "", 25, 100, 800);
    TH1F* h1D_minRLeadingJet_data = new TH1F("minRLeadingJet", "",10,0,5);
    TH1F* h1D_ptj1_data = new TH1F("ptj1", "", 25, 0, 300);
    TH1F* h1D_ptj2_data = new TH1F("ptj2", "", 25, 0, 300);
    TH1F* h1D_massZj1_data = new TH1F("massZj1", "", 25, 100, 600);
    TH1F* h1D_massZj2_data = new TH1F("massZj2", "", 25, 100, 600);
    TH1F* h1D_ptjj_data = new TH1F("ptjj", "", 25, 0, 400);
    TH1F* h1D_massjj_data = new TH1F("massjj", "", 25, 0, 600);

    TH1F* h1D_mtWeemu_data = new TH1F("mtWeem", "", 40, 0, 200); 
    TH1F* h1D_mtWmumue_data = new TH1F("mtWmme", "", 40, 0, 200); 
    TH1F* h1D_mtWeee_data = new TH1F("mtWeee", "", 40, 0, 200); 
    TH1F* h1D_mtWmumumu_data = new TH1F("mtWmmm", "", 40, 0, 200); 
    TH1F* h1D_lepetae_data = new TH1F("lepetae", "", 50, -3.0, 3.0); 
    TH1F* h1D_lepetamu_data = new TH1F("lepetamu", "", 50, -3.0, 3.0); 

    float luminosity = 19.407;
    float ptCut = 20;
    float zmassCut = 10;
    map<int, vector<int> > runLumi;
    cout << "njetsCut will exclude evts with njets<" << njetsCut << endl;
    // DATADATA
    {

        clear_seen();

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

                std::vector<LorentzVector> goodEls;
                std::vector<LorentzVector> goodMus;
                std::vector<JetStruct> goodJets;
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


                    JetStruct myJet = {};
                    myJet.jet = pfjets_p4().at(iJet);
                    myJet.pt = pfjets_p4().at(iJet).pt()*pfjets_corL1FastL2L3().at(iJet);
                    myJet.idx = iJet;

                    goodJets.push_back(myJet);
                }



                vector<int> pair = findZPair(goodEls, goodMus);
                vector<LorentzVector> leps;
                float mtW = 0.0, mass = 0.0;
                if( pair[0] != -1 && pair[1] != -1 && pair[2] != -1 && 
                        pair[3] != -1 && pair[4] != -1 ) {
                    if(pair[0] == 0) { // el
                        mass += (goodEls.at(pair[1])+goodEls.at(pair[2])).mass();
                        leps.push_back(goodEls.at(pair[1]));
                        leps.push_back(goodEls.at(pair[2]));
                    } else { // mu
                        mass += (goodMus.at(pair[1])+goodMus.at(pair[2])).mass();
                        leps.push_back(goodMus.at(pair[1]));
                        leps.push_back(goodMus.at(pair[2]));
                    }

                    if(pair[3] == 0) { // W lep is el
                        if( !tightEl().at(goodToP4MapEl[pair[4]]) ) continue;
                        leps.push_back(goodEls.at(pair[4]));

                        mtW = MT(els_p4().at(pair[4]), pfmet_type1cor(), metphi());

                    } else { // W lep is mu
                        if( !tightMu().at(goodToP4MapMu[pair[4]]) ) continue;
                        leps.push_back(goodMus.at(pair[4]));

                        mtW = MT(goodMus.at(pair[4]), pfmet_type1cor(), metphi());

                    }
                } else {
                    cout << "shouldn't end up with pair[i] == -1" << endl;
                }


                // XXX
                if(abs(mass - 91.2) > zmassCut) continue;
                if(goodJets.size() < njetsCut) continue;

                int nbtags = 0;
                float minRLeadingJet = 9999.0;
                std::sort(goodJets.begin(), goodJets.end(), jetCompare); // sort jets in descending pt
                for(unsigned int iJet = 0; iJet < goodJets.size(); iJet++) {

                    float dR = deltaR(goodJets[iJet].jet, goodJets[0].jet);
                    if(iJet != 0 && dR < minRLeadingJet) minRLeadingJet = dR;

                    // float looseB = 0.244, mediumB = 0.679, tightB = 0.898;
                    float btagval = pfjets_combinedSecondaryVertexBJetTag().at(goodJets[iJet].idx);
                    h1D_btagval_data->Fill(btagval);
                    // fill(h1D_btagval_data, btagval);

                    if(btagval > 0.75) {
                        nbtags++;
                    }
                }

                if(nbtags >= btagCut) continue; // FIXME


                h1D_dummy_data->Fill(5.0);


                // We are now in the region of interest


                for(unsigned int iMu = 0; iMu < goodMus.size(); iMu++) {
                    fill(h1D_leppt_data,goodMus.at(iMu).pt());
                    h1D_lepeta_data->Fill(goodMus.at(iMu).eta());
                    h1D_lepetamu_data->Fill(goodMus.at(iMu).eta());
                }
                for(unsigned int iEl = 0; iEl < goodEls.size(); iEl++) {
                    fill(h1D_leppt_data,goodEls.at(iEl).pt());
                    h1D_lepeta_data->Fill(goodEls.at(iEl).eta());
                    h1D_lepetae_data->Fill(goodEls.at(iEl).eta());
                }

                // FIXME
                if(goodMus.size() == 2 && goodEls.size() == 1) { // mumue
                    mtW = MT(goodEls[0], pfmet_type1cor(), metphi());
                    fill(h1D_mtWmumue_data, mtW);

                } else if(goodMus.size() == 1 && goodEls.size() == 2) { //eemu
                    mtW = MT(goodMus[0], pfmet_type1cor(), metphi());
                    fill(h1D_mtWeemu_data, mtW);

                } else { 
                    if(goodMus.size() == 3) { //mumumu
                        fill(h1D_mtWmumumu_data, mtW);
                    } else { // eee
                        fill(h1D_mtWeee_data, mtW);
                    }
                }

                double mt = MT(leps[0]+leps[1]+leps[2], pfmet_type1cor(), metphi());
                float ptZ = (leps[0]+leps[1]).pt();
                float st = ht + leps[0].pt() + leps[1].pt() + leps[2].pt() + pfmet_type1cor();

                fill(h1D_ptZ_data, ptZ);
                fill(h1D_Wleppt_data,leps[2].pt());
                if(minRLeadingJet < 9000) fill(h1D_minRLeadingJet_data, minRLeadingJet);
                fill(h1D_st_data, st);
                fill(h1D_mt_data,mt);
                fill(h1D_mtW_data, mtW);
                fill(h1D_zmass_data,mass); 
                fill(h1D_njets_data,goodJets.size());
                fill(h1D_ht_data,ht);
                fill(h1D_met_data,met);
                fill(h1D_nbtags_data, nbtags);

                nGoodEvents++;

                if(goodJets.size() < 1) continue;

                float ptj1 = goodJets[0].pt;
                float massZj1 = (goodJets[0].jet + leps[0]+leps[1]).mass();

                fill(h1D_ptj1_data, ptj1);
                fill(h1D_massZj1_data, massZj1);

                if(goodJets.size() < 2) continue;

                float ptjj = (goodJets[0].jet + goodJets[1].jet).pt();
                float massjj = (goodJets[0].jet + goodJets[1].jet).mass();
                float ptj2 = goodJets[1].pt;
                float massZj2 = (goodJets[1].jet + leps[0]+leps[1]).mass();

                fill(h1D_ptjj_data, ptjj);
                fill(h1D_massjj_data, massjj);
                fill(h1D_ptj2_data, ptj2);
                fill(h1D_massZj2_data, massZj2);


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


    vector<TH1F*> h1D_dummy_vec;

    vector<TH1F*> h1D_njets_vec;
    vector<TH1F*> h1D_ht_vec;
    vector<TH1F*> h1D_met_vec;
    vector<TH1F*> h1D_mt_vec;
    vector<TH1F*> h1D_mtW_vec;
    vector<TH1F*> h1D_zmass_vec;
    vector<TH1F*> h1D_lepeta_vec;
    vector<TH1F*> h1D_leppt_vec;
    vector<TH1F*> h1D_Wleppt_vec;

    // btag+jet stuff
    vector<TH1F*> h1D_nbtags_vec;
    vector<TH1F*> h1D_btagval_vec;
    vector<TH1F*> h1D_ptZ_vec;
    vector<TH1F*> h1D_st_vec;
    vector<TH1F*> h1D_minRLeadingJet_vec;
    vector<TH1F*> h1D_ptj1_vec;
    vector<TH1F*> h1D_ptj2_vec;
    vector<TH1F*> h1D_massZj1_vec;
    vector<TH1F*> h1D_massZj2_vec;
    vector<TH1F*> h1D_ptjj_vec;
    vector<TH1F*> h1D_massjj_vec;

    // w mt debugging
    vector<TH1F*> h1D_mtWeemu_vec;
    vector<TH1F*> h1D_mtWmumue_vec;
    vector<TH1F*> h1D_mtWeee_vec;
    vector<TH1F*> h1D_mtWmumumu_vec;
    vector<TH1F*> h1D_lepetae_vec;
    vector<TH1F*> h1D_lepetamu_vec;

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


        TH1F* h1D_dummy_file = new TH1F("dummy"+filename, "dummyhisto", 10, 0, 10); 

        TH1F* h1D_njets_file = new TH1F("njets"+filename, "Njets;;Entries", 15, 0, 15); 
        TH1F* h1D_ht_file = new TH1F("ht"+filename, "H_{T};GeV;Entries", 20, 0, 600); 
        TH1F* h1D_met_file = new TH1F("met"+filename, "#slash{E}_{T};GeV;Entries", 20, 0, 300); 
        TH1F* h1D_mt_file = new TH1F("mt"+filename, "M_{T};GeV;Entries", 20, 0, 400); 
        TH1F* h1D_mtW_file = new TH1F("mtW"+filename, "W lep M_{T};GeV;Entries", 40, 0, 200); 
        TH1F* h1D_zmass_file = new TH1F("zmass"+filename, "Z Mass;GeV;Entries", 30, 70, 120); 
        TH1F* h1D_lepeta_file = new TH1F("lepeta"+filename, "lepton #eta;#eta;Entries", 50, -3.0, 3.0); 
        TH1F* h1D_leppt_file = new TH1F("leppt"+filename, "lepton p_{T};p_{T} [GeV];Entries", 50, 0, 150); 
        TH1F* h1D_Wleppt_file = new TH1F("Wleppt"+filename, "W lepton p_{T};p_{T} [GeV];Entries", 30, 0, 150); 

        TH1F* h1D_nbtags_file = new TH1F("nbtags"+filename, "N btagged jets;;Entries", 15, 0, 15); 
        TH1F* h1D_btagval_file = new TH1F("btagval"+filename, "Value of csv bjet tag;;Entries", 20, 0, 1); 
        TH1F* h1D_ptZ_file = new TH1F("ptZ"+filename, "Z p_{T};[GeV];Entries", 25, 0, 400);
        TH1F* h1D_st_file = new TH1F("st"+filename, "S_{T}=H_{T}+#Sigma p_{T,leps}+#slash{E}_{T};S_{T} [GeV];Entries", 25, 100, 800);
        TH1F* h1D_minRLeadingJet_file = new TH1F("minRLeadingJet"+filename, "Minimum dR between jet 1 and another jet",10,0,5);
        TH1F* h1D_ptj1_file = new TH1F("ptj1"+filename, "j1 p_{T};[GeV];Entries", 25, 0, 300);
        TH1F* h1D_ptj2_file = new TH1F("ptj2"+filename, "j2 p_{T};[GeV];Entries", 25, 0, 300);
        TH1F* h1D_massZj1_file = new TH1F("massZj1"+filename, "Zj1 mass;[GeV];Entries", 25, 100, 600);
        TH1F* h1D_massZj2_file = new TH1F("massZj2"+filename, "Zj2 mass;[GeV];Entries", 25, 100, 600);
        TH1F* h1D_ptjj_file = new TH1F("ptjj"+filename, "j1j2 p_{T};[GeV];Entries", 25, 0, 400);
        TH1F* h1D_massjj_file = new TH1F("massjj"+filename, "j1j2 mass;[GeV];Entries", 25, 0, 600);

        TH1F* h1D_mtWeemu_file = new TH1F("mtWeem"+filename, "ee#mu W lep M_{T};GeV;Entries", 40, 0, 200); 
        TH1F* h1D_mtWmumue_file = new TH1F("mtWmme"+filename, "#mu#mue W lep M_{T};GeV;Entries", 40, 0, 200); 
        TH1F* h1D_mtWeee_file = new TH1F("mtWeee"+filename, "eee W lep M_{T};GeV;Entries", 40, 0, 200); 
        TH1F* h1D_mtWmumumu_file = new TH1F("mtWmmm"+filename, "#mu#mu#mu W lep M_{T};GeV;Entries", 40, 0, 200); 
        TH1F* h1D_lepetae_file = new TH1F("lepetae"+filename, "e #eta;#eta;Entries", 50, -3.0, 3.0); 
        TH1F* h1D_lepetamu_file = new TH1F("lepetamu"+filename, "#mu #eta;#eta;Entries", 50, -3.0, 3.0); 

        h1D_dummy_vec.push_back(h1D_dummy_file); 

        h1D_njets_vec.push_back(h1D_njets_file); 
        h1D_ht_vec.push_back(h1D_ht_file); 
        h1D_met_vec.push_back(h1D_met_file); 
        h1D_mt_vec.push_back(h1D_mt_file); 
        h1D_mtW_vec.push_back(h1D_mtW_file); 
        h1D_zmass_vec.push_back(h1D_zmass_file); 
        h1D_lepeta_vec.push_back(h1D_lepeta_file);
        h1D_leppt_vec.push_back(h1D_leppt_file); 
        h1D_Wleppt_vec.push_back(h1D_Wleppt_file); 

        h1D_nbtags_vec.push_back(h1D_nbtags_file);
        h1D_btagval_vec.push_back(h1D_btagval_file);
        h1D_ptZ_vec.push_back(h1D_ptZ_file);
        h1D_st_vec.push_back(h1D_st_file);
        h1D_minRLeadingJet_vec.push_back(h1D_minRLeadingJet_file);
        h1D_ptj1_vec.push_back(h1D_ptj1_file);
        h1D_ptj2_vec.push_back(h1D_ptj2_file);
        h1D_massZj1_vec.push_back(h1D_massZj1_file);
        h1D_massZj2_vec.push_back(h1D_massZj2_file);
        h1D_ptjj_vec.push_back(h1D_ptjj_file);
        h1D_massjj_vec.push_back(h1D_massjj_file);

        h1D_mtWeemu_vec.push_back(h1D_mtWeemu_file);
        h1D_mtWmumue_vec.push_back(h1D_mtWmumue_file);
        h1D_mtWeee_vec.push_back(h1D_mtWeee_file);
        h1D_mtWmumumu_vec.push_back(h1D_mtWmumumu_file);
        h1D_lepetae_vec.push_back(h1D_lepetae_file);
        h1D_lepetamu_vec.push_back(h1D_lepetamu_file);

        // Loop over Events in current file
        unsigned int nEventsTree = tree->GetEntriesFast();
        for( unsigned int event = 0; event < nEventsTree; ++event) {

            // if(event > 30000) break;

            // Get Event Content
            cms2.GetEntry(event);
            nEventsSoFar++;

            if(event == 0) {
                std::cout << " evt_scale1fb(): " << evt_scale1fb() << " filename: " << filename << std::endl;
            }
            // Progress
            CMS2::progress( nEventsSoFar, nEventsTotal );

            std::vector<LorentzVector> goodEls;
            std::vector<LorentzVector> goodMus;
            std::vector<JetStruct> goodJets;
            std::map<int, int> goodToP4MapEl; // map indices in good{Els,Mus} to {els,mus}_p4 indices
            std::map<int, int> goodToP4MapMu;
            float ht = 0, met = 0;

            float scale = evt_scale1fb()  * luminosity;

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


                JetStruct myJet = {};
                myJet.jet = pfjets_p4().at(iJet);
                myJet.pt = pfjets_p4().at(iJet).pt()*pfjets_corL1FastL2L3().at(iJet);
                myJet.idx = iJet;

                goodJets.push_back(myJet);
            }



            vector<int> pair = findZPair(goodEls, goodMus);
            vector<LorentzVector> leps;
            float mtW = 0.0, mass = 0.0;
            if( pair[0] != -1 && pair[1] != -1 && pair[2] != -1 && 
                    pair[3] != -1 && pair[4] != -1 ) {
                if(pair[0] == 0) { // el
                    mass += (goodEls.at(pair[1])+goodEls.at(pair[2])).mass();
                    leps.push_back(goodEls.at(pair[1]));
                    leps.push_back(goodEls.at(pair[2]));
                } else { // mu
                    mass += (goodMus.at(pair[1])+goodMus.at(pair[2])).mass();
                    leps.push_back(goodMus.at(pair[1]));
                    leps.push_back(goodMus.at(pair[2]));
                }

                if(pair[3] == 0) { // W lep is el
                    if( !tightEl().at(goodToP4MapEl[pair[4]]) ) continue;
                    leps.push_back(goodEls.at(pair[4]));

                    mtW = MT(els_p4().at(pair[4]), pfmet_type1cor(), metphi());

                } else { // W lep is mu
                    if( !tightMu().at(goodToP4MapMu[pair[4]]) ) continue;
                    leps.push_back(goodMus.at(pair[4]));

                    mtW = MT(goodMus.at(pair[4]), pfmet_type1cor(), metphi());

                }
            } else {
                cout << "shouldn't end up with pair[i] == -1" << endl;
            }


            // XXX
            if(abs(mass - 91.2) > zmassCut) continue;
            if(goodJets.size() < njetsCut) continue;


            int nbtags = 0;
            float minRLeadingJet = 9999.0;
            std::sort(goodJets.begin(), goodJets.end(), jetCompare); // sort jets in descending pt
            for(unsigned int iJet = 0; iJet < goodJets.size(); iJet++) {

                float dR = deltaR(goodJets[iJet].jet, goodJets[0].jet);
                if(iJet != 0 && dR < minRLeadingJet) minRLeadingJet = dR;

                // float looseB = 0.244, mediumB = 0.679, tightB = 0.898;
                float btagval = pfjets_combinedSecondaryVertexBJetTag().at(goodJets[iJet].idx);
                h1D_btagval_file->Fill(btagval, scale);
                // fill(h1D_btagval_file, btagval, scale);

                if(btagval > 0.75) {
                    nbtags++;
                }
            }

            if(nbtags >= btagCut) continue; // FIXME


            h1D_dummy_file->Fill(5.0,scale);

            // We are now in the region of interest


            for(unsigned int iMu = 0; iMu < goodMus.size(); iMu++) {
                fill(h1D_leppt_file,goodMus.at(iMu).pt(), scale);
                h1D_lepeta_file->Fill(goodMus.at(iMu).eta(), scale);
                h1D_lepetamu_file->Fill(goodMus.at(iMu).eta(), scale);
            }
            for(unsigned int iEl = 0; iEl < goodEls.size(); iEl++) {
                fill(h1D_leppt_file,goodEls.at(iEl).pt(), scale);
                h1D_lepeta_file->Fill(goodEls.at(iEl).eta(), scale);
                h1D_lepetae_file->Fill(goodEls.at(iEl).eta(), scale);
            }

            // FIXME
            if(goodMus.size() == 2 && goodEls.size() == 1) { // mumue
                mtW = MT(goodEls[0], pfmet_type1cor(), metphi());
                fill(h1D_mtWmumue_file, mtW, scale);

            } else if(goodMus.size() == 1 && goodEls.size() == 2) { //eemu
                mtW = MT(goodMus[0], pfmet_type1cor(), metphi());
                fill(h1D_mtWeemu_file, mtW, scale);

            } else { 
                if(goodMus.size() == 3) { //mumumu
                    fill(h1D_mtWmumumu_file, mtW, scale);
                } else { // eee
                    fill(h1D_mtWeee_file, mtW, scale);
                }
            }

            double mt = MT(leps[0]+leps[1]+leps[2], pfmet_type1cor(), metphi());
            float ptZ = (leps[0]+leps[1]).pt();
            float st = ht + leps[0].pt() + leps[1].pt() + leps[2].pt() + pfmet_type1cor();

            fill(h1D_ptZ_file, ptZ, scale);
            fill(h1D_Wleppt_file,leps[2].pt(), scale);
            if(minRLeadingJet < 9000) fill(h1D_minRLeadingJet_file, minRLeadingJet, scale);
            fill(h1D_st_file, st, scale);
            fill(h1D_mt_file,mt, scale);
            fill(h1D_mtW_file, mtW, scale);
            fill(h1D_zmass_file,mass, scale); 
            fill(h1D_njets_file,goodJets.size(), scale);
            fill(h1D_ht_file,ht, scale);
            fill(h1D_met_file,met, scale);
            fill(h1D_nbtags_file, nbtags, scale);
            fill(error,0.0,scale);

            nGoodEvents++;
            nGoodEventsWeighted+=scale*1.0;

            if(goodJets.size() < 1) continue;

            float ptj1 = goodJets[0].pt;
            float massZj1 = (goodJets[0].jet + leps[0]+leps[1]).mass();

            fill(h1D_ptj1_file, ptj1, scale);
            fill(h1D_massZj1_file, massZj1, scale);

            if(goodJets.size() < 2) continue;

            float ptjj = (goodJets[0].jet + goodJets[1].jet).pt();
            float massjj = (goodJets[0].jet + goodJets[1].jet).mass();
            float ptj2 = goodJets[1].pt;
            float massZj2 = (goodJets[1].jet + leps[0]+leps[1]).mass();

            fill(h1D_ptjj_file, ptjj, scale);
            fill(h1D_massjj_file, massjj, scale);
            fill(h1D_ptj2_file, ptj2, scale);
            fill(h1D_massZj2_file, massZj2, scale);

        }//event loop

    }//file loop MCMC

    std::cout << " nGoodEvents: " << nGoodEvents << " nEventsTotal: " << nEventsTotal << std::endl;
    std::cout << " nGoodEventsWeighted to 19.4 1/fb: " << nGoodEventsWeighted << std::endl;

    std::cout << " error->GetBinContent(1): " << error->GetBinContent(1) << " error->GetBinError(1): " << error->GetBinError(1) << std::endl;

    TString prefix("plots");
    prefix += tag;
    prefix += "/";

    stringstream ss; ss << njetsCut;
    std::string common = " --luminosity 19.4 --percentages --scaletodata --label njets #geq "+ss.str();
    // std::string common = " --luminosity 19.4 --percentages --label njets>="+ss.str();

    float mcScale = -1.0;
    drawStacked(h1D_dummy_data, h1D_dummy_vec,prefix+"h1D_dummy.pdf",""+common, -1.0, &mcScale);
    
    drawStacked(h1D_njets_data, h1D_njets_vec, prefix+"h1D_njets.pdf", "--centerlabel"+common, manualScale);
    drawStacked(h1D_ht_data, h1D_ht_vec,prefix+"h1D_ht.pdf","--logscale --binsize --centerlabel"+common, manualScale);
    drawStacked(h1D_leppt_data, h1D_leppt_vec,prefix+"h1D_leppt.pdf","--logscale --centerlabel"+common, manualScale);
    drawStacked(h1D_Wleppt_data, h1D_Wleppt_vec,prefix+"h1D_Wleppt.pdf","--logscale --centerlabel"+common, manualScale);
    drawStacked(h1D_met_data, h1D_met_vec,prefix+"h1D_met.pdf","--binsize --centerlabel"+common, manualScale);
    drawStacked(h1D_mt_data, h1D_mt_vec,prefix+"h1D_mt.pdf","--binsize --centerlabel"+common, manualScale);
    drawStacked(h1D_mtW_data, h1D_mtW_vec,prefix+"h1D_mtW.pdf","--binsize "+common, manualScale);
    drawStacked(h1D_lepeta_data, h1D_lepeta_vec,prefix+"h1D_lepeta.pdf",""+common, manualScale);
    drawStacked(h1D_zmass_data, h1D_zmass_vec,prefix+"h1D_zmass.pdf","--binsize"+common, manualScale);

    drawStacked(h1D_mtWeemu_data, h1D_mtWeemu_vec,prefix+"h1D_mtWeemu.pdf","--binsize "+common, manualScale);
    drawStacked(h1D_mtWmumue_data, h1D_mtWmumue_vec,prefix+"h1D_mtWmumue.pdf","--binsize "+common, manualScale);
    drawStacked(h1D_mtWmumumu_data, h1D_mtWmumumu_vec,prefix+"h1D_mtWmumumu.pdf","--binsize "+common, manualScale);
    drawStacked(h1D_mtWeee_data, h1D_mtWeee_vec,prefix+"h1D_mtWeee.pdf","--binsize "+common, manualScale);
    drawStacked(h1D_lepetae_data, h1D_lepetae_vec,prefix+"h1D_lepetae.pdf",""+common, manualScale);
    drawStacked(h1D_lepetamu_data, h1D_lepetamu_vec,prefix+"h1D_lepetamu.pdf",""+common, manualScale);

    drawStacked(h1D_btagval_data, h1D_btagval_vec, prefix+"h1D_btagval.pdf", "--centerlabel "+common, manualScale);
    drawStacked(h1D_nbtags_data, h1D_nbtags_vec, prefix+"h1D_nbtags.pdf", "--centerlabel --logscale"+common, manualScale);
    
    drawStacked(h1D_ptZ_data, h1D_ptZ_vec, prefix+"h1D_ptZ.pdf","--centerlabel "+common, manualScale);
    drawStacked(h1D_st_data, h1D_st_vec, prefix+"h1D_st.pdf","--centerlabel "+common, manualScale);
    drawStacked(h1D_minRLeadingJet_data, h1D_minRLeadingJet_vec, prefix+"h1D_minRLeadingJet.pdf","--centerlabel "+common, manualScale);
    drawStacked(h1D_ptj1_data, h1D_ptj1_vec, prefix+"h1D_ptj1.pdf","--centerlabel "+common, manualScale);
    drawStacked(h1D_ptj2_data, h1D_ptj2_vec, prefix+"h1D_ptj2.pdf","--centerlabel "+common, manualScale);
    drawStacked(h1D_massZj1_data, h1D_massZj1_vec, prefix+"h1D_massZj1.pdf","--centerlabel "+common, manualScale);
    drawStacked(h1D_massZj2_data, h1D_massZj2_vec, prefix+"h1D_massZj2.pdf","--centerlabel "+common, manualScale);
    drawStacked(h1D_ptjj_data, h1D_ptjj_vec, prefix+"h1D_ptjj.pdf","--centerlabel "+common, manualScale);
    drawStacked(h1D_massjj_data, h1D_massjj_vec, prefix+"h1D_massjj.pdf","--centerlabel "+common, manualScale);

    return mcScale;
}
