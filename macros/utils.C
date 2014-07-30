#include "TH1.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TString.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TPavesText.h"
#include "Math/VectorUtil.h"

#include <math.h>
#include <stdlib.h>
#include <vector>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std;

///////////////////////////
////// counter stuff //////
///////////////////////////
TH1D * evtCounter = new TH1D("","",1,0,1); 
map<TString, double> evtMap;
void initCounter() {
    evtCounter = new TH1D("","",1,0,1);
    evtMap.clear();
}
void addToCounter(TString filename, double weight) {
    evtCounter->Fill(0.5, weight);
    if(evtMap.find(filename) == evtMap.end() ) evtMap[filename] = weight;
    else evtMap[filename] += weight;
}
void printCounter() {
    cout << string(30, '-') << endl << "Counter totals: " << endl;
    for(map<TString,double>::iterator it = evtMap.begin(); it != evtMap.end(); it++)
        cout << "\t" << it->first << "\t" << it->second << endl;

    cout << "Total: " << evtCounter->GetBinContent(1) << " pm " << evtCounter->GetBinError(1) << endl;
    cout << string(30, '-') << endl;
}


/////////////////////////////
////// root math stuff //////
/////////////////////////////
double deltaR(const LorentzVector& p1, const LorentzVector& p2) {
    return ROOT::Math::VectorUtil::DeltaR(p1,p2);
}

double deltaPhi(float phi1, float phi2) {
    double dPhi = phi1 - phi2;
    while(dPhi > M_PI) dPhi -= 2.0*M_PI;
    while(dPhi <= -M_PI) dPhi += 2.0*M_PI;
    return fabs(dPhi);
}

double MT(LorentzVector p1, float met, float metphi) {
    double dPhi = deltaPhi(p1.Phi(), metphi);
    // 43.61 of pdg.lbl.gov/2013/reviews/rpp2012-rev-kinematics.pdf
    return sqrt( 2.0 * p1.Pt() * met * (1.0-cos(dPhi)) );
}

//////////////////////////////////////
////// analysis utilities stuff //////
//////////////////////////////////////
vector<int> findZPair(std::vector<LorentzVector> goodEls, std::vector<LorentzVector> goodMus) {
    /* Given a vector of electrons and a vector of muons, identify which 2 came from a Z
     * and which came from a W.
     *
     * Returns a 5-tuple (v) with this information:
     * v[0] : 0 if Z lept #1 is e, or 1 if mu
     * v[1] : index in input vector for Z lept #1
     * v[2] : index in input vector for Z lept #2
     * v[3] : 0 if W lept #1 is e, or 1 if mu
     * v[4] : index in input vector for W lept
     */
    int nEl = goodEls.size();
    int nMu = goodMus.size();
    float ZMass = 91.1876;
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
            return v;
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
            return v;
        }
    }
    cout << "shouldn't end up at end of zpair function" << endl;
    return v;
}

///////////////////////////////
////// drawing utilities //////
///////////////////////////////
void drawHist(TH1F* hist, TString filename) {
    TCanvas* c0 = new TCanvas();
    hist->Draw();
    c0->Print(filename);
}

bool integralCompare(TH1F* h1, TH1F* h2) {
    return h1->Integral(0,h1->GetNbinsX()) < h2->Integral(0,h2->GetNbinsX());
}

double maxY(TH1F* data) {
    double maxYval = 0.0;
    for(int ib = 0; ib < data->GetNbinsX(); ib++) {
        double val = data->GetBinContent(ib) + data->GetBinError(ib);
        if(val > maxYval) maxYval = val;
    }
    return maxYval;
}

void fill(TH1F* hist, double value, double scale=1.0) {
    double maximum = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());
    hist->Fill(min(maximum,value),scale);
}

void myStyle() {
    // took these from alex since they look nice
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.04);

    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetTitleOffset(1.20, "X");
    gStyle->SetTitleOffset(1.10, "Y");

    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");
}

void drawLabel(float x1, float y1, TString s, float size=0.04, int align=13) {
    TLatex * label = new TLatex(x1,y1,s);
    label->SetNDC();
    label->SetTextAlign(align);
    label->SetTextSize(size);
    label->Draw();
}

vector<int> getColors() {
    vector<int> colors;

    // int icol = 1000;
    // TColor* col;


    // col = new TColor(icol,0.650980392,0.807843137,0.890196078); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.121568627,0.470588235,0.705882353); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.698039216,0.874509804,0.541176471); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.2,0.62745098,0.17254902); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.984313725,0.603921569,0.6); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.890196078,0.101960784,0.109803922); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.992156863,0.749019608,0.435294118); colors.push_back(icol); icol++;
    // col = new TColor(icol,1,0.498039216,0); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.792156863,0.698039216,0.839215686); colors.push_back(icol); icol++;

    // col = new TColor(icol,0.945,0.404,0.271); colors.push_back(icol); icol++;
    // col = new TColor(icol,1.0,0.776,0.365); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.482,0.784,0.643); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.298,0.765,0.851); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.576,0.392,0.553); colors.push_back(icol); icol++;
    // col = new TColor(icol,0.251,0.251,0.251); colors.push_back(icol); icol++;


    // stole these from Alex and added a few more random colors
    colors.push_back(kMagenta-5);
    colors.push_back(kCyan-3);
    colors.push_back(kOrange-2);
    colors.push_back(kRed-7);
    colors.push_back(kGreen-2);
    colors.push_back(kYellow-7);
    colors.push_back(kBlue-7);
    colors.push_back(kMagenta-9);
    colors.push_back(kCyan-10);
    colors.push_back(kAzure+1);
    colors.push_back(kOrange-8);
    colors.push_back(kGray);

    return colors;
}

int drawStacked(TH1F* data, vector <TH1F*> hists, TString filename, TString options = "", vector<string> titles = vector<string>()) {

    if(hists.size() < 1) return 1;

    myStyle();

    // default option variables
    TString title = hists[0]->GetTitle();
    TString xlabel = hists[0]->GetXaxis()->GetTitle();
    TString ylabel = hists[0]->GetYaxis()->GetTitle();
    vector<TString> labels;
    bool logScale = false;
    bool percentages = false;
    bool haveData = (data->GetEntries() > 0);
    bool centerLabel = false;
    bool scaleToData = false;
    bool reorderStack = false;
    bool printBins = false;
    double luminosity = 0.0;

    TPMERegexp re("--"), reSpace(" ");
    re.Split(options);
    for(int im = 0; im < re.NMatches(); im++) {
        TString param = re[im].Strip(TString::kBoth);
        if(param.Length() < 1) continue;

        reSpace.Split(param,2);
        TString key = reSpace[0], val = reSpace[1];

        // change option variables
        if(key == "title") title = val;
        if(key == "xlabel") xlabel = val;
        if(key == "ylabel") ylabel = val;
        if(key == "logscale") logScale = true;
        if(key == "percentages") percentages = true;
        if(key == "centerlabel") centerLabel = true;
        if(key == "scaletodata") scaleToData = true;
        if(key == "reorderstack") reorderStack = true;
        if(key == "printbins") printBins = true;
        if(key == "luminosity") luminosity = val.Atof();
        if(key == "label") labels.push_back(val);

        if(key == "binsize") {
            float binWidth = hists[0]->GetBinWidth(0);
            if(binWidth < 1) ylabel += Form(" / %.2f GeV", hists[0]->GetBinWidth(0));
            else ylabel += Form(" / %.0f GeV", hists[0]->GetBinWidth(0));
        }
    }

    bool customTitles = titles.size() == hists.size();

    THStack* stack = new THStack("stack", hists[0]->GetTitle());
    TCanvas* c0 = new TCanvas();
    TPad *pTop, *pBot;
    float padYScale = 1.0;

    if(haveData) {
        pTop = new  TPad("pTop","top pad", 0.0, 1.0, 1,0.175);
        padYScale = (1.0-0.175);
        pBot = new  TPad("pBot","bottom pad", 0.0, 0.193, 1.0,0.00);
        pTop->Draw();
        pBot->Draw();
    } else {
        pTop = new  TPad("pTop","top pad", 0.0, 1.0, 1,0);
        pBot = new  TPad("pBot","bottom pad", 0.0, 0.0, 0.0, 0.0);
        pTop->Draw();
    }
    pTop->cd();


    TLegend* leg = new TLegend(0.7,0.55,0.90,0.88);
    stack->Draw();

    std::vector<int> colors = getColors();

    // attach custom titles to histos before we sort so they don't
    // get jumbled!
    if(customTitles) {
        for(unsigned int ih = 0; ih < hists.size(); ih++) {
            TString temp = titles[ih];
            hists[ih]->SetTitle(temp);
        }
    }

    // attach colors to histos so scheme is consistent
    for(unsigned int ih = 0; ih < hists.size(); ih++) {
        if(ih < colors.size()) hists[ih]->SetFillColor(colors[ih]);
    }

    // bigger histograms on top for log scale, but smaller for linear scale
    if(logScale || reorderStack) std::sort(hists.begin(), hists.end(), integralCompare);
    else std::sort(hists.rbegin(), hists.rend(), integralCompare);

    double integral = 0;
    for(unsigned int ih = 0; ih < hists.size(); ih++) {
        TString name(hists[ih]->GetName());
        TRegexp re1(".*baby_");
        TRegexp re2(".root");
        TString cleanName = name;
        cleanName(re1) = "";
        cleanName(re2) = "";
        // fill colors: http://root.cern.ch/root/html/TAttFill.html
        // if(ih < colors.size())
        // hists[ih]->SetFillColor(colors[ih]);

        hists[ih]->SetLineColor(kBlack);

        if(customTitles)
            cleanName = hists[ih]->GetTitle();

        leg->AddEntry(hists[ih],cleanName,"f");
        stack->Add(hists[ih]);

        integral += hists[ih]->Integral(0,hists[ih]->GetNbinsX()+1);
    }

    if(haveData && scaleToData) {
        float scale = 1.0*data->Integral(0,data->GetNbinsX()+1)/integral;
        // if(overrideScale > 0) scale = overrideScale;
        // if(getScale != NULL) *getScale = scale; // return scale to user if he wants it

        // cout << "scaling mc to data by " << scale << endl;

        for(unsigned int ih = 0; ih < hists.size(); ih++) {
            hists[ih]->Scale(scale);
        }
        integral *= scale;

    }

    //set option variables
    stack->SetTitle(title);
    stack->GetXaxis()->SetTitle(xlabel);
    stack->GetYaxis()->SetTitle(ylabel);


    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->Draw();

    float labelX = 0.17;
    float labelYOffset = 0.015, labelDY = 0.04;
    if(centerLabel) labelX = 0.42;

    drawLabel( labelX,0.89-labelYOffset, Form("%.1f events (MC)", integral) );
    labelYOffset += labelDY;

    if(haveData) {
        drawLabel( labelX,0.89-labelYOffset, Form("%.0f events (Data)", data->Integral(0,data->GetNbinsX()+1) ) );
        labelYOffset += labelDY;
    }

    if(luminosity > 0) {
        drawLabel( labelX,0.89-labelYOffset, Form("%.2f fb^{-1}", luminosity) );
        labelYOffset += labelDY;
    }

    if(labels.size() > 0) {
        for(unsigned int ilabel = 0; ilabel < labels.size(); ilabel++) {
            drawLabel( labelX,0.89-labelYOffset, labels.at(ilabel) );
            labelYOffset += labelDY;
        }
    }

    if(percentages) { 
        float dy = (leg->GetY2()-leg->GetY1())/leg->GetNRows();
        float x1 = leg->GetX1()+0.041;
        float y2 = leg->GetY2()-0.30*dy;

        if(haveData) {
            if(hists.size() > 8) dy *= 1.1;
            else if(hists.size() > 7) dy *= 1.067;
            else if(hists.size() > 6) dy *= 1.04;
            else if(hists.size() > 5) dy *= 1.02;
        }

        for(unsigned int irow = 0; irow < hists.size(); irow++) {
            float percentage = 100.0*hists[irow]->Integral(0,hists[irow]->GetNbinsX()+1)/integral;
            drawLabel(x1,y2-irow*dy*padYScale,Form("%.0f%%",percentage),0.023,33);
        }
    }


    if(haveData) {

        stack->SetMaximum(max(maxY(data),stack->GetMaximum()));

        data->UseCurrentStyle();
        data->SetMarkerStyle(20);
        data->SetLineColor(kBlack);
        data->SetMarkerColor(kBlack);
        data->SetMarkerSize(0.7);
        data->Draw("E1 SAME");
        leg->AddEntry(data, "Data","p");

        pBot->cd();

        // +filename to make the "Replacing histogram, possible memory leak" warning go away
        TH1F* comparison = new TH1F("comparison"+filename,"comparison",data->GetNbinsX(),data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());
        TH1F* mcSum = new TH1F("mcSum"+filename,"dummy",data->GetNbinsX(),data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());






        // data/MC
        data->Sumw2(); // important
        comparison->Add(data);
        for(unsigned int ih = 0; ih < hists.size(); ih++) {
            mcSum->Add(hists[ih]);
        }
        comparison->Divide(mcSum);

        // XXX both methods give the same vals,errs, so I will use the above
        // Alex's method: (agreeing results proves that root's doing the sum of squares correctly :) )
        // for(int ib = 1; ib < comparison->GetNbinsX()+1; ib++) {
        //     if(data->GetBinContent(ib) < 1) continue;
        //     float mcVal = 0.0;
        //     float mcErr2 = 0.0;
        //     for(unsigned int ih = 0; ih < hists.size(); ih++) {
        //         mcSum->Add(hists[ih]);
        //         mcVal += hists[ih]->GetBinContent(ib);
        //         mcErr2 += pow(hists[ih]->GetBinError(ib), 2);
        //     }
        //     float dataVal = data->GetBinContent(ib);
        //     float dataErr = data->GetBinError(ib);
        //     float mcErr = sqrt(mcErr2);
        //     float val = dataVal/mcVal;
        //     comparison->SetBinContent(ib, val);
        //     comparison->SetBinError(ib, val*sqrt(pow(mcErr/mcVal,2) + pow(dataErr/dataVal,2)));
        // }

        if(printBins) {
            for(int ib = 0; ib < data->GetNbinsX(); ib++) {
                std::cout << "\tbin: " << ib << "\tdata: " << data->GetBinContent(ib) << "\tmcsum: " << mcSum->GetBinContent(ib) << std::endl;
            }
        }


        // // FIXME
        // float sum = 0.0;
        // int ndof = 0;
        // for(int ib = 0; ib < data->GetNbinsX(); ib++) {
        //     if(data->GetBinContent(ib) < 1) continue;
        //     float weightedResidual = fabs(mcSum->GetBinContent(ib) - data->GetBinContent(ib)) / data->GetBinError(ib);
        //     sum += pow(weightedResidual*weightedResidual, 2);
        //     ndof++;
        // }
        // std::cout << "Probability of data, MC agreement: " << 100.0*TMath::Prob(sum/ndof,ndof) << "%" << std::endl;
        // // FIXME




        gStyle->SetOptStat(0);
        comparison->SetTitle("");
        comparison->GetYaxis()->SetRangeUser(0.0, 2.0);
        comparison->GetYaxis()->SetLabelSize(0.13);
        comparison->GetXaxis()->SetLabelSize(0.13);
        comparison->GetYaxis()->SetTitleSize(0.08);
        comparison->SetMarkerStyle(20);
        comparison->SetLineColor(kBlack);
        comparison->SetLineWidth(1);
        comparison->SetMarkerColor(kBlack);
        comparison->SetMarkerSize(0.7);
        comparison->GetYaxis()->SetNdivisions(505);

        comparison->Draw("E1");

        // draw red line at y=1
        TLine *lineY1 = new TLine(data->GetXaxis()->GetXmin(),1,data->GetXaxis()->GetXmax(),1);
        lineY1->SetLineColor(kRed);
        lineY1->Draw();

        TLatex * label = new TLatex(0.075,0.16,"Data/MC");
        label->SetTextFont(42);
        label->SetNDC();
        label->SetTextAngle(90);
        label->SetTextAlign(13);
        label->SetTextSize(0.20);

        label->Draw();

    }

    if(logScale) pTop->SetLogy(1);
    c0->Print(filename);
    if(logScale) pTop->SetLogy(0);

    return 0;
}

