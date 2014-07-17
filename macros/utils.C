#include "TH1.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TString.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TPavesText.h"
#include "Math/VectorUtil.h"

#include <math.h>
#include <stdlib.h>
#include <vector>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std;

double deltaR(const LorentzVector& p1, const LorentzVector& p2) {
    return ROOT::Math::VectorUtil::DeltaR(p1,p2);
}

double deltaPhi(LorentzVector p1, LorentzVector p2) {
    double dPhi = fabs(p1.Phi() - p2.Phi());
    return dPhi > M_PI ? dPhi : dPhi - M_PI;
}

double deltaPhi(LorentzVector p1, float phi2) {
    double dPhi = fabs(p1.Phi() - phi2);
    return dPhi > M_PI ? dPhi : dPhi - M_PI;
}

double MT(LorentzVector p1, float met, float metphi) {
    double dPhi = deltaPhi(p1, metphi);
    // 43.61 of pdg.lbl.gov/2013/reviews/rpp2012-rev-kinematics.pdf
    return sqrt( 2.0 * p1.Pt() * met * (1.0-cos(dPhi)) );
}

void drawHist(TH1F* hist, TString filename) {
    TCanvas* c0 = new TCanvas();
    hist->Draw();
    c0->Print(filename);
}

bool integralCompare(TH1F* h1, TH1F* h2) {
    return h1->Integral(0,h1->GetNbinsX()) < h2->Integral(0,h2->GetNbinsX());
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
    colors.push_back(kMagenta-5);
    colors.push_back(kCyan-3);
    colors.push_back(kOrange-2);
    colors.push_back(kRed-7);
    colors.push_back(kGreen-2);
    colors.push_back(kYellow-7);
    colors.push_back(kBlue-7);
    colors.push_back(kMagenta-4);
    colors.push_back(kRed+1);
    colors.push_back(kAzure+1);
    colors.push_back(kOrange-8);
    colors.push_back(kGray);
    return colors;
}

int drawStacked(TH1F* data, vector <TH1F*> hists, TString filename, std::string options = "", vector<string> titles = vector<string>()) {

    if(hists.size() < 1) return 1;

    myStyle();

    // default option variables
    TString title = hists[0]->GetTitle();
    TString xlabel = hists[0]->GetXaxis()->GetTitle();
    TString ylabel = hists[0]->GetYaxis()->GetTitle();
    bool logScale = false;
    bool percentages = false;
    bool haveData = (data->GetEntries() > 0);
    // bool haveData = false;
    double luminosity = 0.0;

    TString s(options);
    TPMERegexp re("--"), reSpace(" ");
    re.Split(s);
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
        if(key == "luminosity") luminosity = val.Atof();

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


    //c0->Divide(1,2);
    //c0->cd(1);

    //  xy origin starts at bot left 
    // x1  y1  x2  y2
    //TLegend* leg = new TLegend(0.7,0.55,0.90,0.88,"","NDC");
    TLegend* leg = new TLegend(0.7,0.55,0.90,0.88);
    // NDC methods don't work later on unless I set them (?)
    // even if I use NDC as a draw option
    //leg->SetX1NDC(leg->GetX1()); leg->SetY1NDC(leg->GetY1());
    //leg->SetX2NDC(leg->GetX2()); leg->SetY2NDC(leg->GetY2());
    stack->Draw();

    std::vector<int> colors = getColors();

    // "attach" a color to each histogram so they get sorted too
    // this provides visual consistency among the graphs
    // maybe make this optional, though
    // remove the SetFillColor call in the second hists loop

    // for(unsigned int ih = 0; ih < hists.size(); ih++) {
    //     hists[ih]->SetFillColor(colors[ih]);
    // }

    // bigger histograms on top for log scale, but smaller for linear scale
    if(logScale) std::sort(hists.begin(), hists.end(), integralCompare);
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
        if(ih < colors.size())
            hists[ih]->SetFillColor(colors[ih]);

        hists[ih]->SetLineColor(kBlack);

        if(customTitles)
            cleanName = titles[ih];

        leg->AddEntry(hists[ih],cleanName,"f");
        stack->Add(hists[ih]);

        integral += hists[ih]->Integral(0,hists[ih]->GetNbinsX()+1);
    }

    //set option variables
    stack->SetTitle(title);
    stack->GetXaxis()->SetTitle(xlabel);
    stack->GetYaxis()->SetTitle(ylabel);


    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->Draw();

    float labelYOffset = 0.0, labelDY = 0.04;

    drawLabel( 0.17,0.89-labelYOffset, Form("%.1f events (MC)", integral) );
    labelYOffset += labelDY;

    if(haveData) {
        drawLabel( 0.17,0.89-labelYOffset, Form("%.1f events (Data)", data->Integral(0,data->GetNbinsX()+1) ) );
        labelYOffset += labelDY;
    }

    if(luminosity > 0) {
        drawLabel( 0.17,0.89-labelYOffset, Form("%.2f fb^{-1}", luminosity) );
        labelYOffset += labelDY;
    }

    if(percentages) { 
        float dy = (leg->GetY2()-leg->GetY1())/leg->GetNRows();
        float x1 = leg->GetX1()+0.041;
        float y2 = leg->GetY2()-0.30*dy;

        if(hists.size() > 8) dy *= 1.1;

        for(unsigned int irow = 0; irow < hists.size(); irow++) {
            float percentage = 100.0*hists[irow]->Integral(0,hists[irow]->GetNbinsX()+1)/integral;
            drawLabel(x1,y2-irow*dy*padYScale,Form("%.0f%%",percentage),0.023,33);
        }
    }


    if(haveData) {

        data->UseCurrentStyle();
        data->SetMarkerStyle(20);
        data->SetLineColor(kBlack);
        data->SetMarkerColor(kBlack);
        data->SetMarkerSize(0.7);
        data->Draw("E1 SAME");
        leg->AddEntry(data, "Data","p");

        pBot->cd();

        // +filename to make the "Replacing histogram, possible memory leak" warning go away
        TH1F* dummy = new TH1F("dummy"+filename,"dummy",data->GetNbinsX(),data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());
        TH1F* dummyDenominator = new TH1F("dummyDenominator"+filename,"dummy",data->GetNbinsX(),data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());

        dummy->Sumw2();
        data->Sumw2();

        gStyle->SetOptStat(0);

        dummy->SetTitle("");
        dummy->GetYaxis()->SetLabelSize(0.13);
        dummy->GetXaxis()->SetLabelSize(0.13);
        // dummy->GetYaxis()->SetTitle("#frac{data}{MC}");
        //dummy->GetYaxis()->SetTitle("data/MC");
        dummy->GetYaxis()->SetTitleSize(0.08);

        TLatex * label = new TLatex(0.075,0.37,"data/MC");
        label->SetNDC();
        label->SetTextAngle(90);
        label->SetTextAlign(13);
        label->SetTextSize(0.13);


        dummy->Add(data);
        for(unsigned int ih = 0; ih < hists.size(); ih++) {
            dummyDenominator->Add(hists[ih]);
        }
        dummy->Divide(dummyDenominator);

        dummy->SetMarkerStyle(20);
        dummy->SetLineColor(kBlack);
        dummy->SetLineWidth(1);
        dummy->SetMarkerColor(kBlack);
        dummy->SetMarkerSize(0.7);

        //dummy->GetYaxis()->SetRangeUser(0.6,1.5);
        dummy->Draw("E1");

        TLine *lineY1 = new TLine(data->GetXaxis()->GetXmin(),1,data->GetXaxis()->GetXmax(),1);
        lineY1->SetLineColor(kRed);
        lineY1->Draw();

        label->Draw();

    }

    if(logScale) pTop->SetLogy(1);
    c0->Print(filename);
    if(logScale) pTop->SetLogy(0);

    return 0;
}

