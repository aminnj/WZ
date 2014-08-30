#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TAttFill.h"
#include "TString.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TPavesText.h"
#include "Math/VectorUtil.h"

#include "TROOT.h" 

#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>

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
void addToCounter(TString name, double weight=1.0) {
    evtCounter->Fill(0.5, weight);
    if(evtMap.find(name) == evtMap.end() ) evtMap[name] = weight;
    else evtMap[name] += weight;
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
    double dPhi = fabs( phi1 - phi2 );
    if( dPhi > M_PI ) dPhi = 2.0*M_PI - dPhi;
    return dPhi;

    // double dPhi = phi1 - phi2;
    // while(dPhi > M_PI) dPhi -= 2.0*M_PI;
    // while(dPhi <= -M_PI) dPhi += 2.0*M_PI;
    // return fabs(dPhi);
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
     * and which came from a W. (note, sum of sizes of els, mus vectors must be 3!)
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

////////////////////////////
////// root utilities //////
////////////////////////////
double getIntegral(TH1F* h, bool includeOverflow = true) {
    return h->Integral(0,h->GetNbinsX()+includeOverflow);
}

double getIntegralBetween(TH1F* h, float xmin, float xmax) {
    // http://root.cern.ch/root/roottalk/roottalk03/2211.html
    TAxis *axis = h->GetXaxis();
    int bmin = axis->FindBin(xmin); // built in safeguard to fix
    int bmax = axis->FindBin(xmax); // going out of histogram's range :)
    double integral = h->Integral(bmin,bmax);
    integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/ axis->GetBinWidth(bmin);
    integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/ axis->GetBinWidth(bmax);
    return integral;
}

double getFractionBetween(TH1F* h, float xmin, float xmax) {
    double integral = getIntegral(h);
    if(integral < 0.001) return 0.0;
    return getIntegralBetween(h, xmin, xmax) / getIntegral(h);
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
void fill2D(TH2F* hist, double valuex, double valuey, double scale=1.0) {
    double maximumx = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());
    double maximumy = hist->GetYaxis()->GetBinCenter(hist->GetNbinsY());
    hist->Fill(min(maximumx,valuex), min(maximumy,valuey),scale);
}

////////////////////////////////
////// plotting utilities //////
////////////////////////////////
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

void myStyle2D() {
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.10);
    gStyle->SetPadLeftMargin(0.10);
    gStyle->SetPadRightMargin(0.15);

    gStyle->SetOptStat(0);
    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleSize(0.18, "XYZ");

    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");
}

void setGradient() {
    // http://ultrahigh.org/2007/08/making-pretty-root-color-palettes/
    const int NRGBs = 5;
    const int NCont = 255;
    double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    double red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    double green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    double blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void drawLabel(float x1, float y1, TString s, float size=0.04, int align=13) {
    TLatex * label = new TLatex(x1,y1,s);
    label->SetNDC();
    label->SetTextAlign(align);
    label->SetTextSize(size);
    label->Draw();
}

// void drawHist(TH1F* hist, TString filename) {
//     myStyle();
//     gStyle->SetOptStat(0);
//     TCanvas* c0 = new TCanvas();
//     hist->Draw();
//     c0->Print(filename);
// }

vector<int> getColors() {
    vector<int> colors;

    // stole these from Alex and added a few more random colors
    colors.push_back(kMagenta-5);
    colors.push_back(kCyan-3);
    colors.push_back(kOrange-2);
    colors.push_back(kRed-7);
    colors.push_back(kGreen-2);
    colors.push_back(kYellow-3); // originally kYellow-7
    colors.push_back(kBlue-7);
    colors.push_back(kMagenta-9);
    colors.push_back(kCyan-10);
    colors.push_back(kAzure+1);
    colors.push_back(kOrange-8);
    colors.push_back(kGray);

    return colors;
}

int drawStacked(TH1F* data, vector <TH1F*> hists, TString filename, TString options = "") {

    if(hists.size() < 1) return 1;

    myStyle();
    gStyle->SetOptStat(0);

    // default option variables
    TString title = hists[0]->GetTitle();
    TString xlabel = hists[0]->GetXaxis()->GetTitle();
    TString ylabel = hists[0]->GetYaxis()->GetTitle();
    TString drawOptions = "";
    vector<TString> labels;
    vector<TString> titles;
    bool logScale = false;
    bool percentages = false;
    bool haveData = (data->GetEntries() > 0);
    // bool customTitles = titles.size() == hists.size();
    bool centerLabel = false;
    bool scaleToData = false;
    bool reorderStack = false;
    bool printBins = false;
    bool keepOrder = false;
    bool normalize = false;
    bool noFill = false;
    bool noLegend = false;
    double luminosity = 0.0;
    // double transparency = 1.0;


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
        if(key == "keeporder") keepOrder = true;
        if(key == "normalize") normalize = true;
        if(key == "nofill") noFill = true;
        if(key == "nolegend") noLegend = true;
        if(key == "luminosity") luminosity = val.Atof();
        // if(key == "transparency") transparency = val.Atof();
        if(key == "label") labels.push_back(val);
        if(key == "nostack") drawOptions += "nostack";

        if(key == "titles") {
            TPMERegexp rePipe("\\|");
            rePipe.Split(val);
            // cout << val << endl;
            // cout << "###" << rePipe.NMatches() << endl;
            for(int it = 0; it < rePipe.NMatches(); it++) {
                TString ti = rePipe[it].Strip(TString::kBoth);
                // cout << ti << endl;
                titles.push_back(ti);
                // reSpace.Split(param,2);
            }
        }

        if(key == "binsize") {
            float binWidth = hists[0]->GetBinWidth(0);
            if(binWidth < 1) ylabel += Form(" / %.2f GeV", hists[0]->GetBinWidth(0));
            else ylabel += Form(" / %.0f GeV", hists[0]->GetBinWidth(0));
        }
    }


    THStack* stack = new THStack("stack", hists[0]->GetTitle());
    TLegend* leg = new TLegend(0.7,0.55,0.90,0.88);
    TCanvas* c0 = new TCanvas();

    // handle pads
    TPad *pTop, *pBot;
    if(haveData) {
        pTop = new  TPad("pTop","top pad", 0.0, 1.0, 1,0.175);
        pBot = new  TPad("pBot","bottom pad", 0.0, 0.193, 1.0,0.00);
        pTop->Draw();
        pBot->Draw();
    } else {
        pTop = new  TPad("pTop","top pad", 0.0, 1.0, 1,0);
        pBot = new  TPad("pBot","bottom pad", 0.0, 0.0, 0.0, 0.0);
        pTop->Draw();
    }
    pTop->cd();

    // stack->Draw(drawOptions); // If some stuff starts screwing up (segfault) put this back in

    // attach custom titles to histos before we sort so they don't get jumbled!
    // if(customTitles) {
    for(unsigned int ih = 0; ih < hists.size(); ih++) {
        TString temp = ih < titles.size() ? titles[ih] : " ";
        hists[ih]->SetTitle(temp);
    }
    // }

    std::vector<int> colors = getColors();
    // attach colors to histos so scheme is consistent
    for(unsigned int ih = 0; ih < hists.size(); ih++) {


        if(ih < colors.size()) hists[ih]->SetFillColor(colors[ih]);

        if(drawOptions == "nostack") {
            hists[ih]->SetLineWidth(hists[ih]->GetLineWidth()*2);

            hists[ih]->SetFillStyle(3144);

            // // XXX
            // // In root 5.30/18 and beyond, can just do
            // // h1->SetFillColorAlpha(kRed, 0.5);
            // // BAM! Done. So convoluted right now...
            // TColor *col = gROOT->GetColor(colors[ih]);
            // float r = 0, g = 0, b = 0;
            // col->GetRGB(r,g,b);
            // int icol = 1001+ih;
            // TColor *colTrans = new TColor(icol, r, g, b, "", transparency);
            // std::cout << r << " " << g << " " << b << std::endl;
            // hists[ih]->SetFillColor(icol);
            // // XXX

        }
        hists[ih]->SetLineColor(kBlack);

        if(noFill) {
            if(ih < colors.size()) hists[ih]->SetLineColor(colors[ih]);
            hists[ih]->SetFillStyle(0);
            hists[ih]->SetLineWidth(hists[ih]->GetLineWidth()*1.5);
        }


    }
    stack->Draw(drawOptions);

    // bigger histograms on top for log scale, but smaller for linear scale
    if(!keepOrder) {
        if(logScale || reorderStack) std::sort(hists.begin(), hists.end(), integralCompare);
        else std::sort(hists.rbegin(), hists.rend(), integralCompare);
    }

    double MCintegral = 0;
    for(unsigned int ih = 0; ih < hists.size(); ih++) {
        TString name(hists[ih]->GetName());
        TString cleanName = name;
        TRegexp re1(".*baby_"), re2(".root");
        cleanName(re1) = "";
        cleanName(re2) = "";


        if(titles.size() > 0)
            cleanName = hists[ih]->GetTitle();

        leg->AddEntry(hists[ih],cleanName,"f");
        stack->Add(hists[ih]);

        MCintegral += hists[ih]->Integral(0,hists[ih]->GetNbinsX()+1);

        if(normalize) {
            float scale = getIntegral(hists[ih]);
            hists[ih]->Scale(1.0/scale);
        }
    }

    stack->SetTitle(title);
    stack->GetXaxis()->SetTitle(xlabel);
    stack->GetYaxis()->SetTitle(ylabel);



    if(haveData) {
        leg->AddEntry(data, "Data","lpe");
        stack->SetMaximum(max(maxY(data),stack->GetMaximum()));

        if(scaleToData) {
            float scale = 1.0*data->Integral(0,data->GetNbinsX()+1)/MCintegral;

            for(unsigned int ih = 0; ih < hists.size(); ih++) {
                hists[ih]->Scale(scale);
            }
            MCintegral *= scale;

        }
    }

    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    if(!noLegend) leg->Draw();

    // labels
    float labelX = 0.17;
    float labelYOffset = 0.015, labelDY = 0.04;
    if(centerLabel) labelX = 0.42;

    drawLabel( labelX,0.89-labelYOffset, Form("%.1f events (MC)", MCintegral) );
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
        float dy = (leg->GetY2()-leg->GetY1())/(leg->GetNRows());
        float x1 = leg->GetX1()+0.041;
        float y2 = leg->GetY2()-0.30*dy;

        for(unsigned int irow = 0; irow < hists.size(); irow++) {
            float percentage = 100.0*hists[irow]->Integral(0,hists[irow]->GetNbinsX()+1)/MCintegral;
            drawLabel(x1,(y2-irow*dy),Form("%.0f%%",percentage),0.023,33);
        }
    }

    if(haveData) {

        data->UseCurrentStyle();
        data->SetMarkerStyle(20);
        data->SetLineColor(kBlack);
        data->SetMarkerColor(kBlack);
        data->SetMarkerSize(0.7);
        data->Draw("E1 SAME");

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

        if(printBins) {
            std::cout << "Printing bin contents for " << title << endl;
            for(int ib = 0; ib < data->GetNbinsX(); ib++) {
                std::cout << "\tbin: " << ib << "\tdata: " << data->GetBinContent(ib) << "\tmcsum: " << mcSum->GetBinContent(ib) << std::endl;
            }
        }

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

int drawHist(TH1F* hist, TString filename, TString options = "") {
    TH1F * dog = new TH1F("","",1,0,1);

    vector<TH1F*> temp;
    temp.push_back(hist);

    drawStacked(dog, temp, filename, options);
    return 0;

}

int drawHist2D(TH2F* hist, TString filename, TString options = "") {

    myStyle2D();
    setGradient();

    TString title = hist->GetTitle();
    TString xlabel = hist->GetXaxis()->GetTitle();
    TString ylabel = hist->GetYaxis()->GetTitle();
    TString drawOptions = "colz";
    bool logScale = false;

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
        if(key == "drawoptions") drawOptions = val;
        if(key == "showstats") gStyle->SetOptStat("ne");
    }

    TCanvas* c0 = new TCanvas();
    hist->SetTitle(title);
    hist->GetXaxis()->SetTitle(xlabel);
    hist->GetYaxis()->SetTitle(ylabel);

    hist->Draw(drawOptions);

    c0->SetLogz(logScale);
    c0->Print(filename);

    return 0;
}

int drawGraph(vector<vector<float> > xvecs, vector<vector<float> > yvecs, TString filename, TString options = "") {

    myStyle();

    bool logScale = false;
    bool centerLabel = false;
    TString title = "";
    TString xlabel = "";
    TString ylabel = "";
    TString legendPosition = "";
    vector<TString> labels;
    vector<TString> titles;

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
        if(key == "centerlabel") centerLabel = true;
        if(key == "label") labels.push_back(val);
        if(key == "legendposition") legendPosition = val;

        if(key == "titles") {
            TPMERegexp rePipe("\\|");
            rePipe.Split(val);
            for(int it = 0; it < rePipe.NMatches(); it++) {
                TString ti = rePipe[it].Strip(TString::kBoth);
                titles.push_back(ti);
            }
        }

    }

    TMultiGraph *mg = new TMultiGraph();
    TCanvas* c0 = new TCanvas();

    float y1 = 0.35, y2 = 0.60;
    if(legendPosition == "bottom") y1 = 0.15, y2 = 0.35;
    if(legendPosition == "top") y1 = 0.63, y2 = 0.88;
    TLegend* leg = new TLegend(0.7,y1,0.90,y2);

    mg->Draw("ACP");

    std::vector<int> colors = getColors();
    for(unsigned int ig = 0; ig < xvecs.size(); ig++) {

        TGraph* graph = new TGraph(xvecs[ig].size(), &xvecs[ig][0], &yvecs[ig][0]);
        graph->SetMarkerColor(colors[ig]);
        graph->SetLineColor(colors[ig]);
        graph->SetMarkerStyle(kFullDotLarge);
        graph->SetMarkerSize(1.5);
        graph->SetLineWidth(2);
        mg->Add(graph);

        leg->AddEntry(graph,ig < titles.size() ? titles[ig] : " ","lp");

    }

    mg->SetTitle(title);
    mg->GetXaxis()->SetTitle(xlabel);
    mg->GetYaxis()->SetTitle(ylabel);

    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->Draw();

    // labels
    float labelX = 0.17;
    float labelYOffset = 0.015, labelDY = 0.04;
    if(centerLabel) labelX = 0.42;

    if(labels.size() > 0) {
        for(unsigned int ilabel = 0; ilabel < labels.size(); ilabel++) {
            drawLabel( labelX,0.89-labelYOffset, labels.at(ilabel) );
            labelYOffset += labelDY;
        }
    }

    if(logScale) c0->SetLogy(1);
    c0->Print(filename);
    if(logScale) c0->SetLogy(0);

    return 0;
}

