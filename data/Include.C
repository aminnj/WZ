// C++ includes
#include <string>
#include <set>
#include <vector>
#include <iomanip>

// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

using namespace std;

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int evt_run, evt_event, evt_lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (evt_run != other.evt_run)
    return evt_run < other.evt_run;
  if (evt_event != other.evt_event)
    return evt_event < other.evt_event;
  if(evt_lumi != other.evt_lumi)
    return evt_lumi < other.evt_lumi;
  return false;
}

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (evt_run != other.evt_run)
    return false;
  if (evt_event != other.evt_event)
    return false;
  return true;
}

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

