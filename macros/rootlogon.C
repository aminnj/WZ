{
  #include <stdexcept>

    //gStyle->SetOptStat("emroui");
    //gStyle->SetOptStat("");
    gStyle->SetCanvasColor(0);
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameBorderMode(0);


    string msg("* Welcome to ROOT v");
    msg += gROOT->GetVersion();
    msg += " *";
    string ast(msg.size(), '*');
    cout << ast << endl << msg << endl << ast << endl;

    if (gSystem->Getenv("CMSSW_BASE")) {
        cout<<"loading libFWCoreFWLite.so"<<endl;
        gSystem->Load("libFWCoreFWLite.so");
        AutoLibraryLoader::enable();
        gSystem->Load("libDataFormatsFWLite.so");
        gSystem->Load("libDataFormatsPatCandidates.so");
    } 
    //else if(TString(gSystem->Getenv("ROOTSYS"))=="/code/osgcode/UCSD_root/root_v5.28.00"){
    //   cout<<"loading libminiFWLite from /home/users/cwelke/macros/libMiniFWLiteROOTVer5_28_00.so"<<endl;
    //       gSystem->Load("/home/users/cwelke/macros/libMiniFWLiteROOTVer5_28_00.so");
    //   // AutoLibraryLoader::enable();
    // } 
    else {
        cout<<"loading libminiFWLite from /home/users/cgeorge/macros/MiniFWLite/libMiniFWLite.so"<<endl;
        //gSystem->Load("/home/users/jgran/macros/MiniFWLite/libMiniFWLite.so");
        gSystem->Load("/home/users/cgeorge/macros/MiniFWLite/libMiniFWLite.so");
        // AutoLibraryLoader::enable();
        //       cout<<"Warning: CMSSW_BASE not set. Could not load libFWCoreFWLite.so. ROOTSYS not recognized."<<endl;
        //        cout<<"ROOTSYS = "<<gSystem->Getenv("ROOTSYS")<<endl;
        // AutoLibraryLoader::enable();
    }

Root.ErrorIgnoreLevel: Warning;

}
