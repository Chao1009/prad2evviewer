#include "gain_factor.h"

#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TString.h>
#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

std::string GetChannelName(const std::string &histName);

int main(int argc, char *argv[])
{   
    std::string run_number = "";
    std::string fileDir = ".";
    int opt;
    while ((opt = getopt(argc, argv, "r:s:e:d:")) != -1) {
        switch (opt) {
            case 'r': run_number    = optarg; break;
            case 'd': fileDir       = optarg; break;
        }
    }
    
    TFile* f = new TFile(Form("%s/prad_%s_LMS.root", fileDir.c_str(), run_number.c_str()), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open ROOT file\n";
        return 1;
    }

    std::vector<TH1F*> hists;
    std::vector<prad2::FitResult> results;
    std::vector<std::string> name;
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        if (std::string(key->GetClassName()) != "TH1F") continue;
        TH1F *h = (TH1F*)key->ReadObj();
        if (!h) continue;
        
        results.push_back(prad2::gain_hist_fitter(h, 0.1f));
        name.push_back(GetChannelName(h->GetName()));
        hists.push_back(h);
    }
    printf("Loaded %zu TH1F histograms\n", hists.size());

    TFile* outRoot = new TFile(Form("%s/prad_%s_LMS_fitted.root", fileDir.c_str(), run_number.c_str()), "RECREATE");
    outRoot->cd();
    
    for (auto h : hists) h->Write();
    
    outRoot->Close();
    f->Close();
    
    //format: first 3 line reference channel name, alpha peak position, alpha sigma, alpha fit chi2/ndf, lms peak position, lms sigma, lms fit chi2/ndf
    //format: the rest: HyCal module name, lms peak, lms sigma, lms fit chi2/ndf, and three gain factors using 3 reference PMT
    
    std::ofstream outDatFile;
    outDatFile.open(Form("%s/prad_%s_LMS.dat", fileDir.c_str(), run_number.c_str()));
    
    
    for (unsigned int i=0; i<3; i++){
        outDatFile<<std::setw(9)<<Form("LMS%d", i+1)
                  <<std::setw(15)<<results[i].mean<<std::setw(15)<<results[i].sigma<<std::setw(15)<<results[i].chi2pndf
                  <<std::setw(15)<<results[i+3].mean<<std::setw(15)<<results[i+3].sigma<<std::setw(15)<<results[i+3].chi2pndf<<std::endl;
    }
    
    for (unsigned int i = 6; i < hists.size(); i++){
        float factor[3] = {0., 0., 0.};
        for (int j = 0; j<3; j++){
            if (results[j].mean > 1. && results[j+3].mean > 1.)
                factor[j] = results[i].mean * results[j].mean / results[j+3].mean;
            if (!std::isfinite(factor[j]))
                factor[j] = 0.f;
        }
        outDatFile<<std::setw(9)<<name[i]
                  <<std::setw(15)<<results[i].mean<<std::setw(15)<<results[i].sigma<<std::setw(15)<<results[i].chi2pndf
                  <<std::setw(15)<<factor[0]<<std::setw(15)<<factor[1]<<std::setw(15)<<factor[2]<<std::endl;
    }
    outDatFile.close();

    return 0;
}

std::string GetChannelName(const std::string &histName)
{
    auto pos = histName.find('_');
    if (pos == std::string::npos)
        return histName;
    return histName.substr(0, pos);
}
