/*
 * Copyright (C) 2025 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2025
 *
 */
 
#include "DigitizationRunner.h"
#include "TrackProcessor.h"
#include "Globals.h"
#include "Utils.h"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <limits.h>
#include <utility>
#include <stdexcept>
#include <regex>
#include <algorithm>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "PMTHits.h"
#include "PMTSignal.h"

using namespace std;

DigitizationRunner::DigitizationRunner(const std::string& configFile_,
                                   const std::string& inputDir_,
                                   const std::string& outputDir_)
    : configFile(configFile_),
      inputDir(inputDir_),
      outputDir(outputDir_)
{
    
    
    if (!config.loadConfig(configFile)) {
        cerr << "Failed to load configuration file: " << configFile << endl;
        exit(EXIT_FAILURE);
    }
    

    config.validateAxisMappings();

    

    runCount = config.getInt("start_run_number");
    

    // DEBUG
    //config.printConfig();

    initializeGlobals();
    prepareCameraSettings();
    checkDimensionConsistency();
    initSourceDir();
    
    num_threads = config.getInt("Parallel_threads");

    // DEBUG
    // cout<<"Initialized globals:"<<endl;
    // cout<<"GEM1_gain = "<< GEM1_gain <<" - GEM2_gain = " << GEM2_gain<<" - GEM3_gain = " << GEM3_gain << endl;
    // cout<<"extraction_eff_GEM1 = "<< extraction_eff_GEM1 << " - extraction_eff_GEM2 = "
    //     << extraction_eff_GEM2 << " - extraction_eff_GEM3 = " << extraction_eff_GEM3 << endl;
    // cout<<"omega = "<< omega <<endl;
    // cout<<"x_pix, y_pix = "<<x_pix<<", "<<y_pix<<endl;
    // cout<<"optcounts_per_photon = "<<optcounts_per_photon<<endl;
    // cout<<"y_sensor_size = "<<y_sensor_size<<endl;
    // cout<<"readout_time = "<<readout_time<<endl;
    // cout<<"num_threads = "<<num_threads<<endl;
    
}

void DigitizationRunner::run() {
    auto t0 = std::chrono::steady_clock::now();

    setSeed(3);
    
    processRootFiles();

    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> dur = t1 - t0;
    cout << "Total digitization time: " << dur.count() << " seconds" << endl;
}

void DigitizationRunner::runPedsOnly() {
    auto t0 = std::chrono::steady_clock::now();

    setSeed(3);

    generateHistogramsFromDigi();

    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> dur = t1 - t0;
    cout << "Total joinPeds time: " << dur.count() << " seconds" << endl;
}

void DigitizationRunner::initializeGlobals() {
    double c_G = config.getDouble("c_G");
    double alpha_G = config.getDouble("alpha_G");

    GEM1_gain = c_G * exp(alpha_G * config.getDouble("GEM1_HV"));
    GEM2_gain = c_G * exp(alpha_G * config.getDouble("GEM2_HV"));
    GEM3_gain = c_G * exp(alpha_G * config.getDouble("GEM3_HV"));

    extraction_eff_GEM1 = 0.87319885 * exp(-0.002 * config.getDouble("GEM1_HV"));
    extraction_eff_GEM2 = 0.87319885 * exp(-0.002 * config.getDouble("GEM2_HV"));
    extraction_eff_GEM3 = 0.87319885 * exp(-0.002 * config.getDouble("GEM3_HV"));
}

void DigitizationRunner::prepareCameraSettings() {
    std::string camera = config.get("Camera_type");

    y_pix = 2304;
    if (camera == "Fusion") {
        x_pix = 2304;
        optcounts_per_photon = 1./0.21; // (Fusion sheet) counts/photon = 1. / e–/count (0.21)
        camera_quantum_efficiency = 0.8;
        camera_electron_rms = 0.7; // electrons in ultra quiet scan (UQS) [camera datasheet]
        y_sensor_size = 14.976;    // mm
        readout_time  = 184.4;     // ms in in ultra quiet scan (UQS)
    } else if (camera == "Quest1" || camera == "Quest2") {
        x_pix = 4096;
        optcounts_per_photon = 1./0.107; // (Quest1 and Quest2 sheet) counts/photon = 1. / e–/count (0.107)
        camera_quantum_efficiency = 0.85;
        camera_electron_rms = (camera == "Quest1") ? 0.27 : 0.30; // electrons in ultra quiet scan (UQS) [camera datasheet]
        y_sensor_size = 10.598;      // mm
        readout_time = (camera == "Quest1") ? 199.0 : 39.0; // ms in in ultra quiet scan (UQS)
    }

    double y_dim = config.getDouble("y_dim");
    double demag = y_dim / y_sensor_size;
    double aperture = config.getDouble("camera_aperture");
    omega = 1. / pow(4. * (demag + 1) * aperture, 2);

}

void DigitizationRunner::checkDimensionConsistency() {
    double x_dim = config.getDouble("x_dim");
    double y_dim = config.getDouble("y_dim");
    double testvalue = (x_dim / y_dim) - (double(x_pix) / y_pix);

    if (testvalue > 1e-1) {
        cerr << "Using Quest image dimension with Fusion simulation! Digitization FAILED!" << endl;
        exit(EXIT_FAILURE);
    } else if (testvalue < -1e-1) {
        cerr << "Using Fusion image dimension with Quest simulation! Digitization FAILED!" << endl;
        exit(EXIT_FAILURE);
    }
}

bool DigitizationRunner::isValidInputFile(const std::string& filename) const {
    bool fn_size = (filename.size() >= 5);
    if (!fn_size) {
        return false;
    } else if (fn_size && filename.substr(filename.size() - 5) == ".root") {
        return true;
    }
    return false;
}

bool DigitizationRunner::isValidDigiFile(const std::string& filename) const {
    // Match: starts with "digi_Run", followed by 5 or more digits, ending in ".root"
    string base = std::filesystem::path(filename).filename().string();
    std::regex pattern(R"(digi_Run\d{5,}\.root)");

    return regex_match(base, pattern);
}

int DigitizationRunner::extractDigiRunNumber(const std::string& filename) {
    string base = std::filesystem::path(filename).filename().string();
    std::smatch match;
    std::regex pattern(R"(digi_Run(\d{5,})\.root)");

    if (std::regex_match(base, match, pattern)) {
        return std::stoi(match[1]);  // match[1] contains the run number
    } else {
        throw std::invalid_argument("Invalid filename format: " + filename);
    }
}

void DigitizationRunner::setSeed(int seed) {
    if(config.getBool("fixed_seed")) gRandom->SetSeed(seed);
    return;
}

void DigitizationRunner::initSourceDir() {
    vector<string> path_to_config=Utils::splitString(configFile,'/');
    string parte= "";
    for(unsigned int i=0;i<path_to_config.size()-2;i++) parte=parte+path_to_config[i]+ "/";
    char buffer[PATH_MAX];
    char* ret = getcwd(buffer,sizeof(buffer));
    if(!ret) {
        throw runtime_error("DigitizationRunner::initSourceDir: Failed to get current directory");
    }
    string currentPath(buffer);

    SOURCE_DIR= currentPath+"/"+parte;
}

void DigitizationRunner::SaveValues(shared_ptr<TFile>& outfile) {
    outfile->cd();
    outfile->mkdir("param_dir");
    outfile->cd("param_dir");

    // [Fixme:] it's better to add specific member functions in the ConfigManager class
    map<string,string> options = config.getOptions();

    for(auto const& [key, val] : options)
    {
        // DEBUG
        // cout<<key<<": "<<val<<endl;
        
        if(key!="tag"       && key !="Vig_Map" &&
           key!="bckg_path" && key !="ped_cloud_dir" &&
           key!="noiserun"  && key !="bckg_name" &&
           key !="Camera_type" && key!="detector" &&
           key!="MC_xaxis"  && key !="MC_yaxis" &&
           key!="MC_zaxis"  && key!="digitizers" &&
           key!="pmt_noise_cloud_dir" 
           )
        {
            TH1F h(string(key).c_str(),"",1,0,1);
            
            double value;
            if(val == "True")       value = 1.;
            else if(val == "False") value = 0.;
            else                   value  = stod(val);
            
            h.SetBinContent(1,value);
            h.Write();
        }
    }
    
    if(options["bckg"]=="True") {
        std::vector<int> pedruns;
        
        stringstream ssruns(options["noiserun"]);
        string run;
        vector<string> seglist;
        while(getline(ssruns, run, ';')) {
            int runi = stoi(run);
            pedruns.push_back(runi);
        }
        
        int npeds = (int)pedruns.size();
        TH1F h("pedestal_runs","",npeds,0,npeds);
        for (int i = 0; i < npeds; i++) {
            h.SetBinContent(i+1, pedruns[i]);
        }
        h.Write();
        
    } else {
        TH1F h("pedestal_runs","",1,0,1);
        h.SetBinContent(1, -1);
        h.Write();
    }
    
    
    outfile->cd();
    return;
}

void DigitizationRunner::fillVigmap(TH2F& VignMap) {

    string vignfilename = Form("%sVignettingMap/%s", SOURCE_DIR.c_str(), config.get("Vig_Map").c_str());
    cout<<"Opening "<<vignfilename<<"..."<<endl;
    auto VignFile = unique_ptr<TFile> {TFile::Open(vignfilename.c_str())};
        
    VignMap = (*(TH2F*)VignFile->Get("normmap"));
        
    VignMap.Smooth();
        
    VignFile->Close();

    int x_vign=VignMap.GetNbinsX();
    int y_vign=VignMap.GetNbinsY();

    if(x_vign!=y_vign) {
        if(x_pix==y_pix) {
            cerr<<"You are using a Quest vignette map with a Fusion simulation! Digitization FAILED!\n";
            exit(0);
        }
    } else {
        if(x_pix!=y_pix) {
            cerr<<"You are using a Fusion vignette map with a Quest simulation! Digitization FAILED!\n";
            exit(0);
        }
    }
    return;
}

int DigitizationRunner::is_NR(vector<int> pdgID_hits, int pdg) {
    int ret = -999;
    for(unsigned int i = 0; i<pdgID_hits.size(); i++) {
        if(pdgID_hits[i] == 2212 || pdgID_hits[i] > pdg) {
            ret = pdgID_hits[i];
            break;
        }
    }
    return ret;
}


void DigitizationRunner::AddBckg(std::vector<std::vector<int>>& background) {
        
    string tmpfolder = config.get("bckg_path");
    string tmpname   = config.get("bckg_name");
    
    if(! filesystem::exists(tmpfolder)){
        //DEBUG
        cout<<"Creating tmpfolder..."<<
        system(("mkdir " + tmpfolder).c_str() );
    }
        
    if(! filesystem::exists(tmpfolder+tmpname)) {
        
        vector<int> pedruns;
        
        stringstream ssruns(config.get("noiserun"));
        string run;
        vector<string> seglist;
        while(getline(ssruns, run, ';')) {
            int runi = stoi(run);
            pedruns.push_back(runi);
        }
        
        cygnolib::joinMidasPedestals(pedruns, config.get("ped_cloud_dir"), tmpfolder, tmpname);
        
    }
        
    auto fin = unique_ptr<TFile> {TFile::Open((tmpfolder+tmpname).c_str())};

    // Computing number of images in root file
    int flength = 0;
    for (auto&& keyAsObj : *(fin->GetListOfKeys())){
        if(false) cout<<keyAsObj<<endl; //just to avoid a warning in the compilation
        flength++;
    }
    // Check if pedestal root file has events
    if(flength==0) {
        throw runtime_error("AddBckg: pedestal root file with no events. Please rm "+tmpfolder+tmpname+" before re-running.");
    }
        
    int pic_index = 0;
    if(config.getInt("random_ped")==-1) pic_index = gRandom->Integer(flength);
    else pic_index = config.getInt("random_ped");
    
    cout<<"Using pic # "<<pic_index<<" out of "<<flength<<" as a pedestal..."<<endl;
    TH2I* pic = fin->Get<TH2I>(Form("pic_%d", pic_index));

    int x_ped=pic->GetNbinsX();
    int y_ped=pic->GetNbinsY();
    //Check pedestal and simulation have same camera settings
    if(x_ped!=y_ped) {
        if(x_pix==y_pix) {
            cerr<<"You are using a Quest pedestal map with a Fusion simulation! Digitization FAILED!\n";
            exit(0);
        }
    } else {
        if(x_pix!=y_pix) {
            cerr<<"You are using a Fusion pedestal map with a Quest simulation! Digitization FAILED!\n";
            exit(0);
        }
    }

    for(int i = 0; i<pic->GetNbinsX(); i++) {
        for (int j =0; j<pic->GetNbinsY(); j++) {
            background[i][j] = pic->GetBinContent(i+1,j+1);
            //cout<<background[i][j]<<endl; // DEBUG
        }
    }
    
    return;
}

void DigitizationRunner::addPedsToParamDir(shared_ptr<TFile>& outfile) {
    outfile->cd();
    outfile->cd("param_dir");

    string histname = "pedestal_runs";
    // TH1F* hist = dynamic_cast<TH1F*>(file->Get(histname.c_str()));
    // if (!hist) {
    //     std::cerr << "Histogram not found: " << "pedestal_runs" << std::endl;
    //     outfile->Close();
    //     return;
    // }
    // Remove old histogram from directory
    outfile->Delete((histname + ";*").c_str());  // Delete all versions of the histogram

    // Get pedrun numbers from config file
    std::vector<int> pedruns;
    stringstream ssruns(config.get("noiserun"));
    string run;
    vector<string> seglist;
    while(getline(ssruns, run, ';')) {
        int runi = stoi(run);
        pedruns.push_back(runi);
    }
    int npeds = static_cast<int>(pedruns.size());

    TH1F hnew("pedestal_runs","",npeds,0,npeds);
    for (int i = 0; i < npeds; i++) {
        hnew.SetBinContent(i+1, pedruns[i]);
    }
    hnew.Write();
    outfile->cd();

    return;
}


pair<char, double> DigitizationRunner::getAxisMapping(const string& axis_label) {
    double sign = 1.0;
    char axis = axis_label[0];
    if (axis_label[0] == '-') {
        sign = -1.0;
        axis = axis_label[1];
    }
    return make_pair(axis, sign);
}

void DigitizationRunner::FillRedpix(const std::vector<std::vector<double>>& image,
                                   std::vector<uint16_t>* redpix_ix,
                                   std::vector<uint16_t>* redpix_iy,
                                   std::vector<uint16_t>* redpix_iz)
{
    redpix_ix->clear();
    redpix_iy->clear();
    redpix_iz->clear();

    const int nRows = image.size();
    const int nCols = image[0].size();

    for (int ix = 0; ix < nRows; ++ix) {
        for (int iy = 0; iy < nCols; ++iy) {
            double content = image[ix][iy];
            if (content != 0.0) {
                // N.B. the repix definition of x and y in the reco is physical, so they
                // must be swapped here
                redpix_ix->push_back(static_cast<uint16_t>(ix));
                redpix_iy->push_back(static_cast<uint16_t>(iy));
                redpix_iz->push_back(static_cast<uint16_t>(content));
                
            }
        }
    }
    return;
}


// ================================================================================================

void DigitizationRunner::processRootFiles() {
    
    const string infolder  = Utils::resolvePath(inputDir);
    const string outfolder = Utils::resolvePath(outputDir);

    // TrackProcessor initialization
    TrackProcessor processTrack(config);

    // Probably not needed anymore
    vector<int> eventnumber;
    vector<int> particle_type;
    vector<float> energy_ini;
    vector<float> theta_ini;
    vector<float> phi_ini;

    //std::string inDir = ConfigManager::resolvePath(inputDir);
    for (const auto& fentry : filesystem::directory_iterator(infolder)) {
        const string filename = fentry.path().string();
        if (!isValidInputFile(filename)) continue;
        
        cout << "Processing: " << filename << endl;
        
        auto f         = unique_ptr<TFile> {TFile::Open(filename.c_str())};
        auto inputtree = (TTree*)f->Get("nTuple");
        
        int firstentry = config.getInt("start_event");
        if(firstentry < 0) {
            cerr<<"DigitizationRunner::processRootFiles(): negative 'start_event' is not allowed. "
                <<"Please, check the input config file."<<endl;
            exit(EXIT_FAILURE);
        }
        int MC_events = inputtree->GetEntries();
        int max_events = config.getInt("events");

        int lastentry = -1;
        if(max_events == -1) {
            lastentry = MC_events - 1;
        } else if(max_events != -1 && max_events > (MC_events - firstentry)) {
            cout<<"WARNING: DigitizationRunner::processRootFiles: number of events "
                <<"to digitize greater than number of events available in input file "
                <<"after the first entry n = "<<firstentry<<". Only the latest "
                <<"remaining "<<(MC_events - firstentry)<<" events after event = "<<firstentry
                <<" will be digitized."<<endl;
            lastentry = MC_events - 1;
        } else {
            lastentry = max_events + firstentry - 1;
        }

        if (firstentry>MC_events) throw runtime_error("Error: First entry is larger than dimension of MC file, exiting!");
        cout << "Processing entries from "<<firstentry<<" to "<<lastentry<<"."<<endl;
            
        // output saved in outfolder/filename/
        stringstream ssfilename(filename.c_str());
        string ssfilename_el, filename_tmp;
        while(getline(ssfilename, ssfilename_el, '/')) {
            filename_tmp = ssfilename_el;
        }
        auto delimFN = filename_tmp.find(".root");
        string basefilename   = filename_tmp.substr(0, delimFN);
        string fnameoutfolder = outfolder + "/" + basefilename;
        if(! filesystem::exists(fnameoutfolder) && !config.getBool("queue")){
            int ret = system(("mkdir -p " + fnameoutfolder).c_str() );
            if(ret!=0) {
                cerr << "In DigitizationRunner::processRootFiles: Failed to create outdir: " << fnameoutfolder << endl;
                exit(EXIT_FAILURE);
            }
        }
            
        // Input file branches
        Int_t eventnumber;
        Int_t numhits;
        Double_t energyDep;
        Double_t energyDep_NR;
        vector<int>    *pdgID_hits = 0;
        vector<double> *tracklen_hits = 0;
        vector<double> *px_particle = 0;
        vector<double> *py_particle = 0;
        vector<double> *pz_particle = 0;
        vector<double> *energyDep_hits = 0;
        vector<double> *energyDep_hits_NRQF = 0;
        vector<double> *x_hits = 0;
        vector<double> *y_hits = 0;
        vector<double> *z_hits = 0;
            
        inputtree->SetBranchAddress("eventnumber", &eventnumber);
        inputtree->SetBranchAddress("numhits", &numhits);
        
        inputtree->SetBranchAddress("energyDep",    &energyDep);
        inputtree->SetBranchAddress("energyDep_NR", &energyDep_NR);
        inputtree->SetBranchAddress("pdgID_hits",    &pdgID_hits);
        inputtree->SetBranchAddress("tracklen_hits", &tracklen_hits);
        inputtree->SetBranchAddress("px_particle", &px_particle);
        inputtree->SetBranchAddress("py_particle", &py_particle);
        inputtree->SetBranchAddress("pz_particle", &pz_particle);
            
        inputtree->SetBranchAddress("energyDep_hits", &energyDep_hits);
        inputtree->SetBranchAddress("energyDep_hits_NRQF", &energyDep_hits_NRQF);
        inputtree->SetBranchAddress("x_hits", &x_hits);
        inputtree->SetBranchAddress("y_hits", &y_hits);
        inputtree->SetBranchAddress("z_hits", &z_hits);

        // Fill Vignetting map if needed
        TH2F VignMap;
        if(config.getBool("Vignetting")) {
            fillVigmap(VignMap);
        }
            
        //DEBUG
        //cout<<"DEBUG: "<<VignMap.GetBinContent(100,100)<<endl;

        bool fileStillToDigitize = true;
        int digipart = 0;
        
        while(fileStillToDigitize) {

            // Name of output file
            // standard:       histograms_RunRRRRR.root (R run number)
            // redpix_ouput:         digi_RunRRRRR.root (R run number)
            if(config.getBool("queue")) fnameoutfolder="./";
            string fileoutname = "";
            if(config.getBool("redpix_output")) {
                fileoutname= Form("%s/digi_Run%05d.root",
                                  fnameoutfolder.c_str(),
                                  runCount);
            } else {
                fileoutname= Form("%s/histograms_Run%05d.root",
                                  fnameoutfolder.c_str(),
                                  runCount);
            }
    
    
            auto outfile = shared_ptr<TFile> {TFile::Open(fileoutname.c_str(),
                                                              "RECREATE") };
    
            // Saving the parameters of digitization in output file
            SaveValues(outfile);
            //Output file branches
            Int_t eventnumber_out   = -999;
            Int_t particle_type_out = -999;
            Float_t energy = -999;
            Float_t theta  = -999;
            Float_t phi    = -999;
            Float_t track_length_3D = -1;
            Float_t x_vertex = -1;
            Float_t y_vertex = -1;
            Float_t z_vertex = -1;
            Float_t x_vertex_end = -1;
            Float_t y_vertex_end = -1;
            Float_t z_vertex_end = -1;
            Float_t x_min = -1;
            Float_t x_max = -1;
            Float_t y_min = -1;
            Float_t y_max = -1;
            Float_t z_min = -1;
            Float_t z_max = -1;
            Double_t placement_dx_mm = 0;
            Double_t placement_dy_mm = 0;
            Double_t placement_dz_mm = 0;
            Float_t px = -1;
            Float_t py = -1;
            Float_t pz = -1;
            Float_t proj_track_2D = -1;
            Int_t nhits_og = -1;
            Int_t N_photons = -1;
            Int_t N_primaries_reaching_GEMs = -1;
                
            Int_t row_cut = -1;
            Int_t N_photons_cut = -1;
            Float_t cut_energy = -1;
            Float_t x_min_cut = -1;
            Float_t x_max_cut = -1;
            Float_t y_min_cut = -1;
            Float_t y_max_cut = -1;
            Float_t z_min_cut = -1;
            Float_t z_max_cut = -1;
            Float_t proj_track_2D_cut = -1;

            //First 3 hits for directional algorithms if requested
            const bool save_mc_first3_hits = config.getBool("save_mc_first3_hits");
            std::vector<double> mc_first3_x;
            std::vector<double> mc_first3_y;
            std::vector<double> mc_first3_z;
            

            int nRedpix;
                
            // Smart pointer declarations
            std::unique_ptr<std::vector<uint16_t>> redpix_ix = std::make_unique<std::vector<uint16_t>>();
            std::unique_ptr<std::vector<uint16_t>> redpix_iy = std::make_unique<std::vector<uint16_t>>();
            std::unique_ptr<std::vector<uint16_t>> redpix_iz = std::make_unique<std::vector<uint16_t>>();
            
             // PMT Waveforms (Fast and Slow) ---
            //std::unique_ptr<std::vector<double>> wav_fast_pmt1 = std::make_unique<std::vector<double>>();
            std::vector<double> wav_fast_pmt1;
            std::vector<double> wav_fast_pmt2;
            std::vector<double> wav_fast_pmt3;
            std::vector<double> wav_fast_pmt4;

            std::vector<double> wav_slow_pmt1;
            std::vector<double> wav_slow_pmt2;
            std::vector<double> wav_slow_pmt3;
            std::vector<double> wav_slow_pmt4;
            

            // ROOT TTree
            auto outtree = std::make_unique<TTree>("event_info", "event_info");
                
            outtree->Branch("eventnumber", &eventnumber_out, "eventnumber/I");
            outtree->Branch("particle_type", &particle_type_out, "particle_type/I");
            outtree->Branch("energy", &energy, "energy/F");
            outtree->Branch("theta", &theta, "theta/F");
            outtree->Branch("phi", &phi, "phi/F");
            outtree->Branch("track_length_3D", &track_length_3D, "track_length_3D/F");
            outtree->Branch("proj_track_2D", &proj_track_2D, "proj_track_2D/F");
            outtree->Branch("x_vertex", &x_vertex, "x_vertex/F");
            outtree->Branch("y_vertex", &y_vertex, "y_vertex/F");
            outtree->Branch("z_vertex", &z_vertex, "z_vertex/F");
            outtree->Branch("x_vertex_end", &x_vertex_end, "x_vertex_end/F");
            outtree->Branch("y_vertex_end", &y_vertex_end, "y_vertex_end/F");
            outtree->Branch("z_vertex_end", &z_vertex_end, "z_vertex_end/F");
            outtree->Branch("x_min", &x_min, "x_min/F");
            outtree->Branch("x_max", &x_max, "x_max/F");
            outtree->Branch("y_min", &y_min, "y_min/F");
            outtree->Branch("y_max", &y_max, "y_max/F");
            outtree->Branch("z_min", &z_min, "z_min/F");
            outtree->Branch("z_max", &z_max, "z_max/F");
            outtree->Branch("placement_dx_mm", &placement_dx_mm, "placement_dx_mm/D");
            outtree->Branch("placement_dy_mm", &placement_dy_mm, "placement_dy_mm/D");
            outtree->Branch("placement_dz_mm", &placement_dz_mm, "placement_dz_mm/D");
            outtree->Branch("N_photons", &N_photons, "N_photons/I");
            outtree->Branch("N_primaries_reaching_GEMs", &N_primaries_reaching_GEMs, "N_primaries_reaching_GEMs/I");
            outtree->Branch("px", &px, "px/F");
            outtree->Branch("py", &py, "py/F");
            outtree->Branch("pz", &pz, "pz/F");
            outtree->Branch("nhits_og", &nhits_og, "nhits_og/I");
            if (config.getBool("exposure_time_effect")) {
                outtree->Branch("N_photons_cut", &N_photons_cut, "N_photons_cut/I");
                outtree->Branch("row_cut", &row_cut, "row_cut/I");
                outtree->Branch("cut_energy", &cut_energy, "cut_energy/F");
                outtree->Branch("x_min_cut", &x_min_cut, "x_min_cut/F");
                outtree->Branch("x_max_cut", &x_max_cut, "x_max_cut/F");
                outtree->Branch("y_min_cut", &y_min_cut, "y_min_cut/F");
                outtree->Branch("y_max_cut", &y_max_cut, "y_max_cut/F");
                outtree->Branch("z_min_cut", &z_min_cut, "z_min_cut/F");
                outtree->Branch("z_max_cut", &z_max_cut, "z_max_cut/F");
                outtree->Branch("proj_track_2D_cut", &proj_track_2D_cut, "proj_track_2D_cut/F");
            }
            // Create branches for redpixes
            outtree->Branch("nRedpix", &nRedpix, "nRedpix/I");
            outtree->Branch("redpix_ix", redpix_ix.get());
            outtree->Branch("redpix_iy", redpix_iy.get());
            outtree->Branch("redpix_iz", redpix_iz.get());
            if (config.getBool("pmt_mode")) {
                // Create branches for PMT waveforms
                //outtree->Branch("wav_fast_pmt1", wav_fast_pmt1.get());
                outtree->Branch("wav_fast_pmt1", &wav_fast_pmt1);
                outtree->Branch("wav_fast_pmt2", &wav_fast_pmt2);
                outtree->Branch("wav_fast_pmt3", &wav_fast_pmt3);
                outtree->Branch("wav_fast_pmt4", &wav_fast_pmt4);

                outtree->Branch("wav_slow_pmt1", &wav_slow_pmt1);
                outtree->Branch("wav_slow_pmt2", &wav_slow_pmt2);
                outtree->Branch("wav_slow_pmt3", &wav_slow_pmt3);
                outtree->Branch("wav_slow_pmt4", &wav_slow_pmt4);
            }

            if (save_mc_first3_hits) {
                outtree->Branch("mc_first3_x", &mc_first3_x);
                outtree->Branch("mc_first3_y", &mc_first3_y);
                outtree->Branch("mc_first3_z", &mc_first3_z);
            }

            int start = firstentry + digipart * NMAX_EVENTS;
            int stop  = start + NMAX_EVENTS-1;
            if(stop >= lastentry) stop = lastentry;
            for(int entry = start; entry <= stop; entry++) { // RUNNING ON ENTRIES
                //cout<<"DEBUG: Digitizing entry "<<entry<<", digipart = "<<digipart<<"..."<<endl;

                redpix_ix->clear();
                redpix_iy->clear();
                redpix_iz->clear();
                
                if (config.getBool("pmt_mode")) {
                    wav_fast_pmt1.clear(); wav_fast_pmt2.clear(); wav_fast_pmt3.clear(); wav_fast_pmt4.clear();
                    wav_slow_pmt1.clear(); wav_slow_pmt2.clear(); wav_slow_pmt3.clear(); wav_slow_pmt4.clear();
                }

                inputtree->GetEntry(entry);

                if (save_mc_first3_hits) {
                    //we record the MC truth for the first 3 hits
                    mc_first3_x.clear();
                    mc_first3_y.clear();
                    mc_first3_z.clear();
                
                    const std::size_t ntruth = std::min({
                        static_cast<std::size_t>(3),
                        x_hits->size(),
                        y_hits->size(),
                        z_hits->size()
                    });
                
                    mc_first3_x.assign(x_hits->begin(), x_hits->begin() + ntruth);
                    mc_first3_y.assign(y_hits->begin(), y_hits->begin() + ntruth);
                    mc_first3_z.assign(z_hits->begin(), z_hits->begin() + ntruth);
                }

                
                
                //DEBUG
                //if(options["NR"]=="True") cout<<particle_type<<endl;
                //if (entry==0) {
                //    cout<<numhits<<" - "<<pdgID_hits->size()<<endl;
                //    for(unsigned int i=0; i<pdgID_hits->size(); i++) {
                //        cout<<"---"<< (*pdgID_hits)[i] <<endl;
                //    }
                //}
                
                cout<<"\nEntry "<<entry<<", "<<entry+1-firstentry<< " / "<<lastentry+1-firstentry<<endl;
                    
                cout<<"Energy "<<energyDep    <<" keV"<<endl;
    
                bool NR_flag=false;
                // this would be the energy of the primary particle - equal to
                // deposited energy only if it is completely contained in the sensitive volume
                energy = energyDep;
                if (energyDep_NR>0){
                    particle_type_out = is_NR(*pdgID_hits, int(1.e9));
                    NR_flag = true;
                } 
                else particle_type_out = (*pdgID_hits)[0];  // this will tell us if the deposit was
                // started by a photon or an electron
    
                if(!config.getBool("NR")   && NR_flag==true ) continue;
                    
                //initialize array values - to save info also if the track is skipped (background only)
                row_cut         = -1;
                eventnumber_out = eventnumber;
                    
                    
                // DEBUG
                //cout<<"energyDep_NR = "<<energyDep_NR<<endl;
                //cout<<"particle_type_out = "<<particle_type_out<<endl;
            
                //particle_type_out = -999;
                //energy = -1;
                cut_energy      = -1;
                theta           =  0;
                phi             =  0;
                track_length_3D = -1;
                proj_track_2D   = -1;
                x_vertex        = -1;
                y_vertex        = -1;
                z_vertex        = -1;
                x_vertex_end    = -1;
                y_vertex_end    = -1;
                z_vertex_end    = -1;
                x_min           = -1;
                x_max           = -1;
                y_min           = -1;
                y_max           = -1;
                z_min           = -1;
                z_max           = -1;
                N_photons       =  0;
                N_primaries_reaching_GEMs = -1;
                x_min_cut       = -1;
                x_max_cut       = -1;
                y_min_cut       = -1;
                y_max_cut       = -1;
                z_min_cut       = -1;
                z_max_cut       = -1;
                N_photons_cut   =  0;
                proj_track_2D_cut = -1;
                px    = 0;
                py    = 0;
                pz    = 0;
                nhits_og = numhits;
                placement_dx_mm = placement_dy_mm = placement_dz_mm = 0;
    
    
                vector<vector<int>> background(x_pix,
                                                vector<int>(y_pix, 0)); 
                
                if(config.getBool("bckg")) {
                    AddBckg(background);
                }
      
                if (energy < config.getDouble("ion_pot")){
                    energy = 0;
                    TH2I final_image(Form("pic_run%d_ev%d", runCount, entry-start), "",
                                        x_pix, -0.5, x_pix -0.5,
                                        y_pix, -0.5, y_pix -0.5);
                    
                    for(unsigned int xx =0; xx < background.size(); xx++) {
                        for(unsigned int yy =0; yy < background[0].size(); yy++) {
                            final_image.SetBinContent(xx+1, yy+1, background[xx][yy]);
                        }
                    }

                    // Make sure nRedpix matches
                    nRedpix = redpix_ix->size();
                    
                    outtree->Fill();
                    outfile->cd();
                    if(!config.getBool("redpix_output")) {
                        final_image.Write();
                    }
                    
                    continue;
                }
                    
                vector<double> x_hits_tr;
                vector<double> y_hits_tr;
                vector<double> z_hits_tr;
                
                ///// Coordinate transformation
                // From MC input to digitization reference frame
                
                // Axis mapping from input MC frame to digitization frame
                map<char, std::pair<char, double>> axis_map;
                axis_map['x'] = getAxisMapping(config.get("MC_xaxis"));  // MC x axis  = digit. axis, sign
                axis_map['y'] = getAxisMapping(config.get("MC_yaxis"));  // MC y axis  = digit. axis, sign
                axis_map['z'] = getAxisMapping(config.get("MC_zaxis"));  // MC z axis  = digit. axis, sign
                
                // Map from MC axis name to the corresponding input hit vector
                map<char, const vector<double>*> input_axes;
                input_axes['x'] = x_hits;
                input_axes['y'] = y_hits;
                input_axes['z'] = z_hits;

                // Map from digitization axis name to the output hit vector
                map<char, vector<double>*> output_axes;
                output_axes['x'] = &x_hits_tr;
                output_axes['y'] = &y_hits_tr;
                output_axes['z'] = &z_hits_tr;

                // Perform transformation for each axis in digitization space
                const char axes[3] = {'x', 'y', 'z'};
                for (int i = 0; i < 3; ++i) {
                    char digit_axis = axes[i];
                    char mc_axis = 'n'; // n stands for still not defined
                    
                    // Retrieve the mapped MC axis and sign
                    for (const auto& pair : axis_map) {
                        const char map_mc_axis  = pair.first;
                        char map_digit_axis     = pair.second.first;
                        if (map_digit_axis == digit_axis) {
                            mc_axis = map_mc_axis;
                        }
                    }
                    pair<char, double> mapping = axis_map[mc_axis];
                    double sign = mapping.second;
                    
                    // Get input vector (MC) and output vector (digitization)
                    const vector<double>* input_vec = input_axes[mc_axis];
                    vector<double>* output_vec = output_axes[digit_axis];
                
                    // Get translation offset from config
                    string offset_key = string(1, digit_axis) + "_offset";
                    double offset = config.getDouble(offset_key);

                    // Get extra translation from config
                    string extra_key = string(1, digit_axis) + "_extra";
                    double extra = config.getDouble(extra_key);
                
                    // Transform all hits along this axis
                    transform(input_vec->begin(), input_vec->end(),
                                back_inserter(*output_vec),
                                [sign, offset, extra](double val) {
                                    return sign * val + offset + extra;
                                });
                }
                
                // NOTE: in Geant longitunal TPC axis is towards the GEMs, for the digi it's
                // assume it's towards the cathode
                
                    
                // DEBUG
                //if(entry == 0) {
                //    for(int ihit=0; ihit < numhits; ihit++) {
                //       cout<<x_hits_tr[ihit]<<",";
                //       cout<<y_hits_tr[ihit]<<",";
                //       cout<<z_hits_tr[ihit]<<"\n";
                //    }
                //}
                
                // Definition of energy hits
                vector<double> energy_hits;
                
                if (NR_flag) {
                    energy_hits = (*energyDep_hits_NRQF);
                } else {
                    energy_hits = (*energyDep_hits);
                }
                    
                // Random placement: one rigid translation per event, in digitization mm.
                // Keep Z first to preserve the random sequence of existing Z-only configs.
                map<char, double*> placement_shifts = {
                    {'x', &placement_dx_mm}, {'y', &placement_dy_mm}, {'z', &placement_dz_mm}
                };
                for (char axis : {'z', 'x', 'y'}) {
                    const string key = string("rand") + char(axis - 'a' + 'A');
                    const bool has_min = !config.get(key + "_min").empty();
                    const bool has_max = !config.get(key + "_max").empty();
                    if (has_min != has_max)
                        throw std::invalid_argument(key + " requires both _min and _max (mm)");
                    const double range = has_min || config.get(key + "_range").empty() ? 0. : config.getDouble(key + "_range");
                    const double low = has_min ? config.getDouble(key + "_min") : -range/2.;
                    const double high = has_max ? config.getDouble(key + "_max") : range/2.;
                    if (!std::isfinite(low) || !std::isfinite(high) || low > high)
                        throw std::invalid_argument(key + " bounds must be finite with min <= max (mm)");
                    const double shift = low == high ? low : has_min
                        ? low + gRandom->Uniform() * (high - low)
                        : (gRandom->Uniform() - 0.5) * range;
                    *placement_shifts[axis] = shift;
                    auto& hits = *output_axes[axis];
                    transform(hits.begin(), hits.end(), hits.begin(),
                              [shift](double a) { return a + shift; });
                }
                
                //Compute length and extremes of the complete track
                proj_track_2D = processTrack.GetTrackVariable("proj_track_2D", x_hits_tr, y_hits_tr, z_hits_tr);
                x_vertex      = processTrack.GetTrackVariable("x_vertex", x_hits_tr, y_hits_tr, z_hits_tr);
                y_vertex      = processTrack.GetTrackVariable("y_vertex", x_hits_tr, y_hits_tr, z_hits_tr);
                z_vertex      = processTrack.GetTrackVariable("z_vertex", x_hits_tr, y_hits_tr, z_hits_tr);
                x_vertex_end  = processTrack.GetTrackVariable("x_vertex_end", x_hits_tr, y_hits_tr, z_hits_tr);
                y_vertex_end  = processTrack.GetTrackVariable("y_vertex_end", x_hits_tr, y_hits_tr, z_hits_tr);
                z_vertex_end  = processTrack.GetTrackVariable("z_vertex_end", x_hits_tr, y_hits_tr, z_hits_tr);
                x_min = processTrack.GetTrackVariable("x_min", x_hits_tr, y_hits_tr, z_hits_tr);
                x_max = processTrack.GetTrackVariable("x_max", x_hits_tr, y_hits_tr, z_hits_tr);
                y_min = processTrack.GetTrackVariable("y_min", x_hits_tr, y_hits_tr, z_hits_tr);
                y_max = processTrack.GetTrackVariable("y_max", x_hits_tr, y_hits_tr, z_hits_tr);
                z_min = processTrack.GetTrackVariable("z_min", x_hits_tr, y_hits_tr, z_hits_tr);
                z_max = processTrack.GetTrackVariable("z_max", x_hits_tr, y_hits_tr, z_hits_tr);
                //DEBUG
                //cout<<"proj_track_2D = "<<Form("%.10f", proj_track_2D)<<endl;
                //cout<<"x_vertex = "<<x_vertex<<" ### y_vertex = "<<y_vertex<<" ### z_vertex = "<<z_vertex<<endl;
                //cout<<"x_vertex_end = "<<x_vertex_end<<" ### y_vertex_end = "<<y_vertex_end<<" ### z_vertex_end = "<<z_vertex_end<<endl;
                //cout<<" x_min = "<<x_min<<" x_max = "<<x_max<<" y_min = "<<y_min<<" y_max = "<<y_max<<" z_min = "<<z_min<<" z_max = "<<z_max<<endl;
                    
                vector<vector<double>> array2d_Nph(x_pix, vector<double>(y_pix, 0.0));
                
                auto ta = std::chrono::steady_clock::now();
                // with saturation
                std::vector<PMTVoxel> pmt_voxels; 

                // Select response deposits without changing the complete-track metadata.
                vector<double> x_active, y_active, z_active, energy_active;
                for (size_t i = 0; i < z_hits_tr.size(); ++i) {
                    if (z_hits_tr[i] >= 0.) {
                        x_active.push_back(x_hits_tr[i]);
                        y_active.push_back(y_hits_tr[i]);
                        z_active.push_back(z_hits_tr[i]);
                        energy_active.push_back(energy_hits[i]);
                    }
                }

                cout<<"Starting compute_cmos_with_saturation with size = "<<x_hits_tr.size()<<"..."<<endl;
                if(!processTrack.computeWithSaturation(x_active,
                                                        y_active,
                                                        z_active,
                                                        energy_active,
                                                        VignMap,
                                                        energy,
                                                        NR_flag,
                                                        array2d_Nph,
                                                        pmt_voxels
                                                        )) {
                    std::cerr<<"Warning: DigitizationRunner::processRootFiles: skipping this track because of error from TrackProcessor::computeWithSaturation."<<std::endl;

                    // The reported energy is -1 for those tracks to be recognizable
                    energy = -1;
                    TH2I final_image(Form("pic_run%d_ev%d", runCount, entry-start), "",
                                        x_pix, -0.5, x_pix -0.5,
                                        y_pix, -0.5, y_pix -0.5);
                    
                    for(unsigned int xx =0; xx < background.size(); xx++) {
                        for(unsigned int yy =0; yy < background[0].size(); yy++) {
                            final_image.SetBinContent(xx+1, yy+1, background[xx][yy]);
                        }
                    }

                    // Make sure nRedpix matches
                    nRedpix = redpix_ix->size();
                    
                    outtree->Fill();
                    outfile->cd();
                    if(!config.getBool("redpix_output")) {
                        final_image.Write();
                    }
                    
                    continue;
                    
                }

                cout<<"DEBUG: after compute"<<endl<<flush;
                auto tb = std::chrono::steady_clock::now();
                std::chrono::duration<double> durtmp=tb-ta;
                cout << "Time taken in seconds to compute_cmos_with_saturation is: " << durtmp.count() << endl;
                
                // ========================================================================== //
                //                               PMT SIMULATION                               //
                // ========================================================================== //
                
                if (config.getBool("pmt_mode") && !pmt_voxels.empty()) {
                    
                    auto pmt_start = std::chrono::steady_clock::now();
                    std::cout << "[DEBUG] Starting PMT simulation ..." << std::endl;
                    
                    const std::map<std::string, std::string> pmt_options = config.getOptions();
                    
                    // 1. Properties Setup
                    const double drift_field = config.getDouble("drift_field"); 
                    const double v_drift_mm_ns = Utils::compute_drift_velocity(drift_field) / 100.0;
                    const double qe_factor = 0.0136;        //average QE of PMT 
                
                    // 2. Data Preparation                
                    double min_z = 1.0e9;
                    bool event_detected = false;
                    for (const auto& voxel : pmt_voxels) {
                         if (voxel.n_photons * qe_factor >= 0.5) { 
                            if (voxel.z < min_z) { min_z = voxel.z; event_detected = true; }
                         }
                    }
                    const double time_shift_reference = event_detected ? min_z : 0.0;
                
                    std::vector<double> x_coords, y_coords, arrival_times;
                    std::vector<int> effective_photons;
                    
                    for (const auto& voxel : pmt_voxels) {
                        const int n_photons_eff = std::round(voxel.n_photons * qe_factor);
                        if (n_photons_eff > 0) {
                            x_coords.push_back(voxel.x);
                            y_coords.push_back(voxel.y);
                            effective_photons.push_back(n_photons_eff);
                            
                            const double z_corrected = voxel.z - time_shift_reference;
                            arrival_times.push_back(std::abs(z_corrected) / v_drift_mm_ns);
                        }
                    }
                
                    // 3. Light Propagation and Signal Digitization
                    if (!x_coords.empty()) {
                        
                        // Photon Propagation Step
                        auto prop_start = std::chrono::steady_clock::now();
                        PhotonPropagation propagation(x_coords, y_coords, effective_photons, arrival_times, pmt_options);
                        auto pmt_hits_map = propagation.pmt_hits();
                        auto prop_end = std::chrono::steady_clock::now();
                        std::chrono::duration<double> prop_time = prop_end - prop_start;
                        std::cout << "[DEBUG] Photon propagation (pmt_hits) took: " << prop_time.count() << " s" << std::endl;
                
                        // Signal Simulation Step
                        auto sig_start = std::chrono::steady_clock::now();
                        SignalSimulation simulator(pmt_hits_map, pmt_options, SOURCE_DIR);
                        std::map<std::string, std::vector<double>> fast_signals, slow_signals;
                        simulator.simulated_signals(fast_signals, slow_signals);
                        auto sig_end = std::chrono::steady_clock::now();
                        std::chrono::duration<double> sig_time = sig_end - sig_start;
                        std::cout << "[DEBUG] Signal simulation (pmt_sim) took: " << sig_time.count() << " s" << std::endl;
                
                        // 4. Output Assignment
                        //*wav_fast_pmt1 = fast_signals["pmt_1"];
                        wav_fast_pmt1 = fast_signals["pmt_1"];
                        wav_fast_pmt2 = fast_signals["pmt_2"];
                        wav_fast_pmt3 = fast_signals["pmt_3"];
                        wav_fast_pmt4 = fast_signals["pmt_4"];
                        wav_slow_pmt1 = slow_signals["pmt_1"];
                        wav_slow_pmt2 = slow_signals["pmt_2"];
                        wav_slow_pmt3 = slow_signals["pmt_3"];
                        wav_slow_pmt4 = slow_signals["pmt_4"];
                    }
                
                    auto pmt_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> total_pmt_time = pmt_end - pmt_start;
                    std::cout << "[DEBUG] TOTAL PMT Simulation execution time: " << total_pmt_time.count() << " s" << std::endl;
                }
                
                // ========================================================================== //


                // Integral of the track - if opt.exposure_effect, it's computed anyway after the cut on the original hits
                N_photons = accumulate(array2d_Nph.cbegin(), array2d_Nph.cend(), 0, [](auto sum, const auto& row) {
                    return accumulate(row.cbegin(), row.cend(), sum);
                });
                // DEBUG
                cout<<"N_photons = "<<N_photons<<endl;

                N_primaries_reaching_GEMs = processTrack.GetN_primaries_reaching_GEMs();
                
                // Compute always redpix, before possible track cut by exposure of sensor
                FillRedpix(array2d_Nph, redpix_ix.get(), redpix_iy.get(), redpix_iz.get());
                
                if (config.getBool("exposure_time_effect")) {
                    // Cut image and hits
                    row_cut = processTrack.ApplyExposureCut(array2d_Nph,
                                                            x_hits_tr,
                                                            y_hits_tr,
                                                            z_hits_tr,
                                                            energy_hits);

                    // Integral after camera exposure cut
                    N_photons_cut = accumulate(array2d_Nph.cbegin(),
                                                array2d_Nph.cend(), 0, [](auto sum, const auto& row) {
                        return accumulate(row.cbegin(), row.cend(), sum);
                    });

                    // Energy after camera exposure cut
                    cut_energy = accumulate(energy_hits.begin(), energy_hits.end(), 0.0);

                }

                // Build final TH2I image
                TH2I final_image(Form("pic_run%d_ev%d", runCount, entry-start), "",
                                    x_pix, -0.5, x_pix -0.5,
                                    y_pix, -0.5, y_pix -0.5);
                for(unsigned int xx =0; xx < array2d_Nph.size(); xx++) {
                    for(unsigned int yy =0; yy < array2d_Nph[0].size(); yy++) {
                        int binc = background[xx][yy]+(int)array2d_Nph[xx][yy];
                        final_image.SetBinContent(xx+1, yy+1, binc);
                    }
                }
                
                // Compute variables (after the cut) to be saved in the tree
                if (config.getBool("exposure_time_effect")) {
                    proj_track_2D_cut = processTrack.GetTrackVariable("proj_track_2D", x_hits_tr, y_hits_tr, z_hits_tr);
                    x_min_cut = processTrack.GetTrackVariable("x_min", x_hits_tr, y_hits_tr, z_hits_tr);
                    x_max_cut = processTrack.GetTrackVariable("x_max", x_hits_tr, y_hits_tr, z_hits_tr);
                    y_min_cut = processTrack.GetTrackVariable("y_min", x_hits_tr, y_hits_tr, z_hits_tr);
                    y_max_cut = processTrack.GetTrackVariable("y_max", x_hits_tr, y_hits_tr, z_hits_tr);
                    z_min_cut = processTrack.GetTrackVariable("z_min", x_hits_tr, y_hits_tr, z_hits_tr);
                    z_max_cut = processTrack.GetTrackVariable("z_max", x_hits_tr, y_hits_tr, z_hits_tr);
                }
                //DEBUG
                //cout<<"proj_track_2D_cut = "<<Form("%.10f", proj_track_2D_cut)<<endl;
                //cout<<" x_min_cut = "<<x_min_cut<<" x_max_cut = "<<x_max_cut<<" y_min_cut = "<<y_minv<<" y_max_cut = "<<y_max_cut<<" z_min_cut = "<<z_min_cut<<" z_max_cut = "<<z_max_cut<<endl;

                // if there are at least 2 hits compute theta and phi
                if(x_hits_tr.size() > 1) {
                    phi   = atan2( y_hits_tr[1] - y_hits_tr[0],
                                    x_hits_tr[1] - x_hits_tr[0]
                                    );
                    theta = acos((z_hits_tr[1] - z_hits_tr[0]) /
                                    sqrt(
                                        (x_hits_tr[1] - x_hits_tr[0])*(x_hits_tr[1] - x_hits_tr[0])+
                                        (y_hits_tr[1] - y_hits_tr[0])*(y_hits_tr[1] - y_hits_tr[0])+
                                        (z_hits_tr[1] - z_hits_tr[0])*(z_hits_tr[1] - z_hits_tr[0]))
                                    );
                } else {
                    phi   = -999;
                    theta = -999;
                }
                
            
                auto tc = std::chrono::steady_clock::now();
                std::chrono::duration<double> durtmpc=tc-ta;
                cout << "Total time taken in seconds to digitize track is: " << durtmpc.count() << endl;
                    
                // Make sure nRedpix matches
                nRedpix = redpix_ix->size();
                
                outtree->Fill();
                outfile->cd();
                // DEBUG
                if(!config.getBool("redpix_output")) {
                    final_image.Write();
                }
            }
            //outfile->cd("event_info");
            outtree->Write();
                
            cout<<Form("COMPLETED RUN %d",runCount)<<endl;

            // Carefully closing the ROOT file
            // N.B. this is strictly needed if redpix pointers are defines as smart pointers
            outtree->ResetBranchAddresses(); // prevents ROOT from accessing deleted memory
            outtree.reset(); // must destroy before redpix vectors are gone, i.e. before ending of current scope
            outfile->Close();
            
            digipart ++;
            runCount ++;
            if(stop == lastentry) fileStillToDigitize = false;
        }
        f->Close();
    }
}

// ================================================================================================
// ================================================================================================
// ================================================================================================
// ================================================================================================

void DigitizationRunner::generateHistogramsFromDigi() {
    const string infolder  = Utils::resolvePath(inputDir);
    const string outfolder = Utils::resolvePath(outputDir);

    //std::string inDir = ConfigManager::resolvePath(inputDir);
    for (const auto& fentry : filesystem::directory_iterator(infolder)) {
        const string filename = fentry.path().string();
        if (!isValidDigiFile(filename)) continue;
        int runNumber = extractDigiRunNumber(filename);
        // Generate the output filename
        string outfilename= Form("%s/histograms_Run%05d.root", outfolder.c_str(),runNumber);
        cout << "Generating output filename: " << outfilename << " from input filename: "<< filename << std::endl;

        // Copy the input file to the output file
        filesystem::copy_file(filename, outfilename, filesystem::copy_options::overwrite_existing);
        // DEBUG
        // std::cout << "Copied " << filename << " -> " << outfilename << std::endl;

        // Open the newly created file to add TH2I images
        auto f         = shared_ptr<TFile> {TFile::Open(outfilename.c_str(), "UPDATE")};
        auto inputtree = (TTree*)f->Get("event_info");

        // Add pedestal run number
        addPedsToParamDir(f);

        // Setup redpix readout
        vector<uint16_t> *redpix_ix = 0;
        vector<uint16_t> *redpix_iy = 0;
        vector<uint16_t> *redpix_iz = 0;
        int nRedpix = -1;
        inputtree->SetBranchAddress("redpix_ix", &redpix_ix);
        inputtree->SetBranchAddress("redpix_iy", &redpix_iy);
        inputtree->SetBranchAddress("redpix_iz", &redpix_iz);
        inputtree->SetBranchAddress("nRedpix", &nRedpix);

        // Loop over the events in the digi_RunXXXXX.root file
        int MC_events = inputtree->GetEntries();
        for(int entry = 0; entry < MC_events; entry ++) {
            inputtree->GetEntry(entry);
            cout<< "Producing picture for event "<<entry<<" out of "<< MC_events<<endl;

            // Get random pedestal
            vector<vector<int>> background(x_pix, vector<int>(y_pix, 0));

            AddBckg(background);

            // TH2I image creation
            TH2I final_image(Form("pic_run%d_ev%d", runNumber, entry), "",
                             x_pix, -0.5, x_pix -0.5,
                             y_pix, -0.5, y_pix -0.5);
            for(int idx = 0; idx < nRedpix; idx++) {
                int xx = static_cast<int>((*redpix_ix)[idx]);
                int yy = static_cast<int>((*redpix_iy)[idx]);
                int binc = static_cast<int>((*redpix_iz)[idx]);
                final_image.SetBinContent(xx+1, yy+1, binc);
            }

            // Add pedestal to final image
            for(int xx = 0; xx < x_pix; xx++) {
                for(int yy = 0; yy < y_pix; yy++) {
                    int flatbin = final_image.GetBin(xx + 1, yy + 1);
                    final_image.AddBinContent(flatbin, static_cast<double>(background[xx][yy]));
                }
            }

            f->cd();
            final_image.Write();
        }
        f->Close();
    }

    return;
}
