//CALL: ./progname.exe Configfile -I inputdir -O outputdir

#include <iostream>
#include <unistd.h>
#include <limits.h>
#include <sstream>
#include <fstream>
#include <chrono>
#include <map>
#include <random>
#include <string>
#include <execution>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "s3.h"
#include "cygnolib.h"

using namespace std;


// Global and constant value;
double GEM1_gain;
double GEM2_gain;
double GEM3_gain;
double extraction_eff_GEM1;
double extraction_eff_GEM2;
double extraction_eff_GEM3;
double omega;
int x_ini = 0;
int y_ini = 0;
int z_ini = 0;

vector<string> split(const string& , char );
void ReadConfig(string name, map<string,string>& options);
void ReadG4Isotopes(string name, map<string,string>& dict_isotopes);
void ReadIonlist(string name, vector<vector<string>>& ionlist);
void SaveValues(map<string,string>& options, shared_ptr<TFile>& outfile);

vector<double> NelGEM1(const vector<double>& N_ioniz_el);
vector<double> NelGEM2(const vector<double>& energyDep,const vector<double>& z_hit, map<string,string>& options);

void AddBckg(map<string,string>& options, vector<vector<int>>& background);

bool is_NR(vector<int> pdgID_hits, int pdg);
double angle_between(vector<double>& v1, vector<double>& v2);
vector<double> crossProduct(vector<double>& a, vector<double>& b);
vector<double> rotateByAngleAndAxis(vector<double>& vec, double angle, vector<double>& axis);

void compute_cmos_with_saturation(vector<double>& x_hits_tr,
                                  vector<double>& y_hits_tr,
                                  vector<double>& z_hits_tr,
                                  vector<double>& energy_hits,
                                  map<string,string>& options,
                                  vector<vector<double>>& array2d_Nph);

void compute_cmos_without_saturation(vector<double>& x_hits_tr,
                                     vector<double>& y_hits_tr,
                                     vector<double>& z_hits_tr,
                                     vector<double>& energy_hits,
                                     map<string,string>& options,
                                     vector<vector<double>>& array2d_Nph);

void cloud_smearing3D(vector<double>& x_hits_tr,
                      vector<double>& y_hits_tr,
                      vector<double>& z_hits_tr,
                      vector<double>& energy_hits,
                      map<string,string>& options,
                      vector<double>& S3D_x,
                      vector<double>& S3D_y,
                      vector<double>& S3D_z);

void ph_smearing2D(vector<double>& x_hits_tr,
                   vector<double>& y_hits_tr,
                   vector<double>& z_hits_tr,
                   vector<double>& energy_hits,
                   map<string,string>& options,
                   vector<double>& S2D_x,
                   vector<double>& S2D_y) ;

vector<double> compute_sigma(const double diff_const, const double diff_coeff, const vector<double>& dz);
vector<double> smear(const vector<double>& axis_hit, const vector<double>& axis_sigma, const vector<double>& nel);

vector<double> arange(double start, double stop, double step);
double round_up_to_even(const double f);

//void Nph_saturation(const TH3I& h3d, map<string,string>& options, vector<vector<double>>& hout);
double Nph_saturation(int nel, double A, double beta);

void TrackVignetting(vector<vector<double>>& array2d_Nph, int xpix, int ypix, const TH2F & VignMap);

// Old approach
//string rootlocation(string tag, int run);   // inconsistency for the MC-old tag!!!
//string root_filename(string tag, int run);

int main(int argc, char** argv)
{
    string outfolder;
    string infolder;
   
    if(argc<2) {cerr<<"No Configfile given!!\nSuggested use: ./progname.exe Configfile -I inputdir -O outputdir"; exit(EXIT_FAILURE);}
    string nome=argv[1];
    
    vector<string> path_to_config=split(nome,'/');
    string parte= "";
    for(int i=0;i<path_to_config.size()-2;i++) parte=parte+path_to_config[i]+ "/";
    char buffer[PATH_MAX];
    getcwd(buffer,sizeof(buffer));
    string currentPath(buffer);
    const string SOURCE_DIR= currentPath+"/"+parte;
   
    if(argc<3)
    {
        infolder=SOURCE_DIR+"/input/";
        outfolder=SOURCE_DIR+"/OutDir/";
    }
    else
    {
        if(argc!=6) {cerr<<"Wrong parameter inputs given!!\nCorrect use: ./progname.exe Configfile -I inputdir -O outputdir"; exit(EXIT_FAILURE);}
        string parseop=argv[2];
        string parseop2=argv[4];
        //cout<< parseop << "   " << parseop2<< endl;
        if(parseop!="-I" && parseop2!="-I") {cerr<<"Wrong options name!!\nSuggested use: ./progname.exe Configfile -I inputdir -O outputdir"; exit(EXIT_FAILURE);}
        if(parseop!="-O" && parseop2!="-O") {cerr<<"Wrong options name!!\nSuggested use: ./progname.exe Configfile -I inputdir -O outputdir"; exit(EXIT_FAILURE);}
	
        if(parseop=="-I") {infolder=argv[3]; outfolder=argv[5];}
        else {infolder=argv[5]; outfolder=argv[3];}
    }

    // DEBUG
    // cout<<"Input Folder: "<<infolder<<endl;
    // cout<<"Output Folder: "<<outfolder<<endl;
    
    map<string,string> options;
    ReadConfig(nome,options);						//Function to be checked

    GEM1_gain = 0.0347*exp((0.0209)*stod(options["GEM1_HV"]));
    GEM2_gain = 0.0347*exp((0.0209)*stod(options["GEM2_HV"]));
    GEM3_gain = 0.0347*exp((0.0209)*stod(options["GEM3_HV"]));

    /*
    GEM1_gain = 0.0347*exp((0.025)*stod(options["GEM1_HV"]));
    GEM2_gain = 0.0347*exp((0.025)*stod(options["GEM2_HV"]));
    GEM3_gain = 0.0347*exp((0.025)*stod(options["GEM3_HV"]));
    */

    cout<< "GEM1_gain = " << GEM1_gain;
    cout<< "\nGEM2_gain = " << GEM2_gain;
    cout<< "\nGEM3_gain = " << GEM3_gain << endl;
    
    //dividing Fernando's to Francesco&Karolina's single GEM gain measurement
    extraction_eff_GEM1 = 0.87319885*exp(-0.0020000000*stod(options["GEM1_HV"]));
    extraction_eff_GEM2 = 0.87319885*exp(-0.0020000000*stod(options["GEM2_HV"]));
    extraction_eff_GEM3 = 0.87319885*exp(-0.0020000000*stod(options["GEM3_HV"]));
    cout<< "extraction eff GEM1 = " << extraction_eff_GEM1;
    cout<< "\nextraction eff GEM2 = " << extraction_eff_GEM2;
    cout<< "\nextraction eff GEM3 = " << extraction_eff_GEM3 << endl;
	
    double y_dim=stod(options["y_dim"]);
    double demag=y_dim/stod(options["sensor_size"]);
    double aperture=stod(options["camera_aperture"]);
    omega=1./pow(4.*(demag+1)*aperture,2);
	
    //Code execution
    int run_count = 1;
    auto t0 = std::chrono::steady_clock::now();
    
    if(options["fixed_seed"]=="True" || options["fixed_seed"]=="true") gRandom->SetSeed(10);
    
    vector<int> eventnumber;
    vector<int> particle_type;
    vector<float> energy_ini;
    vector<float> theta_ini;
    vector<float> phi_ini;
    
    if(! filesystem::exists(outfolder)){
        //DEBUG
        cout<<"Creating oufolder..."<<
        system(("mkdir -p" + outfolder).c_str() );
    }
    
    map<string, string> dict_isotopes;
    if(options["GEANT4isotopes"]=="True") ReadG4Isotopes("../Z_A_isotopes.csv", dict_isotopes);
    
    
    
    string ending=".root";
    for(const auto& fentry : filesystem::directory_iterator(infolder))
    {
        string filename=fentry.path();
        cout<<filename<<endl;
        bool ends= false;
        if (ending.size() <= filename.size())    ends=equal(ending.rbegin(), ending.rend(), filename.rbegin());
        
        //DEBUG
        if(!ends) cout<<"Not ending with "<<ending<<endl;
        else cout<<filename<<" is ending with "<<ending<<"!"<<endl;
        
        //if(filename.ends_with(".root"))   //c++20 fix it wtf
        //{
        
        if(ends) {
            
            //DEBUG
            if(filename.find("HeCF4gas_AmBe_part") != string::npos) {
            //if(filename.find("LIME_CADshield") != string::npos) {
                continue;
            }
            
            vector<vector<string>> SRIM_events;
            
            if(options["NR"]=="True" && options["NR_list"]!="") {
                auto delim1 = filename.find("part");
                auto delim2 = filename.find(".root");
                if(delim1==string::npos) throw runtime_error("Cannot determine the 'part' of the file.\n");
                if(delim2==string::npos) throw runtime_error("Input file is not a root file.\n");
                auto part= filename.substr(delim1+4, delim2-delim1-4);
                cout<<"Using NR list from "<<options["NR_list"]<<"_part"<<part<<".py"<<endl;
                
                ReadIonlist(Form("%s/%s_part%s.py", infolder.c_str(), options["NR_list"].c_str(), part.c_str()),
                            SRIM_events);
                
            }
            
            
            
            
            auto f         = unique_ptr<TFile> {TFile::Open(filename.c_str())};
            auto inputtree = (TTree*)f->Get("nTuple");
            
            int max_events = inputtree->GetEntries();
            int totev = (stod(options["events"])==-1) ? max_events : stod(options["events"]);
            totev = min(totev, max_events);
            
            int firstentry = stoi(options["start_event"]);
            if (firstentry>totev) throw runtime_error("Error: First entry is larger than last entry, exiting!");
            cout << "Processing entries from "<<firstentry<<" to "<<totev<<"."<<endl;
            
            
            // output saved in outfolder/filename/
            
            stringstream ssfilename(filename.c_str());
            string ssfilename_el, filename_tmp;
            while(getline(ssfilename, ssfilename_el, '/')) {
                filename_tmp = ssfilename_el;
            }
            auto delimFN = filename_tmp.find(".root");
            string basefilename   = filename_tmp.substr(0, delimFN);
            string fnameoutfolder = outfolder + "/" + basefilename;
            if(! filesystem::exists(fnameoutfolder)){
                system(("mkdir " + fnameoutfolder).c_str() );
            }
            
            // standard: name of output file = histograms_RunRRRRR.root (R run number)
            string fileoutname= Form("%s/histograms_Run%05d.root",
                                     fnameoutfolder.c_str(),
                                     run_count);
            
            // for radioisotope simulation: histograms_RunZZAAANN.root (Z=atomic number, A=mass numbe$
            // NOTE: this is a 7 digit run number, while reconstruction currently looks for 5
            string isot_numb = "0000000";
            if(options["GEANT4isotopes"]=="True") {
                cout<<"GEANT4isotopes option is active."<<endl;
                
                stringstream ssinfile(basefilename);
                string tmpstr;
                int counter = 0;
                while(getline(ssinfile, tmpstr, '_')) {
                    if(counter == 1) {
                        isot_numb = dict_isotopes[tmpstr];
                        // DEBUG
                        //isot_numb = "00000";
                    }
                    counter++;
                }
                
                auto delimBFN = basefilename.find("part");
                if(delimBFN==string::npos) throw runtime_error("Cannot determine the 'part' of the file.\n");
                auto part = basefilename.substr(delimBFN+4, basefilename.size() - delimBFN - 4);
                      
                if(filename.find("part")!= string::npos) {
                    fileoutname = Form("%s/histogram_Runs%05d%02d.root",
                                       fnameoutfolder.c_str(),
                                       stoi(isot_numb),
                                       stoi(part)
                                       );
                    isot_numb = Form("%05d%02d", stoi(isot_numb), stoi(part));
                } else {
                    fileoutname = Form("%s/histogram_Runs%05d00.root",
                                       fnameoutfolder.c_str(),
                                       stoi(isot_numb));
                    
                    isot_numb = Form("%05d00", stoi(isot_numb));
                }
                //DEBUG
                //cout<<"DEBUG: fileoutname = "<<fileoutname<<endl;
                //cout<<"DEBUG: isot_numb = "<<isot_numb<<endl;
            }
            
            if(options["start_event"]!="0"){
                cout<<"out folder "<<fnameoutfolder<<endl;
                int newpart = (int)(stoi(options["start_event"])/500);
                int oldpart = stoi(isot_numb);
                int partnum = oldpart + newpart;
                fileoutname = Form("%s/histograms_Run%07d.root",
                                   fnameoutfolder.c_str(),
                                   partnum);
            }
            
            
            
            
            auto outfile = shared_ptr<TFile> {TFile::Open(fileoutname.c_str(),
                                                          "RECREATE") };
            outfile->mkdir("event_info");
            SaveValues(options,outfile);
            
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
            Float_t px = -1;
            Float_t py = -1;
            Float_t pz = -1;
            Float_t proj_track_2D = -1;
            Int_t nhits_og = -1;
            Int_t N_photons = -1;
            
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
            
            auto outtree = new TTree("event_info", "event_info");
            
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
            outtree->Branch("N_photons", &N_photons, "N_photons/I");
            outtree->Branch("px", &px, "px/F");
            outtree->Branch("py", &py, "py/F");
            outtree->Branch("pz", &pz, "pz/F");
            outtree->Branch("nhits_og", &nhits_og, "nhits_og/I");
            if (options["exposure_time_effect"]=="True") {
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
            
            // Input file branches
            Int_t eventnumber;
            Int_t numhits;
            Double_t energyDep;
            Double_t energyDep_NR;
            Float_t  ekin_particle;
            Int_t particle_type;
            vector<int>    *pdgID_hits = 0;
            vector<double> *tracklen_hits = 0;
            vector<double> *px_particle = 0;
            vector<double> *py_particle = 0;
            vector<double> *pz_particle = 0;
            vector<double> *energyDep_hits = 0;
            vector<double> *x_hits = 0;
            vector<double> *y_hits = 0;
            vector<double> *z_hits = 0;
            
            
            inputtree->SetBranchAddress("eventnumber", &eventnumber);
            inputtree->SetBranchAddress("numhits", &numhits);
            
            if(options["NR"]=="False") { // TO BE CHECK ON REAL SRIM SIMULATIONS
                inputtree->SetBranchAddress("energyDep",    &energyDep);
                inputtree->SetBranchAddress("energyDep_NR", &energyDep_NR);
                inputtree->SetBranchAddress("pdgID_hits",    &pdgID_hits);
                inputtree->SetBranchAddress("tracklen_hits", &tracklen_hits);
                inputtree->SetBranchAddress("px_particle", &px_particle);
                inputtree->SetBranchAddress("py_particle", &py_particle);
                inputtree->SetBranchAddress("pz_particle", &pz_particle);
            }
            
            inputtree->SetBranchAddress("energyDep_hits", &energyDep_hits);
            inputtree->SetBranchAddress("x_hits", &x_hits);
            inputtree->SetBranchAddress("y_hits", &y_hits);
            inputtree->SetBranchAddress("z_hits", &z_hits);
            
            if(options["NR"]=="True") { // TO BE CHECK ON REAL SRIM SIMULATIONS
                inputtree->SetBranchAddress("particle_type", &particle_type);
                inputtree->SetBranchAddress("ekin_particle", &ekin_particle);
            }
            
            TH2F VignMap;
            if(options["Vignetting"]=="True") {
                string vignfilename = Form("%sVignettingMap/%s",SOURCE_DIR.c_str(), options["Vig_Map"].c_str());
                cout<<"Opening "<<vignfilename<<"..."<<endl;
                auto VignFile = unique_ptr<TFile> {TFile::Open(vignfilename.c_str())};
                
                VignMap = (*(TH2F*)VignFile->Get("normmap_lime"));
                
                VignMap.Smooth();
                
                VignFile->Close();
            }
            
            //DEBUG
            //cout<<"DEBUG: "<<VignMap->GetBinContent(0,0)<<endl;
            
            for(int entry=0; entry<totev; entry++) {  // RUNNING ON ENTRIES
                
                inputtree->GetEntry(entry);
                
                //DEBUG
                //if(options["NR"]=="True") cout<<particle_type<<endl;
                //if (entry==0) {
                //    cout<<numhits<<" - "<<pdgID_hits->size()<<endl;
                //    for(unsigned int i=0; i<pdgID_hits->size(); i++) {
                //        cout<<"---"<< (*pdgID_hits)[i] <<endl;
                //    }
                //}
                
                cout<<"Entry "<<entry<<" of "<<totev<<endl;
                
                if (options["NR"]=="True"){
                    cout<<"Energy "<<ekin_particle<<" keV"<<endl;
                } else {
                    cout<<"Energy "<<energyDep    <<" keV"<<endl;
                }
                
                if(options["NR"]=="False" && energyDep>900    ) continue;
                if(options["NR"]=="True"  && ekin_particle>900) continue;
                
                //initialize array values - to save info also if the track is skipped (background only)
                row_cut         = -1;
                eventnumber_out = eventnumber;
                
                if(options["NR"] == "True") {
                    energy            = ekin_particle;
                    particle_type_out = particle_type;
                } else {
                    // this would be the energy of the primary particle - equal to
                    // deposited energy only if it is completely contained in the sensitive volume
                    // energy = ekin_particle * 1000;
                    energy = energyDep;
                    if (energyDep_NR>0) particle_type_out = is_NR(*pdgID_hits, int(1.e9));
                    else particle_type_out = (*pdgID_hits)[0];  // this will tell us if the deposit was
                    // started by a photon or an electron
                }
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
                
                
                vector<vector<int>> background(stoi(options["x_pix"]),
                                               vector<int>(stoi(options["x_pix"]), 0));
                AddBckg(options, background);
                
                if (energy < stod(options["ion_pot"])){
                    energy = 0;
                    TH2I final_image(Form("pic_run%d_ev%d", run_count, entry), "",
                                     stoi(options["x_pix"]), -0.5, stoi(options["x_pix"]) -0.5,
                                     stoi(options["y_pix"]), -0.5, stoi(options["y_pix"]) -0.5);
                    
                    for(unsigned int xx =0; xx < background.size(); xx++) {
                        for(unsigned int yy =0; yy < background[0].size(); yy++) {
                            final_image.SetBinContent(yy, xx, background[xx][yy]);
                        }
                    }
                    
                    outtree->Fill();
                    outfile->cd();
                    final_image.Write();
                    
                    continue;
                }
                
                vector<double> x_hits_tr;
                vector<double> y_hits_tr;
                vector<double> z_hits_tr;
                

                if (options["NR"]=="True") {
                    // x_hits_tr = np.array(tree.x_hits) + opt.x_offset
                    // y_hits_tr = np.array(tree.y_hits) + opt.y_offset
                    // z_hits_tr = np.array(tree.z_hits) + opt.z_offset
                    vector<double> v1 = {1.,0.,0.};
                    vector<double> v2 = {stod(SRIM_events[entry][3])-stod(SRIM_events[entry][2]),
                                         stod(SRIM_events[entry][5])-stod(SRIM_events[entry][4]),
                                         stod(SRIM_events[entry][7])-stod(SRIM_events[entry][6]),
                                        };
                    
                    double        angle = angle_between(v1, v2);
                    vector<double> axis =  crossProduct(v1, v2);
                    
                    
                    double norm = sqrt(inner_product(axis.begin(), axis.end(), axis.begin(), 0.0));
                    vector<double> uaxis = {axis[0]/norm, axis[1]/norm, axis[2]/norm};
                    // DEBUG
                    //std::cout<<angle<<endl;
                    //std::cout<<"-"<<axis[0]<<","<<axis[1]<<","<<axis[2]<<endl;
                    
                    
                    for(int ihit=0; ihit < numhits; ihit++) {
                        vector<double> tmpvec = {(*x_hits)[ihit], (*y_hits)[ihit], (*z_hits)[ihit]};
                        vector<double> rotvec = rotateByAngleAndAxis(tmpvec, angle, uaxis);
                        
                        x_hits_tr.push_back(rotvec[0]+stod(SRIM_events[entry][2])+stod(options["x_offset"]));
                        y_hits_tr.push_back(rotvec[1]+stod(SRIM_events[entry][4])+stod(options["y_offset"]));
                        z_hits_tr.push_back(rotvec[2]+stod(SRIM_events[entry][6])+stod(options["z_offset"]));
                    }
                    
                } else {
                    transform(z_hits->begin(),
                              z_hits->end(),
                              back_inserter(x_hits_tr),
                              [&] (double a) {return -a + stod(options["x_offset"]);});
                    transform(x_hits->begin(),
                              x_hits->end(),
                              back_inserter(y_hits_tr),
                              [&] (double a) {return a + stod(options["y_offset"]);});
                    transform(y_hits->begin(),
                              y_hits->end(),
                              back_inserter(z_hits_tr),
                              [&] (double a) {return -a - stod(options["z_offset"]);});
                    
                    // FIXME: [Check which is the z axis orientation]
                    
                }
                
                // DEBUG
                //for(int ihit=0; ihit < numhits; ihit++) {
                //    cout<<x_hits_tr[ihit]<<",";
                //    cout<<y_hits_tr[ihit]<<",";
                //    cout<<z_hits_tr[ihit]<<"\n";
                //}
                
                vector<double> energy_hits = (*energyDep_hits);
                
                // add random Z to tracks
                if (stod(options["randZ_range"]) != 0) {
                    double rand = (gRandom->Uniform() - 0.5) * stod(options["randZ_range"]);
                    //DEBUG
                    //cout<<"rand = "<<rand<<endl;
                    transform(z_hits_tr.begin(), z_hits_tr.end(), z_hits_tr.begin(),
                              [&] (double a) {return a + rand;}
                              );
                    
                }
                
                //Compute length and extremes of the track before the cut
                proj_track_2D = 0;
                for(int ihit=0; ihit < numhits-1; ihit++){
                    proj_track_2D += sqrt((x_hits_tr[ihit+1]-x_hits_tr[ihit])*(x_hits_tr[ihit+1]-x_hits_tr[ihit])+
                                          (y_hits_tr[ihit+1]-y_hits_tr[ihit])*(y_hits_tr[ihit+1]-y_hits_tr[ihit])
                                          );
                }
                // DEBUG
                //cout<<"proj_track_2D = "<<Form("%.10f", proj_track_2D)<<endl;
                
                
                x_vertex = (x_hits_tr[0] + 0.5 * stod(options["x_dim"]) )*stod(options["x_pix"])/stod(options["x_dim"]); //in pixels
                y_vertex = (y_hits_tr[0] + 0.5 * stod(options["y_dim"]) )*stod(options["y_pix"])/stod(options["y_dim"]); //in pixels
                z_vertex = abs(z_hits_tr[0]+stod(options["z_gem"])); //distance from GEMs in mm
                // DEBUG
                //cout<<"x_vertex = "<<x_vertex<<" ### y_vertex = "<<y_vertex<<" ### z_vertex = "<<z_vertex<<endl;
                
                x_vertex_end = (x_hits_tr[numhits-1] + 0.5 * stod(options["x_dim"])) * stod(options["x_pix"]) / stod(options["x_dim"]); //in pixels
                y_vertex_end = (y_hits_tr[numhits-1] + 0.5 * stod(options["y_dim"])) * stod(options["y_pix"]) / stod(options["y_dim"]); //in pixels
                z_vertex_end = abs(z_hits_tr[numhits-1]+stod(options["z_gem"])); //distance from GEMs in mm
                //DEBUG
                //cout<<"x_vertex_end = "<<x_vertex_end<<" ### y_vertex_end = "<<y_vertex_end<<" ### z_vertex_end = "<<z_vertex_end<<endl;
                
                x_min = (*min_element(x_hits_tr.begin(), x_hits_tr.end()) + 0.5 * stod(options["x_dim"])) * stod(options["x_pix"]) / stod(options["x_dim"]);
                x_max = (*max_element(x_hits_tr.begin(), x_hits_tr.end()) + 0.5 * stod(options["x_dim"])) * stod(options["x_pix"]) / stod(options["x_dim"]);
                y_min = (*min_element(y_hits_tr.begin(), y_hits_tr.end()) + 0.5 * stod(options["y_dim"])) * stod(options["y_pix"]) / stod(options["y_dim"]);
                y_max = (*max_element(y_hits_tr.begin(), y_hits_tr.end()) + 0.5 * stod(options["y_dim"])) * stod(options["y_pix"]) / stod(options["y_dim"]);
                z_min = min(abs(*max_element(z_hits_tr.begin(),
                                             z_hits_tr.end()) + stod(options["z_gem"])),
                            abs(*min_element(z_hits_tr.begin(),
                                             z_hits_tr.end()) + stod(options["z_gem"])));
                z_max = max(abs(*max_element(z_hits_tr.begin(),
                                             z_hits_tr.end()) + stod(options["z_gem"])),
                            abs(*min_element(z_hits_tr.begin(),
                                             z_hits_tr.end()) + stod(options["z_gem"])));
                //DEBUG
                //cout<<" x_min = "<<x_min<<" x_max = "<<x_max<<" y_min = "<<y_min<<" y_max = "<<y_max<<" z_min = "<<z_min<<" z_max = "<<z_max<<endl;
                
                
                //CUT TRACKS due to exposure of camera
                double randcut = gRandom->Uniform(stod(options["exposure_time"])+stod(options["readout_time"]));
                //randcut = 390.0;
                
                if (options["exposure_time_effect"] == "True") {
                    if (randcut<stod(options["readout_time"])) {
                        
                        double y_cut_tmp = stod(options["y_dim"]) * (0.5 - randcut/stod(options["readout_time"]))-3.;
                        
                        // Removing elements from x_hits_tr
                        x_hits_tr.erase(std::remove_if(x_hits_tr.begin(), x_hits_tr.end(), [&](const double& x) {
                            return y_hits_tr[&x-&*x_hits_tr.begin()] < y_cut_tmp;
                        }), x_hits_tr.end());
                        // Removing elements from z_hits_tr
                        z_hits_tr.erase(std::remove_if(z_hits_tr.begin(), z_hits_tr.end(), [&](const double& z) {
                            return y_hits_tr[&z-&*z_hits_tr.begin()] < y_cut_tmp;
                        }), z_hits_tr.end());
                        // Removing elements from energy_hits
                        energy_hits.erase(std::remove_if(energy_hits.begin(), energy_hits.end(), [&](const double& e) {
                            return y_hits_tr[&e-&*energy_hits.begin()] < y_cut_tmp;
                        }), energy_hits.end());
                        
                        // Removing elements from y_hits_tr [must be done after the previous ones]
                        y_hits_tr.erase(std::remove_if(y_hits_tr.begin(), y_hits_tr.end(),[&](const double& y) {
                            return y < y_cut_tmp;
                        }), y_hits_tr.end());
                        
                        row_cut = stoi(options["y_pix"]) - (int)(randcut * stod(options["y_pix"]) / stod(options["readout_time"]));
                        
                        //DEBUG
                        //cout<<"y_cut_tmp = "<<y_cut_tmp<<endl;
                        //cout<<"row_cut = "<<row_cut<<endl;
                        //cout<<"sizes = ["<< x_hits_tr.size()<<","<<y_hits_tr.size()<<","
                        //    <<z_hits_tr.size()<<","<<energy_hits.size()<<"]"<<endl;
                        
                    } else if (randcut>stod(options["exposure_time"])) {
                        double y_cut_tmp = stod(options["y_dim"]) * (0.5 - (randcut - stod(options["exposure_time"])) / stod(options["readout_time"]))-3.;
                        
                        // Removing elements from x_hits_tr
                        x_hits_tr.erase(std::remove_if(x_hits_tr.begin(), x_hits_tr.end(), [&](const double& x) {
                            return y_hits_tr[&x-&*x_hits_tr.begin()] > y_cut_tmp;
                        }), x_hits_tr.end());
                        // Removing elements from z_hits_tr
                        z_hits_tr.erase(std::remove_if(z_hits_tr.begin(), z_hits_tr.end(), [&](const double& z) {
                            return y_hits_tr[&z-&*z_hits_tr.begin()] > y_cut_tmp;
                        }), z_hits_tr.end());
                        // Removing elements from energy_hits
                        energy_hits.erase(std::remove_if(energy_hits.begin(), energy_hits.end(), [&](const double& e) {
                            return y_hits_tr[&e-&*energy_hits.begin()] > y_cut_tmp;
                        }), energy_hits.end());
                        
                        // Removing elements from y_hits_tr [must be done after the previous ones]
                        y_hits_tr.erase(std::remove_if(y_hits_tr.begin(), y_hits_tr.end(),[&](const double& y) {
                            return y > y_cut_tmp;
                        }), y_hits_tr.end());
                        
                        row_cut = stoi(options["y_pix"]) - (int)((randcut-stod(options["readout_time"])) * stod(options["y_pix"]) / stod(options["readout_time"]));
                        
                        //DEBUG
                        //cout<<"y_cut_tmp = "<<y_cut_tmp<<endl;
                        //cout<<"row_cut = "<<row_cut<<endl;
                        //cout<<"sizes = ["<< x_hits_tr.size()<<","<<y_hits_tr.size()<<","
                        //    <<z_hits_tr.size()<<","<<energy_hits.size()<<"]"<<endl;
                        
                    }
                    if(x_hits_tr.size()==0){
                        cut_energy = 0;
                        cout<<"The track was completely cut"<<endl;
                        TH2I final_image(Form("pic_run%d_ev%d", run_count, entry), "",
                                         stoi(options["x_pix"]), -0.5, stoi(options["x_pix"]) -0.5,
                                         stoi(options["y_pix"]), -0.5, stoi(options["y_pix"]) -0.5);
                        
                        for(unsigned int xx =0; xx < background.size(); xx++) {
                            for(unsigned int yy =0; yy < background[0].size(); yy++) {
                                final_image.SetBinContent(yy, xx, background[xx][yy]);
                            }
                        }
                        
                        outtree->Fill();
                        outfile->cd();
                        final_image.Write();
                        
                        continue;
                    }
                }
                
                vector<vector<double>> array2d_Nph(stoi(options["y_pix"]),
                                                   vector<double>(stoi(options["x_pix"]), 0.0));
                
                auto ta = std::chrono::steady_clock::now();
                // with saturation
                if(options["saturation"]=="True") {
                    cout<<"Starting compute_cmos_with_saturation with size = "<<x_hits_tr.size()<<"..."<<endl;
                    compute_cmos_with_saturation(x_hits_tr,
                                                 y_hits_tr,
                                                 z_hits_tr,
                                                 energy_hits,
                                                 options,
                                                 array2d_Nph
                                                 );
                } else {// no saturation
                    compute_cmos_without_saturation(x_hits_tr,
                                                    y_hits_tr,
                                                    z_hits_tr,
                                                    energy_hits,
                                                    options,
                                                    array2d_Nph
                                                    );
                }
                
                // Integral of the track - if opt.exposure_effect, it's computed anyway after the cut on the original hits (to save time we digitize only the part that will be visible)
                N_photons = accumulate(array2d_Nph.cbegin(), array2d_Nph.cend(), 0, [](auto sum, const auto& row) {
                    return accumulate(row.cbegin(), row.cend(), sum);
                });
                // DEBUG
                cout<<"N_photons = "<<N_photons<<endl;
                
                if(options["Vignetting"]=="True") {
                    TrackVignetting(array2d_Nph,
                                    stoi(options["y_pix"]),
                                    stoi(options["x_pix"]),
                                    VignMap);
                }
                
                
                if (options["exposure_time_effect"]=="True") { //cut the track post-smearing
                    if (randcut<stod(options["readout_time"])) {
                        for(unsigned int xx=0; xx < array2d_Nph.size(); xx++) {
                            for(int yy=0; yy < row_cut; yy++) {
                                array2d_Nph[xx][yy] = 0.0;
                            }
                        }
                    } else if(randcut> stod(options["exposure_time"]) ) {
                        for(unsigned int xx=0; xx < array2d_Nph.size(); xx++) {
                            for(int yy=row_cut; yy < (int)array2d_Nph[0].size(); yy++) {
                                array2d_Nph[xx][yy] = 0.0;
                            }
                        }
                    }
                    
                    //integral after camera exposure cut
                    N_photons_cut = accumulate(array2d_Nph.cbegin(),
                                               array2d_Nph.cend(), 0, [](auto sum, const auto& row) {
                        return accumulate(row.cbegin(), row.cend(), sum);
                    });
                }
                
                TH2I final_image(Form("pic_run%d_ev%d", run_count, entry), "",
                                 stoi(options["x_pix"]), -0.5, stoi(options["x_pix"]) -0.5,
                                 stoi(options["y_pix"]), -0.5, stoi(options["y_pix"]) -0.5);
                
                for(unsigned int xx =0; xx < array2d_Nph[0].size(); xx++) {
                    for(unsigned int yy =0; yy < array2d_Nph.size(); yy++) {
                        // DEBUG
                        //int binc = (int)array2d_Nph[yy][array2d_Nph[0].size()-xx];
                        int binc = background[xx][yy]+(int)array2d_Nph[yy][array2d_Nph[0].size()-xx];
                        final_image.SetBinContent(yy, xx, binc);
                    }
                }
                
                //Cut again the hits to save the effective length and energy which is visible in the final image,
                //and compute the number of photons post-cut
                if (options["exposure_time_effect"]=="True") {
                    if (randcut<stod(options["readout_time"])) {
                        double y_cut_tmp = stod(options["y_dim"]) * (0.5 - randcut/stod(options["readout_time"]));
                        
                        // Removing elements from x_hits_tr
                        x_hits_tr.erase(std::remove_if(x_hits_tr.begin(), x_hits_tr.end(), [&](const double& x) {
                            return y_hits_tr[&x-&*x_hits_tr.begin()] < y_cut_tmp;
                        }), x_hits_tr.end());
                        // Removing elements from z_hits_tr
                        z_hits_tr.erase(std::remove_if(z_hits_tr.begin(), z_hits_tr.end(), [&](const double& z) {
                            return y_hits_tr[&z-&*z_hits_tr.begin()] < y_cut_tmp;
                        }), z_hits_tr.end());
                        // Removing elements from energy_hits
                        energy_hits.erase(std::remove_if(energy_hits.begin(), energy_hits.end(), [&](const double& e) {
                            return y_hits_tr[&e-&*energy_hits.begin()] < y_cut_tmp;
                        }), energy_hits.end());
                        
                        // Removing elements from y_hits_tr [must be done after the previous ones]
                        y_hits_tr.erase(std::remove_if(y_hits_tr.begin(), y_hits_tr.end(),[&](const double& y) {
                            return y < y_cut_tmp;
                        }), y_hits_tr.end());
                        
                    } else if (randcut> stod(options["exposure_time"]) ) {
                        double y_cut_tmp = stod(options["y_dim"]) * (0.5 - (randcut - stod(options["exposure_time"])) / stod(options["readout_time"]));
                        
                        // Removing elements from x_hits_tr
                        x_hits_tr.erase(std::remove_if(x_hits_tr.begin(), x_hits_tr.end(), [&](const double& x) {
                            return y_hits_tr[&x-&*x_hits_tr.begin()] > y_cut_tmp;
                        }), x_hits_tr.end());
                        // Removing elements from z_hits_tr
                        z_hits_tr.erase(std::remove_if(z_hits_tr.begin(), z_hits_tr.end(), [&](const double& z) {
                            return y_hits_tr[&z-&*z_hits_tr.begin()] > y_cut_tmp;
                        }), z_hits_tr.end());
                        // Removing elements from energy_hits
                        energy_hits.erase(std::remove_if(energy_hits.begin(), energy_hits.end(), [&](const double& e) {
                            return y_hits_tr[&e-&*energy_hits.begin()] > y_cut_tmp;
                        }), energy_hits.end());
                        
                        // Removing elements from y_hits_tr [must be done after the previous ones]
                        y_hits_tr.erase(std::remove_if(y_hits_tr.begin(), y_hits_tr.end(),[&](const double& y) {
                            return y > y_cut_tmp;
                        }), y_hits_tr.end());
                        
                    }
                    
                    cut_energy = accumulate(energy_hits.begin(), energy_hits.end(), 0.0);
                    
                    if(x_hits_tr.size()==0){
                        cut_energy = 0;
                        cout<<"The track was completely cut after smearing"<<endl;
                        
                        TH2I final_image_cut(Form("pic_run%d_ev%d", run_count, entry), "",
                                         stoi(options["x_pix"]), -0.5, stoi(options["x_pix"]) -0.5,
                                         stoi(options["y_pix"]), -0.5, stoi(options["y_pix"]) -0.5);
                        
                        for(unsigned int xx =0; xx < background.size(); xx++) {
                            for(unsigned int yy =0; yy < background[0].size(); yy++) {
                                final_image_cut.SetBinContent(yy, xx, background[xx][yy]);
                            }
                        }
                        
                        outtree->Fill();
                        outfile->cd();
                        final_image_cut.Write();
                        
                        continue;
                    }
                    
                }
                
                
                // Compute variables to be saved in the tree
                proj_track_2D_cut = 0;
                for(int ihit=0; ihit < numhits-1; ihit++){
                    proj_track_2D_cut += sqrt((x_hits_tr[ihit+1]-x_hits_tr[ihit])*(x_hits_tr[ihit+1]-x_hits_tr[ihit])+
                                              (y_hits_tr[ihit+1]-y_hits_tr[ihit])*(y_hits_tr[ihit+1]-y_hits_tr[ihit])
                                             );
                }
                // DEBUG
                //cout<<"proj_track_2D_cut = "<<Form("%.10f", proj_track_2D_cut)<<endl;
                
                
                x_min_cut = (*min_element(x_hits_tr.begin(), x_hits_tr.end()) + 0.5 * stod(options["x_dim"])) * stod(options["x_pix"]) / stod(options["x_dim"]);
                x_max_cut = (*max_element(x_hits_tr.begin(), x_hits_tr.end()) + 0.5 * stod(options["x_dim"])) * stod(options["x_pix"]) / stod(options["x_dim"]);
                y_min_cut = (*min_element(y_hits_tr.begin(), y_hits_tr.end()) + 0.5 * stod(options["y_dim"])) * stod(options["y_pix"]) / stod(options["y_dim"]);
                y_max_cut = (*max_element(y_hits_tr.begin(), y_hits_tr.end()) + 0.5 * stod(options["y_dim"])) * stod(options["y_pix"]) / stod(options["y_dim"]);
                z_min_cut = min(abs(*max_element(z_hits_tr.begin(),
                                                 z_hits_tr.end()) + stod(options["z_gem"])),
                                abs(*min_element(z_hits_tr.begin(),
                                                 z_hits_tr.end()) + stod(options["z_gem"])));
                z_max_cut = max(abs(*max_element(z_hits_tr.begin(),
                                                 z_hits_tr.end()) + stod(options["z_gem"])),
                                abs(*min_element(z_hits_tr.begin(),
                                                 z_hits_tr.end()) + stod(options["z_gem"])));
                //DEBUG
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
                
                if(options["NR"]=="False") {
                    track_length_3D = accumulate(tracklen_hits->begin(), tracklen_hits->end(), 0.0);
                    px              = (*px_particle)[0];
                    py              = (*py_particle)[0];
                    pz              = (*pz_particle)[0];
                }
                
                auto tb = std::chrono::steady_clock::now();
                std::chrono::duration<double> durtmp=tb-ta;
                cout << "Time taken in seconds to compute_cmos_with_saturation is: " << durtmp.count() << endl;
                
                
                outtree->Fill();
                outfile->cd();
                
                // DEBUG
                final_image.Write();
                
            }
            outfile->cd("event_info");
            outtree->Write();
            
            cout<<Form("COMPLETED RUN %d",run_count)<<endl;
            
            f->Close();
            outfile->Close();
            
            run_count ++;
        }
        
    }
	
    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> dur=t1-t0;
    cout << "Time taken in seconds is: " << dur.count() << endl;
    return 0;
    
}

// ====== FUNCTIONS DEFINITION ======
vector<string> split(const string& str, char delimiter)
{
    vector<string> tokens;
    istringstream stream(str);
    string token;

    while(getline(stream,token,delimiter)) tokens.push_back(token);

    return tokens;
}


void ReadConfig(string name, map<string,string>& options)
{
    ifstream config(name.c_str());
	
    string line;
    while(getline(config,line))
    {
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return isspace(c);}),line.end());
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return c=='\'';}),line.end());
        if(line[0] == '#' || line.empty() || line[0] == '{' || line[0] == '}') continue;
		
        auto delim1= line.find(":");
        auto delim2= line.find(",");
        if(delim2==string::npos) delim2=line.size();
        auto index= line.substr(0,delim1);
        auto val= line.substr(delim1+1,min(delim2,line.size())-delim1-1 );
        options[index]=val;
    }
	
}

void ReadG4Isotopes(string name, map<string,string>& dict_isotopes) {
    ifstream config(name.c_str());
    
    string line;
    while(getline(config, line))
    {
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return isspace(c);}),line.end());
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return c=='\'';}),line.end());
        if(line[0] == '#' || line.empty() || line[0] == '{' || line[0] == '}') continue;
        
        stringstream to_split(line.c_str());
        string element;
        
        int counter = 0;
        string index;
        string val;
        while(getline(to_split, element, ',')) {
            if(counter == 0) index = element;
            else if (counter == 1) val = element;
            else if (counter == 2) val += Form("%03d", stoi(element));
            counter++;
        }
        dict_isotopes[index]=val;
        
        
        // DEBUG
        // cout<<"---"<<index<<": "<<val<<endl;
    }
    
    return;
}


void ReadIonlist(string name, vector<vector<string>>& ionlist) {
    
    cout<<"Opening and parsing "<<name<<" ..."<<endl;
    ifstream config(name.c_str());
    
    string line;
    while(getline(config, line))
    {
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return isspace(c);}),line.end());
        if(line.find("ionlist")!=string::npos) continue;
        stringstream to_split(line.c_str());
        
        string element;
        while(getline(to_split, element, '[')) {
            if(!element.empty()){
                element.erase(remove_if(element.begin(),element.end(),[=] (char c){return c==']';}), element.end());
                stringstream content(element.c_str());

                string var;
                vector<string> row;
                while(getline(content, var, ',')) {
                    if(!var.empty()){
                        row.push_back(var);
                    }
                }
                ionlist.push_back(row);
            }
        }
        break;
    }
    return;
}

vector<double> NelGEM1(const vector<double>& N_ioniz_el) {
    
    vector<double> n_tot_el(N_ioniz_el.size(), 0);
    
    if (N_ioniz_el.size() == 1) { // in case there is only one hit (very low energy)
        //for(int j = 0; j<(int)N_ioniz_el[0]; j++) {
            //double nsec = gRandom->Exp(GEM1_gain) * extraction_eff_GEM1;
            //n_tot_el[0] += nsec;
        n_tot_el[0] = N_ioniz_el[0];
        //}
    } else {
        for(unsigned int i = 0; i<N_ioniz_el.size(); i++) {
            for(int j = 0; j<(int)round(N_ioniz_el[i]); j++) {
                double nsec = gRandom->Exp(GEM1_gain) * extraction_eff_GEM1;
                n_tot_el[i] += nsec;
            }
        }
    }
    //DEBUG
    //for(unsigned int i=0; i<n_tot_el.size(); i++) {
    //    cout<<n_tot_el[i]<<"\n";
    //}
    
    return n_tot_el;
}


vector<double> NelGEM2(const vector<double>& energyDep, const vector<double>& z_hit, map<string,string>& options) {
    vector<double> n_ioniz_el_ini;
    double opt_pot=stod(options["ion_pot"]);
    transform(energyDep.begin(),energyDep.end(),back_inserter(n_ioniz_el_ini), [&] (double a) { return a/opt_pot;});
    
    vector<double> drift_l;
    int opt_gem=stod(options["z_gem"]);
    transform(z_hit.begin(),z_hit.end(),back_inserter(drift_l), [&] (double a) { return abs(a-opt_gem);});
    
    vector<double> n_ioniz_el_mean(n_ioniz_el_ini.size(), 0.0);
    
    double optabsorption_l=stod(options["absorption_l"]);
    for(unsigned int i=0;i<n_ioniz_el_mean.size();i++) n_ioniz_el_mean[i]=abs(n_ioniz_el_ini[i]*exp(-drift_l[i]/optabsorption_l));
    
    vector<double> n_ioniz_el(n_ioniz_el_ini.size(), 0);
    transform(n_ioniz_el_mean.begin(), n_ioniz_el_mean.end(), n_ioniz_el.begin(), [&] (double a) {
        return gRandom->Poisson(a);
    });
    
    // total number of secondary electrons considering the gain in the 2nd GEM foil
    vector<double> n_tot_el = NelGEM1(n_ioniz_el);
    transform(n_tot_el.begin(), n_tot_el.end(), n_tot_el.begin(), [&] (double a) {
        return round(a * GEM2_gain * extraction_eff_GEM2);
    });
    
    return n_tot_el;
    
}

void SaveValues(map<string,string>& options, shared_ptr<TFile>& outfile)
{
    outfile->cd();
    outfile->mkdir("param_dir");
    outfile->cd("param_dir");
    for(auto const& [key, val] : options)
    {
        // DEBUG
        // cout<<key<<": "<<val<<endl;
        
        if(key!="tag"       && key !="Vig_Map" &&
           key!="bckg_path" && key !="ped_cloud_dir" &&
           key!="NR_list"
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
    outfile->cd();
    return;
}

void AddBckg(map<string,string>& options, vector<vector<int>>& background) {
    
    string tmpfolder = options["bckg_path"];
    
    if(! filesystem::exists(tmpfolder)){
        //DEBUG
        cout<<"Creating tmpfolder..."<<
        system(("mkdir " + tmpfolder).c_str() );
    }
    
    
    if(options["bckg"]=="True") {
        
        string tmpname = Form( "%s/run%05d.mid.gz", tmpfolder.c_str(), stoi(options["noiserun"]));
        
        
        if(! filesystem::exists(tmpname)) {
            // DEBUG
            cout << "Looking for file: "<<tmpname<<endl;
            bool verbose = false;
            bool cloud   = true;
            
            int run = stoi(options["noiserun"]);
            
            //Download or find midas file
            string filename = s3::cache_file(s3::mid_file(run, options["ped_cloud_dir"], cloud, verbose),
                                             options["bckg_path"],
                                             cloud,
                                             options["ped_cloud_dir"],
                                             verbose);
        }
        
        int pic_index = gRandom->Integer(100);
        // DEBUG
        pic_index = 0;
        
        // DEBUG
        cout<<"Using pic # "<<pic_index<<" as a pedestal..."<<endl;
        
        //reading data from midas file
        cout<<"Opening midas file "<<tmpname<<" ..."<<endl;
        TMReaderInterface* reader = cygnolib::OpenMidasFile(tmpname);
        bool reading = true;
        
        int counter = 0;
        bool found = false;
        while (reading) {
            
            TMidasEvent event = TMidasEvent();
            reading = TMReadEvent(reader, &event);
            if (!reading) {
                // DEBUG
                // std::cout<<"EOF reached."<<std::endl;
                break;
            }
            
            bool cam_found = cygnolib::FindBankByName(event, "CAM0");
            if(cam_found) {
                if(counter == pic_index) {
                    cygnolib::Picture pic=cygnolib::daq_cam2pic(event, "fusion");
                    
                    //TH2I rootpic(Form("pic_run%d_ev%d", run_count, entry), "" , x_pix, -0.5, x_pix-0.5, y_pix, -0.5, y_pix-0.5);
                    
                    vector<vector<uint16_t>> vecpic = pic.GetFrame();
                    for(unsigned int i = 0; i<vecpic.size(); i++) {
                        for (unsigned int j =0; j<vecpic[0].size(); j++) {
                            //rootpic.SetBinContent(j, i, vecpic[i][j]);
                            background[i][j] = vecpic[i][j];
                        }
                    }
                    //background=rootpic;
                    
                    found = true;
                    break;
                }
                
                counter++;
            }
            
        }
        if(!found) {
            cerr<<"AddBckg: Error: Cannot find pic # "<<pic_index<<" in pedestal run."<<endl;
            exit(EXIT_FAILURE);
        }
        
    }
    
    return;
}


bool is_NR(vector<int> pdgID_hits, int pdg) {
    int ret = -999;
    for(unsigned int i = 0; i<pdgID_hits.size(); i++) {
        if(pdgID_hits[i] > pdg) {
            ret = pdgID_hits[i];
            break;
        }
    }
    return ret;
}

// Returns the angle in radians between vectors 'v1' and 'v2'
// >>> angle_between((1, 0, 0), (0, 1, 0))
// 1.5707963267948966
// >>> angle_between((1, 0, 0), (1, 0, 0))
// 0.0
// >>> angle_between((1, 0, 0), (-1, 0, 0))
// 3.141592653589793

double angle_between(vector<double>& v1, vector<double>& v2) {
    if (v2.size() != 3 || v1.size() != 3) {
        throw std::invalid_argument("angle_between: Both input vectors must have exactly 3 elements.");
    }
    
    double dot    = inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
    double lenSq1 = inner_product(v1.begin(), v1.end(), v1.begin(), 0.0);
    double lenSq2 = inner_product(v2.begin(), v2.end(), v2.begin(), 0.0);
    
    //DEBUG
    //cout<<"~"<<dot<<","<<lenSq1<<","<<lenSq2<<endl;
    
    double angle = acos(dot/sqrt(lenSq1 * lenSq2));
    return angle;
}

vector<double> crossProduct(vector<double>& a, vector<double>& b) {
    if (a.size() != 3 || b.size() != 3) {
        throw std::invalid_argument("crossProduct: Both input vectors must have exactly 3 elements.");
    }
    vector<double> result(3);
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
}


vector<double> rotateByAngleAndAxis(vector<double>& vec, double angle, vector<double>& axis) {
    if (vec.size() != 3 || axis.size() != 3) {
        throw std::invalid_argument("rotateByAngleAndAxis: Both input vectors must have exactly 3 elements.");
    }
    
    // v_rot = (costheta)v + (sintheta)(axis x v) + (1-cos(theta)) (axis dot v) axis
    vector<double> result(3);
    
    vector<double> axisXvec = crossProduct(axis, vec);
    double axisDOTvec       = inner_product(axis.begin(), axis.end(), vec.begin(), 0.0);
    
    result[0] = cos(angle) * vec[0] + sin(angle) * axisXvec[0] + (1.-cos(angle))*axisDOTvec*axis[0];
    result[1] = cos(angle) * vec[1] + sin(angle) * axisXvec[1] + (1.-cos(angle))*axisDOTvec*axis[1];
    result[2] = cos(angle) * vec[2] + sin(angle) * axisXvec[2] + (1.-cos(angle))*axisDOTvec*axis[2];
    
    return result;
}

void compute_cmos_with_saturation(vector<double>& x_hits_tr,
                                  vector<double>& y_hits_tr,
                                  vector<double>& z_hits_tr,
                                  vector<double>& energy_hits,
                                  map<string,string>& options,
                                  vector<vector<double>>& array2d_Nph
                                  ) {
    
    // vectorized smearing
    vector<double> S3D_x;
    vector<double> S3D_y;
    vector<double> S3D_z;
    cloud_smearing3D(x_hits_tr, y_hits_tr, z_hits_tr, energy_hits, options, S3D_x, S3D_y, S3D_z);
    
    // if there are no electrons on GEM3, just use empty image
    if (S3D_x.size() == 0) return;
    // if there are electrons on GEM3, apply saturation effect
    else {
        // numpy histo is faster than ROOT histo
        //histo_cloud_entries = np.array([S3D_x, S3D_y, S3D_z]).transpose();
        //histo_cloud_entries = histo_cloud_entries[histo_cloud_entries[:, 2].argsort()]
        double xmin = (*min_element(S3D_x.begin(), S3D_x.end()));
        double xmax = (*max_element(S3D_x.begin(), S3D_x.end()));
        double ymin = (*min_element(S3D_y.begin(), S3D_y.end()));
        double ymax = (*max_element(S3D_y.begin(), S3D_y.end()));
        double zmin = (*min_element(S3D_z.begin(), S3D_z.end()));
        double zmax = (*max_element(S3D_z.begin(), S3D_z.end()));
        
        double deltaX = abs(xmax-xmin);
        double deltaY = abs(ymax-ymin);
        
        // FIXME: create a function for the saturation loop
        // FIXME: find best value of maxvolume. 1e8 might not me the best one
        long int max_3Dhisto_volume=(long int)1e+8;  // (volume in number of voxels) that's around 0.5*1.6 GB of RAM
        
        double xbin_dim = stod(options["x_vox_dim"]); //opt.x_dim/opt.x_pix
        double ybin_dim = stod(options["y_vox_dim"]); //opt.y_dim/opt.y_pix
        double zbin_dim = stod(options["z_vox_dim"]);
        
        
        long int x_n_bin = round_up_to_even((xmax-xmin)/xbin_dim);
        long int y_n_bin = round_up_to_even((ymax-ymin)/ybin_dim);
        long int z_n_bin_MAX = round_up_to_even((zmax-zmin)/zbin_dim);
        
        // If voxel regions <= 2 -> use the standard algorithm, otherwise use the map algorithm
        bool map_algorithm = true;
        if( x_n_bin * y_n_bin * z_n_bin_MAX < 2*max_3Dhisto_volume) map_algorithm = false;

        long int z_n_bin;
        if (map_algorithm) z_n_bin = z_n_bin_MAX;
        
        // DEBUG
        cout<<"size = "<<S3D_z.size()<<endl;
        
        double optA                 = stod(options["A"]);
        double optbeta              = stod(options["beta"]);
        double optphotons_per_el    = stod(options["photons_per_el"]);
        double optcounts_per_photon = stod(options["counts_per_photon"]);
        
        vector<vector<int>> hout(x_n_bin, vector<int>(y_n_bin, 0.0));
        
        if(map_algorithm) {
            map<long int, int> hcmap;
            
            vector<size_t> indices(S3D_z.size());
            // Fill indices with 0, 1, 2, ..., numbers.size() - 1
            iota(indices.begin(), indices.end(), 0);
            
            // THIS IS THE COMPUTATIONALLY EXPENSIVE PART
            //long int L = x_n_bin-1;
            long int M = y_n_bin;
            long int N = z_n_bin;
            
            for_each(execution::par,indices.begin(), indices.end(), [&](int ihit) {
                long int xx = floor((S3D_x[ihit]- xmin)/ xbin_dim);
                long int yy = floor((S3D_y[ihit]- ymin)/ ybin_dim);
                long int zz = floor((S3D_z[ihit]- zmin)/ zbin_dim);
                long int map_index = xx * (M * N) + yy * N + zz;
                
                if (hcmap.find(map_index) == hcmap.end()) hcmap[map_index] = 1;
                else hcmap[map_index] += 1;
                
            });
            
            for(auto const& [key, val] : hcmap) {
                int zz = key % N ;
                int yy = ((key - zz) / N) % M;
                int xx = (key - yy * N - zz) / N / M;
                hout[xx][yy]+=Nph_saturation(val, optA, optbeta);
                
            }
            
        } else {
            double deltaZ=max(2*stod(options["z_vox_dim"]),
                              stod(options["z_vox_dim"]) * max_3Dhisto_volume /
                                (deltaX / stod(options["x_vox_dim"])) /
                                (deltaY / stod(options["y_vox_dim"]))
                              );
            
            vector<double> split_vals = arange(zmin, zmax, deltaZ);
            if(split_vals[split_vals.size()-1] < zmax) split_vals.push_back(zmax+deltaZ/10.);
        
        
        
            for(unsigned int i=0; i < split_vals.size()-1; i++) {
                
                // Vector to store all indices
                vector<size_t> allindices(S3D_z.size());
                // Fill indices with 0, 1, 2, ..., numbers.size() - 1
                iota(allindices.begin(), allindices.end(), 0);
                
                // Getting the indices where split_vals[i] <= S3D_z < split_vals[i+1]
                vector<size_t> indices;
                copy_if(allindices.begin(), allindices.end(), back_inserter(indices), [&](size_t n) {
                    return ( S3D_z[n] >= split_vals[i] && S3D_z[n] < split_vals[i+1]);
                });
                
                if(indices.size()==0) continue;
                
                z_n_bin  = max(2.0, round_up_to_even((split_vals[i+1]-split_vals[i])/zbin_dim));
                
                vector<vector<vector<double>>> hc(x_n_bin+1,
                                                  vector<vector<double>>(y_n_bin+1,
                                                                         vector<double>(z_n_bin+1, 0.0)));
                
                
                //cout<<"DEBUG "<<hc.size()<<","<<hc[0].size()<<","<<hc[0][0].size()<<endl<<flush;
                cout<<"Amplifying voxel region z=["<<split_vals[i]<<","<<split_vals[i+1]<<"] "<<i<<"/"<<split_vals.size()-1-1<<endl;
                
                // THIS IS THE COMPUTATIONALLY EXPENSIVE PART
                for_each(execution::par,indices.begin(), indices.end(), [&](int ihit) {
                    int xx = floor((S3D_x[ihit]- xmin)/ xbin_dim);
                    int yy = floor((S3D_y[ihit]- ymin)/ ybin_dim);
                    int zz = floor((S3D_z[ihit]-split_vals[i])/ zbin_dim);
                    hc[xx][yy][zz] += 1.;
                });
                
                long int LL = (x_n_bin)*(y_n_bin)*(z_n_bin);
                long int not_empty=0;
                
                // Applying GEM3 amplification
                for(int xx = 0; xx<x_n_bin-1; xx++){
                    for(int yy=0; yy<y_n_bin-1; yy++) {
                        for(int zz=0; zz<z_n_bin-1; zz++){
                            if(hc[xx][yy][zz] != 0.) {
                                not_empty++;
                                hout[xx][yy]+=Nph_saturation(hc[xx][yy][zz], optA, optbeta);
                            }
                        }
                    }
                }
                
                cout<<"Sparse-ness of voxel region: "<< not_empty * 100./LL<<" %"<<endl;
            }
            
        }
        
        // Applying camera response + Poisson smearing
        for_each(hout.begin(), hout.end(),[&](std::vector<int>& v)  {
            transform (v.begin(), v.end(), v.begin(), [&] (int elem){
                return gRandom->Poisson(elem *
                                        omega *
                                        optphotons_per_el *
                                        optcounts_per_photon);
            });
        });
        
        
        // Padding
        // FIXME: Write a function padding()
        //Define a translation vector
        
        int x_center_cloud=(int)round(((xmax+xmin)/2.)/stod(options["x_vox_dim"]));
        int y_center_cloud=(int)round(((ymax+ymin)/2.)/stod(options["y_vox_dim"]));
        //cout<<"x_center_cloud "<<x_center_cloud<<endl;
        //cout<<"y_center_cloud "<<y_center_cloud<<endl;
        vector<int> translation = {x_center_cloud, y_center_cloud};
        // Calculate the center position of the original array in the padded array
        vector<int> center = {(int)(stod(options["x_pix"])/2.)+translation[0],
                              (int)(stod(options["y_pix"])/2.)+translation[1]
                             };
        // cout<<"Center: "<<center[0]<<", "<<center[1]<<endl;
        int x_start = max(0, center[0] -    (int)hout.size()/2);
        int y_start = max(0, center[1] - (int)hout[0].size()/2);
        int x_end   = min(stoi(options["x_pix"]), x_start + (int)hout.size());
        int y_end   = min(stoi(options["y_pix"]), y_start + (int)hout[0].size());
        // cout<<"PADDING ["<<x_start<<":"<<x_end<<","<<y_start<<":"<<y_end<<"]"<<endl;
        for(int xx=x_start; xx<x_end; xx++){
            for(int yy=y_start; yy<y_end; yy++){
                array2d_Nph[xx][yy]=hout[xx-x_start][yy-y_start];
            }
        }
        
    }
    return;
}

void compute_cmos_without_saturation(vector<double>& x_hits_tr,
                                     vector<double>& y_hits_tr,
                                     vector<double>& z_hits_tr,
                                     vector<double>& energy_hits,
                                     map<string,string>& options,
                                     vector<vector<double>>& array2d_Nph) {
    
    vector<vector<double>> signal(stoi(options["x_pix"]),
                                  vector<double>(stoi(options["y_pix"]), 0.0));
    
    vector<double> S2D_x;
    vector<double> S2D_y;
    ph_smearing2D(x_hits_tr,
                  y_hits_tr,
                  z_hits_tr,
                  energy_hits,
                  options,
                  S2D_x,
                  S2D_y
                  );
    
    // Vector to store all indices
    vector<size_t> indices(S2D_x.size());
    // Fill indices with 0, 1, 2, ..., numbers.size() - 1
    iota(indices.begin(), indices.end(), 0);
    
    double optx_dim = stod(options["x_dim"]);
    double optx_pix = stod(options["x_pix"]);
    double opty_dim = stod(options["y_dim"]);
    double opty_pix = stod(options["y_pix"]);
    
    // THIS IS THE COMPUTATIONALLY EXPENSIVE PART
    for_each(indices.begin(), indices.end(), [&](int ihit) {
        int xx = floor((0.5 * optx_dim + S2D_x[ihit]) * optx_pix / optx_dim);
        int yy = floor((0.5 * opty_dim + S2D_y[ihit]) * opty_pix / opty_dim);
        signal[xx][yy] += 1.;
    });
    
    // DEBUG
    //double ntot =0.;
    //for(unsigned int xx = 0; xx<signal.size(); xx++) {
    //    for(unsigned int yy = 0; yy<signal[0].size(); yy++) {
    //        ntot+=signal[xx][yy];
    //    }
    //}
    //cout<<"Tot num of sensor counts after GEM3 without saturation: "<<ntot<<endl;
    
    array2d_Nph = signal;
    
    return;
}

void cloud_smearing3D(vector<double>& x_hits_tr,
                      vector<double>& y_hits_tr,
                      vector<double>& z_hits_tr,
                      vector<double>& energy_hits,
                      map<string,string>& options,
                      vector<double>& S3D_x,
                      vector<double>& S3D_y,
                      vector<double>& S3D_z) {
    
    vector<double> nel = NelGEM2(energy_hits, z_hits_tr, options);
    //DEBUG
    //for(unsigned int i=0; i<nel.size(); i++) {
    //    cout<<nel[i]<<"\n";
    //}
    
    vector<double> dz;
    int opt_gem=stod(options["z_gem"]);
    transform(z_hits_tr.begin(),z_hits_tr.end(),back_inserter(dz), [&] (double a) { return abs(a-opt_gem);});
    
    vector<double> sigma_x = compute_sigma(stod(options["diff_const_sigma0T"]), stod(options["diff_coeff_T"]), dz);
    vector<double> sigma_y = compute_sigma(stod(options["diff_const_sigma0T"]), stod(options["diff_coeff_T"]), dz);
    vector<double> sigma_z = compute_sigma(stod(options["diff_const_sigma0L"]), stod(options["diff_coeff_L"]), dz);
    
    S3D_x = smear(x_hits_tr, sigma_x, nel);
    S3D_y = smear(y_hits_tr, sigma_y, nel);
    S3D_z = smear(z_hits_tr, sigma_z, nel);
    
    // DEBUG
    //for(unsigned int i=0; i<S3D_x.size(); i++) {
    //    cout<<S3D_x[i]<<endl;
    //}
    
    return;
}

void ph_smearing2D(vector<double>& x_hits_tr,
                   vector<double>& y_hits_tr,
                   vector<double>& z_hits_tr,
                   vector<double>& energy_hits,
                   map<string,string>& options,
                   vector<double>& S2D_x,
                   vector<double>& S2D_y) {
    
    // Electrons in GEM2
    vector<double> nel = NelGEM2(energy_hits, z_hits_tr, options);
    
    double optphotons_per_el    = stod(options["photons_per_el"]);
    double optcounts_per_photon = stod(options["counts_per_photon"]);
    double optA                 = stod(options["A"]);
    // Photons in GEM3 (the factor A is added to be able to compare saturated and non-saturated results)
    vector<double> nph;
    transform(nel.begin(), nel.end(), back_inserter(nph), [&](double nel_i) {
        return nel_i * optA * GEM3_gain * omega * optphotons_per_el * optcounts_per_photon;
    });
    
    vector<double> dz;
    int opt_gem=stod(options["z_gem"]);
    transform(z_hits_tr.begin(),z_hits_tr.end(),back_inserter(dz), [&] (double a) { return abs(a-opt_gem);});
    
    vector<double> sigma_xy = compute_sigma(stod(options["diff_const_sigma0T"]), stod(options["diff_coeff_T"]), dz);
    
    S2D_x = smear(x_hits_tr, sigma_xy, nph);
    S2D_y = smear(y_hits_tr, sigma_xy, nph);
    
    return;
}

vector<double> compute_sigma(const double diff_const, const double diff_coeff, const vector<double>& dz) {
    vector<double> sigmas;
    transform(dz.begin(), dz.end(), back_inserter(sigmas), [&] (double a) {
        return sqrt(diff_const + diff_coeff * a /10.0);
    });
    return sigmas;
}
    
vector<double> smear(const vector<double>& axis_hit, const vector<double>& axis_sigma, const vector<double>& nel) {
    
    long int nelsum = accumulate(nel.begin(), nel.end(), (long int)0);
    
    vector<double> X(nelsum, 0.0);
    
    // Create a vector of indices where each index i is repeated nel[i] times
    vector<long int> indices(nelsum);
    vector<long int> positions(axis_hit.size() + 1, 0);
    
    // Compute cumulative sum of nel to determine positions
    partial_sum(nel.begin(), nel.end(), positions.begin() + 1);
    
    // Fill the indices vector
    for_each(positions.begin(), positions.end() - 1, [&, i = 0](long int pos) mutable {
        fill(indices.begin() + pos, indices.begin() + positions[i + 1], i);
        ++i;
    });
    
    // Fill X with Gaussian-distributed values based on axis_hit and axis_sigma
    transform(indices.begin(), indices.end(), X.begin(), [&](int i) {
        return gRandom->Gaus(axis_hit[i], axis_sigma[i]);
    });
    
    return X;
}

vector<double> arange(double start, double stop, double step) {
    
    int length = (stop - start) / step;
    vector<double> result(length+1);
    double value = start;
    generate(result.begin(), result.end(), [&value, step]() mutable {
        double current = value;
        value += step;
        return current;
    });
    return result;
}

double round_up_to_even(const double f) {
    return ceil(f / 2.0) * 2.0;
}

double Nph_saturation(int nel, double A, double beta) {
    return nel * A * GEM3_gain / (1.0 + beta * GEM3_gain * nel);
}


void TrackVignetting(vector<vector<double>>& arrTr, int xpix, int ypix, const TH2F & VignMap) {
    
    for(int xx = 0; xx < xpix; xx++) {
        for(int yy = 0; yy < ypix; yy++) {
            if(arrTr[xx][yy] != 0) {
                arrTr[xx][yy]=round(arrTr[xx][yy] *
                                    VignMap.GetBinContent(VignMap.GetXaxis()->FindBin(yy),
                                                          VignMap.GetYaxis()->FindBin(xx)
                                                          )
                                    );
            }
        }
    }
    
    return;
}

// Old approach:
//string rootlocation(string tag, int run){
//    string sel;
//    if (tag == "Data") {
//        if ((run>=936) && (run<=1601)) {
//            sel = "Data/LTD/Data_Camera/ROOT";
//        } else if ((run>=1632) && (run<4505)) {
//            sel = "Data/LAB";
//        } else if ((run>=4470) && (run<10000)) {
//            sel = "LAB";
//        } else {
//            cerr<<"rootlocation: Error: Data taken with another DAQ or not yet uploaded to the cloud"<<endl;
//            exit(EXIT_FAILURE);
//        }
//    } else if (tag == "DataMango") {
//        sel = (run<3242) ? "Data/Man" : "MAN";
//    } else if (tag ==  "MC") {
//        sel = "Simulation";
//        cerr<<"rootlocation: Error: automatic download for Simulated data not implemented yet"<<endl;
//        exit(EXIT_FAILURE);
//    }
//
//    return sel;
//}

//string root_filename(string tag, int run) {
//
//    string sel = rootlocation(tag,run);
//
//    string BASE_URL = "https://s3.cloud.infn.it/v1/AUTH_2ebf769785574195bde2ff418deac08a/";
//    string bucket;
//    if(tag.find("MC") != std::string::npos) {
//        bucket = (tag=="MC-old") ? "cygnus" : "cygno-sim";
//    } else if (tag == "Data") {
//        bucket = (run<4505) ? "cygnus" : "cygno-data";
//    } else if (tag == "DataMango") {
//        bucket = (run<3242) ? "cygnus" : "cygno-data";
//    }
//
//    BASE_URL = Form("%s%s/", BASE_URL.c_str(), bucket.c_str());
//    string file_root = Form("%s/histogr_Run%05d.root", sel.c_str(), run);
//
//    return Form("%s%s", BASE_URL.c_str(), file_root.c_str());
//}
//
