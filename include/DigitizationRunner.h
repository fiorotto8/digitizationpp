/*
 * Copyright (C) 2025 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2025
 *
 */

#pragma once

#include "ConfigManager.h"
#include "Globals.h"
#include "cygnolib.h"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <utility>
#include "TH2F.h"

// Forward declarations
class TFile;

/**
 * @class DigitizationRunner
 * @author Stefano Piacentini
 * @brief Manages the digitization process from MC simulation input to ROOT image output.
 *
 * @details This class initializes the configuration file, sets up camera parameters, processes ROOT files,
 * and handles digitization including background addition, and application of vignetting maps.
 */
class DigitizationRunner {
public:
    /**
     * @brief Constructor.
     * @param[in] configFile Path to the configuration file.
     * @param[in] inputDir Path to the directory containing input ROOT files.
     * @param[in] outputDir Path to the directory where output will be written.
     */
    DigitizationRunner(const std::string& configFile, const std::string& inputDir, const std::string& outputDir);

    /**
     * @brief Executes the full digitization pipeline.
     */
    void run();

    /**
     * @brief Executes the addition of random pedestals to input "digi_RunXXXXX.root" file.
     */
     void runPedsOnly();

private:

    /**
     * @brief Initializes global variables like GEM gains and extraction efficiencies.
     */
     void initializeGlobals();

     /**
      * @brief Prepares camera parameters such as pixel size, optics, and sensor dimensions.
      */
     void prepareCameraSettings();
 
     /**
      * @brief Validates consistency between simulation dimensions and camera dimensions.
      */
     void checkDimensionConsistency();
 
     /**
      * @brief Loops through ROOT files in the input directory and processes them.
      */
     void processRootFiles();
 
     /**
      * @brief Checks whether a given file is a valid input ROOT file.
      * @param[in] filename Name of the file to check.
      * @return True if valid, false otherwise.
      */
     bool isValidInputFile(const std::string& filename) const;
    
     /**
     * @brief Checks if the input filename matches the expected digi file format.
     *
     * Validates whether the provided filename follows the naming convention:
     * `digi_RunXXXXX.root`, where `XXXXX` is a number with 5 or more digits.
     *
     * @param[in] filename The name of the file to validate.
     * @return true if the filename matches the expected format; false otherwise.
     *
     * @note This check uses a regular expression internally.
     */
     bool isValidDigiFile(const std::string& filename) const;

     /**
     * @brief Extracts the run number from a digi file name.
     *
     * Parses a filename of the form `digi_RunXXXXX.root` and extracts the numeric
     * run number (where `XXXXX` is a number with 5 or more digits).
     *
     * @param[in] filename The digi file name to parse.
     * @return The extracted run number as an integer.
     *
     * @throws std::invalid_argument if the filename does not match the expected format.
     *
     * @see isValidDigiFile() to validate the filename before calling this.
     */
     int extractDigiRunNumber(const std::string& filename);

     /**
      * @brief Sets the random seed for reproducibility.
      * @param[in] seed Seed value (default is 10).
      */
     void setSeed(int seed = 10);
 
     /**
      * @brief Initializes the source directory based on configFile variable.
      */
     void initSourceDir();
 
     /**
      * @brief Saves global and run parameters into the output ROOT file.
      * @param[in] outfile Shared pointer to the output TFile.
      */
     void SaveValues(std::shared_ptr<TFile>& outfile);
 
     /**
      * @brief Fills the vignetting map histogram.
      * @param[in] VignMap Reference to the 2D histogram containing the vignetting map.
      */
     void fillVigmap(TH2F& VignMap);
 
     /**
      * @brief Determines whether an event is a nuclear recoil (NR).
      * @param[in] pdgID_hits List of PDG codes for all hits.
      * @param[in] pdg PDG code to compare against.
      * @return PDG code of the first matching recoil, or -999 if none.
      */
     int is_NR(std::vector<int> pdgID_hits, int pdg);
 
     /**
      * @brief Adds pedestal to an image, according the configfile options.
      * @param[in] background Reference to 2D image matrix that will contain the pedestal.
      */
     void AddBckg(std::vector<std::vector<int>>& background);

     /**
     * @brief Adds pedestal histograms to the "param_dir" directory of the output ROOT file.
     *
     * This function assumes that the output ROOT file (`outfile`) already contains a
     * directory named "param_dir", and that this directory holds existing TH1 histograms.
     * It adds pedestal run numbers to `pedestal_runs` TH1F.
     *
     *
     * @param[in] outfile A shared pointer to the output ROOT file where the histograms
     *                    will be modified or augmented.
     *
     * @note This function modifies the ROOT file in-place, and the changes are saved
     *       directly into `outfile`.
     */
     void addPedsToParamDir(std::shared_ptr<TFile>& outfile);

     /**
     * @brief Parses a MC axis label string and returns the corresponding digititazion axis and sign.
     *
     * This function interprets a string like "x", "-y", or "z" to determine which
     * digitization axis corresponds to a MC axis, and whether that axis is inverted.
     *
     * For example:
     * - "x"   → returns ('x', 1.0)
     * - "-z"  → returns ('z', -1.0)
     *
     * @param axis_label A string representing the axis mapping from the config file.
     *                   Expected format is one of: "x", "y", "z", "-x", "-y", or "-z".
     *
     * @return A pair where:
     *         - first: the digitization axis as a char ('x', 'y', or 'z')
     *         - second: a double representing the sign (1.0 or -1.0)
     */
     std::pair<char, double> getAxisMapping(const std::string& axis_label);


     /**
     * @brief Fills redpix data from a 2D image.
     *
     * Scans all bins of the input vector<vector<double>> and identifies bins with non-zero content.
     * For each such bin, the method stores:
     * - the X bin index in `redpix_ix`
     * - the Y bin index in `redpix_iy`
     * - the bin content in `redpix_iz`
     *
     * All output vectors are cleared at the beginning of the method.
     *
     * @param image     A vector<vector<double>> representing the image to analyze.
     * @param redpix_ix Pointer to a vector to be filled with X bin indices of non-zero pixels.
     * @param redpix_iy Pointer to a vector to be filled with Y bin indices of non-zero pixels.
     * @param redpix_iz Pointer to a vector to be filled with bin contents (e.g., charge or intensity).
     */
     void FillRedpix(const std::vector<std::vector<double>>& image,
                     std::vector<uint16_t>* redpix_ix,
                     std::vector<uint16_t>* redpix_iy,
                     std::vector<uint16_t>* redpix_iz);

    /**
     * @brief Generates TH2I from a digi ROOT file.
     *
     * This function processes a ROOT file named `digi_RunXXXXX.root`, where `XXXXX`
     * represents the run number, and creates a new output ROOT file called
     * `histograms_RunXXXXX.root`. The output file contains TH2F histograms that
     * visualize the tracks recorded in the input digi file.
     *
     * @note The function assumes the input file follows the naming convention
     *       `digi_RunXXXXX.root` and outputs histograms accordingly.
     */
     void generateHistogramsFromDigi();

 
     ConfigManager config;           ///< Configuration manager
     std::string configFile;         ///< Path to the config file
     std::string inputDir;           ///< Path to input directory
     std::string outputDir;          ///< Path to output directory
     int runCount = 0;               ///< Current run number
     int NMAX_EVENTS = 200;          ///< Max number of events to process per output file
     std::string SOURCE_DIR = "";    ///< Path to event source directory
 
};
