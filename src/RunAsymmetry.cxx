#include "RunAsymmetry.h"

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"
#include "Rtypes.h"

#include "JetCorrector.h"
#include "JetSelection.h"
#include "JSON_handler.h"

#include "BranchMapping.h"
#include "EventStructs.h"
#include "JetStruct.h"
#include "Binning.h"
#include "Dijet.h"
#include "DijetHistograms.h"

#include "AnalysisConfig.h"

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

enum class RunMode { MC, NonTriggered, Triggered };

static constexpr Int_t kNRefMax = 200;
static constexpr float kVzCut = 15.0f;

Float_t findCentralityMC( Float_t HFSum ){
    std::vector<Float_t> energy_bounds = {0, 4.10985, 4.77625, 5.26157, 5.68513, 6.07021, 6.42387, 6.75839, 7.07004, 7.37778, 7.68672, 7.98486, 8.28043, 8.57545, 8.86176, 9.14369, 9.42725, 9.71304, 10.0147, 10.3105, 10.6073, 10.9038, 11.2042, 11.5095, 11.8193, 12.1275, 12.446, 12.7678, 13.095, 13.4164, 13.7444, 14.0892, 14.4311, 14.784, 15.1463, 15.5155, 15.8793, 16.2598, 16.6363, 17.0308, 17.4318, 17.8363, 18.2507, 18.6695, 19.1158, 19.5598, 20.0042, 20.4468, 20.9093, 21.3787, 21.8542, 22.3545, 22.8449, 23.3585, 23.8655, 24.3914, 24.928, 25.4831, 26.0269, 26.6125, 27.1776, 27.7521, 28.3604, 28.9708, 29.5734, 30.2092, 30.8553, 31.5125, 32.1607, 32.8557, 33.5435, 34.248, 34.9785, 35.7214, 36.4689, 37.2423, 38.0336, 38.8336, 39.651, 40.464, 41.3046, 42.1642, 43.0373, 43.9275, 44.8478, 45.784, 46.6997, 47.6568, 48.6178, 49.6245, 50.6359, 51.6727, 52.7001, 53.7518, 54.8302, 55.9117, 57.0425, 58.1931, 59.375, 60.5682, 61.7735, 62.9594, 64.1904, 65.4397, 66.7374, 68.0604, 69.3719, 70.7033, 72.0891, 73.4636, 74.8995, 76.3443, 77.8348, 79.3369, 80.8528, 82.3963, 83.9709, 85.5451, 87.1865, 88.8456, 90.4798, 92.1981, 93.9828, 95.7577, 97.5642, 99.4073, 101.278, 103.22, 105.137, 107.099, 109.138, 111.174, 113.263, 115.385, 117.522, 119.691, 121.983, 124.244, 126.547, 128.864, 131.327, 133.831, 136.317, 138.873, 141.404, 144.05, 146.784, 149.55, 152.275, 155.033, 157.909, 160.833, 163.726, 166.732, 169.794, 172.864, 176.065, 179.353, 182.606, 185.968, 189.345, 192.925, 196.385, 199.963, 203.548, 207.266, 211.009, 214.894, 218.785, 222.834, 226.924, 231.082, 235.36, 239.681, 244.136, 248.668, 253.374, 258.162, 263.07, 268.116, 273.287, 278.555, 284.024, 289.552, 295.42, 301.42, 307.779, 314.293, 321.209, 328.506, 336.377, 344.46, 353.27, 362.966, 373.42, 385.572, 399.896, 416.711, 439.198, 473.479, 752.978};
    for( size_t i{0}; i < energy_bounds.size()-1; ++i ){
        if( HFSum >= energy_bounds.at(i) && HFSum < energy_bounds.at(i+1) ){
            return energy_bounds.size()-1-i;
        }
    }
    // if didn't find bound, return 0
    return 0;
}

Float_t findCentrality( Float_t HFSum ){
    // taken from https://indico.cern.ch/event/1663009/contributions/6991992/attachments/3237906/5774499/o_cent_12mar26.pdf
    std::vector<Float_t> energy_bounds = {0,0.91197,1.82394,2.73591,3.64788,4.55985,5.47182,6.38379,7.29576,8.20773,9.1197,10.0317,10.6667,10.7669,11.1084,11.4511,11.7967,12.1453,12.4985,12.8586,13.2279,13.6047,13.9872,14.378,14.7734,15.1775,15.5875,16.0048,16.4344,16.8729,17.3249,17.7878,18.2572,18.7342,19.2173,19.7139,20.221,20.7368,21.2654,21.7955,22.3333,22.8819,23.446,24.0097,24.5881,25.179,25.7856,26.3914,27.0037,27.6247,28.2619,28.9063,29.5473,30.2005,30.8754,31.5574,32.2518,32.9558,33.6616,34.3925,35.1279,35.8711,36.618,37.3905,38.1709,38.9568,39.76,40.5757,41.3886,42.2135,43.0613,43.921,44.8002,45.6699,46.5569,47.4536,48.3673,49.2971,50.2437,51.1962,52.1656,53.1554,54.1508,55.1587,56.1668,57.1942,58.2508,59.3168,60.406,61.5058,62.6151,63.7513,64.8938,66.0488,67.2111,68.4189,69.6168,70.8422,72.0903,73.3433,74.6371,75.9344,77.2465,78.5776,79.9278,81.3019,82.7002,84.105,85.5389,86.9959,88.4836,89.9969,91.5244,93.0648,94.6373,96.2243,97.8198,99.454,101.105,102.798,104.466,106.186,107.925,109.692,111.487,113.265,115.092,116.936,118.856,120.819,122.768,124.744,126.732,128.739,130.809,132.907,135.033,137.173,139.34,141.543,143.777,146.028,148.298,150.674,153.07,155.464,157.912,160.369,162.846,165.358,167.907,170.503,173.094,175.76,178.495,181.231,183.991,186.784,189.674,192.575,195.507,198.523,201.592,204.653,207.796,210.983,214.245,217.554,220.881,224.256,227.741,231.249,234.786,238.444,242.096,245.842,249.627,253.523,257.577,261.697,265.88,270.16,274.601,279.162,283.888,288.745,293.767,298.968,304.526,310.308,316.408,322.893,329.956,337.457,345.746,354.995,365.685,378.512,395.16,420.836,526.238};
    for( size_t i{0}; i < energy_bounds.size()-1; ++i ){
        if( HFSum >= energy_bounds.at(i) && HFSum < energy_bounds.at(i+1) ){
            return energy_bounds.size()-1-i;
        }
    }
    // if didn't find bound, return 0
    return 0;
}

std::vector<TString> input_files( TString input_file ){
    std::vector<TString> output;
    // if input is a root file, add it as the single element for the chain
    if( input_file.Contains( ".root" ) ){
        output.push_back( input_file );
        return output;
    }
    // otherwise, open file and add each line
    std::ifstream file( input_file );
    if( !file ){
        // return empty vector if no file
        std::cout << "File not read!" << std::endl;
        return output;
    }
    std::string str;
    while( std::getline( file, str ) )
    {
        // if line is empty, skip
        if( str.empty() ){ continue; }
        output.push_back(str);
    }
    return output;
}

void runAsymmetry( TString input, TString output, TString modeFlag, Long64_t maxEvents ){
    std::vector<TString> in_list = input_files( input );
    const AnalysisConfig& cfg = Config();
    PrintConfigSummary( cfg );
    RunMode mode = RunMode::Triggered;
    if( modeFlag == "mc" ) mode = RunMode::MC;
    else if( modeFlag == "non-triggered" ) mode = RunMode::NonTriggered;
    else if( modeFlag == "triggered" ) mode = RunMode::Triggered;
    else {
        throw std::invalid_argument( Form( "ERROR: invalid runAsymmetry mode '%s'. Expected mc, non-triggered, or triggered.", modeFlag.Data() ) );
    }
    // checking configuration
    const size_t nCones = cfg.coneLabels.size();
    if( cfg.jecFilesPerCone.size() != nCones || cfg.jetTreePaths.size() != nCones ){
        std::cerr << "ERROR: cone labels, JEC files, and jet tree paths must be the same length\n";
        return;
    }
    // ?
    size_t trigConeIdx = nCones;
    for( size_t c = 0; c < nCones; c++ ){
        if( cfg.coneLabels[c] == cfg.trigCone ){ trigConeIdx = c; break; }
    }
    if( trigConeIdx == nCones ){
        std::cerr << "ERROR: trigger cone '" << cfg.trigCone << "' not found in cone labels\n";
        return;
    }

    // JetCorrector — chain per-cone L2Residual files onto the base JEC, but
    // only for data (Triggered/NonTriggered); MC must never see residual
    // corrections, so this check is the only thing that decides that, not
    // what's in the TOML.
    std::vector<JetCorrector> jecs;
    jecs.reserve( nCones );
    for( size_t c = 0; c < nCones; c++ ){
        std::vector<std::string> chain = cfg.jecFilesPerCone[c];
        if( mode != RunMode::MC && !cfg.residualFilesPerCone.empty() ){
            for( const auto& f : cfg.residualFilesPerCone[c] ){ chain.push_back( f ); }
        }
        jecs.emplace_back( chain );
    }
    // loading jet ID, jet veto map, and golden json
    JetSelect js( cfg.vetoMapPath, cfg.vetoMapHist );
    std::unique_ptr<JSON_handler> dcs;
    if( mode != RunMode::MC ){ dcs = std::make_unique<JSON_handler>( cfg.jsonPath.Data() ); }

    // structures
    std::vector<JetStruct<kNRefMax>> jets( nCones );
    EventStruct event;
    FilterTriggerStruct filters;
    FilterTriggerStruct triggers;
    // input
    // TFile* fi = TFile::Open( input, "read" );
    // if( !fi || fi->IsZombie() ){
    //     std::cerr << "Cannot open " << input << "\n";
    //     return;
    // }

    // ttrees: jet collections, event info, filter, trigger
    // const size_t kEvtIdx = nCones;
    // const size_t kSkimIdx = nCones + 1;
    // const size_t kTrigIdx = nCones + 2;
    std::vector<TChain*> trees( nCones, nullptr );

    for( size_t c = 0; c < nCones; c++ ){
        trees[c] = new TChain( cfg.jetTreePaths[c] );
    }
    TChain *eventChain = new TChain( cfg.hiTreePath );
    TChain *skimChain  = new TChain( cfg.skimTreePath );
    TChain *trigChain = new TChain( cfg.trigTreePath );
    for( auto iter : in_list ){
        TFile *cf = new TFile(iter, "READ");
        for( size_t c = 0; c < trees.size(); c++ ){
            TTree *check_tree = ( TTree* )cf->Get( trees[c]->GetName() );
            if( !check_tree ){ std::cerr << "Missing jet tree " << trees[c]->GetName() << " in " << input << "\n"; return; }
            trees[c]->Add(iter);
        }
        TTree *check_event = ( TTree* )cf->Get( cfg.hiTreePath );
        if( !check_event ){ std::cerr << "Missing HiTree in " << input << "\n"; return; }
        eventChain->Add( iter );
        TTree *check_skim = ( TTree* )cf->Get( cfg.skimTreePath );
        if( !check_skim ){ std::cerr << "Missing skim tree in " << input << "\n"; return; }
        skimChain->Add( iter );
        if( mode == RunMode::Triggered or mode == RunMode::NonTriggered ){
            TTree *check_trig = ( TTree* )cf->Get( cfg.trigTreePath );
            if( !check_trig ){ std::cerr << "Missing HLT tree in " << input << "\n(check cfg/2024ppRef.toml)\n"; return; }
            trigChain->Add( iter );
        }
        cf->Close();
    }

    // branch mapping
    const bool isMC = ( mode == RunMode::MC );
    SetBranches( eventChain, event.BranchMap( isMC ) );
    for( size_t c = 0; c < nCones; c++ ){ SetBranches( trees[c], jets[c].BranchMap( isMC ) ); }
    if( mode != RunMode::MC ){ SetBranches( skimChain, filters.BranchMap( cfg.filterBranches ) ); }
    if( mode == RunMode::Triggered or mode == RunMode::NonTriggered ){ SetBranches( trigChain, triggers.BranchMap( cfg.hltTriggerBranches ) ); }

    // event histograms
    TH1D* hvz_all = new TH1D( "hvz_all", "all events;v_{z} (cm);N", 40, -20, 20 );
    TH1D* hvz = new TH1D( "hvz", "after vz+filter;v_{z} (cm);N", 40, -20, 20 );
    TH1I* hfilt = new TH1I( "hfilt", "filter combination;filter;N", 2, 0, 2 );
    // TH1I* hfilt_all = new TH1I( "hfilt_all", "filter combination all;filter;N", 2, 0, 2 );
    TH1D* hcentr = new TH1D( "hcentr", "centrality;filter;N", 100, 0, 100 );
    TH1D* hselcentr = new TH1D( "sel_hcentr", "selected centrality;filter;N", 100, 0, 100 );
    TH1I* htriggers = new TH1I( "htrigger", "HLT trigger;bit;N", 2, 0, 2 );
    TH1I* hjetID = new TH1I( "hjetID", "Jet ID; pass;N", 2, 0, 2);
    TH2D* residual_vs_genpt = new TH2D("residual_vs_genpt", "p_{T} residual; gen p_{T}; (p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}", 100, 30, 600, 50, -1, 2);
    TH2D* residual_vs_eta = new TH2D("residual_vs_eta", "p_{T} residual; #eta; (p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}", 100, -5, 5, 50, -1, 2);
    

    TH1D *hpt = new TH1D( "lead_sublead_pt", "lead + sublead pt;p_{T} (GeV); N", 100, 40, 1000 );


    TH1D *reco_selections = new TH1D("evt_selections", "Event satisfying selections; event selections", 10, 0, 10);

    std::vector<TH1D*> dijet_avgpt_histos = {};
    for(size_t i{1}; i < cfg.ptavgEdges.size(); ++i ){
        TH1D* avgpt_dijet = new TH1D( Form("h_dijet_avgpt_%.f_%.f", cfg.ptavgEdges.at(i-1), cfg.ptavgEdges.at(i)), Form("Selected dijet p_{T} %.f-%.f; p_{T}^{avg}; N", cfg.ptavgEdges.at(i-1), cfg.ptavgEdges.at(i)), 200, cfg.ptavgEdges.at(i-1), cfg.ptavgEdges.at(i) );
        dijet_avgpt_histos.push_back(avgpt_dijet);
    }
    

    // histograms
    BinningConfig bins;
    std::vector<ConeHistograms> cones( nCones );
    for( size_t c = 0; c < nCones; c++ ){ cones[c].Init( cfg.coneLabels[c], bins, isMC ); }

    // corrected pT buffer: corrPt[cone][jet]
    std::vector<std::vector<float>> corrPt( nCones, std::vector<float>( kNRefMax, 0.0f ) );

    // sorted leading-jet indices
    std::vector<SortedJets> sorted( nCones );

    // event loop
    const Long64_t nEvents = eventChain->GetEntries();
    const Long64_t nLoop = ( maxEvents < 0 ) ? nEvents : std::min( maxEvents, nEvents );

    for( Long64_t i = 0; i < nLoop; i++ ){
        eventChain->GetEntry( i );
        skimChain->GetEntry( i );
        trigChain->GetEntry( i );
        float weight = event.w;
        // trigger selections
        Bool_t toBeScaled = 1;
        Int_t passTrigger = 0;
        for( size_t t{0}; t < cfg.hltTriggerBranches.size(); ++t ){
            if( triggers.filterTriggers[t] == 1 ){
                passTrigger = 1;
                // if any of the triggers apart from the last one what gets prescaled fires, don't apply PS
                if( t < cfg.hltTriggerBranches.size()-1 ){
                    toBeScaled = 0;
                }
            }
        }

        // if we satisfy only the lowest trigger in the combination, scale by the input PS for the trigger
        if( mode != RunMode::MC and passTrigger and toBeScaled ) weight = event.w*cfg.trigLowestPrescale;

        reco_selections->Fill( "Before any selections", weight );

        htriggers->Fill( passTrigger, weight );

        if( mode != RunMode::MC and !passTrigger ) continue;

        reco_selections->Fill( "Trigger", weight );
        // vz
        hvz_all->Fill( event.vz, weight );
        if( TMath::Abs( event.vz ) > kVzCut ){ continue; }
        hvz->Fill( event.vz, weight );
        reco_selections->Fill("Vz", weight);

        // filters

        Bool_t passFilters = true;
        for( size_t f{0}; f < filters.filterTriggers.size(); ++f ){
            if( filters.filterTriggers[f] == 0 ){ passFilters = false; }
        }
        hfilt->Fill( passFilters, weight );
        if( !passFilters ){
            continue;
        }

        reco_selections->Fill( "Filters", weight );


        // centrality
        Float_t centrality = 0;
        // if it is MC, use the MC centrality calibration...
        if( mode == RunMode::MC ){
            centrality = findCentralityMC( event.hiHF_pf );
            // centrality = event.hiBin;
        }
        // ..., if not MC, use data calibration
        else centrality = findCentrality( event.hiHF_pf );

        hcentr->Fill( centrality, weight );
        // std::cout << "Centrality value: " << centrality << std::endl;

        // *2 because this is the hiBin variable (2x centrality)
        if( !isMC and ( centrality < cfg.centrLower*2. or centrality > cfg.centrHigher*2. ) ){ continue; }
        hselcentr->Fill( centrality, weight );

        reco_selections->Fill( "Centrality", weight );

        // applying golden JSON
        if( mode != RunMode::MC and dcs && !dcs->isGood( event.run, event.lumi ) ){ continue; }
        reco_selections->Fill( "JSON", weight );

        // jet trees (nref)
        for( size_t c = 0; c < nCones; c++ ){ trees[c]->GetEntry( i ); }
        if( jets[trigConeIdx].reco.nref < 2 ){ continue; }
        reco_selections->Fill( "Njets < 2", weight );


        // JEC: fill corrPt[c][j] for every cone and jet
        for( size_t c = 0; c < nCones; c++ ){
            for( int j = 0; j < jets[c].reco.nref; j++ ){
                jecs[c].SetJetPT( jets[c].reco.rawpt[j] );
                jecs[c].SetJetEta( jets[c].reco.eta[j] );
                jecs[c].SetJetPhi( jets[c].reco.phi[j] );
                // std::cout << "For jet with pt=" << jets[c].reco.rawpt[j] << " get factor " << jecs[c].GetCorrectedPT() << std::endl;
                corrPt[c][j] = ( float )jecs[c].GetCorrectedPT();
            }
        }

        // sort
        for( size_t c = 0; c < nCones; c++ ){ sorted[c] = FindLeadingJets( corrPt[c].data(), jets[c].reco.nref ); }
            
        // trigger efficiency cut
        if( mode == RunMode::Triggered ){
            if( corrPt[trigConeIdx][sorted[trigConeIdx].lead] <= cfg.hltJ80Thresh ){ continue; }
        }
        else if( mode == RunMode::NonTriggered ){
            if( corrPt[trigConeIdx][sorted[trigConeIdx].lead] > cfg.hltJ80Thresh ){ continue; }
        }



        // per-cone: jet ID → dijet → |A| → fill
        for( size_t c = 0; c < nCones; c++ ){
            // incl jets: all corrected jets in this cone passing cfg.minJetPt
            for( int j = 0; j < jets[c].reco.nref; j++ ){
                if( corrPt[c][j] >= cfg.minJetPt ){
                    cones[c].FillInclJet( corrPt[c][j], jets[c].reco.eta[j], jets[c].reco.phi[j], weight );
                    if( isMC ){
                        cones[c].FillInclJetResp( corrPt[c][j], jets[c].reco.rawpt[j], jets[c].reco.pt[j],
                                                 jets[c].reco.eta[j], jets[c].ref.pt[j], weight );
                    }
                }
            }

            if( sorted[c].sublead == -1 ){ continue; }

            reco_selections->Fill( "Sorting", weight );

            // jet ID
            auto& pf = jets[c].reco.pf;
            std::vector<Double_t> NHF = {};
            std::vector<Double_t> CHF = {};
            std::vector<Double_t> CEF = {};
            std::vector<Double_t> NEF = {};
            std::vector<Double_t> MUF = {};

            for( Int_t cf{0}; cf < jets[c].reco.nref; ++cf){
                
                NHF.push_back((double)jets[c].reco.pf.neutralSum[cf]/jets[c].reco.rawpt[cf]);
                CHF.push_back((double)jets[c].reco.pf.chargedSum[cf]/jets[c].reco.rawpt[cf]);
                CEF.push_back((double)jets[c].reco.pf.eSum[cf]/jets[c].reco.rawpt[cf]);
                NEF.push_back((double)jets[c].reco.pf.photonSum[cf]/jets[c].reco.rawpt[cf]);
                MUF.push_back((double)jets[c].reco.pf.muSum[cf]/jets[c].reco.rawpt[cf]);
            }
            if( !js.JetID( jets[c].reco.eta[sorted[c].lead], jets[c].reco.phi[sorted[c].lead],
                                CHF[sorted[c].lead], NHF[sorted[c].lead], CEF[sorted[c].lead],
                                NEF[sorted[c].lead], MUF[sorted[c].lead],
                                pf.CHM[sorted[c].lead], pf.NHM[sorted[c].lead], pf.CEM[sorted[c].lead],
                                pf.NEM[sorted[c].lead], pf.MUM[sorted[c].lead] ) ){ hjetID->Fill(int(0), weight); continue; }

            reco_selections->Fill( "JetID 1st jet", weight );

            if( !js.JetID( jets[c].reco.eta[sorted[c].sublead], jets[c].reco.phi[sorted[c].sublead],
                                CHF[sorted[c].sublead], NHF[sorted[c].sublead], CEF[sorted[c].sublead],
                                NEF[sorted[c].sublead], MUF[sorted[c].sublead],
                                pf.CHM[sorted[c].sublead], pf.NHM[sorted[c].sublead], pf.CEM[sorted[c].sublead],
                                pf.NEM[sorted[c].sublead], pf.MUM[sorted[c].sublead] ) ){ hjetID->Fill(int(0), weight); continue; }

            reco_selections->Fill( "JetID 2nd jet", weight );

            if( !js.JetVeto( jets[c].reco.eta[sorted[c].lead], jets[c].reco.phi[sorted[c].lead] ) ){ continue; }
            reco_selections->Fill( "JetVeto 1st jet", weight );

            if( !js.JetVeto( jets[c].reco.eta[sorted[c].sublead], jets[c].reco.phi[sorted[c].sublead] ) ){ continue; }
            reco_selections->Fill( "JetVeto 2nd jet", weight );

            hjetID->Fill(1, weight);
            bool hasThird = ( sorted[c].third != -1 );
            if( hasThird ){
                bool thirdID = js.JetID( jets[c].reco.eta[sorted[c].third], jets[c].reco.phi[sorted[c].third],
                                CHF[sorted[c].third], NHF[sorted[c].third], CEF[sorted[c].third],
                                NEF[sorted[c].third], MUF[sorted[c].third],
                                pf.CHM[sorted[c].third], pf.NHM[sorted[c].third], pf.CEM[sorted[c].third],
                                pf.NEM[sorted[c].third], pf.MUM[sorted[c].third] );
                bool thirdVeto = js.JetVeto( jets[c].reco.eta[sorted[c].third], jets[c].reco.phi[sorted[c].third] );
                if ( thirdID and thirdVeto ) hasThird = true;
                else hasThird = false;
            }

            // the below code replaces the bit above in case we are processing the datasets which contain a different convention for the jetID variables,
            // see also the JetStruct.h for the appropriate branches 
            // jet ID
            // auto& pf = jets[c].reco.pf;
            // if( !js.JetID( jets[c].reco.eta[sorted[c].lead], jets[c].reco.phi[sorted[c].lead],
            //                     pf.CHF[sorted[c].lead], pf.NHF[sorted[c].lead], pf.CEF[sorted[c].lead],
            //                     pf.NEF[sorted[c].lead], pf.MUF[sorted[c].lead],
            //                     pf.CHM[sorted[c].lead], pf.NHM[sorted[c].lead], pf.CEM[sorted[c].lead],
            //                     pf.NEM[sorted[c].lead], pf.MUM[sorted[c].lead] ) ){ continue; }

            // if( !js.JetID( jets[c].reco.eta[sorted[c].sublead], jets[c].reco.phi[sorted[c].sublead],
            //                     pf.CHF[sorted[c].sublead], pf.NHF[sorted[c].sublead], pf.CEF[sorted[c].sublead],
            //                     pf.NEF[sorted[c].sublead], pf.MUF[sorted[c].sublead],
            //                     pf.CHM[sorted[c].sublead], pf.NHM[sorted[c].sublead], pf.CEM[sorted[c].sublead],
            //                     pf.NEM[sorted[c].sublead], pf.MUM[sorted[c].sublead] ) ){ continue; }
            // hjetID->Fill(1, weight);

            // if( !js.JetVeto( jets[c].reco.eta[sorted[c].lead], jets[c].reco.phi[sorted[c].lead] ) ){ continue; }
            // reco_selections->Fill( "JetVeto 1st jet", weight );

            // if( !js.JetVeto( jets[c].reco.eta[sorted[c].sublead], jets[c].reco.phi[sorted[c].sublead] ) ){ continue; }
            // reco_selections->Fill( "JetVeto 2nd jet", weight );
            // bool hasThird = ( sorted[c].third != -1 );
            // if( hasThird ){
            //     bool thirdID = js.JetID( jets[c].reco.eta[sorted[c].third], jets[c].reco.phi[sorted[c].third],
            //                     pf.CHF[sorted[c].third], pf.NHF[sorted[c].third], pf.CEF[sorted[c].third],
            //                     pf.NEF[sorted[c].third], pf.MUF[sorted[c].third],
            //                     pf.CHM[sorted[c].third], pf.NHM[sorted[c].third], pf.CEM[sorted[c].third],
            //                     pf.NEM[sorted[c].third], pf.MUM[sorted[c].third] );
            //     bool thirdVeto = js.JetVeto( jets[c].reco.eta[sorted[c].third], jets[c].reco.phi[sorted[c].third] );
            //     if ( thirdID and thirdVeto ) hasThird = true;
            //     else hasThird = false;
            // }

            if( isMC and ( jets[c].ref.pt[sorted[c].lead] < 15 or jets[c].ref.pt[sorted[c].sublead] < 15 ) ){ continue; };

            reco_selections->Fill( "gen pt veto", weight );
            
            // dijet logic (see include/Dijet.h)
            DijetResult dijet = MakeDijet( sorted[c], hasThird, corrPt[c].data(), jets[c].reco.eta, jets[c].reco.phi,
                                          event.event, cfg.minJetPt, cfg.dphiCut );
            if( !dijet.valid ){ continue; }

            reco_selections->Fill( "Valid dijet", weight );
            // dijet-level |A| acceptance cut
            if( TMath::Abs( dijet.A ) > cfg.maxAbsA ){ continue; }

            reco_selections->Fill( Form("Asymmetry >%.f", cfg.maxAbsA ), weight );

            hpt->Fill( corrPt[c][sorted[c].lead], weight );
            hpt->Fill( corrPt[c][sorted[c].sublead], weight );
            // fill histograms
            for(size_t h{0}; h < dijet_avgpt_histos.size(); ++h){
                if( !dijet_avgpt_histos.at(h)->IsBinOverflow(dijet_avgpt_histos.at(h)->FindBin(dijet.ptavg)) and !dijet_avgpt_histos.at(h)->IsBinUnderflow(dijet_avgpt_histos.at(h)->FindBin(dijet.ptavg)) ){
                    dijet_avgpt_histos.at(h)->Fill(dijet.ptavg, weight);
                    break;
                }
            }
            cones[c].Fill( dijet, corrPt[c].data(), jets[c].reco.eta, jets[c].reco.phi, weight );
            if( isMC ){
                cones[c].FillResp( dijet, corrPt[c].data(), jets[c].reco.rawpt, jets[c].reco.pt,
                                  jets[c].reco.eta, jets[c].ref.pt, weight );
            }
        }
    }

    // output
    { Ssiz_t sl = output.Last( '/' ); if( sl != kNPOS ){ gSystem->mkdir( TString( output( 0, sl ) ), kTRUE ); } }
    TFile* fo = new TFile( output, "recreate" );
    fo->cd();
    hvz_all->Write();
    hvz->Write();
    hcentr->Write();
    hselcentr->Write();
    hjetID->Write();
    hpt->Write();
    reco_selections->Write();
    for(size_t h{0}; h < dijet_avgpt_histos.size(); ++h){
        dijet_avgpt_histos.at(h)->Write();
    }
    if( hfilt ){ hfilt->Write(); }
    if( htriggers ){ htriggers->Write(); }
    for( size_t c = 0; c < nCones; c++ ){
        TDirectory* dir = fo->mkdir( cfg.coneLabels[c].Data() );
        cones[c].Write( dir );
        fo->cd();
    }
    fo->Close();
    // fi->Close();
}
