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
    FiltersStruct filters;
    Int_t hlt_j80 = 0;

    // input
    // TFile* fi = TFile::Open( input, "read" );
    // if( !fi || fi->IsZombie() ){
    //     std::cerr << "Cannot open " << input << "\n";
    //     return;
    // }

    // ttrees: jet collections, event info, filter, trigger
    const size_t kEvtIdx = nCones;
    const size_t kSkimIdx = nCones + 1;
    const size_t kTrigIdx = nCones + 2;
    std::vector<TChain*> trees( nCones + 3, nullptr );

    for( size_t c = 0; c < nCones; c++ ){
        trees[c] = new TChain( cfg.jetTreePaths[c] );
    }
    trees[kEvtIdx] = new TChain( cfg.hiTreePath );
    trees[kSkimIdx] = new TChain( cfg.skimTreePath );
    trees[kTrigIdx] = new TChain( cfg.trigTreePath );

    for( auto iter : in_list ){
        for( size_t c = 0; c < trees.size(); c++ ){
            TFile *cf = new TFile(iter, "READ");
            TTree *check_tree = ( TTree* )cf->Get( trees[c]->GetName() );
            if( !check_tree ){
                std::cerr << "Missing jet tree " << trees[c]->GetName() << " in " << input << "\n";
                return;
            }
            cf->Close();
            trees[c]->Add(iter);
        }
    }
    // trees[kEvtIdx] = ( TTree* )fi->Get( cfg.hiTreePath );
    // if( !trees[kEvtIdx] ){ std::cerr << "Missing HiTree in " << input << "\n"; return; }
    // if( mode != RunMode::MC ){
    //     trees[kSkimIdx] = ( TTree* )fi->Get( cfg.skimTreePath );
    //     if( !trees[kSkimIdx] ){
    //         std::cerr << "Missing skim tree in " << input << "\n(check cfg/2024ppRef.toml)\n";
    //         return;
    //     }
    // }
    // if( mode == RunMode::Triggered ){
    //     trees[kTrigIdx] = ( TTree* )fi->Get( cfg.trigTreePath );
    //     if( !trees[kTrigIdx] ){
    //         std::cerr << "Missing HLT tree in " << input << "\n(check cfg/2024ppRef.toml)\n";
    //         return;
    //     }
    // }

    // branch mapping
    const bool isMC = ( mode == RunMode::MC );
    SetBranches( trees[kEvtIdx], event.BranchMap( isMC ) );
    for( size_t c = 0; c < nCones; c++ ){ SetBranches( trees[c], jets[c].BranchMap( isMC ) ); }
    if( mode != RunMode::MC ){ SetBranches( trees[kSkimIdx], filters.BranchMap( cfg.filterBranch ) ); }
    if( mode == RunMode::Triggered ){ SetBranches( trees[kTrigIdx], {{cfg.hltJ80Branch, &hlt_j80}} ); }

    // event histograms
    TH1D* hvz_all = new TH1D( "hvz_all", "all events;v_{z} (cm);N", 40, -20, 20 );
    TH1D* hvz = new TH1D( "hvz", "after vz+filter;v_{z} (cm);N", 40, -20, 20 );
    TH1I* hfilt = ( mode != RunMode::MC ) ? new TH1I( "hfilt", "ppvF;filter;N", 2, 0, 2 ) : nullptr;
    TH1D* hcentr = new TH1D( "hcentr", "centrality;filter;N", 100, 0, 100 );
    TH1I* h_j80 = ( mode == RunMode::Triggered ) ? new TH1I( "h_hlt_j80", "HLT_AK4PFJet80;bit;N", 2, 0, 2 ) : nullptr;

    // histograms
    BinningConfig bins;
    std::vector<ConeHistograms> cones( nCones );
    for( size_t c = 0; c < nCones; c++ ){ cones[c].Init( cfg.coneLabels[c], bins, isMC ); }

    // corrected pT buffer: corrPt[cone][jet]
    std::vector<std::vector<float>> corrPt( nCones, std::vector<float>( kNRefMax, 0.0f ) );

    // sorted leading-jet indices
    std::vector<SortedJets> sorted( nCones );

    // event loop
    const Long64_t nEvents = trees[kEvtIdx]->GetEntries();
    const Long64_t nLoop = ( maxEvents < 0 ) ? nEvents : std::min( maxEvents, nEvents );

    for( Long64_t i = 0; i < nLoop; i++ ){

        // vz
        trees[kEvtIdx]->GetEntry( i );
        hvz_all->Fill( event.vz );
        if( TMath::Abs( event.vz ) > kVzCut ){ continue; }

        // ppvF filter
        if( mode != RunMode::MC ){
            trees[kSkimIdx]->GetEntry( i );
            hfilt->Fill( filters.ppvF );
            if( filters.ppvF == 0 ){ continue; }
        }
        hvz->Fill( event.vz );

        // trig
        if( mode == RunMode::Triggered ){
            trees[kTrigIdx]->GetEntry( i );
            h_j80->Fill( hlt_j80 );
        }

        // centrality
        Float_t centrality = findCentrality( event.hiHF_pf );
        hcentr->Fill( centrality );
        if( centrality < cfg.centrLower or centrality > cfg.centrHigher ){ continue; }

        // applying golden JSON
        if( dcs && !dcs->isGood( event.run, event.lumi ) ){ continue; }

        // jet trees (nref)
        for( size_t c = 0; c < nCones; c++ ){ trees[c]->GetEntry( i ); }
        if( jets[trigConeIdx].reco.nref < 2 ){ continue; }

        const float weight = event.w;

        // JEC: fill corrPt[c][j] for every cone and jet
        for( size_t c = 0; c < nCones; c++ ){
            for( int j = 0; j < jets[c].reco.nref; j++ ){
                jecs[c].SetJetPT( jets[c].reco.rawpt[j] );
                jecs[c].SetJetEta( jets[c].reco.eta[j] );
                jecs[c].SetJetPhi( jets[c].reco.phi[j] );
                corrPt[c][j] = ( float )jecs[c].GetCorrectedPT();
            }
        }

        // sort
        for( size_t c = 0; c < nCones; c++ ){ sorted[c] = FindLeadingJets( corrPt[c].data(), jets[c].reco.nref ); }

        // trigger efficiency cut
        if( mode == RunMode::Triggered ){
            if( hlt_j80 == 0 ){ continue; }
            if( hlt_j80 == 1 && corrPt[trigConeIdx][sorted[trigConeIdx].lead] <= cfg.hltJ80Thresh ){ continue; }
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

            // jet ID
            // auto& pf = jets[c].reco.pf;
            // if( !js.JetSelection( jets[c].reco.eta[sorted[c].lead], jets[c].reco.phi[sorted[c].lead],
            //                     pf.CHF[sorted[c].lead], pf.NHF[sorted[c].lead], pf.CEF[sorted[c].lead],
            //                     pf.NEF[sorted[c].lead], pf.MUF[sorted[c].lead],
            //                     pf.CHM[sorted[c].lead], pf.NHM[sorted[c].lead], pf.CEM[sorted[c].lead],
            //                     pf.NEM[sorted[c].lead], pf.MUM[sorted[c].lead] ) ){ continue; }

            // if( !js.JetSelection( jets[c].reco.eta[sorted[c].sublead], jets[c].reco.phi[sorted[c].sublead],
            //                     pf.CHF[sorted[c].sublead], pf.NHF[sorted[c].sublead], pf.CEF[sorted[c].sublead],
            //                     pf.NEF[sorted[c].sublead], pf.MUF[sorted[c].sublead],
            //                     pf.CHM[sorted[c].sublead], pf.NHM[sorted[c].sublead], pf.CEM[sorted[c].sublead],
            //                     pf.NEM[sorted[c].sublead], pf.MUM[sorted[c].sublead] ) ){ continue; }

            bool hasThird = ( sorted[c].third != -1 );
            // if( hasThird ){
            //     hasThird = js.JetSelection( jets[c].reco.eta[sorted[c].third], jets[c].reco.phi[sorted[c].third],
            //                     pf.CHF[sorted[c].third], pf.NHF[sorted[c].third], pf.CEF[sorted[c].third],
            //                     pf.NEF[sorted[c].third], pf.MUF[sorted[c].third],
            //                     pf.CHM[sorted[c].third], pf.NHM[sorted[c].third], pf.CEM[sorted[c].third],
            //                     pf.NEM[sorted[c].third], pf.MUM[sorted[c].third] );
            // }

            // dijet logic (see include/Dijet.h)
            DijetResult dijet = MakeDijet( sorted[c], hasThird, corrPt[c].data(), jets[c].reco.eta, jets[c].reco.phi,
                                          event.event, cfg.minJetPt, cfg.dphiCut );
            if( !dijet.valid ){ continue; }

            // dijet-level |A| acceptance cut
            if( TMath::Abs( dijet.A ) > cfg.maxAbsA ){ continue; }

            // fill histograms
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
    if( hfilt ){ hfilt->Write(); }
    if( h_j80 ){ h_j80->Write(); }
    for( size_t c = 0; c < nCones; c++ ){
        TDirectory* dir = fo->mkdir( cfg.coneLabels[c].Data() );
        cones[c].Write( dir );
        fo->cd();
    }
    fo->Close();
    // fi->Close();
}
