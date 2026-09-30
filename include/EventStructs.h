#ifndef EVENTSTRUCTS_H
#define EVENTSTRUCTS_H

#include "TString.h"
#include "Rtypes.h"

#include <vector>
#include <utility>

struct EventStruct {
    // HFSum
    Float_t hiHF_pf = 0.0f;
    // weight is 1 for data and is mapped to a branch for MC
    Float_t w = 1.0f;
    // vertex position
    Float_t vz;
    // run number
    UInt_t run;
    // event number
    ULong64_t event;
    // lumisection
    UInt_t lumi;
    // centrality
    Int_t hiBin;

    // mapping from variables to branches
    std::vector<std::pair<TString, void*>> BranchMap( bool isMC ){
        std::vector<std::pair<TString, void*>> branches = {
            { "vz", &vz },
            { "evt", &event },
            { "hiHF_pf", &hiHF_pf },
            { "hiBin", &hiBin },
        };
        if( isMC ){
            branches.push_back( { "weight", &w } );
        }
        else {
            branches.insert( branches.end(), {
                { "run", &run },
                { "lumi", &lumi }
            } );
        }
        return branches;
    }
};

struct FilterTriggerStruct {
    // Int_t ppvF;
    std::vector<Int_t> filterTriggers = {1, 1, 1, 1};
    std::vector<std::pair<TString, void*>> BranchMap( const std::vector<TString>& inputBranches ){
        std::vector<std::pair<TString, void*>> vect = {};
        for( size_t i{0}; i < inputBranches.size(); ++i ){
            vect.push_back( { inputBranches[i], &filterTriggers[i] } ); 
        }
        return vect;
    }
};

#endif
