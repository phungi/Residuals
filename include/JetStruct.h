#ifndef JETSTRUCT_H
#define JETSTRUCT_H

#include "TString.h"
#include "Rtypes.h"

#include <vector>
#include <utility>

// declaring event stucture
template <Int_t MAXNREF>
struct JetStruct {

    // jet variables //
    // more than any number of jets in an event in ttrees being processed
    static constexpr Int_t maxnref = MAXNREF;

    // reco jets
    struct RecoMomenta {
        // number in event
        Int_t nref;
        // momenta
        Float_t rawpt[MAXNREF];
        Float_t pt[MAXNREF];
        Float_t eta[MAXNREF];
        Float_t phi[MAXNREF];
        struct PFInfo {
            // PF fractions
            Float_t CHF[MAXNREF];
            Float_t NHF[MAXNREF];
            Float_t CEF[MAXNREF];
            Float_t NEF[MAXNREF];
            Float_t MUF[MAXNREF];

            Float_t photonSum[MAXNREF];
            Float_t chargedSum[MAXNREF];
            Float_t eSum[MAXNREF];
            Float_t neutralSum[MAXNREF];
            Float_t muSum[MAXNREF];
            // PF multiplicities
            Int_t CHM[MAXNREF];
            Int_t NHM[MAXNREF];
            Int_t CEM[MAXNREF];
            Int_t NEM[MAXNREF];
            Int_t MUM[MAXNREF];

            // Int_t chargedHadN[MAXNREF];
            // Int_t photonN[MAXNREF];
            // Int_t elN[MAXNREF];
            // Int_t neutralN[MAXNREF];
            // Int_t muN[MAXNREF];
        } pf;
    } reco;

    // ref jets
    struct RefMomenta {
        Float_t pt[MAXNREF];
        Float_t eta[MAXNREF];
        Float_t phi[MAXNREF];
    } ref;

    // gen jets
    struct GenMomenta {
        Float_t pt[MAXNREF];
        Float_t eta[MAXNREF];
        Float_t phi[MAXNREF];
    } gen;

    // mapping variables to branches
    std::vector<std::pair<TString, void*>> BranchMap( bool isMC ){
        std::vector<std::pair<TString, void*>> map = {
            { "nref", &reco.nref },
            { "jtptUncorrected", reco.rawpt },
            // the below code replaces the bit above in case we are processing the datasets which contain a different convention for the jetID variables,
            // see also the RunAsymmetry.cxx for the appropriate branches 
            // { "rawpt", reco.rawpt },
            { "jtpt", reco.pt },
            { "jteta", reco.eta },
            { "jtphi", reco.phi },

            { "chargedN", reco.pf.CHM },
            { "neutralN", reco.pf.NHM },
            { "eN", reco.pf.CEM },
            { "photonN", reco.pf.NEM },
            { "muN", reco.pf.MUM },
            {"photonSum", &reco.pf.photonSum},
            {"chargedSum", &reco.pf.chargedSum},
            {"neutralSum", &reco.pf.neutralSum},
            {"eSum", &reco.pf.eSum},
            {"muSum", &reco.pf.muSum}


            // the below code replaces the bit above in case we are processing the datasets which contain a different convention for the jetID variables,
            // see also the RunAsymmetry.cxx for the appropriate branches 

            // { "jtPfCHF", reco.pf.CHF },
            // { "jtPfNHF", reco.pf.NHF },
            // { "jtPfCEF", reco.pf.CEF },
            // { "jtPfNEF", reco.pf.NEF },
            // { "jtPfMUF", reco.pf.MUF },
            // { "jtPfCHM", reco.pf.CHM },
            // { "jtPfNHM", reco.pf.NHM },
            // { "jtPfCEM", reco.pf.CEM },
            // { "jtPfNEM", reco.pf.NEM },
            // { "jtPfMUM", reco.pf.MUM },
        };

        // including gen and ref jet collections iff MC
        if( isMC ){
            map.insert( map.end(), {
                { "refpt", ref.pt },
                { "refeta", ref.eta },
                { "refphi", ref.phi },
                // { "genpt", gen.pt },
                // { "geneta", gen.eta },
                // { "genphi", gen.phi },
            } );
        }
        return map;
    }
};

#endif
