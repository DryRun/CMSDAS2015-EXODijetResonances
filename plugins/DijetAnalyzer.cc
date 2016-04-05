// -*- C++ -*-
//
// Package:    CMSDAS2015/EXODijetResonances
// Class:      DijetAnalyzer
// 
/**\class DijetAnalyzer DijetAnalyzer.cc CMSDAS2015/EXODijetResonances/plugins/DijetAnalyzer.cc

 Description: [one line class summary]

 Implementation:
	  [Notes on implementation]
*/
//
// Original Author:  Dinko Ferencek
//         Created:  Tue, 06 Jan 2015 03:17:30 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

// For JECs
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "TH1F.h"
#include "TH2F.h"

#include "MyTools/AnalysisTools/interface/ObjectSelector.h"
#include "MyTools/AnalysisTools/interface/EventSelector.h"

//
// class declaration
//

class DijetAnalyzer : public edm::EDAnalyzer {
	public:
		explicit DijetAnalyzer(const edm::ParameterSet&);
		~DijetAnalyzer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		template<typename T>
		static bool compare_JetPt(const T& jet1, const T& jet2) {
		  return ( jet1.pt() > jet2.pt() );
		}

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

		// Selectors
		ObjectSelector<pat::Jet>* jet_selector_;
		EventSelector* event_selector_;

		// ----------member data ---------------------------
		const edm::InputTag jetSrc;
		const edm::InputTag vertexSrc;
		const edm::InputTag rhoSrc;
		const edm::InputTag genEventSrc;
		const double JESbias;
		const double dijetMassCut;

		edm::FileInPath L1corr, L2corr, L3corr;
		JetCorrectorParameters *L1Par;
		JetCorrectorParameters *L2Par;
		JetCorrectorParameters *L3Par;
		FactorizedJetCorrector *JetCorrector;

		// file service
		edm::Service<TFileService> fs;

		// histograms to be filled
		TH1F* hVertexZ;
		TH1F* hJetCorrPt;
		TH1F* hJetRawPt;
		TH1F* hJetEta;
		TH1F* hJetPhi;
		TH1F* hJetCHF;
		TH1F* hJetNHF;
		TH1F* hJetEMF;

		TH1F* hRawDijetMass;
		TH1F* hCorrDijetMass;
		TH1F* hJet1Pt;
		TH1F* hJet1Eta;
		TH1F* hJet1Phi;
		TH1F* hJet1CHF;
		TH1F* hJet1NHF;
		TH1F* hJet1EMF;
		TH1F* hJet2Pt;
		TH1F* hJet2Eta;
		TH1F* hJet2Phi;
		TH1F* hJet2CHF;
		TH1F* hJet2NHF;
		TH1F* hJet2EMF;
		TH1F* hDijetDeltaPhi;
		TH1F* hDijetDeltaEta;

		TH2F* hDijetDeltaPhiNJets;
		TH2F* hDijetEta1Eta2;

		TH1F* hEventCount;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DijetAnalyzer::DijetAnalyzer(const edm::ParameterSet& iConfig) :

  jetSrc(iConfig.getParameter<edm::InputTag>("jetSrc")),
  vertexSrc(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  rhoSrc(iConfig.getParameter<edm::InputTag>("rhoSrc")),
  genEventSrc( iConfig.existsAs<edm::InputTag>("genEventSrc") ? iConfig.getParameter<edm::InputTag>("genEventSrc") : edm::InputTag("generator") ),
  JESbias(iConfig.getParameter<double>("JESbias")),
  dijetMassCut(iConfig.getParameter<double>("dijetMassCut"))

{
	//now do what ever initialization is needed
	const int NBINS=90;

	double BOUNDARIES[NBINS+1] = {    1,    3,    6,   10,   16,   23,   31,   40,   50,
												61,   74,   88,  103,  119,  137,  156,  176,  197,
											  220,  244,  270,  296,  325,  354,  386,  419,  453,
											  489,  526,  565,  606,  649,  693,  740,  788,  838,
											  890,  944, 1000, 1058, 1118, 1181, 1246, 1313, 1383,
											 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132,
											 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147,
											 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,
											 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328,
											 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752,
											 9067 };
	TString cut_name;
	std::vector<TString> cut_descriptors;
	std::vector<double> cut_parameters;

	jet_selector_ = new ObjectSelector<pat::Jet>;
	jet_selector_->Configure();
	// Add cuts
	cut_parameters.clear();
	cut_descriptors.clear();
	cut_name = "MinPt";
	cut_parameters.push_back(25);
	jet_selector_->RegisterCut(cut_name, cut_descriptors, cut_parameters);

	cut_parameters.clear();
	cut_descriptors.clear();
	cut_name = "MaxAbsEta";
	cut_parameters.push_back(2.5);
	jet_selector_->RegisterCut(cut_name, cut_descriptors, cut_parameters);

	event_selector_ = new EventSelector();
	event_selector_->Configure();
	event_selector_->AddJetSelector(jet_selector_);

	cut_name = "VertexMinNdof";
	cut_parameters.clear();
	cut_descriptors.clear();
	cut_parameters.push_back(4);
	event_selector_->RegisterCut(cut_name, cut_descriptors, cut_parameters);

	cut_name = "VertexMaxAbsZ";
	cut_parameters.clear();
	cut_descriptors.clear();
	cut_parameters.push_back(24.);
	event_selector_->RegisterCut(cut_name, cut_descriptors, cut_parameters);

	cut_name = "VertexMaxRho";
	cut_parameters.clear();
	cut_descriptors.clear();
	cut_parameters.push_back(2.);
	event_selector_->RegisterCut(cut_name, cut_descriptors, cut_parameters);

	cut_name = "MinNJets";
	cut_parameters.clear();
	cut_descriptors.clear();
	cut_parameters.push_back(2.);
	event_selector_->RegisterCut(cut_name, cut_descriptors, cut_parameters);

	// setup histograms
	hVertexZ   = fs->make<TH1F>("hVertexZ",   "Z position of the Primary Vertex;z [cm]", 60, -30, 30);
	hVertexZ->Sumw2();
	hVertexZ->SetDefaultSumw2(kTRUE);
	hJetRawPt  = fs->make<TH1F>("hJetRawPt",  "Raw Jet p_{T}; p_{T} [GeV]",           100, 0, 5000);
	hJetCorrPt = fs->make<TH1F>("hJetCorrPt", "Corrected Jet p_{T}; p_{T} [GeV]",     100, 0, 5000);
	hJet1Pt    = fs->make<TH1F>("hJet1Pt",    "Corrected Jet_{1} p_{T}; p_{T} [GeV]", 100, 0, 5000);
	hJet2Pt    = fs->make<TH1F>("hJet2Pt",    "Corrected Jet_{2} p_{T}; p_{T} [GeV]", 100, 0, 5000);

	hJetEta  = fs->make<TH1F>("hJetEta",  "Corrected Jet #eta;#eta",     100, -5, 5);
	hJet1Eta = fs->make<TH1F>("hJet1Eta", "Corrected Jet_{1} #eta;#eta", 100, -5, 5);
	hJet2Eta = fs->make<TH1F>("hJet2Eta", "Corrected Jet_{2} #eta;#eta", 100, -5, 5);

	hJetPhi  = fs->make<TH1F>("hJetPhi",  "Corrected Jet #phi;#phi",     50, -3.1416, 3.1416);
	hJet1Phi = fs->make<TH1F>("hJet1Phi", "Corrected Jet_{1} #phi;#phi", 50, -3.1416, 3.1416);
	hJet2Phi = fs->make<TH1F>("hJet2Phi", "Corrected Jet_{2} #phi;#phi", 50, -3.1416, 3.1416);

	hJetCHF  = fs->make<TH1F>("hJetCHF",  "CH Fraction of Jets;CH Fraction",    100, 0, 1);
	hJet1CHF = fs->make<TH1F>("hJet1CHF", "CH Fraction of Jet_{1};CH Fraction", 100, 0, 1);
	hJet2CHF = fs->make<TH1F>("hJet2CHF", "CH Fraction of Jet_{2};CH Fraction", 100, 0, 1);

	hJetNHF  = fs->make<TH1F>("hJetNHF",  "NH Fraction of Jets;NH Fraction",    100, 0, 1);
	hJet1NHF = fs->make<TH1F>("hJet1NHF", "NH Fraction of Jet_{1};NH Fraction", 100, 0, 1);
	hJet2NHF = fs->make<TH1F>("hJet2NHF", "NH Fraction of Jet_{2};NH Fraction", 100, 0, 1);

	hJetEMF  = fs->make<TH1F>("hJetEMF",  "EM Fraction of Jets;EM Fraction",    100, 0, 1);
	hJet1EMF = fs->make<TH1F>("hJet1EMF", "EM Fraction of Jet_{1};EM Fraction", 100, 0, 1);
	hJet2EMF = fs->make<TH1F>("hJet2EMF", "EM Fraction of Jet_{2};EM Fraction", 100, 0, 1);

	hRawDijetMass  = fs->make<TH1F>("hRawDijetMass",  "Raw Dijet Mass;m_{jj} [GeV]", NBINS, BOUNDARIES);
	hCorrDijetMass = fs->make<TH1F>("hCorrDijetMass", "Corrected Dijet Mass;m_{jj} [GeV]", NBINS, BOUNDARIES);
	hDijetDeltaPhi = fs->make<TH1F>("hDijetDeltaPhi", "Dijet |#Delta#phi|;|#Delta#phi|",  50, 0, 3.1416);
	hDijetDeltaEta = fs->make<TH1F>("hDijetDeltaEta", "Dijet |#Delta#eta|;|#Delta#eta|",  25, 0,    1.5);

	hDijetDeltaPhiNJets = fs->make<TH2F>("hDijetDeltaPhiNJets", "Dijet |#Delta#phi| vs the number of jets;|#Delta#phi|;Number of jets", 50,  0, 3.1416,  7, 0.5, 7.5);
	hDijetEta1Eta2      = fs->make<TH2F>("hDijetEta1Eta2",      "#eta_{1} vs #eta_{2} of dijet events;#eta_{1};#eta_{2}",               50, -5,      5, 50,  -5,   5);

	hEventCount = fs->make<TH1F>("hEventCount", "hEventCount", 4, 0.5, 4.5);

	// For JEC
	L1corr = iConfig.getParameter<edm::FileInPath>("L1corr");
	L2corr = iConfig.getParameter<edm::FileInPath>("L2corr");
	L3corr = iConfig.getParameter<edm::FileInPath>("L3corr");

	L1Par = new JetCorrectorParameters(L1corr.fullPath());
	L2Par = new JetCorrectorParameters(L2corr.fullPath());
	L3Par = new JetCorrectorParameters(L3corr.fullPath());

	std::vector<JetCorrectorParameters> vPar;
	vPar.push_back(*L1Par);
	vPar.push_back(*L2Par);
	vPar.push_back(*L3Par);

	JetCorrector = new FactorizedJetCorrector(vPar);
}


DijetAnalyzer::~DijetAnalyzer()
{
 
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

// for PF jet ID
PFJetIDSelectionFunctor pfjetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE );
PFJetIDSelectionFunctor pfjetIDTight( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT );

pat::strbitset retpf = pfjetIDLoose.getBitTemplate();

//
// member functions
//

// ------------ method called for each event  ------------
void
DijetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	/////// Get weight, if this is MC //////////
	double mWeight = 1.;
	// magic to get the GenEventInfoProduct from EDM
	edm::Handle<GenEventInfoProduct> eventInfo;
	iEvent.getByLabel(genEventSrc, eventInfo);
	if( eventInfo.isValid() )
	  mWeight = eventInfo->weight();

	hEventCount->Fill(1.);          // counts unweighted events
	hEventCount->Fill(2., mWeight); // counts weighted events

	////////////////////////////////////////////
	// Get Primary Vertex Information
	////////////////////////////////////////////

	// magic to get the vertices from EDM
	edm::Handle<std::vector<reco::Vertex> > vertices;
	iEvent.getByLabel(vertexSrc, vertices);
	if( !vertices.isValid() ) return;

	// require in the event that there is at least one reconstructed vertex
	if( !(vertices->size()>0) ) return;

	// pick the first (i.e. the highest sum pt^2) vertex
	const reco::Vertex* theVertex = &(vertices->front());

	// require that the vertex meets certain criteria
	if( !(theVertex->ndof()>=4. &&
			std::abs(theVertex->z())<=24. &&
			std::abs(theVertex->position().rho())<=2.) ) return;

	// fill the vertex histogram
	hVertexZ->Fill(theVertex->z(), mWeight);

	////////////////////////////////////////////
	// Get rho density for JEC
	////////////////////////////////////////////

	// magic to get the rho from EDM
	edm::Handle<double> rho;
	iEvent.getByLabel(rhoSrc,rho);

	////////////////////////////////////////////
	// Get Jet Information
	////////////////////////////////////////////

	// magic to get the jets from EDM
	edm::Handle<std::vector<pat::Jet> > jets;
	iEvent.getByLabel(jetSrc, jets);
	if( !jets.isValid()) return;

	jet_selector_->ClassifyObjects(jets);
	event_selector_->ProcessEvent(&iEvent, &vertices);

	// collection of selected jets
	std::vector<pat::Jet> selectedJets;

	// loop over the jet collection
	for (auto& it_jet : jet_selector_->GetGoodObjectIndices()) {
	//for(std::vector<pat::Jet>::const_iterator j_it = jets->begin(); j_it != jets->end(); ++j_it) {
		pat::Jet this_jet = (*jets)[it_jet];
		this_jet.scaleEnergy(JESbias);

		hJetCorrPt->Fill(this_jet.pt(), mWeight);
		hJetRawPt->Fill(this_jet.correctedJet("Uncorrected").pt(), mWeight);
		hJetEta->Fill(this_jet.eta(), mWeight);
		hJetPhi->Fill(this_jet.phi(), mWeight);
		hJetCHF->Fill(this_jet.chargedHadronEnergyFraction(), mWeight);
		hJetNHF->Fill(this_jet.neutralHadronEnergyFraction(), mWeight);
		hJetEMF->Fill(this_jet.electronEnergyFraction() + this_jet.photonEnergyFraction(), mWeight);

	}

	hEventCount->Fill(3.);          // counts unweighted events after full event selection
	hEventCount->Fill(4., mWeight); // counts weighted events after full event selection

	return;
}


// ------------ method called once each job just before starting event loop  ------------
void 
DijetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DijetAnalyzer::endJob() 
{
	jet_selector_->MakeCutflowHistograms(&(*fs));
	event_selector_->MakeCutflowHistograms(&(*fs));
}

// ------------ method called when starting to processes a run  ------------
/*
void 
DijetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DijetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DijetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DijetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DijetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DijetAnalyzer);
