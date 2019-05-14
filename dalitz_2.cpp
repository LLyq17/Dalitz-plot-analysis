//modified to analysis: D+ -> Ks K+ pi0
//xiexh 2019 04 13

// ROOT stuff
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TMath.h>
#include <TCut.h>

// System stuff
#include <fstream>
#include <sys/time.h>
#include <sys/times.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/utilities/Style.h>

using namespace std;
using namespace GooFit;
Variable fixedf2Mass("f2_mass", 1.2755, 0.01, 1.2, 1.4);
Variable fixedf2Width("f2_width", 0.1867, 0.001, 1e-5, 5e-1);



const fptype _mDs = 1.96828;
const fptype _mDs2      = _mDs * _mDs;
const fptype _mDs2inv   =1. / _mDs2;
const fptype KPlusMass = 0.493677;
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;


// Constants used in more than one PDF component.

Variable motherM("motherM",  _mDs);
Variable dau1M("dau1M", piPlusMass);
Variable dau2M("dau2M", piPlusMass);
Variable dau3M("dau3M", piPlusMass);
Variable massSum("massSum",  _mDs *_mDs + 3 * piPlusMass * piPlusMass); // = 3.53481 
Variable constantOne("constantOne", 1);
Variable constantZero("constantZero", 0);

fptype cpuGetM23(fptype massPZ, fptype massPM) {
    return (_mDs2 + piPlusMass * piPlusMass + piPlusMass * piPlusMass + piPlusMass * piPlusMass - massPZ - massPM);
}

void getToyData(std::string toyFileName, GooFit::Application &app, DataSet &data) {
    toyFileName = app.get_filename(toyFileName, "examples/dalitz");

    auto obs               = data.getObservables();
    Observable m12         = obs.at(0);
    Observable m13         = obs.at(1);
    Observable eventNumber = obs.at(2);

    TH2F dalitzplot("dalitzplot",
                    "Original Data",
                    m12.getNumBins(),
                    m12.getLowerLimit(),
                    m12.getUpperLimit(),
                    m13.getNumBins(),
                    m13.getLowerLimit(),
                    m13.getUpperLimit());
    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);


      std::cout<<"Reading file "<<toyFileName<<std::endl;
      TFile*f = TFile::Open(toyFileName.c_str());
      TTree*t = (TTree*)f->Get("DecayTree");
      assert(t);
      std::cout<<"Entries: "<<t->GetEntries()<<std::endl;
      double m2_12, m2_13;
      t->SetBranchAddress("s12", &m2_12);
      t->SetBranchAddress("s13", &m2_13);
      for (int i=0;i<t->GetEntries();i++){
          t->GetEntry(i);
           m12=m2_12;
          m13=m2_13;
         eventNumber.setValue(data.getNumEvents());
//          eventNumber.setValue()= data.getNumEvents();
          data.addEvent();
          dalitzplot.Fill(m12.getValue(), m13.getValue());
      }
      f->Close();

    GOOFIT_INFO("Read in {} events", data.getNumEvents());

    TCanvas foo;
    dalitzplot.SetStats(false);
   dalitzplot.Draw("colz");
    foo.SaveAs("dalitzplot.png");

}
std::vector<fptype> HH_bin_limits;
std::vector<Variable> pwa_coefs_reals;
std::vector<Variable> pwa_coefs_imags;

char strbuffer[1000];

void name(const string fname = "PiPiPi_PWA_COEFFS_30_SAda_ampphs.txt") {

std::ifstream reader;
reader.open(fname.c_str());
assert(reader.good());
HH_bin_limits.clear();
pwa_coefs_reals.clear();
pwa_coefs_imags.clear();


//double Rmag = 1., Rphs = 1.;

int i=0;
double e1,e2,e3;
while (reader >> e1 >> e2 >> e3 ) {
      HH_bin_limits.push_back(e1);

      sprintf(strbuffer, "pwa_coef_%d_mag", i);
      Variable vr(strbuffer, e2, 0.01, 0, 0);//0, 10*e2);

      sprintf(strbuffer, "pwa_coef_%d_phase", i);
      Variable vi(strbuffer, e3, 0.01, 0, 0);//-2*TMath::Pi(), 2.5*TMath::Pi());

      pwa_coefs_reals.push_back(vr);
      pwa_coefs_imags.push_back(vi);
      i++;
  }
}

void makeToyData(DalitzPlotter &dplotter, UnbinnedDataSet &data) {}

Amp3Body *makeSignalPdf(Observable m12, Observable m13, EventNumber eventNumber, GooPdf *eff = 0) {
    DecayInfo3 dtop0pp;
    dtop0pp.motherMass   = _mDs;
    dtop0pp.daug1Mass    = piPlusMass;
    dtop0pp.daug2Mass    = piPlusMass;
    dtop0pp.daug3Mass    = piPlusMass;
    dtop0pp.meson_radius = 1.5;

  bool fixAmps = false; // Takes ~400x longer
  bool fixAmp =true;
  Variable rho0_mass("rho0_mass", 0.77526, 0.01, 0.4, 1.1);
  Variable rho0_width("rho0_width", 0.1491, 0.01, 0.1, 0.3);

    ResonancePdf *rho_1 = new Resonances::RBW(
        "rho0",
        fixAmps ? Variable("rho0_amp_real1", 0.19) : Variable("rho0_amp_real1", 0.19, 0.001, 0, 0),
        fixAmps ? Variable("rho0_amp_imag1", -1.1) : Variable("rho0_amp_imag1", -1.1, 0.1, 0, 0),
        rho0_mass,
        rho0_width,
        1,
        PAIR_12 ,true);



   Variable sharedMass("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
   Variable shareWidth("rhop_1450_width", 0.400, 0.01, 0.01, 5.0);


    ResonancePdf *rho0_1450_1 = new Resonances::RBW(
        "rho0_1450",
        fixAmps ? Variable("rho0_1450_amp_real1", 1.2) : Variable("rho0_1450_amp_real1",1, 0.01, 0, 0),
        fixAmps ? Variable("rho0_1450_amp_imag1", -4.1) : Variable("rho0_1450_amp_imag1", 0, 0.01, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_12,true);


    ResonancePdf *f2_1270_1 = new Resonances::RBW(
        "f2_1270",
         fixAmp ? Variable("f2_1270_amp_real1", 1.0)  : Variable("f2_1270_amp_real1", 1.0, 0.01, 0, 0),
         fixAmp ? Variable("f2_1270_amp_imag1", 0)    : Variable("f2_1270_amp_imag1", 0, 0.01, 0, 0),
         fixedf2Mass,
         fixedf2Width,
         2,
         PAIR_12,true);



  const fptype scale = 1;
    ResonancePdf *swave1 = new Resonances::Spline(
          "swave",
          fixAmp ? Variable("swave_amp_real1", 1)   : Variable("swave_amp_real1", -1.86396*scale, 0.001, 0, 0),
          fixAmp ? Variable("swave_amp_imag1",0 )    : Variable("swave_amp_imag1", .775892256*scale, 0.1, 0, 0),
          HH_bin_limits,
          pwa_coefs_reals,
          pwa_coefs_imags,
          PAIR_12,
          true);


    dtop0pp.resonances.push_back(rho_1);

   dtop0pp.resonances.push_back(rho0_1450_1);
   dtop0pp.resonances.push_back(f2_1270_1);

    bool fitMasses = false;

    if(!fitMasses) {
        for(vector<ResonancePdf *>::iterator res = dtop0pp.resonances.begin(); res != dtop0pp.resonances.end(); ++res) {
            (*res)->setParameterConstantness(true);
        }
    }

  dtop0pp.resonances.push_back(swave1);


    if(!eff) {
        // By default create a constant efficiency.
        vector<Variable> offsets;
        vector<Observable> observables;
        vector<Variable> coefficients;

        observables.push_back(m12);
        observables.push_back(m13);
        offsets.push_back(constantZero);
        offsets.push_back(constantZero);
        coefficients.push_back(constantOne);
        eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
    }


    return new Amp3Body("signalPDF", m12, m13, eventNumber, dtop0pp,eff);
}


int runToyFit(Amp3Body *signal, UnbinnedDataSet *data) {
    // EXERCISE 1 (real part): Create a PolynomialPdf which models
    // the efficiency you imposed in the preliminary, and use it in constructing
    // the signal PDF.

    // EXERCISE 2: Create a K0 veto function and use it as the efficiency.

    // EXERCISE 3: Make the efficiency a product of the two functions
    // from the previous exercises.
	auto obs               = data->getObservables();
	Observable m12         = obs.at(0);
	Observable m13         = obs.at(1);
	Observable eventNumber = obs.at(2);

	signal->setData(data);
	signal->setDataSize(data->getNumEvents());
	FitManager datapdf(signal);
	datapdf.fit();

	ProdPdf prodpdf{"prodpdf", {signal}};

	DalitzPlotter plotter(&prodpdf, signal);

	//Draw option
	TCanvas *foo = new TCanvas("c1","c1",800,700);
	foo->Divide(2,2);

	TH2F *dalitzplot = plotter.make2D();
	dalitzplot->GetXaxis()->SetTitle("M^{2}(K^{+}#pi^{0}) (GeV^{2}/c^{4})");
	dalitzplot->GetYaxis()->SetTitle("M^{2}(K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");
/*
	fptype m12_tmp,m13_tmp,m23_tmp;
	TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
	m12_pdf_hist.SetStats(false); 
	m12_pdf_hist.SetLineColor(kBlue); 
	m12_pdf_hist.SetLineWidth(2); 
	m12_pdf_hist.GetXaxis()->SetTitle("M^{2}(K^{+}#pi^{0}) (GeV^{2}/c^{4})");

	TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	m13_pdf_hist.SetStats(false); 
	m13_pdf_hist.SetLineColor(kBlue); 
	m13_pdf_hist.SetLineWidth(2); 
	m13_pdf_hist.GetXaxis()->SetTitle("M^{2}(K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");

	TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	m23_pdf_hist.SetStats(false); 
	m23_pdf_hist.SetLineColor(kBlue); 
	m23_pdf_hist.SetLineWidth(2); 
	m23_pdf_hist.GetXaxis()->SetTitle("M^{2}(K^{+}K_{S}^{0}) (GeV^{2}/c^{4})");

	UnbinnedDataSet toyMC({m12, m13, eventNumber});
	plotter.fillDataSetMC(toyMC, 5000000);
	cout << "toyMC size = " << eventNumber.getValue() << endl;
	for (int i = 0; i < eventNumber.getValue(); ++i){
		m12_tmp = toyMC.getValue(m12, i);
		m12_pdf_hist.Fill(m12_tmp);

		m13_tmp = toyMC.getValue(m13, i);
		m13_pdf_hist.Fill(m13_tmp);

		m23_tmp = cpuGetM23(m12_tmp, m13_tmp);
		m23_pdf_hist.Fill(m23_tmp);
	}
*/
	//for comparison
	TChain *tr_data = new TChain("DecayTree");tr_data->Add("FastMC_DsPiPiPi_PWA_11K_30pt_0.root");



	TH1F *m12_data = new TH1F("m12_data","", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
	tr_data->Project("m12_data","s12");
	m12_data->SetMarkerStyle(20);
	m12_data->SetMarkerSize(0.6);
	m12_data->SetLineWidth(2);
	m12_data->SetLineColor(1);


	TH1F *m13_data = new TH1F("m13_data","", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	tr_data->Project("m13_data","s13");
	m13_data->SetMarkerStyle(20);
	m13_data->SetMarkerSize(0.6);
	m13_data->SetLineWidth(2);
	m13_data->SetLineColor(1);

	TH1F *m23_data = new TH1F("m23_data","",m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	tr_data->Project("m23_data","s23");
	m23_data->SetMarkerStyle(20);
	m23_data->SetMarkerSize(0.6);
	m23_data->SetLineWidth(2);
	m23_data->SetLineColor(1);


	foo->cd(1);
	dalitzplot->Draw("colz");

foo->cd(2);
    m12_data->Draw("e");
foo->SaveAs("m12.png");
/*
	foo->cd(2);
//	m12_data->Draw("e");
//	m12_pdf_hist.Scale(m12_data->Integral()/m12_pdf_hist.Integral());
//	m12_pdf_hist.Draw("Hsame");
	m12_pdf_hist.Scale(m12_data->Integral()/m12_pdf_hist.Integral());
	m12_pdf_hist.Draw("H");
	m12_data->Draw("Esame");

	foo->cd(3);
//	m13_pdf_hist.Draw("Hsame");
	m13_pdf_hist.Scale(m13_data->Integral()/m13_pdf_hist.Integral());
	m13_pdf_hist.Draw("H");
	m13_data->Draw("Esame");

	foo->cd(4);
//	m23_data->Draw("e");
//	m23_pdf_hist.Scale(m23_data->Integral()/m23_pdf_hist.Integral());
//	m23_pdf_hist.Draw("Hsame");
	m23_pdf_hist.Scale(m23_data->Integral()/m23_pdf_hist.Integral());
	m23_pdf_hist.Draw("H");
	m23_data->Draw("Esame");

	foo->SaveAs("plots/dalitz_with_projections.C");
*/
//	makeToyDalitzPdfPlots(&prodpdf,data);
	return datapdf;
}

int main(int argc, char **argv) {
	GooFit::Application app("Dalitz example", argc, argv);

        name("PiPiPi_PWA_COEFFS_30_SAda_ampphs.txt");

	std::string filename = "";
	app.add_option("-f,--filename,filename", filename, "File to read in", true)->check(GooFit::ExistingFile);

//    bool make_toy;
//    app.add_flag("-m,--make-toy", make_toy, "Make a toy instead of reading a file in");

	GOOFIT_PARSE(app);

	GooFit::setROOTStyle();

    // Observables setup
//	fptype m12_lower = (KpMass+pi0Mass)*(KpMass+pi0Mass);
//	fptype m12_upper = (_mDp-KsMass)*(_mDp-KsMass);
//
//	fptype m13_lower = (KsMass+pi0Mass)*(KsMass+pi0Mass);
//	fptype m13_upper = (_mDp-KpMass)*(_mDp-KpMass);

	Observable m12("m12", 0, 3.5);
	Observable m13("m13", 0, 3.5);
	EventNumber eventNumber("eventNumber");
	m12.setNumBins(350);
	m13.setNumBins(350);

	// Prepare the data
	UnbinnedDataSet data({m12, m13, eventNumber});

	//Set up efficiency pdf

	// Set up the model
	Amp3Body *signal = makeSignalPdf(m12, m13, eventNumber);

	// A wrapper for plotting without complex number segfault
	ProdPdf prodpdf{"prodpdf", {signal}};

	// Add nice tool for making data or plotting
//	DalitzPlotter dplotter{&prodpdf, signal};

	// Read in data
	getToyData(filename, app, data);

	try {
		return runToyFit(signal, &data);
	} catch(const std::runtime_error &e) {
		std::cerr << e.what() << std::endl;
		return 7;
	}
}

