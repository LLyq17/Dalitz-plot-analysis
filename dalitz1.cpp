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
#include <TMath.h>
#include <TChain.h>
#include <TBranch.h>


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
#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/utilities/Style.h>

using namespace std;
using namespace GooFit;
//TCanvas *foodal;
//TCanvas foo;
//double x;

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


void normalize(TH1 *dat) {
    double integral = 0;

    for(int i = 1; i <= dat->GetNbinsX(); ++i) {
        integral += dat->GetBinContent(i);
    }

    integral = 1.0 / integral;

    for(int i = 1; i <= dat->GetNbinsX(); ++i) {
        dat->SetBinContent(i, integral * dat->GetBinContent(i));
        dat->SetBinError(i, integral * dat->GetBinError(i));
    }
}


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

// TH1 * hx= dalitzplot.ProjectionX();
// TH1 * hy= dalitzplot.ProjectionY();
//    hx->Rebin(8.75);
//    hy->Rebin(8.75);
//    x=hx->Integral();
//  hx->SetOption("lego");
//    hy->Draw("colz");
//  hx->SetLineColor(kRed);
//  hx->Draw("e1");
//  foodal.SaveAs("dalitzplot.png");

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
/*
void makeDalitzPlots(GooPdf *overallSignal, std::string plotdir = "./plots_from_mixfit/") {
    std::string mkplotdir{"mkdir " + plotdir};
    system(mkplotdir.c_str());
    foo->cd();

    TH1F m12_dat_hist("m12_dat_hist", "", m12->getNumBins(), m12->getLowerLimit(), m12->getUpperLimit());
    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(8);
    m12_dat_hist.SetMarkerSize(1.2);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m12_pdf_hist("m12_pdf_hist", "", m12->getNumBins(), m12->getLowerLimit(), m12->getUpperLimit());
    m12_pdf_hist.SetStats(false);
    m12_pdf_hist.SetLineColor(kBlue);
    m12_pdf_hist.SetLineWidth(3);

    TH1F m13_dat_hist("m13_dat_hist", "", m13->getNumBins(), m13->getLowerLimit(), m13->getUpperLimit());
    m13_dat_hist.SetStats(false);
    m13_dat_hist.SetMarkerStyle(8);
    m13_dat_hist.SetMarkerSize(1.2);
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    m13_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m13_pdf_hist("m13_pdf_hist", "", m13->getNumBins(), m13->getLowerLimit(), m13->getUpperLimit());
    m13_pdf_hist.SetStats(false);
    m13_pdf_hist.SetLineColor(kBlue);
    m13_pdf_hist.SetLineWidth(3);

    TH1F m23_dat_hist("m23_dat_hist", "", m13->getNumBins(), m13->getLowerLimit(), m13->getUpperLimit());
    m23_dat_hist.SetStats(false);
    m23_dat_hist.SetMarkerStyle(8);
    m23_dat_hist.SetMarkerSize(1.2);
    m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
    m23_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m23_pdf_hist("m23_pdf_hist", "", m13->getNumBins(), m13->getLowerLimit(), m13->getUpperLimit());
    m23_pdf_hist.SetStats(false);
    m23_pdf_hist.SetLineColor(kBlue);
    m23_pdf_hist.SetLineWidth(3);

    TH2F dalitzpm_dat_hist("dalitzpm_dat_hist",
                           "",
                           m12->getNumBins(),
                           m12->getLowerLimit(),
                           m12->getUpperLimit(),
                           m13->getNumBins(),
                           m13->getLowerLimit(),
                           m13->getUpperLimit());
    dalitzpm_dat_hist.SetStats(false);
    dalitzpm_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    dalitzpm_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    TH2F dalitzpm_pdf_hist("dalitzpm_pdf_hist",
                           "",
                           m12->getNumBins(),
                           m12->getLowerLimit(),
                           m12->getUpperLimit(),
                           m13->getNumBins(),
                           m13->getLowerLimit(),
                           m13->getUpperLimit());
    dalitzpm_pdf_hist.SetStats(false);
    double totalPdf = 0;
    double totalDat = 0;

    for(unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
        double currm12 = data->getValue(*m12, evt);
        m12_dat_hist.Fill(currm12);

        double currm13 = data->getValue(*m13, evt);
        m13_dat_hist.Fill(currm13);

        dalitzpm_dat_hist.Fill(currm12, currm13);

        double currm23 = cpuGetM23(currm12, currm13);
        m23_dat_hist.Fill(currm23);

       totalDat++
       }

        wBkg1->setValue(0);
    const int division = 2;

    for(int half = 0; half < division; ++half) {
        std::vector<Observable> vars;
        vars.push_back(*m12);
        vars.push_back(*m13);
        vars.push_back(*eventNumber);

        UnbinnedDataSet currData(vars);

        int evtCounter = 0;

        for(int i = 0; i < m12->getNumBins(); ++i) {
            m12->setValue(m12->getLowerLimit()
                          + (m12->getUpperLimit() - m12->getLowerLimit()) * (i + 0.5) / m12->getNumBins());

            for(int j = 0; j < m13->getNumBins(); ++j) {
                m13->setValue(m13->getLowerLimit()
                              + (m13->getUpperLimit() - m13->getLowerLimit()) * (j + 0.5) / m13->getNumBins());

                if(!cpuDalitz(m12->getValue(), m13->getValue(), _mD0, piZeroMass, piPlusMass, piPlusMass))
                    continue;

                    eventNumber->setValue(evtCounter);
                    evtCounter++;
                    currData.addEvent();
                }
            }
           std::vector<std::vector<double>> pdfValues = overallSignal->getCompProbsAtDataPoints();

                   for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {

                // Um... these two are switched? Weirdness...
                double currm12 = currData.getValue(*m13, j);
                m12_pdf_hist.Fill(currm12, pdfValues[0][j]);

                double currm13 = currData.getValue(*m12, j);
                m13_pdf_hist.Fill(currm13, pdfValues[0][j]);
                dalitzpm_pdf_hist.Fill(currm12, currm13, pdfValues[0][j]);

                double currm23 = cpuGetM23(currm12, currm13);
                m23_pdf_hist.Fill(currm23, pdfValues[0][j]);


                totalPdf += pdfValues[0][j];

            }

            // If PDF doesn't depend on sigma, don't project from that dimension.
        }
   for(int i = 1; i <= m12->getNumBins(); ++i) {
        m12_pdf_hist.SetBinContent(i, m12_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
    }

    for(int i = 1; i <= m13->getNumBins(); ++i) {
        m13_pdf_hist.SetBinContent(i, m13_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
        m23_pdf_hist.SetBinContent(i, m23_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
    }

    for(int i = 1; i <= m12->getNumBins(); ++i) {
        for(int j = 1; j <= m13->getNumBins(); ++j) {
            dalitzpm_pdf_hist.SetBinContent(i, j, dalitzpm_pdf_hist.GetBinContent(i, j) * totalDat / totalPdf);
            }
         }
    m13_dat_hist.Draw("p");
    m13_pdf_hist.Draw("lsame");
    foo->SaveAs((plotdir + "/m13_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/m13_fit_log.png").c_str());
    foo->SetLogy(false);

    m12_dat_hist.Draw("p");
    m12_pdf_hist.Draw("lsame");
    foo->SaveAs((plotdir + "/m12_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/m12_fit_log.png").c_str());
    foo->SetLogy(false);

    m23_dat_hist.Draw("p");
    m23_pdf_hist.Draw("lsame");
    foo->SaveAs((plotdir + "/m23_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/m23_fit_log.png").c_str());
    foo->SetLogy(false);

    foodal->cd();
    dalitzpm_dat_hist.Draw("colz");
}      

*/     


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

    ProdPdf prodpdf{"prodpdf",{signal}};


    DalitzPlotter plotter(&prodpdf, signal);

    TCanvas foo;
    foo.cd();
//foo.Clear();
    TH2F *dalitzplot = plotter.make2D();
    dalitzplot->Draw("colz");

//    TH1* hx = dalitzplot->ProjectionX();
//     hx->Draw("same");
//TH1* hy = dalitzplot->ProjectionY();
//hy->Rebin(8.75);
//const double sca = x/hy->Integral();
//cout<<sca<<endl;
// hy->Scale(sca);//hy->Rebin(8.75);

//normalize(hy);
//hy->SetLineColor(kGreen+3);
//  hy->Scale(sca);

// hy->Draw("Csame");

//foodal->cd();


    foo.SaveAs("dalitzpdf.png");
TCanvas *foodal = new TCanvas("c1","c1",800,700);



TChain *tr_data = new TChain("DecayTree");
tr_data->Add("FastMC_DsPiPiPi_PWA_11K_30pt_0.root");


TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
tr_data->Project("m12_dat_hist","s12");

//    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(20);
    m12_dat_hist.SetMarkerSize(0.6);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    
//foo.SaveAs("dalitzpdf.png");
foodal->cd();
m12_dat_hist.Draw("e");
foodal->SaveAs("m12_dat_hist.png");

 //for comparison
/*
TCanvas *foo = new TCanvas("c1","c1",800,700);
        foo->Divide(2,2);


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

*/
    return datapdf;
}

int main(int argc, char **argv) {
    GooFit::Application app("Dalitz example", argc, argv);

    name("PiPiPi_PWA_COEFFS_30_SAda_ampphs.txt");

    std::string toyFileName = "";
    app.add_option("-f,--filename,filename",toyFileName, "File to read in", true)->check(GooFit::ExistingFile);

//    bool make_toy;
//    app.add_flag("-m,--make-toy", make_toy, "Make a toy instead of reading a file in");

    GOOFIT_PARSE(app);

    GooFit::setROOTStyle();

    // Observables setup
    Observable m12("m12", 0, 3.5);
    Observable m13("m13", 0, 3.5);
    EventNumber eventNumber("eventNumber");
    m12.setNumBins(350);
    m13.setNumBins(350);

    // Prepare the data
    UnbinnedDataSet data({m12, m13, eventNumber});

    // Set up the model
    Amp3Body *signal = makeSignalPdf(m12, m13, eventNumber);

    // A wrapper for plotting without complex number segfault
    ProdPdf prodpdf{"prodpdf", {signal}};

    // Add nice tool for making data or plotting
    DalitzPlotter dplotter{&prodpdf, signal};

    // Read in data
        getToyData(toyFileName, app, data);
 
    try {
        return runToyFit(signal, &data);
    } catch(const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        return 7;
    }
}
