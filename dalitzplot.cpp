#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TTree.h>

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


void getToyData(std::string toyFileName, GooFit::Application &app, UnbinnedDataSet &data) {
    toyFileName = app.get_filename(toyFileName, "examples/dalitzplot");

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
          data.addEvent();
          dalitzplot.Fill(m12.getValue(), m13.getValue());
      }
       f->Close();
 
    TCanvas foo;
    dalitzplot.SetStats(false);
    dalitzplot.Draw("colz");
    foo.SaveAs("dalitzplot_1.png");

    TCanvas *foodal = new TCanvas("c1","c1",800,700);

    double totalDat = 0;
    TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(8);
    m12_dat_hist.SetMarkerSize(1.2);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    
    TH1F m13_dat_hist("m13_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_dat_hist.SetStats(false);
    m13_dat_hist.SetMarkerStyle(8);
    m13_dat_hist.SetMarkerSize(1.2);
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    m13_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
   
    TH2F dalitzpm_dat_hist("dalitzpm_dat_hist",
                           "",
                           m12.getNumBins(),
                           m12.getLowerLimit(),
                           m12.getUpperLimit(),
                           m13.getNumBins(),
                           m13.getLowerLimit(),
                           m13.getUpperLimit());
    dalitzpm_dat_hist.SetStats(false);
    dalitzpm_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    dalitzpm_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");



 for(unsigned int evt = 0; evt < eventNumber.getValue(); ++evt) {
        double currm12 = data.getValue(m12, evt);
        m12_dat_hist.Fill(currm12);
        double currm13 = data.getValue(m13, evt);
        m13_dat_hist.Fill(currm13);

        dalitzpm_dat_hist.Fill(currm12, currm13);

       totalDat++;
       }
     foodal->cd();
//    m12_dat_hist.Draw("colz");
    dalitzpm_dat_hist.Draw("colz");

    foodal->SaveAs("m13_fit.png");


}

int main(int argc, char **argv) {
    GooFit::Application app("Dalitz example", argc, argv);
    std::string toyFileName = "FastMC_DsPiPiPi_PWA_11K_30pt_0.root";

    app.add_option("-f,--filename,filename", toyFileName, "File to read in", true)->check(GooFit::ExistingFile);

    bool make_toy;
    app.add_flag("-m,--make-toy", make_toy, "Make a toy instead of reading a file in");


    GOOFIT_PARSE(app);

    GooFit::setROOTStyle();

    Observable m12("m12", 0, 3.5);
    Observable m13("m13", 0, 3.5);
    EventNumber eventNumber("eventNumber");
    m12.setNumBins(240);
    m13.setNumBins(240);

    // Prepare the data
    UnbinnedDataSet data({m12, m13, eventNumber});

    getToyData(toyFileName, app, data);

 try {
      app.run();
  } catch (const GooFit::ParseError& e) {
      return app.exit(e);
  }

}
                                                                                                                     
