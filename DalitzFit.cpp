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
#include <TText.h>
#include <TStyle.h>



// System stuff
#include <fstream>
#include <sys/time.h>
#include <sys/times.h>
#include <cassert>
#include <climits>

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

#include <goofit/PDFs/combine/EventWeightedAddPdf.h>

#include <mcbooster/GContainers.h>
#include <goofit/utilities/Uncertain.h>

#include <goofit/fitting/FitManagerMinuit1.h>
#include <goofit/fitting/FitManagerMinuit2.h>

#include <goofit/FunctorWriter.h>



using namespace std;
using namespace GooFit;


TCanvas foodal;
TCanvas *foo;
fptype x;


Observable m12("m12", 0, 3.5);
Observable m13("m13", 0, 3.5);
EventNumber eventNumber("eventNumber");
UnbinnedDataSet data();

//Observable m12();
//Observable m13;
//EventNumber eventNumber;
//UnbinnedDataSet data;


Amp3Body* signaldalitz;
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

bool cpuDalitz(fptype m12, fptype m13, fptype bigM, fptype dm1, fptype dm2, fptype dm3) {
    if(m12 < pow(dm1 + dm2, 2))
        return false; // This m12 cannot exist, it's less than the square of the (1,2) particle mass.

    if(m12 > pow(bigM - dm3, 2))
        return false; // This doesn't work either, there's no room for an at-rest 3 daughter.

    // Calculate energies of 1 and 3 particles in m12 rest frame.
    fptype e1star = 0.5 * (m12 - dm2 * dm2 + dm1 * dm1) / sqrt(m12);
    fptype e3star = 0.5 * (bigM * bigM - m12 - dm3 * dm3) / sqrt(m12);

    // Bounds for m13 at this value of m12.
    fptype minimum
        = pow(e1star + e3star, 2) - pow(sqrt(e1star * e1star - dm1 * dm1) + sqrt(e3star * e3star - dm3 * dm3), 2);

    if(m13 < minimum)
        return false;

    fptype maximum
        = pow(e1star + e3star, 2) - pow(sqrt(e1star * e1star - dm1 * dm1) - sqrt(e3star * e3star - dm3 * dm3), 2);

    if(m13 > maximum)
        return false;

    return true;
}

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

void getToyData(std::string toyFileName, GooFit::Application &app, UnbinnedDataSet &data) {
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
    data =UnbinnedDataSet(vars);
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
    
 
    GOOFIT_INFO("Read in {} events", data.getNumEvents());

    TCanvas fooda;
    dalitzplot.SetStats(false);
    dalitzplot.Draw("colz");
    fooda.SaveAs("dalitzplot.png");

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
void makeDalitzPlots(Amp3Body *overallSignal, std::string plotdir = "./plots_from_mixfit/") {
//    std::string mkplotdir{"mkdir " + plotdir};
//    system(mkplotdir.c_str());
  //  foo->cd();

    TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(8);
    m12_dat_hist.SetMarkerSize(1.2);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_pdf_hist.SetStats(false);
    m12_pdf_hist.SetLineColor(kBlue);
    m12_pdf_hist.SetLineWidth(3);

    TH1F m13_dat_hist("m13_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_dat_hist.SetStats(false);
    m13_dat_hist.SetMarkerStyle(8);
    m13_dat_hist.SetMarkerSize(1.2);
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    m13_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_pdf_hist.SetStats(false);
    m13_pdf_hist.SetLineColor(kBlue);
    m13_pdf_hist.SetLineWidth(3);

    TH1F m23_dat_hist("m23_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_dat_hist.SetStats(false);
    m23_dat_hist.SetMarkerStyle(8);
    m23_dat_hist.SetMarkerSize(1.2);
    m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
    m23_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_pdf_hist.SetStats(false);
    m23_pdf_hist.SetLineColor(kBlue);
    m23_pdf_hist.SetLineWidth(3);

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
    TH2F dalitzpm_pdf_hist("dalitzpm_pdf_hist",
                           "",
                           m12.getNumBins(),
                           m12.getLowerLimit(),
                           m12.getUpperLimit(),
                           m13.getNumBins(),
                           m13.getLowerLimit(),
                           m13.getUpperLimit());
    dalitzpm_pdf_hist.SetStats(false);
    double totalPdf = 0;
    double totalDat = 0;

    for(unsigned int evt = 0; evt < data.getNumEvents(); ++evt) {
        double currm12 = data.getValue(m12, evt);
        m12_dat_hist.Fill(currm12);

        double currm13 = data.getValue(m13, evt);
        m13_dat_hist.Fill(currm13);

        dalitzpm_dat_hist.Fill(currm12, currm13);

        double currm23 = cpuGetM23(currm12, currm13);
        m23_dat_hist.Fill(currm23);

       totalDat++;
       }

  

    
        std::vector<Observable> vars;
        vars.push_back(m12);
        vars.push_back(m13);
        vars.push_back(eventNumber);

        UnbinnedDataSet currData(vars);

        int evtCounter = 0;

        for(int i = 0; i < m12.getNumBins(); ++i) {
            m12.setValue(m12.getLowerLimit()
                          + (m12.getUpperLimit() - m12.getLowerLimit()) * (i + 0.5) / m12.getNumBins());

            for(int j = 0; j < m13.getNumBins(); ++j) {
                m13.setValue(m13.getLowerLimit()
                              + (m13.getUpperLimit() - m13.getLowerLimit()) * (j + 0.5) / m13.getNumBins());

                if(!cpuDalitz(m12.getValue(), m13.getValue(), _mDs, piPlusMass, piPlusMass, piPlusMass))
                    continue;

                    eventNumber.setValue(evtCounter);
                    evtCounter++;
                    currData.addEvent();
                }
            }
            overallSignal->setData(&currData);
//        signaldalitz->setDataSize(currData.getNumEvents());

           std::vector<std::vector<double>> pdfValues = overallSignal->getCompProbsAtDataPoints();

                   for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {

                // Um... these two are switched? Weirdness...
                double currm12 = currData.getValue(m13, j);
                m12_pdf_hist.Fill(currm12, pdfValues[0][j]);

                double currm13 = currData.getValue(m12, j);
                m13_pdf_hist.Fill(currm13, pdfValues[0][j]);
                dalitzpm_pdf_hist.Fill(currm12, currm13, pdfValues[0][j]);

                double currm23 = cpuGetM23(currm12, currm13);
                m23_pdf_hist.Fill(currm23, pdfValues[0][j]);


                totalPdf += pdfValues[0][j];

            }


*/
            // If PDF doesn't depend on sigma, don't project from that dimension.
 /*       
   for(int i = 1; i <= m12.getNumBins(); ++i) {
        m12_pdf_hist.SetBinContent(i, m12_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
    }

    for(int i = 1; i <= m13.getNumBins(); ++i) {
        m13_pdf_hist.SetBinContent(i, m13_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
        m23_pdf_hist.SetBinContent(i, m23_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
    }

    for(int i = 1; i <= m12.getNumBins(); ++i) {
        for(int j = 1; j <= m13.getNumBins(); ++j) {
            dalitzpm_pdf_hist.SetBinContent(i, j, dalitzpm_pdf_hist.GetBinContent(i, j) * totalDat / totalPdf);
            }
         }
*/
/*
std::vector<std::vector<double>> pdfValue = signaldalitz->getResProbsAtDataPoints();

//  signalDalitz->getResProbsAtDataPoints(pdfValues);
  const int nRes = pdfValue.size();
  TH1* m12_pdf_hist_res[nRes];
  TH1* m13_pdf_hist_res[nRes];
  TH1* m23_pdf_hist_res[nRes];
  for (int i=0;i<nRes;i++){
      sprintf(strbuffer, "%s_res%d", m12_pdf_hist.GetName(), i);
      m12_pdf_hist_res[i] = (TH1*)m12_pdf_hist.Clone(strbuffer);
      m12_pdf_hist_res[i]->Reset();
      sprintf(strbuffer, "%s_res%d", m13_pdf_hist.GetName(), i);
      m13_pdf_hist_res[i] = (TH1*)m13_pdf_hist.Clone(strbuffer);
      m13_pdf_hist_res[i]->Reset();
      sprintf(strbuffer, "%s_res%d", m23_pdf_hist.GetName(), i);
      m23_pdf_hist_res[i] = (TH1*)m23_pdf_hist.Clone(strbuffer);
      m23_pdf_hist_res[i]->Reset();
  }
//  cout<<"******************Total:" << pdfValues[0].size()<<endl;  
  for (unsigned int j = 0; j < pdfValue[0].size(); ++j) {
      double currm12 = currData.getValue(*m12, j);
      double currm13 = currData.getValue(*m13, j);
      double currm23 = cpuGetM23(currm12, currm13);
      for (int i=0;i<nRes;i++){
          m12_pdf_hist_res[i]->Fill(currm12, pdfValues[i][j]);
          m13_pdf_hist_res[i]->Fill(currm13, pdfValues[i][j]);
          m23_pdf_hist_res[i]->Fill(currm23, pdfValues[i][j]);
//          cout<<"For Res #"<<i<<", "<<j<<": "<<currm12<<','<<currm13<<": "<<pdfValues[i][j]<<endl;
      }
  }
const int colors[] = {kRed, kMagenta, kGray+2, kGreen+3, kYellow+3, kOrange};
for (int i=0;i<nRes;i++){
        m12_pdf_hist_res[i]->Scale(1.1);
        m12_pdf_hist_res[i]->Rebin(8.75);
        for (int j=0;j<m12_pdf_hist_res[i]->GetNbinsX();j++)
            m12_pdf_hist_res[i]->SetBinError(j+1,0);
//        hres[i]->SetLineStyle(kDashed);
        m12_pdf_hist_res[i]->SetLineColor(colors[i]);
        m12_pdf_hist_res[i]->Draw("csame");
}
    foo->SaveAs((plotdir+"/m12_fit.png").c_str());
*/

/*
foo->Divide(2,2);
   foo->cd(1);
    m13_dat_hist.Draw("p");
    m13_pdf_hist.Draw("lsame");
    foo->SaveAs("m13_fit.png");
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/m13_fit_log.png").c_str());
    foo->SetLogy(false);


    foo->cd(1);
    m13_dat_hist.Draw("p");
    m13_pdf_hist.Draw("lsame");
    foo->SaveAs((plotdir + "/m13_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/m13_fit_log.png").c_str());
    foo->SetLogy(false);
    foo->cd(2);
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

    foodal.cd();
    dalitzpm_dat_hist.Draw("colz");

}      

  */   
void  makeDalitzPlots(UnbinnedDataSet *data){


TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(8);
    m12_dat_hist.SetMarkerSize(1.2);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_pdf_hist.SetStats(false);
    m12_pdf_hist.SetLineColor(kBlue);
    m12_pdf_hist.SetLineWidth(3);

    TH1F m13_dat_hist("m13_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_dat_hist.SetStats(false);
    m13_dat_hist.SetMarkerStyle(8);
    m13_dat_hist.SetMarkerSize(1.2);
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    m13_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_pdf_hist.SetStats(false);
    m13_pdf_hist.SetLineColor(kBlue);
    m13_pdf_hist.SetLineWidth(3);

    TH1F m23_dat_hist("m23_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_dat_hist.SetStats(false);
    m23_dat_hist.SetMarkerStyle(8);
    m23_dat_hist.SetMarkerSize(1.2);
    m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
    m23_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_pdf_hist.SetStats(false);
    m23_pdf_hist.SetLineColor(kBlue);
    m23_pdf_hist.SetLineWidth(3);

   double totalPdf = 0;
    double totalDat = 0;

    for(unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
        double currm12 = data->getValue(m12, evt);
        m12_dat_hist.Fill(currm12);

        double currm13 = data->getValue(m13, evt);
        m13_dat_hist.Fill(currm13);

       }

      foo->Divide(2,2);
   foo->cd(1);
    m13_dat_hist.Draw("p");
   m13_pdf_hist.Draw("lsame");
    foo->SaveAs("m13_fit.png");


}

int runToyFit(Amp3Body *signal, UnbinnedDataSet *data) {
    // EXERCISE 1 (real part): Create a PolynomialPdf which models
    // the efficiency you imposed in the preliminary, and use it in constructing
    // the signal PDF.

    // EXERCISE 2: Create a K0 veto function and use it as the efficiency.

    // EXERCISE 3: Make the efficiency a product of the two functions
    // from the previous exercises.
    signal->setData(data);
    signal->setDataSize(data->getNumEvents());
    FitManager datapdf(signal);

    datapdf.fit();

    ProdPdf prodpdf{"prodpdf",{signal}};

//    prodpdf.setData(data);
    DalitzPlotter plotter(&prodpdf, signal);
//plotter.fillDataSetMC(data, 1000000);
//    TCanvas foo;
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
    foodal.SaveAs("dalitzpdf.png");
//  signaldalitz=signal;
//   makeDalitzPlots(signal,"plots_from_toy");
//    return datapdf;
//Observable m12("m12", 0, 3.5);
//Observable m13("m13", 0, 3.5);
//EventNumber eventNumber("eventNumber");
//UnbinnedDataSet ToyMC({m12, m13, eventNumber});
//plotter.fillDataSetMC(ToyMC, 1000000);
/*

    TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(8);
    m12_dat_hist.SetMarkerSize(1.2);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_pdf_hist.SetStats(false);
    m12_pdf_hist.SetLineColor(kBlue);
    m12_pdf_hist.SetLineWidth(3);

    TH1F m13_dat_hist("m13_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_dat_hist.SetStats(false);
    m13_dat_hist.SetMarkerStyle(8);
    m13_dat_hist.SetMarkerSize(1.2);
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    m13_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_pdf_hist.SetStats(false);
    m13_pdf_hist.SetLineColor(kBlue);
    m13_pdf_hist.SetLineWidth(3);

    TH1F m23_dat_hist("m23_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_dat_hist.SetStats(false);
    m23_dat_hist.SetMarkerStyle(8);
    m23_dat_hist.SetMarkerSize(1.2);
    m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
    m23_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_pdf_hist.SetStats(false);
    m23_pdf_hist.SetLineColor(kBlue);
    m23_pdf_hist.SetLineWidth(3);


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
    TH2F dalitzpm_pdf_hist("dalitzpm_pdf_hist",
                           "",
                           m12.getNumBins(),
                           m12.getLowerLimit(),
                           m12.getUpperLimit(),
                           m13.getNumBins(),
                           m13.getLowerLimit(),
                           m13.getUpperLimit());
    dalitzpm_pdf_hist.SetStats(false);
    double totalPdf = 0;
    double totalDat = 0;
  

    for(unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
        double currm12 = data->getValue(m12, evt);
        m12_dat_hist.Fill(currm12);

        double currm13 = data->getValue(m13, evt);
        m13_dat_hist.Fill(currm13);

        dalitzpm_dat_hist.Fill(currm12, currm13);

        double currm23 = cpuGetM23(currm12, currm13);
        m23_dat_hist.Fill(currm23);

       totalDat++;
       }

*/
/*
   std::vector<Observable> vars;
        vars.push_back(m12);
        vars.push_back(m13);
        vars.push_back(eventNumber);

        UnbinnedDataSet currData(vars);

        int evtCounter = 0;

        for(int i = 0; i < m12.getNumBins(); ++i) {
            m12.setValue(m12.getLowerLimit()
                          + (m12.getUpperLimit() - m12.getLowerLimit()) * (i + 0.5) / m12.getNumBins());

            for(int j = 0; j < m13.getNumBins(); ++j) {
                m13.setValue(m13.getLowerLimit()
                              + (m13.getUpperLimit() - m13.getLowerLimit()) * (j + 0.5) / m13.getNumBins());

                if(!cpuDalitz(m12.getValue(), m13.getValue(), _mDs, piPlusMass, piPlusMass, piPlusMass))
                    continue;

                    eventNumber.setValue(evtCounter);
                    evtCounter++;
                    currData.addEvent();
                }
            }
*/
//   foo->cd();
//    m13_dat_hist.Draw("p");
//   m13_pdf_hist.Draw("lsame");
//    foo->SaveAs("m13_fit.png");
//    foo->SetLogy(true);
//    foo->SaveAs((plotdir + "/m13_fit_log.png").c_str());
//    foo->SetLogy(false);
// makeDalitzPlots(data);
    return datapdf;

}

int main(int argc, char **argv) {
    GooFit::Application app("Dalitz example", argc, argv);

    name("PiPiPi_PWA_COEFFS_30_SAda_ampphs.txt");

    std::string toyFileName = "";
    app.add_option("-f,--filename,filename",toyFileName, "File to read in", true)->check(GooFit::ExistingFile);

    bool make_toy;
    app.add_flag("-m,--make-toy", make_toy, "Make a toy instead of reading a file in");

    GOOFIT_PARSE(app);

    GooFit::setROOTStyle();

    // Observables setup
/*    Observable m12("m12", 0, 3.5);
    Observable m13("m13", 0, 3.5);
    EventNumber eventNumber("eventNumber");
    m12.setNumBins(350);
    m13.setNumBins(350);
*/

m12=Observable("m12", 0, 3.5);
m13=Observable("m13", 0, 3.5);
eventNumber= EventNumber("eventNumber");
m12.setNumBins(350);
m13.setNumBins(350);

    // Prepare the data
 //   UnbinnedDataSet data({m12, m13, eventNumber});
//data=UnbinnedDataSet({m12, m13, eventNumber});

    // Set up the model
    Amp3Body *signal = makeSignalPdf(m12, m13, eventNumber);

    // A wrapper for plotting without complex number segfault
//    ProdPdf prodpdf{"prodpdf", {signal}};

    // Add nice tool for making data or plotting
//    DalitzPlotter dplotter{&prodpdf, signal};

    // Read in data
 //   dplotter.fillDataSetMC(data, 1000000);

        getToyData(toyFileName, app,data);
//   makeDalitzPlots(&data);
 
//    try {
        return runToyFit(signal, &data);
//    } catch(const std::runtime_error &e) {
//        std::cerr << e.what() << std::endl;
//        return 7;
//    }
}
