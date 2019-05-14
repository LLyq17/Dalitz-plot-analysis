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
#include <TH1.h>



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
#include <goofit/PDFs/physics/Amp3Body_IS.h>
#include <goofit/PDFs/physics/Amp3Body_TD.h>

#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/utilities/Style.h>

using namespace std;
using namespace GooFit;
TCanvas fooo;
TCanvas foodal;
TCanvas *foo = new TCanvas("c1","c1",800,700);
char strbuffer[1000];

Observable m12("m12", 0, 3.5);
Observable m13("m13", 0, 3.5);
EventNumber eventNumber("eventNumber");
UnbinnedDataSet *data     = nullptr;

Amp3Body *incsum1      = nullptr;
Amp3Body *incsum2      = nullptr;
Amp3Body *incsum3      = nullptr;
Amp3Body *incsum4      = nullptr;
std::vector<PdfBase *> comps;
//Amp3Body* signaldalitz= nullptr;
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


void getToyData(std::string toyFileName, GooFit::Application &app) {
    toyFileName = app.get_filename(toyFileName, "examples/dalitz");

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
//    data =new UnbinnedDataSet(vars);
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
         eventNumber.setValue(data->getNumEvents());
          data->addEvent();
          dalitzplot.Fill(m12.getValue(), m13.getValue());
      }
    
 
    GOOFIT_INFO("Read in {} events", data->getNumEvents());

    TCanvas fooda;
    dalitzplot.SetStats(false);
    dalitzplot.Draw("colz");
    fooda.SaveAs("dalitzplot.png");


}



std::vector<fptype> HH_bin_limits;
std::vector<Variable> pwa_coefs_reals;
std::vector<Variable> pwa_coefs_imags;



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

fpcomplex coefs_1;
ResonancePdf *rho_1=nullptr;
ResonancePdf *rho0_1450_1=nullptr;
ResonancePdf *f2_1270_1=nullptr;
ResonancePdf *swave1=nullptr;

//incsum1 =new Amp3Body_IS("incsum1", *m12, *m13, *eventNumber, special_rho_decay, bkg2_rho_mods);

Amp3Body *ResonanceSinglePdf(GooPdf *eff = 0) {

 DecayInfo3 dtop0pp1;
    dtop0pp1.motherMass   = _mDs;
    dtop0pp1.daug1Mass    = piPlusMass;
    dtop0pp1.daug2Mass    = piPlusMass;
    dtop0pp1.daug3Mass    = piPlusMass;
    dtop0pp1.meson_radius = 1.5;

    bool fixAmps = false; // Takes ~400x longer
  bool fixAmp =true;
  Variable rho0_mass("rho0_mass", 0.77526, 0.01, 0.4, 1.1);
  Variable rho0_width("rho0_width", 0.1491, 0.01, 0.1, 0.3);


   ResonancePdf* rho_11 = new Resonances::RBW(
        "rho0",
        Variable("rho0_amp_real1", 0.19, 0.001, 0, 0),
        Variable("rho0_amp_imag1", -1.1, 0.1, 0, 0),
        rho0_mass,
        rho0_width,
        1,
        PAIR_12 ,true);
   

    dtop0pp1.resonances.push_back(rho_11); 
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



    return new Amp3Body("signalPDF", m12, m13, eventNumber, dtop0pp1,eff);


}

Amp3Body *ResonanceSignalPdf(int bkg,GooPdf *eff = 0) {

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


   ResonancePdf* rho_11 = new Resonances::RBW(
        "rho0",
        rho_1->get_amp_real(),
        rho_1->get_amp_img(),
        rho0_mass,
        rho0_width,
        1,
        PAIR_12 ,true);



   Variable sharedMass("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
   Variable shareWidth("rhop_1450_width", 0.400, 0.01, 0.01, 5.0);


   ResonancePdf *rho0_1450_11 = new Resonances::RBW(
        "rho0_1450",
        fixAmps ? Variable("rho0_1450_amp_real1", 1.2) : Variable("rho0_1450_amp_real1",1, 0.01, 0, 0),
        fixAmps ? Variable("rho0_1450_amp_imag1", -4.1) : Variable("rho0_1450_amp_imag1", 0, 0.01, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_12,true);


   ResonancePdf *f2_1270_11 = new Resonances::RBW(
        "f2_1270",
         fixAmp ? Variable("f2_1270_amp_real1", 1.0)  : Variable("f2_1270_amp_real1", 1.0, 0.01, 0, 0),
         fixAmp ? Variable("f2_1270_amp_imag1", 0)    : Variable("f2_1270_amp_imag1", 0, 0.01, 0, 0),
         fixedf2Mass,
         fixedf2Width,
         2,
         PAIR_12,true);



    const fptype scale = 1;
   ResonancePdf *swave11 = new Resonances::Spline(
          "swave",
          fixAmp ? Variable("swave_amp_real1", 1)   : Variable("swave_amp_real1", -1.86396*scale, 0.001, 0, 0),
          fixAmp ? Variable("swave_amp_imag1",0 )    : Variable("swave_amp_imag1", .775892256*scale, 0.1, 0, 0),
          HH_bin_limits,
          pwa_coefs_reals,
          pwa_coefs_imags,
          PAIR_12,
          true);


   switch(bkg) {
    case 4:
        dtop0pp.resonances.push_back(rho_11);

        break;

    case 3:
        dtop0pp.resonances.push_back(rho0_1450_11);
        break;
    
    case 2:
        dtop0pp.resonances.push_back(swave11);
        break;

    case 1:
    default:
         dtop0pp.resonances.push_back(f2_1270_11);
        break;
    }

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

Amp3Body *makeSignalPdf(GooPdf *eff = 0) {
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
  

    rho_1 = new Resonances::RBW(
        "rho0",
        fixAmps ? Variable("rho0_amp_real1", 0.19) : Variable("rho0_amp_real1", -0.063, 0.001, 0, 0),
        fixAmps ? Variable("rho0_amp_imag1", -1.1) : Variable("rho0_amp_imag1", -1.42, 0.1, 0, 0),
        rho0_mass,
        rho0_width,
        1,
        PAIR_12 ,true);
   

 
   Variable sharedMass("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
   Variable shareWidth("rhop_1450_width", 0.400, 0.01, 0.01, 5.0);


   rho0_1450_1 = new Resonances::RBW(
        "rho0_1450",
        fixAmps ? Variable("rho0_1450_amp_real1", 1.2) : Variable("rho0_1450_amp_real1",0.78, 0.01, 0, 0),
        fixAmps ? Variable("rho0_1450_amp_imag1", -4.1) : Variable("rho0_1450_amp_imag1", 0.66, 0.01, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_12,true);


   f2_1270_1 = new Resonances::RBW(
        "f2_1270",
         fixAmp ? Variable("f2_1270_amp_real1", 1.0)  : Variable("f2_1270_amp_real1", 1.0, 0.01, 0, 0),
         fixAmp ? Variable("f2_1270_amp_imag1", 0)    : Variable("f2_1270_amp_imag1", 0, 0.01, 0, 0),
         fixedf2Mass,
         fixedf2Width,
         2,
         PAIR_12,true);



    const fptype scale = 1;
   swave1 = new Resonances::Spline(
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

//cout<<rho0_amp_imag1.value<<endl;
//    cout<<rho0_amp_real1.value<<endl; 

    return new Amp3Body("signalPDF", m12, m13, eventNumber, dtop0pp,eff);
}

void drawFitPlotsProjection(TH1F m13_dat_hist, TH1F m13_pdf_hist, TH1** m13_pdf_hist_res , int nRes, char pairname[10],std::string plotdir = "./plots_from_mixfit/"){


    foo->cd(); 
    foo->Clear();
    const int colors[] = {kRed, kMagenta, kGray+2, kGreen+3, kYellow+3, kOrange};

    m13_dat_hist.Rebin(8.75);
    m13_pdf_hist.Rebin(8.75);
//    m13_dat_hist.Draw("ep");
    m13_pdf_hist.Scale(m13_dat_hist.Integral()/m13_pdf_hist.Integral());
    for (int i=0;i<m13_pdf_hist.GetNbinsX();i++)
        m13_pdf_hist.SetBinError(i+1,0);

    m13_pdf_hist.Draw("c");
    m13_dat_hist.Draw("e1same");

    for (int i=0;i<nRes;i++){
       m13_pdf_hist_res[i]->Scale(m13_dat_hist.Integral()/m13_pdf_hist.Integral());
       m13_pdf_hist_res[i]->Rebin(8.75);

       for (int j=0;j<m13_pdf_hist_res[i]->GetNbinsX();j++)
            m13_pdf_hist_res[i]->SetBinError(j+1,0);
//        hres[i]->SetLineStyle(kDashed);
            m13_pdf_hist_res[i]->SetLineColor(colors[i]);
            m13_pdf_hist_res[i]->Draw("csame");
    }

    sprintf(strbuffer, "%s/%s_pwa_fit.png", plotdir.c_str(), pairname);
    foo->SaveAs(strbuffer);
    sprintf(strbuffer, "%s/%s_pwa_fit.pdf", plotdir.c_str(), pairname);
    foo->SaveAs(strbuffer);

//    foo->SaveAs((plotdir + "/m13_fit.png").c_str());

//    foo->SaveAs((plotdir + "/m13_fit_log.png").c_str());

/*
    foo->cd();
    foo->Clear();
    m12_dat_hist.Rebin(8.75);
    m12_pdf_hist.Rebin(8.75);
    m12_pdf_hist.Scale(m12_dat_hist.Integral()/m12_pdf_hist.Integral());
    for (int i=0;i<m12_pdf_hist.GetNbinsX();i++)
        m12_pdf_hist.SetBinError(i+1,0);

    m12_pdf_hist.Draw("c");
    m12_dat_hist.Draw("e1same");

    for (int i=0;i<nRes;i++){
       m12_pdf_hist_res[i]->Scale(1.1);
       m12_pdf_hist_res[i]->Rebin(8.75);

       for (int j=0;j<m12_pdf_hist_res[i]->GetNbinsX();j++)
            m12_pdf_hist_res[i]->SetBinError(j+1,0);
//        hres[i]->SetLineStyle(kDashed);
            m12_pdf_hist_res[i]->SetLineColor(colors[i]);
            m12_pdf_hist_res[i]->Draw("csame");
    }

    foo->SaveAs((plotdir + "/m12_fit.png").c_str());

    foo->SaveAs((plotdir + "/m12_fit_log.png").c_str());


    foo->cd();
    foo->Clear();
    m23_dat_hist.Rebin(8.75);
    m23_pdf_hist.Rebin(8.75);
    m23_pdf_hist.Scale(m23_dat_hist.Integral()/m23_pdf_hist.Integral());
    for (int i=0;i<m23_pdf_hist.GetNbinsX();i++)
        m23_pdf_hist.SetBinError(i+1,0);

    m23_pdf_hist.Draw("c");
    m23_dat_hist.Draw("e1same");

    for (int i=0;i<nRes;i++){
       m12_pdf_hist_res[i]->Scale(1.1);
       m12_pdf_hist_res[i]->Rebin(8.75);

       for (int j=0;j<m12_pdf_hist_res[i]->GetNbinsX();j++)
            m12_pdf_hist_res[i]->SetBinError(j+1,0);
//        hres[i]->SetLineStyle(kDashed);
            m12_pdf_hist_res[i]->SetLineColor(colors[i]);
            m12_pdf_hist_res[i]->Draw("csame");
    }


    foo->SaveAs((plotdir + "/m23_fit.png").c_str());
//    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/m23_fit_log.png").c_str());
//    foo->SetLogy(false);

*/
}

void makeDalitzPlots(GooPdf *overallSignal,Amp3Body *signaldalitz, std::string plotdir = "./plots_from_mixfit/") {

    TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(20);
    m12_dat_hist.SetMarkerSize(0.6);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_pdf_hist.SetStats(false);
    m12_pdf_hist.SetLineColor(kBlue);
    m12_pdf_hist.SetLineWidth(2);

    TH1F m13_dat_hist("m13_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_dat_hist.SetStats(false);
    m13_dat_hist.SetMarkerStyle(20);
    m13_dat_hist.SetMarkerSize(0.6);
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    m13_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_pdf_hist.SetStats(false);
    m13_pdf_hist.SetLineColor(kBlue);
    m13_pdf_hist.SetLineWidth(2);

    TH1F m23_dat_hist("m23_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_dat_hist.SetStats(false);
    m23_dat_hist.SetMarkerStyle(20);
    m23_dat_hist.SetMarkerSize(0.6);
    m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
    m23_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_pdf_hist.SetStats(false);
    m23_pdf_hist.SetLineColor(kBlue);
    m23_pdf_hist.SetLineWidth(2);

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
     signaldalitz->setDataSize(currData.getNumEvents());

     std::vector<std::vector<double>> pdfValues = overallSignal->getCompProbsAtDataPoints();

     for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {

     // Um... these two are switched? Weirdness...
           double currm12 = currData.getValue(m13, j);
           m12_pdf_hist.Fill(currm12, pdfValues[1][j]);

           double currm13 = currData.getValue(m12, j);
           m13_pdf_hist.Fill(currm13, pdfValues[0][j]);
           dalitzpm_pdf_hist.Fill(currm12, currm13, pdfValues[0][j]);

           double currm23 = cpuGetM23(currm12, currm13);
           m23_pdf_hist.Fill(currm23, pdfValues[0][j]);


           totalPdf += pdfValues[0][j];

      }



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

//    std::vector<std::vector<double>> pdfValue = signaldalitz->getResProbsAtDataPoints();
     
      
//    const int nRes = pdfValue.size();
//    cout<<"nRes"<<":"<<nRes<<endl;
//    cout<< pdfValue[0].size()<<endl;
      
/*      
 //   TH1* m12_pdf_hist_res[nRes]={nullptr,nullptr,nullptr,nullptr};
//    TH1* m13_pdf_hist_res[nRes]={nullptr,nullptr,nullptr,nullptr};
//    TH1* m23_pdf_hist_res[nRes]={nullptr,nullptr,nullptr,nullptr};
   TH1 *m12_pdf_hist_res[nRes];
   TH1 *m13_pdf_hist_res[nRes];
   TH1 *m23_pdf_hist_res[nRes];
  
 
   for (int i=0;i<nRes;i++){ 
      m12_pdf_hist_res[i]=new TH1F("m12_pdf_hist_res[i]", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
      m13_pdf_hist_res[i]=new TH1F("m13_pdf_hist_res[i]", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
      m23_pdf_hist_res[i]=new TH1F("m23_pdf_hist_res[i]", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
      
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
       double currm12 = currData.getValue(m12, j);
       double currm13 = currData.getValue(m13, j);
       double currm23 = cpuGetM23(currm12, currm13);

      for (int i=0;i<nRes;i++){
          m12_pdf_hist_res[i]->Fill(currm12, pdfValue[i][j]);
          m13_pdf_hist_res[i]->Fill(currm13, pdfValue[i][j]);
          m23_pdf_hist_res[i]->Fill(currm23, pdfValue[i][j]);
          cout<<"For Res #"<<i<<", "<<j<<": "<<currm12<<','<<currm13<<": "<<pdfValues[i][j]<<endl;
        }
    }
  
*/


//    TH1F* m12_pdf_res=new TH1F("m12_pdf_res", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    
 
 
    
    foo->cd();
    m12_dat_hist.Rebin(8.75);
    m12_pdf_hist.Rebin(8.75);
   
//    m13_dat_hist.Draw("ep");
    
    m12_pdf_hist.Scale(m12_dat_hist.Integral()/m12_pdf_hist.Integral());
    for (int i=0;i<m12_pdf_hist.GetNbinsX();i++)
          m12_pdf_hist.SetBinError(i+1,0);

//    for (int j=0;j<m12_pdf_res->GetNbinsX();j++)
//            m12_pdf_res->SetBinError(j+1,0);


    m12_pdf_hist.Draw("c");
    m12_dat_hist.Draw("e1same");
    foo->SaveAs((plotdir + "/m12_fit.png").c_str());

    foo->SaveAs((plotdir + "/m12_fit_log.png").c_str());

    foo->cd();
    foo->Clear();
    m13_dat_hist.Rebin(8.75);
    m13_pdf_hist.Rebin(8.75);

    m13_pdf_hist.Scale(m13_dat_hist.Integral()/m13_pdf_hist.Integral());
    for (int i=0;i<m13_pdf_hist.GetNbinsX();i++)
        m13_pdf_hist.SetBinError(i+1,0);

    m13_pdf_hist.Draw("c");
    m13_dat_hist.Draw("e1same");
    foo->SaveAs((plotdir + "/m13_fit.png").c_str());

    foo->SaveAs((plotdir + "/m13_fit_log.png").c_str());


    foo->cd();
    foo->Clear();
    m23_dat_hist.Rebin(8.75);
    m23_pdf_hist.Rebin(8.75);
    m23_pdf_hist.Scale(m23_dat_hist.Integral()/m23_pdf_hist.Integral());
    for (int i=0;i<m23_pdf_hist.GetNbinsX();i++)
        m23_pdf_hist.SetBinError(i+1,0);

    m23_pdf_hist.Draw("c");
    m23_dat_hist.Draw("e1same");
    foo->SaveAs((plotdir + "/m23_fit.png").c_str());
//    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/m23_fit_log.png").c_str());
//    foo->SetLogy(false);

  
    char pairname[10]="m12";
//    drawFitPlotsProjection(m12_dat_hist, m12_pdf_hist, m12_pdf_hist_res , nRes, pairname,"plots_from_toy");
    foo->cd();
    foo->Clear();
    dalitzpm_dat_hist.Draw("colz");
    foo->SaveAs((plotdir + "/dalitzpm_dat_hist.png").c_str());

}      

     
void  getDalitzdataPlots(){
            
         
    TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(20);
    m12_dat_hist.SetMarkerSize(0.6);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_pdf_hist.SetStats(false);
    m12_pdf_hist.SetLineColor(kBlue);
    m12_pdf_hist.SetLineWidth(2);
      
    TH1F m13_dat_hist("m13_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_dat_hist.SetStats(false);
    m13_dat_hist.SetMarkerStyle(20);
    m13_dat_hist.SetMarkerSize(0.6);
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    m13_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_pdf_hist.SetStats(false);
    m13_pdf_hist.SetLineColor(kBlue);
    m13_pdf_hist.SetLineWidth(2);
       
    TH1F m23_dat_hist("m23_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_dat_hist.SetStats(false);
    m23_dat_hist.SetMarkerStyle(20);
    m23_dat_hist.SetMarkerSize(0.6);
    m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
    m23_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_pdf_hist.SetStats(false);
    m23_pdf_hist.SetLineColor(kBlue);
    m23_pdf_hist.SetLineWidth(2);

    double totalPdf = 0;
    double totalDat = 0;
    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);
//    data =new UnbinnedDataSet(vars);

    for(unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
        double currm12 = data->getValue(m12, evt);
        m12_dat_hist.Fill(currm12);

        double currm13 = data->getValue(m13, evt);
        m13_dat_hist.Fill(currm13);

       }

   foo->Divide(2,2);
   foo->cd(1);
   foo->Clear();
   m13_dat_hist.Draw("p");
   foo->SaveAs("m13_data.png");


}
void draw_1(){


 incsum1=ResonanceSignalPdf(1);

incsum1->setData(data);
incsum1->setDataSize(data->getNumEvents());
 ProdPdf prodpdf1{"prodpdf1", {incsum1}};

// DalitzPlotter plotter_1(&prodpdf1, incsum1);

TH1F m12_dat_hist_1("m12_dat_hist_1", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_dat_hist_1.SetStats(false);
    m12_dat_hist_1.SetMarkerStyle(20);
    m12_dat_hist_1.SetMarkerSize(0.6);
    m12_dat_hist_1.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist_1.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m12_pdf_hist_1("m12_pdf_hist_1", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_pdf_hist_1.SetStats(false);
    m12_pdf_hist_1.SetLineColor(kBlue);
    m12_pdf_hist_1.SetLineWidth(2);

//   UnbinnedDataSet toyMC({m12, m13, eventNumber});
//   plotter_1.fillDataSetMC(toyMC, 5000000);

//    fptype m12_tmp,m13_tmp,m23_tmp;
//    cout << "toyMC size = " << eventNumber.getValue() << endl;
//        for (int i = 0; i < eventNumber.getValue(); ++i){
//                m12_tmp = toyMC.getValue(m12, i);
//               m12_pdf_hist_1.Fill(m12_tmp);

           //     m13_tmp = toyMC.getValue(m13, i);
           //     m13_pdf_hist.Fill(m13_tmp);

           //     m23_tmp = cpuGetM23(m12_tmp, m13_tmp);
           //     m23_pdf_hist.Fill(m23_tmp);
//        }
    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);

    UnbinnedDataSet currData(vars);

    int evtCounter = 0;
    double totalPdf = 0;
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

     prodpdf1.setData(&currData);
//     signaldalitz->setDataSize(currData.getNumEvents());

     std::vector<std::vector<double>> pdfValues = prodpdf1.getCompProbsAtDataPoints();

     for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {

     // Um... these two are switched? Weirdness...
           double currm12 = currData.getValue(m13, j);
           m12_pdf_hist_1.Fill(currm12, pdfValues[1][j]);
}

   fooo.cd();
m12_pdf_hist_1.Draw("csame");
  
   fooo.SaveAs("pdf.png");



}
void  drawFitPlots(){
//    incsum1=ResonanceSignalPdf(1);
//    incsum2=ResonanceSignalPdf(2);
//    incsum3=ResonanceSignalPdf(3);
//    incsum4=ResonanceSignalPdf(4);

//    ProdPdf prodpdf1{"prodpdf1", {incsum1}};

//   DalitzPlotter plotter_1(&prodpdf1, incsum1);
//    signal_1->setData(data);
//    signal_1->setDataSize(data->getNumEvents());

    Amp3Body *signal_1 = ResonanceSinglePdf();
    signal_1->setData(data);
    signal_1->setDataSize(data->getNumEvents());
    ProdPdf prodpdf{"prodpdf", {signal_1}};
    TCanvas doo;
    doo.cd();
    DalitzPlotter plotter1(&prodpdf, signal_1);
    TH2F* dalitzplot1 =new TH2F("dalitzplot1",
                    "Original Data",
                    m12.getNumBins(),
                    m12.getLowerLimit(),
                    m12.getUpperLimit(),
                    m13.getNumBins(),
                    m13.getLowerLimit(),
                    m13.getUpperLimit());

    dalitzplot1 = plotter1.make2D();
    dalitzplot1->Draw("colz");

    doo.SaveAs("pdf_1.png");
    TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(20);
    m12_dat_hist.SetMarkerSize(0.6);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_pdf_hist.SetStats(false);
    m12_pdf_hist.SetLineColor(kBlue);
    m12_pdf_hist.SetLineWidth(2);

    TH1F m13_dat_hist("m13_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_dat_hist.SetStats(false);
    m13_dat_hist.SetMarkerStyle(20);
    m13_dat_hist.SetMarkerSize(0.6);
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    m13_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_pdf_hist.SetStats(false);
    m13_pdf_hist.SetLineColor(kBlue);
    m13_pdf_hist.SetLineWidth(2);

    TH1F m23_dat_hist("m23_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_dat_hist.SetStats(false);
    m23_dat_hist.SetMarkerStyle(20);
    m23_dat_hist.SetMarkerSize(0.6);
    m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
    m23_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_pdf_hist.SetStats(false);
    m23_pdf_hist.SetLineColor(kBlue);
    m23_pdf_hist.SetLineWidth(2);

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

    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);

    UnbinnedDataSet currData(vars);

    int evtCounter = 0;
    double totalPdf = 0;
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

     prodpdf.setData(&currData);
//     signaldalitz->setDataSize(currData.getNumEvents());

     std::vector<std::vector<double>> pdfValues = prodpdf.getCompProbsAtDataPoints();

     for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {

     // Um... these two are switched? Weirdness...
           double currm12 = currData.getValue(m13, j);
           m12_pdf_hist.Fill(currm12, pdfValues[1][j]);

           double currm13 = currData.getValue(m12, j);
           m13_pdf_hist.Fill(currm13, pdfValues[0][j]);
           dalitzpm_pdf_hist.Fill(currm12, currm13, pdfValues[0][j]);

           double currm23 = cpuGetM23(currm12, currm13);
           m23_pdf_hist.Fill(currm23, pdfValues[0][j]);


           totalPdf += pdfValues[0][j];

      }


/*
    UnbinnedDataSet toyMC({m12, m13, eventNumber});
   plotter1.fillDataSetMC(toyMC, 5000000);

    fptype m12_tmp,m13_tmp,m23_tmp;
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
//   TCanvas fooo;
//   fooo.cd();      
//   m12_pdf_hist.Draw("c");
//   fooo.SaveAs("PDF.png");
//draw_1();

}

int runToyFit(Amp3Body *signal) {
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

//    FitManager datapdf(&prodpdf);
//    datapdf.fit();


//    prodpdf.setData(data);
    DalitzPlotter plotter(&prodpdf, signal);
//  TH1D* funct_hist = prodpdf.plotToROOT(m12);         
//  funct_hist->Draw("colz");  
    foodal.cd();
    TH2F *dalitzplot = plotter.make2D();
    dalitzplot->Draw("colz");
           
    foodal.SaveAs("dalitzpdf.png");
//    signaldalitz=signal;
//    const int nSamples = 1000;
//    cout<<rho0_amp_imag1.value<<endl;
//    cout<<rho0_amp_real1.value<<endl; 
    std::vector<std::vector<fptype>> ff =signal->fit_fractions();
    size_t n_res=ff.size();
    cout<<"n_res"<<":"<<n_res<<endl;
/*
    float mean[n_res];
    float rms[n_res];
    vector <float> fractions[n_res];

    for (int i=0;i<n_res;i++) mean[i] = rms[i] = 0;
    for (int i=0;i<nSamples;i++){
    
      datapdf.fit();
      for (int j=0;j<n_res; j++) {
          std::vector<std::vector<fptype>> _ff=signal->fit_fractions();
          fractions[j].push_back(_ff[i][i]);
          mean[j] += _ff[j][j];
          rms[j] += _ff[j][j]*_ff[j][j];
         }
     }

    double errfracList[n_res];
    
    TH1F* hFracs[n_res];
    for (int i=0;i<n_res;i++) {
      mean[i] /= nSamples;
      rms[i] = sqrt(rms[i]/nSamples-mean[i]*mean[i]);
      sprintf(strbuffer, "hfrac_res%d", i);
      hFracs[i] = new TH1F(strbuffer, "", 100, mean[i]-4*rms[i], mean[i]+4*rms[i]);
      for (int j=0;j<nSamples;j++)
          hFracs[i]->Fill(fractions[i][j]);
      
      errfracList[i] = hFracs[i]->GetRMS();//gaus1->GetParameter(2);
     }
*/
    for(size_t i = 0; i < n_res; i++)
//        for(size_t j = 0; j < n_res; j++)
        cout<<"Integral contribution for res # " << i << ": "<<ff[i][i]<<endl;

 

//    std::vector<Observable> vars;
//    vars.push_back(m12);
//    vars.push_back(m13);
//    vars.push_back(eventNumber);
//    UnbinnedDataSet grid(vars);
  
//   prodpdf.setData(&grid);
   
// Evaluate PDF and components at set data values
//   std::vector<std::vector<fptype>> pdfVals
//    = prodpdf.getCompProbsAtDataPoints();

//   incsum1=ResonanceSignalPdf(1);
//   incsum2=ResonanceSignalPdf(2);
//   incsum3=ResonanceSignalPdf(3);
//   incsum4=ResonanceSignalPdf(4);
  

 
    makeDalitzPlots(&prodpdf,signal,"plots_from_toy");
          
            
//    getDalitzdataPlots();
cout<<rho_1->get_amp_real()<<endl;
cout<<rho_1->get_amp_img()<<endl;
//drawFitPlots();
//draw_1();
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
//    Observable m12("m12", 0, 3.5);
//    Observable m13("m13", 0, 3.5);
//    EventNumber eventNumber("eventNumber");
//    m12.setNumBins(350);
//    m13.setNumBins(350);
        

    m12=Observable("m12", 0, 3.5);
    m13=Observable("m13", 0, 3.5);
    eventNumber= EventNumber("eventNumber");
    m12.setNumBins(350);
    m13.setNumBins(350);

    // Prepare the data
//   UnbinnedDataSet data({m12, m13, eventNumber});

    data= new UnbinnedDataSet({m12, m13, eventNumber});

    // Set up the model
    Amp3Body *signal = makeSignalPdf();

    // A wrapper for plotting without complex number segfault
    ProdPdf prodpdf{"prodpdf", {signal}};

    // Add nice tool for making data or plotting
    DalitzPlotter dplotter{&prodpdf, signal};

    // Read in data

    getToyData(toyFileName, app);
    

 
    try {
      return runToyFit(signal);
   } catch(const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        return 7;
    }


}

