{
ifstream reader;
reader.open("PiPiPi_PWA_COEFFS_30_SAda_ampphs.txt");
float e1,e2,e3,e4;
int count1 = 0;
float x[1000];
float ym[1000];
float yp[1000];
while(reader>>e1>>e2>>e3){
    x[count1] = e1;
    ym[count1] = e2;
    yp[count1] =  e3;
    count1 ++;
}
cout<<"Lines found: "<<count1<<endl;
reader.close();
//reader.open("new_test.log");
//reader.open("new_test_2.log");
//reader.open("toyfit_t1.log");
//reader.open("pureswave_1.log");
//reader.open("toysigfit_n350_1.log");
//reader.open("toysigfit_n350_r0.log");
//reader.open("new_test_600x600.log");
//reader.open("test_1M.log");
char buff[1024];
const char *ms = "MIGRAD MINIMIZATION HAS CONVERGED.";
const char *mn = "pwa_coef_";
const int mnlen = strlen(mn);
int i,j; char nstr[50];
const int NF = 100;
float ym2[NF][50];
float yp2[NF][50];
float eym2[NF][50];
float eyp2[NF][50]; 
int foundarray[NF];
TNtuple *ntp = new TNtuple("ntp", "", "fcn");

float minfcn = 1e8;
float fcnarry[NF];
for (int iF=0;iF<NF;iF++){
//    cout<<"iF: "<<iF<<'\t'<<Form("toylogs/toysigfit_n350_r%d.log",iF+1)<<endl;
float yr = 0, eyr = 0;
bool found = 0;
foundarray[iF] = 0;
//ifstream *readerp = new ifstream(Form("toylogs/toysigfit_n1400_r%d.log",iF+1));
//ifstream *readerp = new ifstream(Form("toylogs/toysigfit_n2800_r%d.log",iF+1));
//ifstream *readerp = new ifstream(Form("toylogs/toysigfit_n700_r%d.log",iF+1));
//ifstream *readerp = new ifstream(Form("toylogs/toysigfit_n2000000_r%d.log",iF+1));
//_paraEff_histBkg_
//ifstream *readerp = new ifstream(Form("toylogs/datafit_paraEff_histBkg_c_n2000000_r%d.log",iF+1));
ifstream *readerp = new ifstream(Form("toylogs/toysigfit_n5000000_toy%d_defstart.log",iF));
//ifstream *readerp = new ifstream(Form("toylogs/datafit_paraEff_histBkg_mc100000_r%d.log",iF+1));
//ifstream *readerp = new ifstream(Form("toylogs/datafit_paraEff_histBkg_n1200_r%d.log",iF+1));
//ifstream *readerp = new ifstream(Form("toylogs/toysigfit_mcb100000_r%d.log",iF+1));
//ifstream *readerp = new ifstream(Form("toylogs/toysigfit_mcb200000_r%d.log",iF+1));
//ifstream *readerp = new ifstream(Form("toylogs/toysigfit_mcb500000_r%d.log",iF+1));
bool fcnfound=0;
readerp->seekg (0, readerp->beg);
while(!(readerp->eof())){
    readerp->getline(buff,2047);
if (strstr(buff, ms)) found = 1;
else {
    if (found){
        foundarray[iF] = 1;
        float fcnv = 0.;
        if (!fcnfound){
        const char *ptr = strstr(buff,"FCN=");
        if (ptr){
        sscanf(ptr+4, "%f FROM MIGRAD", &fcnv);
        ntp->Fill(fcnv); fcnfound = 1; fcnarry[iF] = fcnv; 
        if (minfcn> fcnv) minfcn = fcnv; 
        }
        }
        int ne = sscanf(buff, "%d pwa_coef_%d_%s %f %f %f %f", &i, &j, nstr, &e1, &e2, &e3, &e4);
        if (ne>1){
            float err = e2;
            if (!strstr(nstr, "phase")) {
                yr = e1; eyr = err;
            }
            else{
                ym2[iF][j] = yr;
                yp2[iF][j] = e1 ;
                cout<<iF<<'\t'<<j<<'\t'<<yp2[iF][j]<<'\t'<<yp[j]<<endl;
                if (yp2[iF][j]- yp[j] > 1.5*TMath::Pi() ) yp2[iF][j] -= 2*TMath::Pi();
                else if (yp2[iF][j]- yp[j] < -1.5*TMath::Pi() ) yp2[iF][j] += 2*TMath::Pi();
//                if (yp2[j]<0) yp2[j] += 2*TMath::Pi();
//               if (j>0&&yp2[iF][j]+TMath::PiOver2()<yp2[iF][j-1]) yp2[iF][j] += 2*TMath::Pi();
//                float a = ym2[iF][j];
                eym2[iF][j] = eyr;
                eyp2[iF][j] = err;
                if (j+1== count1) break;
            }
        }
}
}
}
cout<<"Here "<<iF<<endl;
readerp->close(); delete readerp;}
TCanvas* cv = new TCanvas("cv", "",1000,500);
cv->Divide(2);
cv->cd(1);
TGraph* gr = new TGraph(count1, x, ym);
gr->SetMarkerStyle(20);
gr->SetMarkerColor(kRed);
gr->SetTitle("S-Wave Magnitude");
gr->GetXaxis()->SetTitle("m^{2}(#pi^{+}#pi^{-}) [GeV^{2}/c^{4}]");
gr->Draw("AP");
TGraphErrors* gre[NF];
for (int iF=0;iF<NF;iF++){
if (!foundarray[iF]) continue;
//TGraphErrors* gre = new TGraphErrors(count1, x, ym2, 0, eym2);
gre[iF] = new TGraphErrors(count1, x, ym2[iF], 0, eym2[iF]);
gre[iF]->SetMarkerStyle(4);
int mycolor = fcnarry[iF] >= minfcn+10?kGray+2:kBlue;
mycolor = kBlue;
gre[iF]->SetMarkerColor(mycolor);
gre[iF]->SetLineColor(mycolor);
//gre->SetLineColor(kBlue);
/*gre->SetTitle("S-Wave Magnitude");
gre->Draw("AP");
gre->GetXaxis()->SetTitle("m^{2}(K^{+}K^{-}) [GeV^{2}/c^{4}]");*/
//if ( fcnarry[iF] < minfcn+1)
    gre[iF]->Draw("Psame");
}
gr->Draw("P");
cv->cd(2);
TGraph* gr2 = new TGraph(count1, x, yp);
gr2->SetMarkerStyle(20);
gr2->SetMarkerColor(kRed);
gr2->SetTitle("S-Wave Phase");
gr2->GetXaxis()->SetTitle("m^{2}(K^{+}K^{-}) [GeV^{2}/c^{4}]");
gr2->Draw("AP");
for (int iF=0;iF<NF;iF++){
if (!foundarray[iF]) continue;
//TGraphErrors* gre2 = new TGraphErrors(count1, x, yp2, 0, eyp2);
TGraphErrors* gre2 = new TGraphErrors(count1, x, yp2[iF], 0, eyp2[iF]);
gre2->SetMarkerStyle(4);
int mycolor = fcnarry[iF] >= minfcn+10?kGray+2:kBlue;
mycolor = kBlue;
gre2->SetMarkerColor(mycolor);
gre2->SetLineColor(mycolor);
/*gre2->SetTitle("S-Wave Phase");
gre2->Draw("AP");
gre2->GetXaxis()->SetTitle("m^{2}(K^{+}K^{-}) [GeV^{2}/c^{4}]");*/
//if ( fcnarry[iF] < minfcn+1)
    gre2->Draw("Psame");}
gr2->Draw("P");
TCanvas* c2 = new TCanvas("c2", "");
ntp->Draw("fcn");
}
