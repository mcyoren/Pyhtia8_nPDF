// main113.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Leif Lonnblad <leif.lonnblad@thep.lu.se>.

// Keywords: heavy ions; charged multiplicity; centrality;
// angantyr;

// This test program will generate Pb-Pb collisions at
// sqrt(S_NN)=2.76TeV using the Angantyr model for Heavy Ion
// collisions. The analysis will divide the event in centrality
// classes using the same observable as was used for p-Pb in the ATLAS
// analysis in arXiv:1508.00848 [hep-ex] (see main112.cc). The
// centrality classes are same as in the ALICE analysis in
// arXiv:1012.1657 [nucl-ex] although the actual observable used is
// not the same. Histograms of multiplicity distributions are measured
// for each centrality percentile.

// Note that heavy ion collisions are computationally quite CPU
// intensive and generating a single event will take around a second
// on a reasonable desktop. To get reasonable statistics, this program
// will take a couple of hours to run.

#include <TH1D.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include "Pythia8/Pythia.h"
#include "LHAPDF/LHAPDF.h"

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;
using namespace std;

vector<double> weights_PDF(vector<LHAPDF::PDF*> lha_npdf, LHAPDF::PDFSet lha_npdf_set, int pdfid, double pdfx, double pdfscale, double pdfpdf);
vector<int> create_dictionary(vector<string> *vector_string_with_conicidences);

int main(int argc, char* argv[]) {

  string str_seed = argv[1];
  string n_seeds = argv[3];
  string projectile = argv[6];
  string target     = argv[7];
  int n_seeds_int = std::stoi(n_seeds);
  string nEv = argv[4];
  int nEvents = std::stoi(nEv);
  string iftrig_str = argv[8];
  int iftrig = stoi(iftrig_str);
    cout << str_seed << "  " << n_seeds <<"  "<<nEv <<endl;
    cout << projectile << "  " << target <<"  "<<argv[5]<<"  "<< argv[2] <<endl;
  string nPDF_length_str = argv[9];
  int nPDF_length = stoi(nPDF_length_str);
  string N_opt_str = argv[10+2*nPDF_length];
  int N_opt = stoi(N_opt_str);
  vector<string> *nPDFA_str = new vector<string>();
  vector<string> *nPDFB_str = new vector<string>();
    for (int i = 0; i < nPDF_length; ++i) {
        nPDFB_str->push_back(argv[9+1+i]);
        nPDFA_str->push_back(argv[9+nPDF_length+1+i]);
        cout<<argv[9+1+i]<<"  "<<argv[9+nPDF_length+1+i]<<";  ";
    }
    cout<<endl;
    vector<int> nPDFA_dict = create_dictionary(nPDFA_str);
    vector<int> nPDFB_dict = create_dictionary(nPDFB_str);
    for(string i : *nPDFA_str) cout<<i<<" ";
    cout<<endl;
    for(string i : *nPDFB_str) cout<<i<<" ";
    cout<<endl;
    for (int i = 0; i < nPDF_length; ++i) {
        cout<<nPDFA_dict[i]<<"  "<<nPDFB_dict[i]<<";  ";
    }
    cout<<endl;
  for (int i = 0; i < N_opt; ++i) {
      cout<<argv[10+2*nPDF_length+1+i]<<endl;
  }


  int params[] = {1+3*nPDF_length, 0, 1+3*nPDF_length, 200, 0, 20};
  //double pt_bins[13] = {1.4, 1.7, 1.9, 2.1, 2.3, 2.6, 2.9, 3.4, 4.0, 4.5, 5., 6.5, 8.5};
  TH2D *pt_kstar = new TH2D("pt_kstar", "pt_kstar", params[0], params[1], params[2], params[3], params[4], params[5]);
  TH2D *pt_phi = new TH2D("pt_phi", "pt_phi", params[0], params[1], params[2], params[3], params[4], params[5]);
  TH2D *pt_pi = new TH2D("pt_pi", "pt_pi", params[0], params[1], params[2], params[3], params[4], params[5]);
  TH1D *hist = new TH1D("hist", "hist", 100, 0, 10);
  TH1D* sumch = new TH1D("sumch", "sumch", 500, 0, 500);
  TH2D* nwhist = new TH2D("nwhist", "nwhist", 500, 0, 500,100,0,100);
  TH2D* ncollhist = new TH2D("ncollhist", "ncollhist", 500, 0, 500,100,0,100);
  TH2D* ncollhist1 = new TH2D("ncollhist1", "ncollhist1", 500, 0, 500,100,0,100);
  TH2D* h_weights_centrality = new TH2D("h_weights_centrality", "h_weights_centrality", params[0], params[1], params[2], 5, 0, 5);
  TH1D* nchhist = new TH1D("nchhist", "nchhist", 500, 0, 500);
  TH1D* ncoll2hist = new TH1D("ncoll2hist", "ncoll2hist", 500, 0, 500);
  TH1D* infohist = new TH1D("infohist", "infohist", 2, 0, 2);
  TH1D* weighted_events = new TH1D("weighted_events", "weighted_events", params[0], params[1], params[2]);
  TH3D* y_charg = new TH3D("y_charg", "y_charg", params[0], params[1], params[2], 100, -5, 5,params[3], params[4], params[5]);
  TH3D* y_phi = new TH3D("y_phi", "y_phi", params[0], params[1], params[2], 100, -5, 5,params[3], params[4], params[5]);
  TH3D* y_pi = new TH3D("y_pi", "y_pi", params[0], params[1], params[2], 100, -5, 5,params[3], params[4], params[5]);

  Pythia pythia;

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + str_seed);

  /// Setup the beams.
  pythia.readString("Beams:idA = " + projectile); // The coper ion
  pythia.readString("Beams:idB = " + target); // The lead ions.
  pythia.readString("Beams:eCM = 200.");
  pythia.readString("Beams:frameType = 1");

  // Initialize the Angantyr model to fit the total and semi-includive
  // cross sections in Pythia within some tolerance.
 // pythia.readString("HeavyIon:SigFitErr = 0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
  // These parameters are typicall suitable for sqrt(S_NN)=5TeV
  //pythia.readString("HeavyIon:SigFitDefPar = 9.82,1.69,0.29,0.0,0.0,0.0,0.0,0.0");
  // A simple genetic algorithm is run for 20 generations to fit the
  // parameters.
  //pythia.readString("HeavyIon:SigFitNGen = 0");

  for (int i = 0; i < N_opt; ++i) {
      cout<<argv[10+2*nPDF_length+1+i]<<endl;
      pythia.readString(argv[10+2*nPDF_length+1+i]);
  }


    string file_name;
    file_name.append("");
    file_name.append(argv[5]);
    file_name.append("_");
    file_name.append(argv[2]);
    file_name.append(".root");

    TFile *outFile = TFile::Open(file_name.c_str(), "RECREATE");


  // Initialise Pythia.
  pythia.init();


  //////////////zalupa start////////////////////////////////
    Info info;
    //string pdfPath;
    //string pdfSet = nPDF_set;//"LHAPDF6:nCTEQ15np_197_79";
    //pdfPath = pythia.settings.word("xmlPath") + "../pdfdata";
    //PDFPtr oldPDF = make_shared<LHAPDF>(2212, "cteq6l", pdfPath, &info);
    //PDFPtr newPDF = make_shared<LHAPDF>(2212, pdfSet, &info);
    //Pythia8::LHAPDF pdfs_nnpdf(2112, "LHAPDF6:cteq6l1/0", &info);
    //Pythia8::LHAGrid1 pdfs_nnpdf(stoi(projectile), "GKG18_DPDF_FitA_NLO_0000.dat", pdfPath, &info);
   // LHAPDF pdfs_nnpdf_lha(stoi(projectile), pdfSet, &info);

    vector<LHAPDF::PDF*> def_npdf = LHAPDF::mkPDFs("cteq6l1");

    vector<vector<LHAPDF::PDF*>> lha_npdf_vec;
    vector<LHAPDF::PDFSet> lha_npdf_set_vec;

    for(string nPDFB : *nPDFB_str)
    {
        lha_npdf_vec.push_back(LHAPDF::mkPDFs(nPDFB));
        lha_npdf_set_vec.push_back(LHAPDF::getPDFSet(nPDFB));
    }

    vector<vector<LHAPDF::PDF*>> lha_npdfA_vec;
    vector<LHAPDF::PDFSet> lha_npdfA_set_vec;

    vector<bool> usepdfA ; string zero = "0"; vector<LHAPDF::PDF*> NULLPDF; LHAPDF::PDFSet NULLPDFSet;

    for(string nPDFA : *nPDFA_str)
    {
        if(nPDFA!=zero)
        {
            usepdfA.push_back(true);
            lha_npdfA_vec.push_back(LHAPDF::mkPDFs(nPDFA));
            lha_npdfA_set_vec.push_back(LHAPDF::getPDFSet(nPDFA));
        }
        else
        {
            usepdfA.push_back(false);
            lha_npdfA_vec.push_back(NULLPDF);
            lha_npdfA_set_vec.push_back(NULLPDFSet);
        }
    }



  // Loop over events.

  double sumw = 0.0;
  double sigmaALL = 0.0;


  for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {
    if (!pythia.next()) continue;

    if ((iEvent+1) % 10000 == 0) cout << "N_Event: " << iEvent << endl;


      //////////////zalupa start////////////////////////////////
      auto pdfid1 = pythia.info.id1pdf();
      auto pdfid2 = pythia.info.id2pdf();
      auto pdfx1 = pythia.info.x1pdf();
      auto pdfx2 = pythia.info.x2pdf();
      auto pdfscale = pythia.info.Q2Ren();
      auto pdfpdf1 = pythia.info.pdf1();
      auto pdfpdf2 = pythia.info.pdf2();
      //auto pdfscaleQ1 = pythia.info.QRen();      //auto pdfscalefac = pythia.info.Q2Fac();      //auto pdfscaleQ1fac = pythia.info.QFac();
     // cout<<pdfid2<<" "<<pdfpdf2<<endl;
      if(stoi(target)==2112)
      {
          if(pdfid2==1) pdfid2*=2;
          else if(pdfid2==2) pdfid2/=2;
      }
      if(stoi(projectile)==2112)
      {
          if(pdfid1==1) pdfid1*=2;
          else if(pdfid1==2) pdfid1/=2;
      }
      //cout<<pdfid2<<" "<<pdfpdf2<<endl;  //cout<<pdfpdf1<<" yolo "<<pdfpdf2<<endl;

      pdfpdf1 = def_npdf[0]->xfxQ2(pdfid1, pdfx1, pdfscale) / pdfx1;
      pdfpdf2 = def_npdf[0]->xfxQ2(pdfid2, pdfx2, pdfscale) / pdfx2;

      //cout<<def_npdf[0]->xfxQ2(pdfid2, pdfx2, pdfscale)<<" "<<def_npdf[0]->xfxQ(pdfid2, pdfx2, pdfscaleQ1)<<" "
      //<<def_npdf[0]->xfxQ2(pdfid2, pdfx2, pdfscalefac)<<" "<<def_npdf[0]->xfxQ(pdfid2, pdfx2, pdfscaleQ1fac)<<endl;
      //cout<<pdfpdf1<<" "<<pdf
      // pdf2<<endl;
      //cout<<pdfid2<<" "<<pdfid1<<" "<<pdfs_nnpdf.xf(pdfid2, pdfx2, pdfscale)<<" "<<def_npdf[0]->xfxQ2(pdfid1, pdfx2, pdfscale)<<endl;
      vector<double> pdf_weights;
      pdf_weights.push_back(1.0);
      for (int i: nPDFB_dict) {
          vector<double> pdf_weights_loc;
          pdf_weights_loc = weights_PDF(lha_npdf_vec[i],lha_npdf_set_vec[i],pdfid2,pdfx2,pdfscale,pdfpdf2);
          for(double pdf : pdf_weights_loc) pdf_weights.push_back(pdf);
      }

      int k = 0;
      for (int i: nPDFA_dict) {
          vector<double> pdf_weights_loc;
          if(usepdfA[i]){
              pdf_weights_loc = weights_PDF(lha_npdfA_vec[i],lha_npdfA_set_vec[i],pdfid1,pdfx1,pdfscale,pdfpdf1);
              for (int j = 0; j < 3; ++j) pdf_weights[1+3*k+j]*=pdf_weights_loc[j];
          }
          k++;
      }


      for (int i = 0; i < 3*nPDF_length+1; ++i) {
          //cout<<pdf_weights[i]<<"   ";
          weighted_events->Fill(i,pdf_weights[i]);
      }
      //cout<<"yolo: "<<pdf_weight_c<<" "<<pdf_weight_up<<" "<<pdf_weight_down<<endl;

      //////////////zalupa ends////////////////////////////////

    double summcharge = 0.0;
    int nch = 0;
    bool fwd = false, bkw = false;
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle &p = pythia.event[i];
      if (p.isFinal() && p.isCharged()) {
        double eta = p.eta();
        if ( fabs(eta) > 3.1 && fabs(eta) < 4. ){
          summcharge += fabs(p.charge());
        }
        if ( fabs(eta) < 0.5 ){
            nch++;
        }
        if (eta > 3.1 && eta < 4.) fwd = true;
        if (eta < -3.1 && eta > -4.) bkw = true;
      }
    }
    if((!fwd || !bkw)&&iftrig==1) continue;

    nchhist->Fill(nch);

    // Keep track of the sum of waights
    double weight = pythia.info.weight();
    sumw += weight;
    double sigma = pythia.info.sigmaGen();
    sigmaALL += sigma;


    // Also fill the number of (absorptively and diffractively)
    // wounded nucleaons.
      int nw = 2, ncoll = 1, ncoll1 = 1;
      if (stoi(target)>2212)
    {
        nw = pythia.info.hiInfo->nAbsTarg() +
              pythia.info.hiInfo->nDiffTarg() +
              pythia.info.hiInfo->nAbsProj() +
              pythia.info.hiInfo->nDiffProj();
            ///ncoll calc by Norbert
        ncoll = pythia.info.hiInfo->nCollNDTot();
        ncoll1 =pythia.info.hiInfo->nAbsProj() +
                pythia.info.hiInfo->nDiffProj() +
                pythia.info.hiInfo->nAbsTarg() +
                pythia.info.hiInfo->nDiffTarg() -
                pythia.info.hiInfo->nCollND() -
                pythia.info.hiInfo->nCollDD();
    }

    sumch->Fill(summcharge, weight);
    nwhist->Fill(summcharge,nw,weight);

    ncollhist->Fill(summcharge,ncoll,weight);
    ncollhist1->Fill(summcharge,ncoll1,weight);
    ncoll2hist->Fill(ncoll1,weight);

    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle &p = pythia.event[i];
      if(fabs(p.eta()) < 0.5)
      {
          for (int j = 0; j < 3*nPDF_length+1; ++j) {
              if (p.isFinal() && p.isCharged()) {
                  pt_kstar->Fill(j, p.pT(), weight*pdf_weights[j]);
                  hist->Fill(p.pT());
              }

              if (p.id() == 111) {
                  pt_pi->Fill(j, p.pT(), weight*pdf_weights[j]);
              }

              if (p.id() == 333) {
                  pt_phi->Fill(j, p.pT(), weight*pdf_weights[j]);
              }
          }
      }
        for (int j = 0; j < 3*nPDF_length+1; ++j) {
            if (p.isFinal() && p.isCharged()) {
                y_charg->Fill(j, p.eta(), p.pT(), weight*pdf_weights[j]);
            }

            if (p.id() == 111) {
                y_pi->Fill(j ,p.eta(), p.pT(), weight*pdf_weights[j]);
            }

            if (p.id() == 333) {
                y_phi->Fill(j, p.eta(), p.pT(), weight*pdf_weights[j]);
            }
        }
    }
  }


  sumch->Write();
  pt_kstar->Write();
  pt_phi->Write();
  pt_pi->Write();
  hist->Write();
  nwhist->Write();
  nchhist->Write();
  ncollhist->Write();
  ncollhist1->Write();
  ncoll2hist->Write();
  weighted_events->Write();
  y_charg->Write();
  y_phi->Write();
  y_pi->Write();


    infohist->SetBinContent(1,sigmaALL/sumw);
  infohist->SetBinContent(2,sumw);

  outFile->Close();

    std::cout << "SIGMA = " << sigmaALL/sumw*1000 << std::endl;
    std::cout << "SUMW = " << sumw << std::endl;


  // And we're done!
  return 0;
}

vector<double> weights_PDF(vector<LHAPDF::PDF*> lha_npdf, LHAPDF::PDFSet lha_npdf_set, int pdfid, double pdfx, double pdfscale, double pdfpdf)
{
    vector<double> pdf_vec;
    for(auto const &pdf:  lha_npdf)
    {
        pdf_vec.push_back(pdf->xfxQ2(pdfid, pdfx, pdfscale));
    }
    //auto xgAll = lha_npdf[0]->xfxQ2(pdfid1, pdfx1, pdfscale);
    auto xgErr = lha_npdf_set.uncertainty(pdf_vec);
//
    auto pdf2_c = xgErr.central / pdfx;
    auto pdf2_up = (xgErr.central + xgErr.errplus) / pdfx;
    auto pdf2_down = (xgErr.central - xgErr.errminus) / pdfx;

    double pdf_weight_c=1.0, pdf_weight_up=1.0, pdf_weight_down = 1.0;
    vector<double> result_w;
    if(pdfpdf>0.000001&&pdf2_c>0.000001)
    {

        pdf_weight_c = ( pdf2_c) / ( pdfpdf);
        pdf_weight_up = ( pdf2_up) / ( pdfpdf);
        pdf_weight_down = ( pdf2_down) / ( pdfpdf);
    }else{
        //pdf_weight_c = 0;
        //pdf_weight_up = 0;
        //pdf_weight_down = 0;
    }
    result_w.push_back(pdf_weight_c);
    result_w.push_back(pdf_weight_up);
    result_w.push_back(pdf_weight_down);

    return result_w;
}

vector<int> create_dictionary(vector<string> *vector_string_with_conicidences)
{
    vector<int> dictionary;
    vector<bool> coincidences;
    for (int i = 0; i < vector_string_with_conicidences->size(); ++i)
    {
        dictionary.push_back(-999);
        coincidences.push_back(false);
    }
    int k = -1;
    for (int i = 0; i < vector_string_with_conicidences->size(); ++i) {
        if(!coincidences[i])
        {
            k++;
            dictionary[i]=k;
        }
        for (int j = i+1; j < vector_string_with_conicidences->size(); ++j) {
            if(vector_string_with_conicidences->at(i)==vector_string_with_conicidences->at(j)&&!coincidences[j])
            {
                dictionary[j]=k;
                coincidences[j]=true;
            }
        }
    }
    for (int i = vector_string_with_conicidences->size()-1; i >=0 ; i--) {
        if(coincidences[i]) vector_string_with_conicidences->erase(vector_string_with_conicidences->begin()+i);
    }
    return dictionary;
}


