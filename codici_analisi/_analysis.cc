// analisi file di output. quello che fa è plottare la distribuzione di probabilità del numero di hits per layer, per evento. l'obiettivo è leggere, per ogni riga pari, il numero di hits per layer; costruire quindi un istogramma 1D per ciascun layer.
// modifica per cluster: consente di catalogare anche i cluster (strip consecutive attive), in un vettorone che ha per elementi (layer, (strips coinvolte), molteplicità del cluster)

#include "_detector_utils.h"  // serve per le robe sul detector
#include "_readfile.cc"  // serve per le funzioni

#include <TFile.h>      // gestire file ROOT
#include <TH1D.h>
#include <TH2D.h>       // istogrammi 1d
#include <TCanvas.h>

// --- ISTOGRAMMI DISTRIBUZIONE DELLE HITS
void HitsDistributions(){
    // definisco (vettore di) istogrammi
    std::vector<TH1D*> h;
    // crea istogrammi
    for (int i = 0; i < layerNames.size(); i++){
        h.push_back(
            new TH1D(
                ("h_" + layerNames[i]).c_str(),
                (layerNames[i] + ";hits;counts").c_str(),
                21, 0, 21)
                    );
    }
    
    // vettore con tutti gli eventi letti dal file
    std::vector<Event> events = readFile(basePath + fileName); // ****** MODIFICARE

    // for loop per riempire istogrammi
    for (Event &event : events) {
        for (int i = 0; i < layerNames.size(); i++) {
                h[i]->Fill(event.hitsPerLayer[i]);
                
            }
        }

    // disegna istogrammi
    TCanvas *c = new TCanvas("c", "Hits distributions", 1200, 1000);
    c->Divide(5,2); // divide il canvas in 2x5 riquadri
    
    for(int i = 0; i < layerNames.size(); i++){
        c->cd(i+1); // seleziona il riquadro del canvas in cui disegnare
        h[i]->Draw("HIST");
    }
    
}

// --- ISTOGRAMMI POSIZIONE (non servono tanto guardi i cluster...)
void HitsPosDistributions(){
    // definisco (vettore di) istogrammi
    std::vector<TH1D*> h;
    // crea istogrammi
    for (int i = 0; i < layerNames.size(); i++){
        h.push_back(
            new TH1D(
                ("h_" + layerNames[i]).c_str(),
                (layerNames[i] + ";position [mm];counts").c_str(),
                360, 0, 360)
                    );
    }
    
    // vettore con tutti gli eventi letti dal file
    std::vector<Event> events = readFile(basePath + fileName); // ****** MODIFICARE

    // for loop per riempire istogrammi
    for (Event &event : events) {
        for (int i = 0; i < layerNames.size(); i++) {
            for (int j = 0; j < event.positionsPerLayer[i].size(); j++){
                h[i]->Fill(event.positionsPerLayer[i][j]);
            }
            }
        }

    // disegna istogrammi
    TCanvas *c = new TCanvas("c", "Hits distributions", 1200, 1000);
    c->Divide(5,2); // divide il canvas in 2x5 riquadri
    
    for(int i = 0; i < layerNames.size(); i++){
        c->cd(i+1); // seleziona il riquadro del canvas in cui disegnare
        h[i]->SetFillColor(kBlue);
        h[i]->SetLineColor(kWhite);
        h[i]->Draw("HIST");
    }
    
}

// --- ISTOGRAMMI DISTRIBUZIONE DEI CLUSTER (strips colpite)
void ClusterMultDistributions(){
    // definisco (vettore di) istogrammi
    std::vector<TH1D*> h;
    // crea istogrammi
    for (int i = 0; i < layerNames.size(); i++){
        h.push_back(
            new TH1D(
                ("h_cluster_" + layerNames[i]).c_str(),
                (layerNames[i] + ";cluster size [strips];counts").c_str(),
                 60, 0, 60)
                    );
    }
    
    // vettore con tutti gli eventi letti dal file
    std::vector<Event> events = readFile(basePath + fileName); // ****** MODIFICARE

    // for loop per riempire istogrammi
    for (Event &event : events) {
        for (int i = 0; i < layerNames.size(); i++) {
            for (const Cluster &cl : event.clustersPerLayer[i]) {
                h[i]->Fill(cl.mult);
            }
        }
    }

    // disegna istogrammi
    TCanvas *c = new TCanvas("c", "Cluster multiplicity distributions", 1200, 1000);
    c->Divide(5,2); // divide il canvas in 2x5 riquadri
    
    for(int i = 0; i < layerNames.size(); i++){
        c->cd(i+1); // seleziona il riquadro del canvas in cui disegnare
       // h[i]->SetFillColor(kBlue);
       // h[i]->SetLineColor(kBlack);
        h[i]->Draw("HIST");
    }
    
}

void ClusterPosDistributions(){
    // definisco (vettore di) istogrammi
    std::vector<TH2D*> h;
    // crea istogrammi
    for (int i = 0; i < layerNames.size(); i++){
        h.push_back(
            new TH2D(
                ("h_cluster_" + layerNames[i]).c_str(),
                (layerNames[i] + ";cluster size [strips];strip;counts").c_str(),
                50, 10, 60,
                1536, 0, 1536
                )
                    );
    }
    
    // vettore con tutti gli eventi letti dal file
    std::vector<Event> events = readFile(basePath + fileName); // ****** MODIFICARE

    // for loop per riempire istogrammi
    for (Event &event : events) {
        for (int i = 0; i < layerNames.size(); i++) {
            for (const Cluster &cl : event.clustersPerLayer[i]) {
                int center = cl.strips[cl.strips.size() / 2]; // calcolo la strip al centro del cluster
                h[i]->Fill(cl.mult,center);
            }
        }
    }

    // disegna istogrammi
    TCanvas *c = new TCanvas("c", "Cluster position distributions", 1200, 1000);
    c->Divide(5,2); // divide il canvas in 2x5 riquadri
    
    for(int i = 0; i < layerNames.size(); i++){
        c->cd(i+1); // seleziona il riquadro del canvas in cui disegnare
        h[i]->Draw("LEGO");
    }
    
}

