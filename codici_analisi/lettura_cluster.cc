// analisi file di output. quello che fa è plottare la distribuzione di probabilità del numero di hits per layer, per evento. l'obiettivo è leggere, per ogni riga pari, il numero di hits per layer; costruire quindi un istogramma 1D per ciascun layer.
// modifica per cluster: consente di catalogare anche i cluster (strip consecutive attive), in un vettorone che ha per elementi (layer, (strips coinvolte), molteplicità del cluster)


#include <TFile.h>      // gestire file ROOT
#include <TH2D.h>       // istogrammi 1d
#include <TCanvas.h>
#include <TTree.h>      // gestire TTree (strutture dati ROOT)
#include <iostream>
#include <fstream>      // leggere / scrivere file di testo
#include <sstream>      // gestire stringhe
#include <vector>
#include <algorithm>    // necessario per std::sort (ordinare le strip colpite nei cluster)


// --- CARTELLA FILE
std::string basePath = "/Users/martina/Desktop/LabIntFond/3_glast/acquisizioni_prova/15-04-30-all-masked/"; // ****** MODIFICARE

// --- DEFINISCO STRUTTURA DEI CLUSTER
struct Cluster {
    int mult; // molteplicità
    std::vector<int> strips; // strip interessate dal cluster
};

// --- DEFINISCO LE STRUTTURE DI HIT E TOT IN EVENT
struct Event{
    int totalHits = 0; // numero totale di hit per evento, inizializzato a 0
    int hitsPerLayer[10] = {0}; // vettore coi numeri di hit per layer
    std::vector<int> stripsPerLayer[10]; // vettore delle strip colpite per layer, inizializzato a 0 di default;
    std::vector<Cluster> clustersPerLayer[10]; // vettore dei cluster per layer, inizializzato a 0 di default;
    int totalTOT = 0;
    int TOTsPerLayer[10] = {0};
};


std::vector<std::string> layerNames = {"X0", "X1", "X2", "X3", "X4", "Y0", "Y1", "Y2", "Y3", "Y4"}; // vettore di stringhe per i layer


// --- FUNZIONE PER LEGGERE I CLUSTER, PER OGNI EVENTO, PER OGNI LAYER
std::vector<Cluster> findClusters(std::vector<int> strips){
    // ordina le strip colpite in ordine crescente
    std::sort(strips.begin(), strips.end());
    // definisce vettore finale di TUTTI i cluster trovati
    std::vector<Cluster> clusters;
    // inizializza un cluster e ci mette la prima strip colpita (poi posso toglierla)
    Cluster currentCluster;
    currentCluster.strips.push_back(strips[0]);
    for(int i = 1; i < strips.size(); i++){
        if( strips[i] - strips[i-1] == 1 ){
            // faccio push back di qualcosa e incremento la molteplicità
            currentCluster.strips.push_back(strips[i]);
        }
        else{ // COMMENTARE LA CONDIZIONE DELL'IF SE VUOI ANCHE MOLTEPLICITA' 1
          //  if(currentCluster.strips.size() > 1){
                currentCluster.mult = currentCluster.strips.size(); // calcola molteplicità
                clusters.push_back(currentCluster); // salvo il cluster corrente nel vettorone dei clusters
           // }
            // reset prima del prossimo cluster!!
            currentCluster = Cluster();
            currentCluster.strips.push_back(strips[i]);
        }
    }
    // questo per l'ultimo cluster
    // COMMENTARE LA CONDIZIONE DELL'IF SE VUOI ANCHE MOLTEPLICITA' 1
 //  if(currentCluster.strips.size() > 1){
        currentCluster.mult = currentCluster.strips.size();
        clusters.push_back(currentCluster);
  //  }
    return clusters;
}


// --- FUNZIONE CHE LEGGE IL FILE E RITORNA UN VETTORE DI EVENTI, UNO PER OGNI RIGA (pari) DI FILE
std::vector<Event> readFile(const std::string& filename){
    // apre file
    std::ifstream infile(filename);
    // salta la prima riga
    std::string line;
    getline(infile, line);

    std::vector<Event> events; // crea vettore di eventi
    
    // loop while su ogni riga: definisce n (primo numerino che mi dice quali righe prendere, solo quelle con n pari) e eventID (secondo numerozzo della riga)
    int n, evtID, totalHits;
    while(infile >> n >> evtID >> totalHits ){
        Event event = Event(); // crea oggetto event
        event.totalHits = totalHits;
            for(int i = 0; i < totalHits; i++){
                // MODIFICA IL VETTORE DI HIT PER LAYER
                // incremento contatore di layer se lo trova in layerNames
                std::string layer;
                int layerEnd, strip, chipID, channelID;
                infile >> layer >> layerEnd >> strip >> chipID >> channelID;
               // std::cout << "Letto layer -> '" << layer << "' e strip ->'" << strip << "'" << std::endl;
                for(int j = 0; j < 10; j++){
                    if(layerNames[j] == layer){
                        event.hitsPerLayer[j]++; // incrementa il contatore degli hit per layer
                        event.stripsPerLayer[j].push_back(strip); // aggiunge strips per layer
                        
                    }
                }
            }
            infile >> event.totalTOT;
            for(int i = 0; i < event.totalTOT; i++){
               // salta i blocchi dei TOT
                std::string layer;
                int layerEnd, TDCcounts, TOT;
                infile >> layer >> layerEnd >> TDCcounts >> TOT;
            }
            if(n % 2 == 0){
                for(int j = 0; j < 10; j++){
                    if(event.stripsPerLayer[j].size() > 0){
                        event.clustersPerLayer[j] = findClusters(event.stripsPerLayer[j]);
                    }
                    for(const Cluster &cl : event.clustersPerLayer[j]){
                        std::cout << "Layer -> '" << layerNames[j] << "' | cluster mult -> " << cl.mult << " | strips -> ";
                        for(int s : cl.strips){
                            std::cout << s << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            events.push_back(event);
            }
    }
    return events;
}



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
    std::vector<Event> events = readFile(basePath + "TkrDataTaking_333001103.lif"); // ****** MODIFICARE

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
                (layerNames[i] + ";strip;counts").c_str(),
                1536, 0, 1536)
                    );
    }
    
    // vettore con tutti gli eventi letti dal file
    std::vector<Event> events = readFile(basePath + "TkrDataTaking_333001103.lif"); // ****** MODIFICARE

    // for loop per riempire istogrammi
    for (Event &event : events) {
        for (int i = 0; i < layerNames.size(); i++) {
            for (int j = 0; j < event.stripsPerLayer[i].size(); j++){
                h[i]->Fill(event.stripsPerLayer[i][j]);
            }
            }
        }

    // disegna istogrammi
    TCanvas *c = new TCanvas("c", "Hits distributions", 1200, 1000);
    c->Divide(5,2); // divide il canvas in 2x5 riquadri
    
    for(int i = 0; i < layerNames.size(); i++){
        c->cd(i+1); // seleziona il riquadro del canvas in cui disegnare
        h[i]->SetFillColor(kBlue);
        h[i]->SetLineColor(kBlue);
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
                7, 0.5, 7.5)
                    );
    }
    
    // vettore con tutti gli eventi letti dal file
    std::vector<Event> events = readFile(basePath + "TkrDataTaking_333001103.lif"); // ****** MODIFICARE

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
                7, 0.5, 7.5,
                1536, 0, 1536
                )
                    );
    }
    
    // vettore con tutti gli eventi letti dal file
    std::vector<Event> events = readFile(basePath + "TkrDataTaking_333001103.lif"); // ****** MODIFICARE

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

