// analisi file di output. quello che fa è plottare la distribuzione di probabilità del numero di hits per layer, per evento. l'obiettivo è leggere, per ogni riga pari, il numero di hits per layer; costruire quindi un istogramma 1D per ciascun layer.
// le cose da modificare sono nelle righe 15 e 84, dopo ****** MODIFICARE

#include <TFile.h>      // gestire file ROOT
#include <TH2D.h>       // istogrammi 1d
#include <TCanvas.h>
#include <TTree.h>      // gestire TTree (strutture dati ROOT)
#include <iostream>
#include <fstream>      // leggere / scrivere file di testo
#include <sstream>      // gestire stringhe
#include <vector>


// --- CARTELLA FILE
std::string basePath = "/Users/martina/Desktop/LabIntFond/3_glast/acquisizioni_prova/16-04-stackedY-30/"; // ****** MODIFICARE

// --- DEFINISCO LE STRUTTURE DI HIT E TOT
struct Event{
    int totalHits = 0; // numero totale di hit per evento, inizializzato a 0
    int hitPerLayer[10] = {0}; // vettore coi numeri di hit per layer
    int totalTOT = 0;
    int TOTPerLayer[10] = {0};
};

std::vector<std::string> layerNames = {"X0", "X1", "X2", "X3", "X4", "Y0", "Y1", "Y2", "Y3", "Y4"}; // vettore di stringhe per i layer

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
                std::cout << "Letto layer -> '" << layer << "'" << std::endl;
                for(int j = 0; j < 10; j++){
                    if(layerNames[j] == layer){
                        event.hitPerLayer[j]++;
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
        if(n % 2 == 0){events.push_back(event);
        }
    }
    return events;
}


// --- ISTOGRAMMI
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
    std::vector<Event> events = readFile(basePath + "TkrDataTaking_333001110.lif"); // ****** MODIFICARE

    // for loop per riempire istogrammi
    for (Event &event : events) {
        for (int i = 0; i < layerNames.size(); i++) {
                h[i]->Fill(event.hitPerLayer[i]);
                
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

