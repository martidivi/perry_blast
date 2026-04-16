// analisi noise occupancy. per ogni acquisizione (soglia fissata), costruisce istogrammi 2d di occupancy vs (soglia, strip). 
// il codice legge i file relativi a ciascun layer con la funzione readLayer, che ritorna il nome del layer e un elemento LayerData (vettore di elementi (strip, occupancy));
// legge la cartella relativa a ogni acquisizione a soglia fissata con la funzione readAcquisition, che ritorna la soglia e un vettore di elementi LayerData;
// infine l'analisi viene fatta con la funzione NoiseOccupancy(), che crea un vettore di 10 istogrammi 2d (uno per ogni layer), li riempie e li disegna.

// per runnare: modificare la cartella (riga 18, dove c'è  ******* MODIFICARE)

#include <TFile.h>      // gestire file ROOT
#include <TH2D.h>       // istogrammi 2d
#include <TCanvas.h>
#include <iostream>
#include <fstream>      // leggere / scrivere file di testo
#include <sstream>      // gestire stringhe
#include <vector>
#include <map>          // dizionari

// --- CARTELLA FILE ******* MODIFICARE
std::string basePath = "/Users/martina/Desktop/LabIntFond/3_glast/noise_occupancy/";

// --- DICHIARAZIONE STRUTTURE
// struttura di ogni file
struct StripOccupancy{
    int strip;
    double occupancy;
};

struct LayerData{
    std::string name; // nome del layer (es. "X0")
    std::vector<StripOccupancy> strips; // un vettore di occupancy delle strip
};

// struttura di ogni cartella
struct Acquisition{
    int threshold;  // soglia, es. 30
    std::vector<LayerData> layers; // vettore coi layer
};

// soglie, runID, nomi di layer da usare per leggere i file
std::vector<int> thresholdValues = {20, 25, 30, 35, 40, 45, 50, 55, 60}; // vettore di interi
std::vector<std::string> runID = {"333001085", "333001086", "333001082", "333001089", "333001083", "333001090", "333001084", "333001091", "333001087"}; // vettore di stringhe
std::vector<std::string> layerNames = {"X0", "X1", "X2", "X3", "X4", "Y0", "Y1", "Y2", "Y3", "Y4"}; // vettore di stringhe anche qua


// --- LETTURA FILE PER SINGOLO LAYER
LayerData readLayer(const std::string& filename, const std::string& layerName){
    // crea layer vuoto
    LayerData layer;
    layer.name = layerName;
    // apre file
    std::ifstream infile(filename);
    // per saltare le prime righe di file
    std::string line;
    getline(infile, line);
    getline(infile, line);
    // crea oggetto "event" con struttura StripOccupancy (strip, occupancy)
    StripOccupancy event;
    // legge righe e le cataloga nel vettore
    while(infile >> event.strip >> event.occupancy){
        // aggiunge l'ultimo oggetto "event" alla fine del vettore
        layer.strips.push_back(event);
    }
    return layer;
}

// --- LETTURA CARTELLA PER OGNI ACQUISIZIONE A SOGLIA FISSATA
Acquisition readAcquisition(int threshold, std::string runID){
    // crea cartella di acquisizione vuota
    Acquisition acq;
    acq.threshold = threshold;

    //
    std::string folder = basePath + "14-04-NO-" + std::to_string(threshold);
    
    for (int i = 0; i < layerNames.size(); i++) {// file tipo: TkrNoiseOccupancy_LayerX2_333001087.tnt
        std::string layername = layerNames[i];
        
        std::string filename = folder + "/TkrNoiseOccupancy_Layer" + layername + "_" + runID + ".tnt";
        
        LayerData layer = readLayer(filename, layername);

        acq.layers.push_back(layer);
    }
    
    return acq;
}


// --- ANALISI
void NoiseOccupancy(){
    // definisco (vettore di) istogrammi
    std::vector<TH2D*> h;
    // crea istogrammi
    for (int i = 0; i < layerNames.size(); i++){
        h.push_back(
            new TH2D(
                ("h_" + layerNames[i]).c_str(),
                (layerNames[i] + ";strip;threshold;occupancy").c_str(),
                1536, 0, 1536,
                9, 17.5, 62.5)
                    );
    }
    
    // vettore con tutte le cartelle di acquisizione
    std::vector<Acquisition> all;
    // per ogni soglia, riempie il vettore
    for (int i = 0; i < thresholdValues.size(); i++) {
            Acquisition acq = readAcquisition(thresholdValues[i], runID[i]);
            all.push_back(acq);
        }
    
    // riempie istogrammi
    for (Acquisition &acq : all) {
        for (int i = 0; i < acq.layers.size(); i++) {
            for (const StripOccupancy &event : acq.layers[i].strips) {
                h[i]->Fill(event.strip, acq.threshold, event.occupancy);
            }
        }
    }
    
    // disegna istogrammi
    TCanvas *c = new TCanvas("c", "Noise Occupancy", 1200, 1000);
    c->Divide(5,2); // divide il canvas in 2x5 riquadri
    
    for(int i = 0; i < layerNames.size(); i++){
        c->cd(i+1); // seleziona il riquadro del canvas in cui disegnare
        h[i]->Draw("LEGO");
    }
    
}
