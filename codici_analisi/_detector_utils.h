// file header per definire cose usate in tutti i file di analisi: path e nome del file da leggere, specifiche del detector, e le struct utilizzate.

#ifndef _DETECTOR_UTILS_H
#define _DETECTOR_UTILS_H

#include <vector>
#include <string>
#include <map>
#include <math.h>
#include "TMath.h"



// --- CARTELLA FILE
inline std::string basePath = "/Users/martina/Desktop/LabIntFond/3_glast/acquisizioni_prova/22-04-all-masked-30/"; // ****** MODIFICARE
inline std::string fileName = "TkrDataTaking_333001116.lif";

// --- DEFINIZIONE LUNGHEZZE E COSE DETECTOR [mm]
inline int EdgeWidth = 1;
inline double StripPitch = .228;
inline double LadderSeparation = .220;
inline double StripLength = 430; // lunghezza strip
inline std::vector<std::string> layerNames = {"Y0", "X0", "X1", "Y1", "Y2", "X2", "X3", "Y3", "Y4", "X4"}; // vettore di stringhe per i layer
inline std::map<std::string, int> layerIndex = {
    {"Y0",0}, {"X0",1}, {"X1",2}, {"Y1",3}, {"Y2",4},
    {"X2",5}, {"X3",6}, {"Y3",7}, {"Y4",8}, {"X4",9}
}; // mappa per gli indici, giusto per non stare a smattare
inline double z_coords[10] = {0, 41, 41+24, 2*41+24, 2*41+2*24, 3*41+2*24, 3*41+3*24, 4*41+3*24, 4*41+4*24, 5*41+4*24};
inline double BoxHeight = 24;
inline double InterBoxGap = 41;

// --- FUNZIONI PER CALCOLARE LA COORDINATE x,y e z IN mm
// Coordinate = EdgeWidth + StripPitch*StripNumber + (LadderSeparation +2*EdgeWidth -StripPitch)*int(StripNumber/384)
inline double getCoordinate(double strip) {
    double coord = EdgeWidth + StripPitch * strip +
                   (LadderSeparation + 2 * EdgeWidth - StripPitch) * (int(strip / 384));
    return coord;
}

inline int GetStripNumber(double coord) {
    for (int ladder = 0; ladder < 4; ++ladder) {
        double start = EdgeWidth + ladder * (384 * StripPitch + (LadderSeparation + 2 * EdgeWidth - StripPitch));
        double end = start + (384 * StripPitch);
        if (coord >= start && coord <= end) return (int)((coord - start) / StripPitch) + (ladder * 384);
    }
    return -1; // returna -1 se non ha funzionato, ad esempio perché la coordinata sta fuori dal layer
}

inline double getZ(int layerIdx){
    return z_coords[layerIdx];
}


// --- DEFINISCO STRUTTURE
// Cluster (più strip consecutive in stesso layer, in stesso evento)
struct Cluster{
    int mult; // molteplicità
    std::vector<int> strips; // strip interessate dal cluster
    double center; // centro del cluster (in unità di strip)
    double center_err; // errore (distr. uniforme nel cluster)
};

// Hit (un punto su un layer: strip singola colpita, o cluster)
struct Hit{
    int layerIdx;
    double pos; // unità in mm
    double err; // unità in mm
};

// Evento (riga pari di file -> total hits [hit block: layer, layer end, strip, chipID, channel ID] |  total time over threshold [TOT block: layer, layer end, TDC counts, TOT in us])
struct Event{
    int totalHits = 0; // numero totale di hit per evento, inizializzato a 0
    int hitsPerLayer[10] = {0}; // vettore coi numeri di hit per layer
    std::vector<int> stripsPerLayer[10]; // vettore delle strip colpite per layer, inizializzato a 0 di default;
    std::vector<int> positionsPerLayer[10]; // vettore per posizioni
    std::vector<Cluster> clustersPerLayer[10]; // vettore dei cluster per layer, inizializzato a 0 di default;
    std::vector<Hit> x_hits; // vettore con le hit sulle x
    std::vector<Hit> y_hits; // vettore con le hit sulle y
    int totalTOT = 0;
    int TOTsPerLayer[10] = {0};
    int evtID;
};

#endif
