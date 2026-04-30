#include "_detector_utils.h"
#include <iostream>
#include <fstream>      // leggere / scrivere file di testo
#include <sstream>      // gestire stringhe
#include <algorithm>    // necessario per std::sort (ordinare le strip colpite nei cluster)

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
                currentCluster.center = (currentCluster.strips.front() + currentCluster.strips.back())/2;
                currentCluster.center_err = currentCluster.mult/sqrt(12);
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
        currentCluster.center = (currentCluster.strips.front() + currentCluster.strips.back())/2;
        currentCluster.center_err = currentCluster.mult/sqrt(12);
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
        event.evtID = evtID;
            for(int i = 0; i < totalHits; i++){
                // MODIFICA IL VETTORE DI HIT PER LAYER
                // incremento contatore di layer se lo trova in layerNames
                std::string layer;
                int layerEnd, strip, chipID, channelID, position;
                infile >> layer >> layerEnd >> strip >> chipID >> channelID;
               // std::cout << "Letto layer -> '" << layer << "' e strip ->'" << strip << "'" << std::endl;
                for(int j = 0; j < 10; j++){
                    if(layerNames[j] == layer){ // STRIP DA MASCHERARE
                        if(!( (layer == "Y2" && strip == 962) || (layer == "X3" && (strip == 193 || strip == 195)) )){
                            event.hitsPerLayer[j]++; // incrementa il contatore degli hit per layer
                            event.stripsPerLayer[j].push_back(strip); // aggiunge strips per layer
                            position = EdgeWidth + StripPitch*strip + (LadderSeparation +2*EdgeWidth -StripPitch)*int(strip/384);
                            event.positionsPerLayer[j].push_back(position);
                        }
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
                for(int j = 0; j < 10; j++){// per ogni layer
                    if(event.stripsPerLayer[j].size() > 0){// se ci sono strip accese in quel layer
                        event.clustersPerLayer[j] = findClusters(event.stripsPerLayer[j]); // trova cluster
                    }
                    for(const Cluster &cl : event.clustersPerLayer[j]){ // per ogni cluster in quel layer
                        // OCCHIO A DECOMMENTARE if(currentCluster.strips.size() > 1)
                        // costruisco l'hit associata
                        Hit hit;
                        hit.layerIdx = j;
                        hit.pos = getCoordinate(cl.center); // unità in mm
                        hit.err = cl.center_err*StripPitch; // unità in mm
                        // riempio x_hits e y_hits
                        if(layerNames[j][0] == 'X') event.x_hits.push_back(hit);
                        else event.y_hits.push_back(hit);

                        /*std::cout << "Layer -> '" << layerNames[j] << "' | cluster mult -> " << cl.mult << " | strips -> ";
                        for(int s : cl.strips){
                            std::cout << s << " ";
                        }
                        std::cout << std::endl;*/
                    }
                }
            events.push_back(event);
            }
    }
    return events;
}
