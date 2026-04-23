// analisi file di output. quello che fa è plottare la distribuzione di probabilità del numero di hits per layer, per evento. l'obiettivo è leggere, per ogni riga pari, il numero di hits per layer; costruire quindi un istogramma 1D per ciascun layer.
// modifica per cluster: consente di catalogare anche i cluster (strip consecutive attive), in un vettorone che ha per elementi (layer, (strips coinvolte), molteplicità del cluster)


#include <TFile.h>      // gestire file ROOT
#include <TH1D.h>
#include <TH2D.h>       // istogrammi 1d
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TView.h>
#include <TTree.h>      // gestire TTree (strutture dati ROOT)
#include <iostream>
#include <fstream>      // leggere / scrivere file di testo
#include <sstream>      // gestire stringhe
#include <vector>
#include <algorithm>    // necessario per std::sort (ordinare le strip colpite nei cluster)
#include <math.h>

// --- CARTELLA FILE
std::string basePath = "/Users/martina/Desktop/LabIntFond/3_glast/acquisizioni_prova/21-04-30-all-masked/"; // ****** MODIFICARE
std::string fileName = "TkrDataTaking_333001114.lif";

// --- DEFINIZIONE LUNGHEZZE E COSE DETECTOR [mm]
int EdgeWidth = 1;
double StripPitch = .228;
double LadderSeparation = .220;
std::vector<std::string> layerNames = {"Y0", "X0", "X1", "Y1", "Y2", "X2", "X3", "Y3", "Y4", "X4"}; // vettore di stringhe per i layer
double z_coords[10] = {0, 25, 75, 100, 150, 175, 225, 250, 300, 325}; // in mm

// --- FUNZIONI PER CALCOLARE LA COORDINATE x,y e z IN mm
// Coordinate = EdgeWidth + StripPitch*StripNumber + (LadderSeparation +2*EdgeWidth -StripPitch)*int(StripNumber/384)
double getCoordinate(double strip) {
    double coord = EdgeWidth + StripPitch * strip +
                   (LadderSeparation + 2 * EdgeWidth - StripPitch) * (int(strip / 384));
    return coord;
}

double getZ(int layerIdx){
    return z_coords[layerIdx];
}

// DISEGNINI
void DrawPlane(double z) {
    double L = 430;
    TPolyLine3D *plane = new TPolyLine3D(5);

    plane->SetPoint(0, 0, 0, z);
    plane->SetPoint(1, L, 0, z);
    plane->SetPoint(2, L, L, z);
    plane->SetPoint(3, 0,  L, z);
    plane->SetPoint(4, 0, 0, z); // chiusura

    plane->SetLineColor(kGray+2);
    plane->SetLineWidth(1);
    plane->Draw();
}
void DrawStrip(double pos, int layerIdx, bool isX) {
    double z = getZ(layerIdx);
    double L = 430; // lunghezza strip

    TPolyLine3D *strip = new TPolyLine3D(2);

    if (isX) {
        // misura x -> strip lungo y
        strip->SetPoint(0, pos, 0, z);
        strip->SetPoint(1, pos, L, z);
    } else {
        // misura y -> strip lungo x
        strip->SetPoint(0, 0, pos, z);
        strip->SetPoint(1, L, pos, z);
    }

    strip->SetLineColor(kRed);
    strip->SetLineWidth(2);
    strip->Draw();
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

// --- FUNZIONE PER DISEGNO 3D ---
void DrawEvent3D(Event& ev, double mx, double qx, double my, double qy, int eventNum) {
    std::cout << "Drawing event " << ev.evtID << std::endl;
    TCanvas *c3d = new TCanvas(Form("c3d_%d", eventNum), Form("evento 3D ID %d", ev.evtID), 600, 600); // Form crea le stringhe
    TView *view = TView::CreateView(1); // crea view 3D (1 dice in prospettiva)
    view->SetRange(-100, -100, -10, 500, 500, 350); // range geometrico (x,y,z)

    // Disegna i punti (Hits)
/*    TPolyMarker3D *pm = new TPolyMarker3D();
    for(Hit &hx : ev.x_hits){ // loop sulle hit sulle X (misura)
        // un punto 3D richiede x,y,z; uso la retta Y(z) per stimare la Y del punto X
        pm->SetNextPoint(hx.pos, my*getZ(hx.layerIdx) + qy, getZ(hx.layerIdx)); // aggiunge un altro punto
    }
    for(Hit &hy : ev.y_hits) { // analogo a prima
        pm->SetNextPoint(mx*getZ(hy.layerIdx) + qx, hy.pos, getZ(hy.layerIdx));
    }
    pm->SetMarkerStyle(20);
    pm->SetMarkerColor(kRed);
    pm->Draw();
*/
    // Disegna le strip
    // strip x
    for (Hit &hx : ev.x_hits) {
        DrawStrip(hx.pos, hx.layerIdx, true);
    }

    //  strip y
    for (Hit &hy : ev.y_hits) {
        DrawStrip(hy.pos, hy.layerIdx, false);
    }
    
    // Disegna la traccia (Retta fittata)
    TPolyLine3D *line = new TPolyLine3D();
    double zMin = -5;
    double zMax = 330;
    line->SetPoint(0, mx*zMin + qx, my*zMin + qy, zMin); // inizio linea
    line->SetPoint(1, mx*zMax + qx, my*zMax + qy, zMax); // fine linea
    line->SetLineColor(kBlue);
    line->SetLineWidth(2);
    line->Draw();
    
    // Disegna i layer
    for(int i = 0; i < 10; i++){
        DrawPlane(getZ(i));
    }
    c3d->Update();
}

// --- FUNZIONE PER RICOSTRUIRE LE TRACCE
void CosmicReconstruction(){
    // apro eventi
    std::vector<Event> events = readFile(basePath + fileName);
    // istogrammi
    TH1D* hAngles = new TH1D("hAngles", "Distribuzione Angolare;Theta [deg];Eventi", 18, 0, 90);
    int drawnCount = 0;
    // ciclo sugli eventi
    for(Event &ev : events){
        // scarto se ho meno di 2 x o di 2 y
        if(ev.x_hits.size() < 2 || ev.y_hits.size() < 2) continue;
        
        // controllo coppie attive
        bool hitX[5] = {false};
        bool hitY[5] = {false};
        
        for(const Hit &h : ev.x_hits){
            int hit = layerNames[h.layerIdx][1] - '0'; // converte una stringa in un numero, tipo '3' in int 3
            hitX[hit] = true;
        }
        for(const Hit &h : ev.y_hits){
            int hit = layerNames[h.layerIdx][1] - '0'; // converte una stringa in un numero, tipo '3' in int 3
            hitY[hit] = true;
        }
        int nCouples = 0; // numero di coppie
        for(int i = 0; i < 5; i++){
            if(hitX[i] && hitY[i]) nCouples++;
        }
        if(nCouples > 2) continue;
        
        // FIT
        // fit su x
        TGraphErrors grX;
        for(Hit &h : ev.x_hits){
            grX.SetPoint(grX.GetN(), getZ(h.layerIdx), h.pos);}
        grX.Fit("pol1", "Q"); // fa il fit lineare
        TF1* fitX = grX.GetFunction("pol1"); // ripiglia la funzione di best fit
        // fit su y
        TGraphErrors grY;
        for(Hit &h : ev.y_hits){
            grY.SetPoint(grY.GetN(), getZ(h.layerIdx), h.pos);}
        grY.Fit("pol1", "Q");
        TF1* fitY = grY.GetFunction("pol1");
        // parametri
        double mx = fitX->GetParameter(1);
        double qx = fitX->GetParameter(0);
        double my = fitY->GetParameter(1);
        double qy = fitY->GetParameter(0);
        double cosTheta = 1.0 / sqrt(mx*mx + my*my + 1.0); // roba geometrica
        double thetaDeg = acos(cosTheta) * 180.0 / M_PI;
        hAngles->Fill(thetaDeg);
        // DISEGNA LE PRIME 5 TRACCE
        if (drawnCount < 2) {
            DrawEvent3D(ev, mx, qx, my, qy, drawnCount);
            drawnCount++;
        }
    }
        TCanvas* c1 = new TCanvas("c1", "Analisi Cosmici", 800, 600);
        hAngles->SetFillColor(kCyan);
        hAngles->Draw();

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
                 55, 5, 60)
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

