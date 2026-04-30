// ricostruzione tracce

#include "_detector_utils.h"  // serve per le robe sul detector
#include "_readfile.cc"  // serve per le funzioni

#include <TH1D.h>
#include <TH2D.h>       // istogrammi 1d
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TView.h>

// DEFINIZIONE STRUTTURE
struct Track {
    double mx, qx;   // parametri fit vista X
    double my, qy;   // parametri fit vista Y
    double chi2x, chi2y;
    double theta;
};

// DEFINIZIONE VETTORE DI EVENTS (usato praticamente in tutte le funzioni)
// apro eventi
std::vector<Event> events = readFile(basePath + fileName);

// DISEGNINI
void DrawPlane(double z) {
    TPolyLine3D *plane = new TPolyLine3D(5);

    plane->SetPoint(0, 0, 0, z);
    plane->SetPoint(1, StripLength, 0, z);
    plane->SetPoint(2, StripLength, StripLength, z);
    plane->SetPoint(3, 0,  StripLength, z);
    plane->SetPoint(4, 0, 0, z); // chiusura

    plane->SetLineColor(kGray);
    plane->SetLineWidth(1);
    plane->Draw();
}

void DrawStrip(double pos, int layerIdx, bool isX) {
    double z = getZ(layerIdx);

    TPolyLine3D *strip = new TPolyLine3D(2);

    if (isX) {
        // misura x -> strip lungo y
        strip->SetPoint(0, pos, 0, z);
        strip->SetPoint(1, pos, StripLength, z);
    } else {
        // misura y -> strip lungo x
        strip->SetPoint(0, 0, pos, z);
        strip->SetPoint(1, StripLength, pos, z);
    }

    strip->SetLineColor(kRed);
    strip->SetLineWidth(2);
    strip->Draw();
}

// --- FUNZIONI
// FUNZIONE PER DISEGNO 3D
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

// FUNZIONE PER FITTARE LE TRACCE -> esce una track
Track fitTrack(std::vector<Hit> x_hits, std::vector<Hit> y_hits){
    Track track;
    // fa il fit e restituisce i parametri
    TGraphErrors grX, grY;
    // fit su X
    for(Hit &h : x_hits){
        grX.SetPoint(grX.GetN(), getZ(h.layerIdx), h.pos);}
    grX.Fit("pol1", "Q"); // fa il fit lineare
    TF1* fitX = grX.GetFunction("pol1"); // ripiglia la funzione di best fit
    // fit su y
    for(Hit &h : y_hits){
        grY.SetPoint(grY.GetN(), getZ(h.layerIdx), h.pos);}
    grY.Fit("pol1", "Q");
    TF1* fitY = grY.GetFunction("pol1");
    // parametri (li salvo nella traccia)
    track.mx = fitX->GetParameter(1);
    track.qx = fitX->GetParameter(0);
    track.my = fitY->GetParameter(1);
    track.qy = fitY->GetParameter(0);
    track.chi2x = ( fitX->GetChisquare() / fitX->GetNDF() );
    track.chi2y = ( fitY->GetChisquare() / fitY->GetNDF() );
    
    double cosTheta = 1.0 / sqrt(track.mx*track.mx + track.my*track.my + 1.0);
    track.theta = acos(cosTheta) * TMath::RadToDeg();
    
    return track;
}


// FUNZIONE PER RICOSTRUIRE LE TRACCE
void CosmicReconstruction(){
    // istogrammi
    TH1D* hAngles = new TH1D("hAngles", "Distribuzione Angolare;Theta [deg];Eventi", 45, 0, 90);
    int drawnCount = 0;
    // vettore di tracce
    std::vector<Track> tracks;
    // per controllare quanti passano la selezione rispetto a quelli totali
    int nTotal = 0;
    int nSelected = 0;
    // ciclo sugli eventi
    for(Event &ev : events){
        nTotal++;
        // scarto se ho meno di 3 x o di 3 y
        if(ev.x_hits.size() < 3 || ev.y_hits.size() < 3) continue;
        // scarto se ho più hit per layer
        bool multiHit = false;
        bool singleCluster = false;
        for(int i = 0; i < 10; i++){
            if(ev.clustersPerLayer[i].size()>1) multiHit = true;
        }
        if(multiHit) continue;
        for(int i = 0; i < 10; i++){
            for(const Cluster &cl : ev.clustersPerLayer[i]){
                if(cl.mult == 1){singleCluster = true;}
            }
        }
      //  if(singleCluster) continue; // scarta i cluster di molteplicità singola
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
   /*    for(int i = 0; i < 5; i++){
            if(hitX[i] && hitY[i]) nCouples++;
        }
        if(nCouples > 2) continue;
    */
        nSelected++;
        // chiamo la funzione di creazione tracce
        Track track = fitTrack(ev.x_hits, ev.y_hits);
        tracks.push_back(track);
        
        hAngles->Fill( track.theta );
        // DISEGNA LE PRIME 5 TRACCE
        if (drawnCount < 5) {
            DrawEvent3D(ev, track.mx, track.qx, track.my, track.qy, drawnCount);
            std::cout << "ch2X = " << track.chi2x << "\n";
            std::cout << "ch2Y = " << track.chi2y << "\n";
            drawnCount++;
        }
    }
    TCanvas* c1 = new TCanvas("c1", "Analisi Cosmici", 800, 600);
    hAngles->SetFillColor(kCyan);
    hAngles->Draw();
    std::cout << "selezionati / totali = " << nSelected << " / " << nTotal << " = " << (double)nSelected/nTotal<< "\n";

}

void LayerEfficiency(){
    // definizione layer di trigger. dentro alle {,}, il primo è il layer di cui voglio studiare l'efficienza e il secondo è uno dei due layer di trigger. si accede con trigger1["layer di cui voglio studiare l'efficienza"] e l'analogo col trigger2
    std::map<std::string, std::string> trigger1 = {
        {"X0", "Y4"}, {"X1", "Y0"}, {"X2", "Y2"}, {"X3", "Y2"}, {"X4", "Y0"},
        {"Y0", "X0"}, {"Y1", "X1"}, {"Y2", "X1"}, {"Y3", "X3"}, {"Y4", "X3"}
    };
    std::map<std::string, std::string> trigger2 = {
        {"X0", "Y0"}, {"X1", "Y1"}, {"X2", "Y3"}, {"X3", "Y3"}, {"X4", "Y4"},
        {"Y0", "X4"}, {"Y1", "X2"}, {"Y2", "X2"}, {"Y3", "X4"}, {"Y4", "X4"}
    };
    std::map<std::string, std::string> trigger3 = { // solo per layer estremali
        {"X4", "Y2"},
        {"Y0", "X2"}
    };
   
    // loop sui layer da testare
    for(int l = 0; l < 10; l++){
        std::string layerName = layerNames[l];
        int nTotal = 0;   // eventi con trigger attivo, inizializzati a 0
        int nHit = 0;     // eventi con trigger attivo e hit sul layer testato, inizializzati a 0
        
        // per ogni evento
        for(Event &ev : events){
            // controlla se i layer di trigger sono colpiti
            bool hit1 = ev.clustersPerLayer[ layerIndex[trigger1[layerName]] ].size() > 0;
            bool hit2 = ev.clustersPerLayer[ layerIndex[trigger2[layerName]] ].size() > 0;
            bool hitLayer = ev.clustersPerLayer[ layerIndex[layerName] ].size() > 0;
            
            if(layerName != "Y0" || layerName != "X4"){
                if( hit1 && hit2 ){ // se sì, incrementa nTotal
                nTotal++;
                }
                // controlla se il layer testato è colpito compatibilmente con la traccia (STRIP MASCHERATE?????)
                if( hit1 && hit2 && hitLayer ){
                    // se sì, incrementa nHit
                    nHit++;
                }
            }
            else{ // per layer estremali aggiungo ulteriore layer di trigger
                bool hit3 = ev.clustersPerLayer[ layerIndex[trigger3[layerName]] ].size() > 0;
                if( hit1 && hit2 && hit3 ){ // se sì, incrementa nTotal
                nTotal++;
                }
                // controlla se il layer testato è colpito compatibilmente con la traccia (STRIP MASCHERATE?????)
                if( hit1 && hit2 && hit3 && hitLayer ){
                    // se sì, incrementa nHit
                    nHit++;
                }

            }
        }
        // efficienza = nHit / nTotal
        double eff = (double)nHit / nTotal;
        std::cout << "nTotal =" << nTotal << " nHit =" << nHit << "\n";
        std::cout << "--------";
        std::cout << "efficienza layer '" << layerNames[l] << "' = " << eff << "\n";
        std::cout << "--------";
        std::cout << "--------";

    }
}

void Alignment(){
    // inizializzo vettori di istogrammi
    std::vector<TH1D*> hResiduals;
    std::vector<TH2D*> hResvsPos;

    for(int l = 0; l < 10; l++){
        hResiduals.push_back(new TH1D(
            ("hRes_" + layerNames[l]).c_str(),
            (layerNames[l] + ";res. [mm];counts").c_str(),
            40, -20, 20));
        hResvsPos.push_back(new TH2D(
            ("hResvsPos_" + layerNames[l]).c_str(),
            (layerNames[l] + ";pos [mm];res. [mm]").c_str(),
            400, 0, 400, 50, -50, 50));
        
    }
    // inizializzo vettore di correzioni, una correzione per ogni layer
    double corrections[10] = {0};
    // itera finché converge
       // for(int i = 0; i < 10; i++){
            for(int l = 0; l < 10; l++){ // per ogni layer, fit senza quel layer e calcola residuo medio
                // esclude i layer di riferimento dall'analisi
             //   if(l == layerIndex["X0"] || l == layerIndex["Y0"]) continue;
                double sumResiduals = 0;
                int nResiduals = 0;
                std::vector<Track> tracks;

                // loop sugli eventi
                for(Event &ev : events){
                    // prendi gli eventi di quando è colpito il layer l
                    if(ev.clustersPerLayer[l].size() == 0) continue;
                    // controlla che sia stato colpito il layer di riferimento per la stessa vista del layer l
              /*      if(layerNames[l][0] == 'X'){
                        if(ev.clustersPerLayer[ layerIndex["X0"] ].size() == 0) continue;
                    }
                    if(layerNames[l][0] == 'Y'){
                        if(ev.clustersPerLayer[ layerIndex["Y0"] ].size() == 0) continue;
                    }
              */
                    // fit senza il layer l
                    // costruisco x_hits e y_hits escludendo il layer l
                    std::vector<Hit> x_hits_wo_l, y_hits_wo_l;
                    for(const Hit &h : ev.x_hits){
                        if(h.layerIdx != l) x_hits_wo_l.push_back(h);
                    }
                    for(const Hit &h : ev.y_hits){
                        if(h.layerIdx != l) y_hits_wo_l.push_back(h);
                    }
                    // escludo tracce con poche hit
                    if(x_hits_wo_l.size() < 4 || y_hits_wo_l.size() < 4) continue;
                    // scarto se ho più hit per layer
                    bool multiHit = false;
                    for(int i = 0; i < 10; i++){
                        if(ev.clustersPerLayer[i].size()>1) multiHit = true;
                    }
                    if(multiHit) continue;
                    
                    Track track = fitTrack(x_hits_wo_l, y_hits_wo_l);
                    tracks.push_back(track);
                    
                    // calcolo del residuo: ( posizione prevista - posizione misurata ) in mm
                    if(layerNames[l][0] == 'X'){
                        // posizione prevista dalla traccia
                        double expectedX = track.mx * getZ(l) + track.qx;
                        // posizione misurata
                        for(const Hit &h : ev.x_hits){
                            if(h.layerIdx == l){
                                double measuredX = h.pos;
                                double errX = h.err;
                                double res = (measuredX - expectedX);
                                double norm_res = (measuredX - expectedX)/errX;
                                sumResiduals += measuredX - expectedX;  // non normalizzato, per la correzione
                                hResiduals[l]->Fill(res);
                                hResvsPos[l]->Fill(measuredX, res);
                                nResiduals++;
                                break;
                            }
                        }
                    }
                    if(layerNames[l][0] == 'Y'){
                        // posizione prevista dalla traccia
                        double expectedY = track.my * getZ(l) + track.qy;
                        // posizione misurata
                        for(const Hit &h : ev.y_hits){
                            if(h.layerIdx == l){
                                double measuredY = h.pos;
                                double errY = h.err;
                                double res = (measuredY - expectedY);
                                double norm_res = (measuredY - expectedY)/errY;
                                sumResiduals += measuredY - expectedY;  // non normalizzato, per la correzione
                                hResiduals[l]->Fill(res);
                                hResvsPos[l]->Fill(measuredY, res);
                                nResiduals++;
                                break;
                            }
                        }
                    }
                }
                //corrections[l] = sumResiduals / nResiduals;
            }
    //    }
    // istogramma residui
    TCanvas *c1 = new TCanvas("cRes", "Residui", 1200, 1000);
    c1->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c1->cd(l+1);
        hResiduals[l]->Draw("HIST");
    }
    
    TCanvas *c2 = new TCanvas("cResvsPos", "Residui vs posizione", 1200, 1000);
    c2->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c2->cd(l+1);
        hResvsPos[l]->Draw("COLZ");
    }
}

