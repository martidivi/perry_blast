// ricostruzione tracce
// --- HOUGH
// per ogni evento (-> loop sugli eventi) faccio un istogramma con tante curve (theta,rho) quanti sono le hit (un istogramma per ciascuna vista). per fare l'istogramma mi serve un vettore di theta (generati unif. in [0,pi]) e rho (calcolati in base a ciascuna hit).
// i punti di intersezione delle curve identificano i parametri (theta,rho) che parametrizzano le best tracce
// analisi separata per le due viste
// per ogni hit (su una vista!!), costruisco un vettore di theta e un vettore di rho.
// alla fine dell'analisi per hit, i vettori di theta e di rho vanno incapsulati in due nuovi vettoroni con cui fillare l'istogramma


#include "_detector_utils.h"  // serve per le robe sul detector
#include "_readfile.cc"  // serve per le funzioni

#include <TH1D.h>
#include <TH2D.h>       // istogrammi 1d
#include "TSpectrum2.h"
#include <TF1.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TView.h>
#include "TRandom3.h" // generazione casuale x hough. ma non serve, è molto lento

std::string views[2] = {"X","Y"};

// DEFINIZIONE STRUTTURE
struct Track {
    double mx, qx;   // parametri fit vista X
    double mx_err, qx_err;
    double my, qy;   // parametri fit vista Y
    double my_err, qy_err;
    double chi2x, chi2y;
    double theta;
    double theta_err;

};

// tracce su X e su Y per un evento (seppur raramente, posso avere un numero diverso di tracce su X e su Y). le tracks sono comunque pensate per essere 3D, quindi per tracksX avrò tutti i parametri Y a 0 e viceversa.
struct EventTracks {
    std::vector<Track> tracksX;
    std::vector<Track> tracksY;
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

    std::vector<TGraph*> g(2);
    g[0] = new TGraph();
    g[1] = new TGraph();

    int nx = 0, ny = 0;

    // X hits -> X vs Z
    for (auto &hx : ev.x_hits) {
        double z = getZ(hx.layerIdx);
        g[0]->SetPoint(nx, hx.pos, z);
        nx++;
    }

    // Y hits -> Y vs Z
    for (auto &hy : ev.y_hits) {
        double z = getZ(hy.layerIdx);
        g[1]->SetPoint(ny, hy.pos, z);
        ny++;
    }

    TCanvas *c3d = new TCanvas(Form("c3d_%d", eventNum),
                               Form("evento 3D ID %d", ev.evtID),
                               1600, 900);

    c3d->Divide(3,1);

    // PAD 1 → 3D
    c3d->cd(1);

    TView *view = TView::CreateView();
    view->SetRange(-100, -100, -10, 500, 500, 350);

    for (Hit &hx : ev.x_hits) DrawStrip(hx.pos, hx.layerIdx, true);
    for (Hit &hy : ev.y_hits) DrawStrip(hy.pos, hy.layerIdx, false);

    TPolyLine3D *line = new TPolyLine3D();

    double zMin = -10;
    double zMax = 310;

    line->SetPoint(0, mx*zMin + qx, my*zMin + qy, zMin);
    line->SetPoint(1, mx*zMax + qx, my*zMax + qy, zMax);

    line->SetLineColor(kBlue);
    line->SetLineWidth(2);
    line->Draw();

    for(int i = 0; i < 10; i++) {
        DrawPlane(getZ(i));
    }

    // PAD 2 → X vs Z
    c3d->cd(2);
    g[0]->SetTitle("X vs Z;X[mm];Z[mm]");
    g[0]->GetYaxis()->SetRangeUser(-10, 310);   // asse Z
    g[0]->GetXaxis()->SetLimits(-10, 360); // asse X
    g[0]->SetMarkerStyle(21);
    g[0]->SetMarkerColor(2);
    g[0]->Draw("AP");

    // PAD 3 → Y vs Z
    c3d->cd(3);
    g[1]->SetTitle("Y vs Z;Y[mm];Z[mm]");
    g[1]->GetYaxis()->SetRangeUser(-10, 310);   // asse Z
    g[1]->GetXaxis()->SetLimits(-10, 360); // asse X
    g[1]->SetMarkerStyle(21);
    g[1]->SetMarkerColor(2);
    g[1]->Draw("AP");

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
    if(!fitX) return Track();  // fit fallito, ritorna traccia vuota
    // fit su y
    for(Hit &h : y_hits){
        grY.SetPoint(grY.GetN(), getZ(h.layerIdx), h.pos);}
    grY.Fit("pol1", "Q");
    TF1* fitY = grY.GetFunction("pol1");
    if(!fitY) return Track();  // fit fallito, ritorna traccia vuota
    // parametri (li salvo nella traccia)
    track.mx = fitX->GetParameter(1);
    track.mx_err = fitX->GetParError(1);
    track.qx = fitX->GetParameter(0);
    track.qx_err = fitX->GetParError(0);
    track.my = fitY->GetParameter(1);
    track.my_err = fitY->GetParError(1);
    track.qy = fitY->GetParameter(0);
    track.qy_err = fitY->GetParError(0);
    track.chi2x = ( fitX->GetChisquare() / fitX->GetNDF() );
    track.chi2y = ( fitY->GetChisquare() / fitY->GetNDF() );
    
    double cosTheta = 1.0 / sqrt(track.mx*track.mx + track.my*track.my + 1.0);
    track.theta = acos(cosTheta) * TMath::RadToDeg();
    
    return track;
}

/*std::vector<Track> fitHoughTransform(std::vector<Hit> x_hits, std::vector<Hit> y_hits, int n){
    // ritorna un vettore di tracce, la useremo per ogni evento -> se il vettore ha più elementi, significa che ho avuto tracce multiple
    std::vector<Track> tracks;
    // fa il fit e restituisce i parametri
    TGraphErrors grX, grY;
    // fit su X
    for(Hit &h : x_hits){
        for(int i = 0; i < n; i++){
            // genera un theta (parametro, chiamo tHeta per distinguerlo dall'angolo polare) in [0,pi]
            TRandom3 *r = new TRandom3();
            double tHeta = r->Uniform(0, M_PI);
            double rho = getZ(h.layerIdx)*cos(tHeta) + h.pos*sin(tHeta);
            
        }
    }
}
*/


// FUNZIONE PER RICOSTRUIRE LE TRACCE
void CosmicReconstruction(){
    // istogrammi
    TH1D* hAngles = new TH1D("hAngles", "Distribuzione Angolare;Theta [deg];Eventi", 45, 0, 90);
    int drawnCount = 0;
    // vettore di tracce
    std::vector<Track> tracks;
    // per controllare quanti passano la selezione rispetto a quelli totali
    int nTotal = 0, nSelected = 0;
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

// FUNZIONE PER TRACCE MULTIPLE CON HOUGH
void fitTrackHough(){
    // ------
    // mi faccio ritornare un vettorone mega gigante
 //   std::vector<EventTracks> allEventsTracks;
    // ------

    // dichiaro istogrammi
    std::vector<TH1D*> hClusterMax;
    TH1D *hClusterMaxX = new TH1D("hClusterMaxX","Max bin content per cluster X;max value;entries",10,0,10);
    TH1D *hClusterMaxY = new TH1D("hClusterMaxY","Max bin content per cluster Y;max value;entries",10,0,10);
    std::vector<TH2D*> hHough;
    for(int i = 0; i < 2; i++){
        hHough.push_back(new TH2D(Form("hHough%s", views[i].c_str()), Form("Spazio dei parametri (theta%s, rho%s);theta [rad];rho", views[i].c_str(), views[i].c_str()),250,0,M_PI, 500,-500,500));
        hClusterMax.push_back(new TH1D(Form("hClusterMax%s", views[i].c_str()),Form("Valori dei massimi %s;max values;counts", views[i].c_str()),10,0,10));
    }
    gStyle->SetOptStat("e"); // opzione per i grafici, mostra solo le entries invece della sbrodolata
    // contatori per disegnare, controllare quanti eventi passano la selezione rispetto a quelli totali
    int drawnCount = 0, nTotal = 0, nSelected = 0;
    // CICLO SUGLI EVENTI
    for(Event &ev : events){
        // incremento il numero totale di eventi
        nTotal++;
        // SELEZIONI: scarto se ho meno di tot x o di tot y colpite, scarto se ho meno di tot hit per ogni layer colpito (cerco tracce multiple), scarto se ho più di tot hit totali
      //  if(ev.x_hits.size() < 6 || ev.y_hits.size() < 6) continue;
        int multiHit = 0;
        for(int i = 0; i < 10; i++){
            if(ev.clustersPerLayer[i].size()>1) multiHit++;
        }
        if(multiHit<4) continue;
        if(ev.totalHits>50) continue;
    
       if(ev.x_hits.size() < 3 || ev.y_hits.size() < 3) continue;
        // scarto se ho più hit per layer
     /*   bool multiHit = false;
        bool singleCluster = false;
        for(int i = 0; i < 10; i++){
            if(ev.clustersPerLayer[i].size()>1) multiHit = true;
        }
        if(multiHit) continue;
        for(int i = 0; i < 10; i++){
            for(const Cluster &cl : ev.clustersPerLayer[i]){
                if(cl.mult == 1){singleCluster = true;}
            }
        }*/
        // incremento il numero di eventi che passano la selezione
        nSelected++;
        // ----
        // creo vettore vuoto di tracce
        EventTracks evTracks;
        // ----
        // per tentare local maximum scan
        double threshold = 3; // soglia di quanto è popolato un bin per definirlo "massimo" (cioè quanti punti richiedo per una traccia)
        std::vector<std::pair<int,int>> peaks[2]; // uno per X, uno per Y
        std::vector<std::vector<std::pair<int,int>>> clusters[2]; // uno per X, uno per Y
        
        for(int i = 0; i < 2; i++){ // loop per X e per Y
            // resetto istogrammi dopo ogni evento nel ciclo
            hHough[i]->Reset();
            // fillo gli istogrammi
            auto& hits = (i == 0) ? ev.x_hits : ev.y_hits;
            for(Hit &h : hits){
                for(int ib = 1; ib <= hHough[i]->GetNbinsX(); ib++){
                    double theta = hHough[i]->GetXaxis()->GetBinCenter(ib);
                    double rho = getZ(h.layerIdx)*cos(theta) + h.pos*sin(theta);
                    hHough[i]->Fill(theta, rho);
                }
            } // qui sono ancora dentro all'evento e ho creato istogrammi con theta e rho
        // nel loop parto da 2 perché poi controllo il bin 2-1, quindi se partissi da 1 andrei a 0 che però come bin non esiste
            for(int ix = 2; ix < hHough[i]->GetNbinsX(); ix++){ // faccio il loop sui bin x dell'istogramma
                for(int iy = 2; iy < hHough[i]->GetNbinsY(); iy++){ // faccio il loop sui bin y dell'istogramma
                    // trovo il valore al centro del bin
                    double center = hHough[i]->GetBinContent(ix, iy);
                    // controllo che sia sopra soglia
                    if(center < threshold) continue;
                    // se sta sopra la soglia, è un massimo locale eheh
                    bool isMaximum = true;
                    // controllo vicini 3x3
                    for(int dx = -1; dx <= 1; dx++){
                        for(int dy = -1; dy <= 1; dy++){
                            // salto il centro
                            if(dx == 0 && dy == 0) continue;
                            double neighbour = hHough[i]->GetBinContent(ix + dx, iy + dy);
                            // se un vicino è maggiore, non è massimo locale
                            if(neighbour > center){
                                isMaximum = false;
                                break;
                            }
                        }
                        if(!isMaximum) break;
                    }
                    // trovato massimo locale
                    if(isMaximum){
                        double theta = hHough[i]->GetXaxis()->GetBinCenter(ix);
                        double rho = hHough[i]->GetYaxis()->GetBinCenter(iy);
                        peaks[i].push_back({ix, iy});
                    }
                }
            } // esco dal loop per la ricerca di massimi locali
            // definisco un vettore di booleani, lungo quanto il vettore dei picchi
            std::vector<bool> used(peaks[i].size(), false);
            // loop sul vettore peaks[i] per identificare come un unico peak un cluster di più peak vicini
            for(int j = 0; j < peaks[i].size(); j++){
                // se è già true, continua
                if(used[j]) continue;

                used[j] = true;
                // lo mette nel cluster
                std::vector<std::pair<int,int>> cluster;
                cluster.push_back(peaks[i][j]);
                // loop su tutti i peak successivi per evitare duplicazioni
                for(int k = j+1; k < peaks[i].size(); k++){
                    if(used[k]) continue;
                    // first = indica theta, second = indica rho
                    int dx = abs(peaks[i][j].first - peaks[i][k].first);
                    int dy = abs(peaks[i][j].second - peaks[i][k].second);
                    // se stanno sufficientemente vicini, allora sono la stessa traccia
                    if(dx <= 2 && dy <= 2){
                        used[k] = true;
                        cluster.push_back(peaks[i][k]);
                    }
                }
                // infine fa push back nel vettorone di clustersX o clustersY
                clusters[i].push_back(cluster);
                // loop sui cluster
                double maxVal = 0;
                double sumW = 0, sumTheta = 0, sumRho = 0;

                for(auto &p : cluster){
                    double w = hHough[i]->GetBinContent(p.first, p.second);
                    // per fillare poi l'istogramma dei massimi
                    if(w > maxVal) maxVal = w;
                    sumTheta += w * hHough[i]->GetXaxis()->GetBinCenter(p.first);
                    sumRho += w * hHough[i]->GetYaxis()->GetBinCenter(p.second);
                    sumW += w;
                }
                double theta_peak = sumTheta / sumW;
                double rho_peak = sumRho / sumW;
                double m = -cos(theta_peak) / sin(theta_peak);
                double q = rho_peak / sin(theta_peak);
            /* -------------
                if(i == 0){ // traccia X
                    Track track;
                    track.mx = m;
                    track.qx = q;
                    evTracks.tracksX.push_back(track);
                }
                else{ // traccia Y
                    Track track;
                    track.my = m;
                    track.qy = q;
                    evTracks.tracksY.push_back(track);
                }
                */

                hClusterMax[i]->Fill(maxVal);
            }
        }
        // aggiunge nel vettorone finale il vettore di tracce per l'evento
        // ------------- allEventsTracks.push_back(evTracks);
        // std::cout << "numero di tracce ricostruite (su X e su Y) = " << allEventsTracks.size() << std::endl; // mi f
        // DISEGNA LE PRIME tot TRACCE, ancora dentro al loop degli eventi
       if (drawnCount < 5) {
           DrawEvent3D(ev, 0, 0, 0, 0, drawnCount);
           std::cout << " x_hits=" << ev.x_hits.size()
                     << " y_hits=" << ev.y_hits.size() << "\n"
                     << " hHoughX entries=" << hHough[0]->GetEntries() << "\n"
                     << " hHoughY entries=" << hHough[1]->GetEntries() << "\n";

           TCanvas *c = new TCanvas(Form("c_%d_hits%d", ev.evtID, ev.totalHits), Form("Trasformata di Hough per evento %d, hits%d", ev.evtID, ev.totalHits), 1600, 900); // Form crea le stringhe

           c->Divide(2,1);
           for(int i = 0; i < 2; i++){
               c->cd(i+1);
               hHough[i]->Draw("COLZ");
           }
           c->Update();
           drawnCount++;
       }
    }
    TCanvas *c2 = new TCanvas("c2","Max values",800,400);
    c2->Divide(2,1);
    for(int i = 0; i < 2; i++){
        c2->cd(i+1);
        hClusterMax[i]->Draw("HIST");
    }
    c2->Update();
    std::cout << "selezionati / totali = " << nSelected << " / " << nTotal << " = " << (double)nSelected/nTotal<< "\n";
    // -------
   // return allEventsTracks;
    // ----
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
    std::vector<TH2D*> hResvsPos_rot;
    std::vector<TH2D*> hResvsPos_theta;
    std::vector<TGraph*> gConvTrasl, gConvDeltaZ, gConvAlpha;

    for(int l = 0; l < 10; l++){
        hResiduals.push_back(new TH1D(
            ("hRes_" + layerNames[l]).c_str(),
            (layerNames[l] + ";res. [mm];counts").c_str(),
            40, -20, 20));
        
        hResvsPos.push_back(new TH2D(
            ("hResvsPos_" + layerNames[l]).c_str(),
            (layerNames[l] + ";pos [mm];res. [mm]").c_str(),
            400, 0, 400, 60, -30, 30));
        
        hResvsPos_rot.push_back(new TH2D(
            ("hResvsPos_rot_" + layerNames[l]).c_str(),
            (layerNames[l] + ";pos [mm];res. [mm]").c_str(),
            400, 0, 400, 60, -30, 30));
        
        hResvsPos_theta.push_back(new TH2D(
            ("hResvsPos_theta_" + layerNames[l]).c_str(),
            (layerNames[l] + ";theta [deg];res. [mm]").c_str(),
             35, 0, 70, 60, -30, 30));
        
        gConvTrasl.push_back(new TGraph());
        gConvDeltaZ.push_back(new TGraph());
        gConvAlpha.push_back(new TGraph());
        gConvTrasl[l]->SetTitle((layerNames[l] + ";iter;trasl [mm]").c_str());
        gConvDeltaZ[l]->SetTitle((layerNames[l] + ";iter;deltaZ [mm]").c_str());
        gConvAlpha[l]->SetTitle((layerNames[l] + ";iter;alpha [deg]").c_str());
            }
    // inizializzo vettori di correzioni, una correzione per ogni layer. più gli errori associati
    double trasl[10] = {0};     // traslazione X/Y dei layer
    double deltaZ[10] = {0};     // correzione lungo z
    double alpha[10]  = {0};     // rotazione rispetto all'asse z
    double trasl_err[10] = {0};     // traslazione X/Y dei layer
    double deltaZ_err[10] = {0};     // correzione lungo z
    double alpha_err[10]  = {0};     // rotazione rispetto all'asse z
    
    int iter = 20; // numero di iterazioni

    for(int i = 0; i < iter; i++){
        std::cout << "\n iterazione " << i << "\n";
            for(int l = 0; l < 10; l++){ // per ogni layer, fit senza quel layer e calcola residuo medio. poi: per ogni layer, fit sulla vista opposta, calcolo residuo e plotto in funzione della vista opposta (alla z del layer)
                // reset istogrammi
                hResiduals[l]->Reset();
                hResvsPos[l]->Reset();
                hResvsPos_rot[l]->Reset();
                hResvsPos_theta[l]->Reset();
                
                double sumResiduals = 0;
                int nResiduals = 0;

                // coefficiente angolare
                std::vector<double> slopes;
                std::vector<double> slopes_err;
                // residui
                std::vector<double> residuals;
                std::vector<double> residuals_err;
                
                std::vector<double> otherView;
                std::vector<double> otherView_err;
                
                std::vector<Track> tracks;
                // layer di riferimento
               // if(l == layerIndex["X4"] || l == layerIndex["Y0"]) continue;

                // loop sugli eventi
                for(Event &ev : events){
                    // prendi gli eventi di quando è colpito il layer l
                    if(ev.clustersPerLayer[l].size() == 0) continue;
                    // escludo tracce con poche hit (se toglie <5 significa che sta richiedendo TUTTI i layer accesi)
                    if(ev.x_hits.size() < 5 || ev.y_hits.size() < 5) continue;
                    // fit senza il layer l
                    // costruisco x_hits e y_hits escludendo il layer l
                    std::vector<Hit> x_hits_wo_l, y_hits_wo_l;

                    for(const Hit &h : ev.x_hits){
                        if(h.layerIdx != l){
                            x_hits_wo_l.push_back(h);
                        }

                    }
                    for(const Hit &h : ev.y_hits){
                        if(h.layerIdx != l){
                            y_hits_wo_l.push_back(h);
                        }
                    }
                    // scarto se ho più hit per layer
                    bool multiHit = false;
                    for(int i = 0; i < 10; i++){
                        if(ev.clustersPerLayer[i].size()>1) multiHit = true;
                    }
                    if(multiHit) continue;
                    
                    Track track = fitTrack(x_hits_wo_l, y_hits_wo_l);
                    tracks.push_back(track);
                    
                    // z corretta
                    double zCorr = getZ(l) + deltaZ[l];
                    
                    // calcolo del residuo: ( posizione prevista - posizione misurata ) in mm
                    // calcolo posizione prevista dalla traccia
                    double expectedX = track.mx * zCorr + track.qx;
                    double expectedX_err = sqrt( pow(zCorr * track.mx_err,2) + pow(track.qx_err,2) );
                    double expectedY = track.my * zCorr + track.qy;
                    double expectedY_err = sqrt( pow(zCorr * track.my_err,2) + pow(track.qy_err,2) );
                    double theta = track.theta;
                    trasl_err[l] = sqrt(pow(trasl_err[l], 2) ); // + pow(sumResiduals_err/nResiduals, 2)
                    deltaZ_err[l] = sqrt(pow(deltaZ_err[l], 2)); //  + pow(fitDeltaZ->GetParError(1), 2)
                    alpha_err[l] = sqrt(pow(alpha_err[l], 2) ); // + pow(fitRot->GetParError(1), 2)
                    
                    if(layerNames[l][0] == 'X'){
                        // posizione misurata
                        for(const Hit &h : ev.x_hits){
                            if(h.layerIdx == l){
                                // correzione traslazione
                                double measuredX = h.pos - trasl[l];
                                // correzione rotazione
                                measuredX -= alpha[l] * expectedY;
                                double errX = h.err;
                                double res = measuredX - expectedX;
                                double res_err = sqrt( pow(errX,2) + pow(expectedX_err,2));
                                sumResiduals += res;  // non normalizzato, per la correzione
                                hResiduals[l]->Fill(res);
                                hResvsPos[l]->Fill(measuredX, res);
                                hResvsPos_rot[l]->Fill(expectedY,res);
                                hResvsPos_theta[l]->Fill(theta,res);
                                
                                slopes.push_back(track.mx);
                                slopes_err.push_back(track.mx_err);
                                residuals.push_back(res);
                                residuals_err.push_back(res_err);
                                otherView.push_back(expectedY);
                                otherView_err.push_back(expectedY_err);
                                nResiduals++;
                                break;
                            }
                        }
                    }
                    if(layerNames[l][0] == 'Y'){
                        // posizione misurata
                        for(const Hit &h : ev.y_hits){
                            if(h.layerIdx == l){
                                // correzione traslazione
                                double measuredY = h.pos - trasl[l];
                                // correzione rotazione
                                measuredY -= alpha[l] * expectedX;
                                double errY = h.err;
                                double res = measuredY - expectedY;
                                double res_err = sqrt( pow(errY,2) + pow(expectedY_err,2));
                                sumResiduals += res;  // non normalizzato, per la correzione
                                hResiduals[l]->Fill(res);
                                hResvsPos[l]->Fill(measuredY, res);
                                hResvsPos_rot[l]->Fill(expectedX,res);
                                hResvsPos_theta[l]->Fill(theta,res);
                                
                                slopes.push_back(track.my);
                                slopes_err.push_back(track.my_err);
                                residuals.push_back(res);
                                residuals_err.push_back(res_err);
                                otherView.push_back(expectedX);
                                otherView_err.push_back(expectedX_err);
                                nResiduals++;
                                break;
                            }
                        }
                    }
                    
                }
                TGraphErrors* gr = new TGraphErrors( residuals.size(), slopes.data(), residuals.data(), slopes_err.data(), residuals_err.data() );
                gr->SetTitle((layerNames[l]).c_str());
                gr->Fit("pol1", "Q");
                TF1* fitDeltaZ = gr->GetFunction("pol1");
                if(fitDeltaZ){
                    deltaZ[l] += fitDeltaZ->GetParameter(1);
                  //  deltaZ[layerIndex["X4"]] = 0;
                  //  deltaZ[layerIndex["Y0"]] = 0;
                    std::cout << "Layer " << layerNames[l] << " -> deltaZ = " << fitDeltaZ->GetParameter(1) << " +/- " << fitDeltaZ->GetParError(1) << " mm\n";
                }
                //grPosZ.push_back(gr);
                
                TGraphErrors* grRot = new TGraphErrors( residuals.size(), otherView.data(), residuals.data(), otherView_err.data(), residuals_err.data() );
                grRot->SetTitle((layerNames[l]).c_str());
                grRot->Fit("pol1", "Q");
                TF1* fitRot = grRot->GetFunction("pol1");
                if(fitRot){
                    alpha[l]  += fitRot->GetParameter(1);
                  //  alpha[layerIndex["X4"]] = 0;
                  //  alpha[layerIndex["Y0"]] = 0;

                    std::cout << "Layer " << layerNames[l] << " -> alpha = " << fitRot->GetParameter(1)* TMath::RadToDeg() << "° +/- " << fitRot->GetParError(1)* TMath::RadToDeg() << "°\n";
                }
                // grRotation.push_back(grRot);
                trasl[l] += sumResiduals / nResiduals;
               // trasl[layerIndex["X4"]] = 0;
              //  trasl[layerIndex["Y0"]] = 0;
                std::cout << "Layer " << layerNames[l] << " -> trasl = " << trasl[l] << " +/- " << "0" << " mm\n";
                
                gConvTrasl[l]->SetPoint(i, i, trasl[l]);
                gConvDeltaZ[l]->SetPoint(i, i, deltaZ[l]);
                gConvAlpha[l]->SetPoint(i, i, alpha[l]);
            }
        }
    // istogramma residui
    TCanvas *c1 = new TCanvas("cRes", "Residui", 1600, 900);
    c1->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c1->cd(l+1);
        hResiduals[l]->Draw("HIST");
    }
    
    TCanvas *c2 = new TCanvas("cResvsPos", "Residui vs posizione", 1600, 900);
    c2->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c2->cd(l+1);
        hResvsPos[l]->Draw("COLZ");
    }
    
    TCanvas *c3 = new TCanvas("cResvsPos_rot", "Residui vs posizione altra vista", 1600, 900);
    c3->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c3->cd(l+1);
        hResvsPos_rot[l]->Draw("COLZ");
    }
    
    TCanvas *c4 = new TCanvas("cResvsPos_theta", "Residui vs angolo polare", 1600, 900);
    c4->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c4->cd(l+1);
        hResvsPos_theta[l]->Draw("COLZ");
    }
    
/*    TCanvas *c5 = new TCanvas("cPosZ", "Residui vs tangente angolo polare", 1600, 900);
    c5->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c5->cd(l+1);
        grPosZ[l]->Draw("AP");
    }

    TCanvas *c6 = new TCanvas("cRot", "Residui vs posizione altra vista", 1600, 900);
    c6->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c6->cd(l+1);
        grRotation[l]->Draw("AP");
    }
 */
    TCanvas *c7 = new TCanvas("cTrasl", "Traslazioni", 1600, 900);
    c7->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c7->cd(l+1);
        gConvTrasl[l]->Draw("ALP");
    }
    
    TCanvas *c8 = new TCanvas("cAlpha", "Rotazioni", 1600, 900);
    c8->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c8->cd(l+1);
        gConvAlpha[l]->Draw("ALP");
    }
    
    TCanvas *c9 = new TCanvas("cDeltaZ", "Delta Z", 1600, 900);
    c9->Divide(5, 2);
    for(int l = 0; l < 10; l++){
        c9->cd(l+1);
        gConvDeltaZ[l]->Draw("ALP");
    }
}

