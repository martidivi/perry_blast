/*   Simulazione di 10 piani di rivelatori di cosmici con strip di silicio.
 
     Instruciton usage:
     To simulate 100 cosmic rays with a polar distribution proportional to Cos(theta)^2.5, root prompt:
     .L CosmicMC.cc++
     CosmicMC(100,3.5)
 */
// Units: length in meters
// Units: time in nseconds
// Axis: up vertical +z, right +x, z x x = y

#include "_detector_utils.h"
#include "_readfile.cc"

#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TView.h>

#include "TRandom3.h"
#include <iostream>

// --- DEFINISCO STRUTTURE
struct Layer {
    std::string name;
    char orientation;
    double z_pos;
    int box_id;
    int numControllers;
};

struct CosmicRay {
    double x0, y0, z0; // punto di generazione
    double vx, vy, vz; // versore della direzione
    double thetaDeg;
};

// --- FUNZIONI DI SUPPORTO ---
CosmicRay GenerateCosmic(TRandom3 &gen) {
    CosmicRay r;
    r.z0 = 172.5;
    r.x0 = gen.Uniform(-100, 500);
    r.y0 = gen.Uniform(-100, 500);
    double u = gen.Uniform(0, 1);
    double theta = acos(pow(u, 1.0/3.0));
    double phi = gen.Uniform(0, 2 * M_PI);
    r.thetaDeg = theta * TMath::RadToDeg();
    r.vx = sin(theta) * cos(phi);
    r.vy = sin(theta) * sin(phi);
    r.vz = cos(theta);
    return r;
}

// --- SIMULAZIONE UNIFICATA ---
void CosmicSim(int nGen = 10000) {
    TRandom3 gen(0);
    std::vector<Layer> detector;
    
    // 1. Inizializzazione Geometria
    double current_z = 0.0;
    for (int i = 0; i < 10; ++i) {
        Layer l;
        l.name = layerNames[i];
        l.orientation = layerNames[i][0];
        l.z_pos = current_z;
        l.box_id = (i == 0) ? 5 : (i == 9) ? 0: 5 - (i+1)/2;
        l.numControllers = (i ==2  || i ==6 ) ? 1 : 2;

        if (i > 0) {
            if (i % 2 != 0) current_z+=InterBoxGap;
            else current_z+= BoxHeight;
        } else current_z = 0.0;
        detector.push_back(l);
    }

    // 2. Setup Grafica
    TCanvas *c1 = new TCanvas("c1", "Primi 100 eventi simulati", 1000, 1000);
    TView *view = TView::CreateView(1);
    view->SetRange(-200, -200,-200, 600, 600, 600);
    for(auto& l : detector) {
        TPolyLine3D *pl = new TPolyLine3D(5);
        pl->SetPoint(0,0,0,l.z_pos); pl->SetPoint(1,400,0,l.z_pos);
        pl->SetPoint(2,400,400,l.z_pos); pl->SetPoint(3,0,400,l.z_pos);
        pl->SetPoint(4,0,0,l.z_pos);
        pl->SetLineColor(kGray); pl->Draw();
    }

    TCanvas *c2 = new TCanvas("c2", "Analisi Statistica", 1000, 500);
    c2->Divide(2,1);
    TH1D *hSorgente = new TH1D("hSorgente", "Sorgente (cos^2);Theta [deg];Conteggi", 90, 0, 90);
    TH1D *hTrigger = new TH1D("hTrigger", "Triggerati (Accettanza);Theta [deg];Conteggi", 90, 0, 90);

    // 3. Loop
    int triggerCount = 0;
    int oneLayerCount = 0; // <--- NUOVO CONTATORE

    for (int i = 0; i < nGen; i++) {
        CosmicRay ray = GenerateCosmic(gen);
        hSorgente->Fill(ray.thetaDeg);

        std::map<int, bool> hit_layer;
        int totalHit = 0;
        bool geomHit = false;

        // 1. Dobbiamo usare un indice i per poter riempire la mappa hit_layer[i]
        for (int i = 0; i < (int)detector.size(); ++i) {
            auto& l = detector[i];
            double t = (l.z_pos - ray.z0) / ray.vz;
            double xi = ray.x0 + t * ray.vx;
            double yi = ray.y0 + t * ray.vy;

            if (xi >= 0 && xi <= 400 && yi >= 0 && yi <= 400) {
                geomHit = true;
                
                double coord = (l.orientation == 'x') ? yi : xi;
                if (GetStripNumber(coord) != -1) {
                    hit_layer[i] = true; // riempie la mappa!!
                    totalHit++;
                }
            }
        }

        if (geomHit) oneLayerCount++;

        bool threeLayers = false;
        // 2. Il ciclo sulle triplette ora trova i dati dentro hit_layer
        for(int l = 0; l < 8; l++) {
           if (hit_layer[l] && hit_layer[l+1] && hit_layer[l+2]) {
               threeLayers = true;
                break;
           }
        }

        bool isTriggered = threeLayers;

//for(int l=0; l<=9; l++)
//          if( (l % 2 == 0 &&  x_hit_box[l.box_id] && y_hit_box[(l.box_id)-1]) || (l % 2 != 0 &&  x_hit_box[l.box_id] && y_hit_box[(l.box_id)+1]))
//         coincidenza = true;
//      bool isTriggered = (coincidenza && totalX >= 2 && totalY >= 2);

        if (isTriggered) {
            triggerCount++;
            hTrigger->Fill(ray.thetaDeg);
        }

        if (i < 100) {
            c1->cd();
            TPolyLine3D *tr = new TPolyLine3D(2);
            tr->SetPoint(0,  ray.x0 - 300*ray.vx, ray.y0 - 300*ray.vy, ray.z0 - 300*ray.vz);
            tr->SetPoint(1, ray.x0 + 300*ray.vx, ray.y0 + 300*ray.vy, ray.z0 + 300*ray.vz);
            tr->SetLineColor(isTriggered ? kRed : kBlue);
            tr->SetLineStyle(isTriggered ? 1 : 3);
            tr->Draw();
        }
    }

    c2->cd(1); hSorgente->Draw();
    c2->cd(2); hTrigger->SetLineColor(kRed); hTrigger->Draw();
    
    std::cout << "--- RISULTATI ---" << std::endl;
    std::cout << "Raggi totali generati (60x60): " << nGen << std::endl;
    std::cout << "Raggi che hanno colpito almeno un layer (40x40): " << oneLayerCount << std::endl;
    std::cout << "Eventi Triggerati: " << triggerCount << std::endl;
    std::cout << "Accettanza Geometrica Assoluta: " << (double)triggerCount/nGen * 100 << "%" << std::endl;
    std::cout << "Efficienza di Trigger (su raggi entranti): " << (double)triggerCount/oneLayerCount * 100 << "%" << std::endl;
}
