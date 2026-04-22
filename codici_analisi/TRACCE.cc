#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TView.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

// --- COSTANTI GEOMETRICHE ---
const double StripPitch = 0.228;
const double EdgeWidth = 1.0;
const double LadderSeparation = 0.200;
const double PI = 3.141592653589793;

double getZ(int layerIdx) {
    double z_coords[10] = {0, 25, 75, 100, 150, 175, 200, 250, 275, 300};
    return z_coords[layerIdx];
}

double getCoordinate(int stripNum) {
    double coord = EdgeWidth + StripPitch * stripNum + 
                   (LadderSeparation + 2 * EdgeWidth - StripPitch) * (int(stripNum / 384));
    return coord;
}

struct Hit {
    int layerIdx;
    double pos;
    double err;
};

struct Event {
    int id;
    std::vector<Hit> x_hits;
    std::vector<Hit> y_hits;
    double totalTime = 0;
};

// --- FUNZIONE PER DISEGNO 3D ---
void DrawEvent3D(Event& ev, double mx, double qx, double my, double qy, int eventNum) {
    TCanvas *c3d = new TCanvas(Form("c3d_%d", eventNum), Form("Evento 3D ID %d", ev.id), 600, 600);
    TView *view = TView::CreateView(1);
    view->SetRange(-50, -50, -10, 400, 400, 350); // Range geometrico dello strumento (mm)

    // Disegna i punti (Hits)
    TPolyMarker3D *pm = new TPolyMarker3D();
    for(auto& hx : ev.x_hits) {
        // Un punto 3D richiede X, Y e Z. Usiamo la retta Y(z) per stimare la Y del punto X
        pm->SetNextPoint(hx.pos, my*getZ(hx.layerIdx) + qy, getZ(hx.layerIdx));
    }
    for(auto& hy : ev.y_hits) {
        pm->SetNextPoint(mx*getZ(hy.layerIdx) + qx, hy.pos, getZ(hy.layerIdx));
    }
    pm->SetMarkerStyle(20);
    pm->SetMarkerColor(kRed);
    pm->Draw();

    // Disegna la traccia (Retta fittata)
    TPolyLine3D *line = new TPolyLine3D(2);
    double zMin = -20; double zMax = 320;
    line->SetPoint(0, mx*zMin + qx, my*zMin + qy, zMin);
    line->SetPoint(1, mx*zMax + qx, my*zMax + qy, zMax);
    line->SetLineColor(kBlue);
    line->SetLineWidth(2);
    line->Draw();
    
    c3d->Update();
}

// [Mantenere la funzione parseLIF identica a quella fornita dall'utente]
std::vector<Event> parseLIF(std::string path) {
    std::ifstream file(path);
    std::vector<Event> events;
    if (!file.is_open()) { std::cerr << "Errore apertura file!" << std::endl; return events; }
    std::string line;
    getline(file, line); 
    int n, evtID, nHits;
    while (file >> n >> evtID >> nHits) {
        Event ev;
        ev.id = evtID;
        std::vector<int> layerStrips[10];
        for (int i = 0; i < nHits; i++) {
            std::string lName; int lEnd, strip, chip, chan;
            file >> lName >> lEnd >> strip >> chip >> chan;
            int idx = -1;
            if(lName == "Y0") idx=0; else if(lName == "X0") idx=1;
            else if(lName == "X1") idx=2; else if(lName == "Y1") idx=3;
            else if(lName == "Y2") idx=4; else if(lName == "X2") idx=5;
            else if(lName == "X3") idx=6; else if(lName == "Y3") idx=7;
            else if(lName == "Y4") idx=8; else if(lName == "X4") idx=9;
            if(idx != -1) layerStrips[idx].push_back(strip);
        }
        for (int i = 0; i < 10; i++) {
            if (layerStrips[i].empty()) continue;
            std::sort(layerStrips[i].begin(), layerStrips[i].end());
            std::vector<std::vector<int>> clusters;
            std::vector<int> current;
            for(int s : layerStrips[i]) {
                if(current.empty() || s == current.back()+1 || s == current.back()-1) current.push_back(s);
                else { clusters.push_back(current); current.clear(); current.push_back(s); }
            }
            clusters.push_back(current);
            for(auto& cl : clusters) {
                double meanStrip = (cl.front() + cl.back()) / 2.0;
                Hit h; h.layerIdx = i; h.pos = getCoordinate(meanStrip);
                h.err = (cl.size() * StripPitch) / sqrt(12);
                if ( i == 1 || i == 2 || i == 5 || i == 6 || i == 9) ev.x_hits.push_back(h);
                else ev.y_hits.push_back(h);
            }
        }
        int nTOT; file >> nTOT;
        for (int i = 0; i < nTOT; i++) {
            std::string ln; int le, tdc; double tus;
            file >> ln >> le >> tdc >> tus;
            ev.totalTime += tus;
        }
        if (n % 2 == 0) events.push_back(ev);
    }
    return events;
}

void CosmicReconstruction() {
    std::string path = R"(/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/21-04-30-all-masked/TkrDataTaking_333001114.lif)";
    std::vector<Event> events = parseLIF(path);
    TH1D* hAngles = new TH1D("hAngles", "Distribuzione Angolare;Theta [deg];Eventi", 18, 0, 90);

    int drawnCount = 0;
    for (auto& ev : events) {
        if (ev.x_hits.size() < 2 || ev.y_hits.size() < 2) continue;

        TGraphErrors grX;
        for(auto& h : ev.x_hits) grX.SetPoint(grX.GetN(), getZ(h.layerIdx), h.pos);
        grX.Fit("pol1", "Q");
        TF1* fitX = grX.GetFunction("pol1");

        TGraphErrors grY;
        for(auto& h : ev.y_hits) grY.SetPoint(grY.GetN(), getZ(h.layerIdx), h.pos);
        grY.Fit("pol1", "Q");
        TF1* fitY = grY.GetFunction("pol1");

        double mx = fitX->GetParameter(1); double qx = fitX->GetParameter(0);
        double my = fitY->GetParameter(1); double qy = fitY->GetParameter(0);

        double cosTheta = 1.0 / sqrt(mx*mx + my*my + 1.0);
        double thetaDeg = acos(cosTheta) * 180.0 / PI;
        hAngles->Fill(thetaDeg);

        // DISEGNA LE PRIME 5 TRACCE
        if (drawnCount < 5) {
            DrawEvent3D(ev, mx, qx, my, qy, drawnCount);
            drawnCount++;
        }
    }

    TCanvas* c1 = new TCanvas("c1", "Analisi Cosmici", 800, 600);
    hAngles->SetFillColor(kCyan-10);
    hAngles->Draw();
}