#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

// --- QUESTA FUNZIONE DISEGNA NEL SUBPLOT ASSEGNATO ---
void EseguiPlotLayer(std::map<double, std::string> files, std::string layerName, int slotIdx, TCanvas *mainCanvas) {
    if (files.empty()) return;

    int nStrips = 1536;
    double vMin = 0.0;
    double vMax = 130.0;
    int nBinsV = 130;

    std::string hName = "h3D_" + layerName;
    TH2F *h3D = new TH2F(hName.c_str(), 
                         (layerName + ";Strip;V [V];Occ").c_str(), 
                         nStrips, 0, nStrips, 
                         nBinsV, vMin, vMax);

    for (auto const& [voltage, path] : files) {
        if (gSystem->AccessPathName(path.c_str())) {
            std::cout << "!!! MANCANTE: " << layerName << " @ " << voltage << "V" << std::endl;
            continue;
        }

        std::ifstream infile(path);
        std::string line;
        while (std::getline(infile, line)) {
            if (line.empty() || line.find("Noise") != std::string::npos || line.find("Strip") != std::string::npos) continue;
            std::stringstream ss(line);
            int strip; double occupancy;
            if (ss >> strip >> occupancy) {
                h3D->SetBinContent(h3D->GetXaxis()->FindBin(strip), h3D->GetYaxis()->FindBin(voltage), occupancy);
            }
        }
        infile.close();
    }

    // Ci spostiamo nel riquadro corretto della griglia
    mainCanvas->cd(slotIdx);
    gPad->SetLogz(); // Scala logaritmica per il singolo subplot
    
    h3D->GetXaxis()->SetTitleOffset(1.2);
    h3D->GetYaxis()->SetTitleOffset(1.2);
    h3D->SetStats(0); // Rimuove box statistiche per pulizia
    
    h3D->Draw("LEGO2Z");
}

void plot_occupancy() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    // Creiamo un'unica Canvas gigante e la dividiamo in 2 colonne e 5 righe
    TCanvas *cMain = new TCanvas("cMain", "Occupancy All Layers", 1200, 1000);
    cMain->Divide(2, 5);

    // --- CARICAMENTO DATI (I tuoi percorsi) ---
    
    // LAYER 1: X0
    std::map<double, std::string> L1;
    L1[3.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-3V/TkrNoiseOccupancy_LayerX0_333001100.tnt";
    L1[7.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-7V/TkrNoiseOccupancy_LayerX0_333001099.tnt";
    L1[15.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-15V/TkrNoiseOccupancy_LayerX0_333001098.tnt";
    L1[40.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-40V/TkrNoiseOccupancy_LayerX0_333001097.tnt";
    L1[60.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-60V/TkrNoiseOccupancy_LayerX0_333001096.tnt";
    L1[90.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-90V/TkrNoiseOccupancy_LayerX0_333001092.tnt";
    L1[120.0] = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30/TkrNoiseOccupancy_LayerX0_333001082.tnt";
    EseguiPlotLayer(L1, "X0", 1, cMain);

    // LAYER 2: X1
    std::map<double, std::string> L2;
    L2[3.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-3V/TkrNoiseOccupancy_LayerX1_333001100.tnt";
    L2[7.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-7V/TkrNoiseOccupancy_LayerX1_333001099.tnt";
    L2[15.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-15V/TkrNoiseOccupancy_LayerX1_333001098.tnt";
    L2[40.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-40V/TkrNoiseOccupancy_LayerX1_333001097.tnt";
    L2[60.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-60V/TkrNoiseOccupancy_LayerX1_333001096.tnt";
    L2[90.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-90V/TkrNoiseOccupancy_LayerX1_333001092.tnt";
    L2[120.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30/TkrNoiseOccupancy_LayerX1_333001082.tnt";
    EseguiPlotLayer(L2, "X1", 2, cMain);

    // LAYER 3: X2
    std::map<double, std::string> L3;
    L3[3.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-3V/TkrNoiseOccupancy_LayerX2_333001100.tnt";
    L3[7.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-7V/TkrNoiseOccupancy_LayerX2_333001099.tnt";
    L3[15.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-15V/TkrNoiseOccupancy_LayerX2_333001098.tnt";
    L3[40.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-40V/TkrNoiseOccupancy_LayerX2_333001097.tnt";
    L3[60.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-60V/TkrNoiseOccupancy_LayerX2_333001096.tnt";
    L3[90.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-90V/TkrNoiseOccupancy_LayerX2_333001092.tnt";
    L3[120.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30/TkrNoiseOccupancy_LayerX2_333001082.tnt";
    EseguiPlotLayer(L3, "X2", 3, cMain);

    // LAYER 4: X3
    std::map<double, std::string> L4;
    L4[3.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-3V/TkrNoiseOccupancy_LayerX3_333001100.tnt";
    L4[7.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-7V/TkrNoiseOccupancy_LayerX3_333001099.tnt";
    L4[15.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-15V/TkrNoiseOccupancy_LayerX3_333001098.tnt";
    L4[40.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-40V/TkrNoiseOccupancy_LayerX3_333001097.tnt";
    L4[60.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-60V/TkrNoiseOccupancy_LayerX3_333001096.tnt";
    L4[90.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-90V/TkrNoiseOccupancy_LayerX3_333001092.tnt";
    L4[120.0]="/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30/TkrNoiseOccupancy_LayerX3_333001082.tnt";
    EseguiPlotLayer(L4, "X3", 4, cMain);

    // LAYER 5: X4
    std::map<double, std::string> L5;
    L5[3.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-3V/TkrNoiseOccupancy_LayerX4_333001100.tnt";
    L5[7.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-7V/TkrNoiseOccupancy_LayerX4_333001099.tnt";
    L5[15.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-15V/TkrNoiseOccupancy_LayerX4_333001098.tnt";
    L5[40.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-40V/TkrNoiseOccupancy_LayerX4_333001097.tnt";
    L5[60.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-60V/TkrNoiseOccupancy_LayerX4_333001096.tnt";
    L5[90.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-90V/TkrNoiseOccupancy_LayerX4_333001092.tnt";
    L5[120.0] = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30/TkrNoiseOccupancy_LayerX4_333001082.tnt";
    // ... segui lo schema sopra per gli altri file di X4 ...
    EseguiPlotLayer(L5, "X4", 5, cMain);

    // LAYER 6: Y0
    std::map<double, std::string> L6;
    L6[3.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-3V/TkrNoiseOccupancy_LayerY0_333001100.tnt";
    L6[7.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-7V/TkrNoiseOccupancy_LayerY0_333001099.tnt";
    L6[15.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-15V/TkrNoiseOccupancy_LayerY0_333001098.tnt";
    L6[40.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-40V/TkrNoiseOccupancy_LayerY0_333001097.tnt";
    L6[60.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-60V/TkrNoiseOccupancy_LayerY0_333001096.tnt";
    L6[90.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-90V/TkrNoiseOccupancy_LayerY0_333001092.tnt";
    L6[120.0] = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30/TkrNoiseOccupancy_LayerY0_333001082.tnt";
    EseguiPlotLayer(L6, "Y0", 6, cMain);

    // LAYER 7: Y1
    std::map<double, std::string> L7;
    L7[3.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-3V/TkrNoiseOccupancy_LayerY1_333001100.tnt";
    L7[7.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-7V/TkrNoiseOccupancy_LayerY1_333001099.tnt";
    L7[15.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-15V/TkrNoiseOccupancy_LayerY1_333001098.tnt";
    L7[40.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-40V/TkrNoiseOccupancy_LayerY1_333001097.tnt";
    L7[60.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-60V/TkrNoiseOccupancy_LayerY1_333001096.tnt";
    L7[90.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-90V/TkrNoiseOccupancy_LayerY1_333001092.tnt";
    L7[120.0] = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30/TkrNoiseOccupancy_LayerY1_333001082.tnt";
    EseguiPlotLayer(L7, "Y1", 7, cMain);

    // LAYER 8: Y2
    std::map<double, std::string> L8;
    L8[3.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-3V/TkrNoiseOccupancy_LayerY2_333001100.tnt";
    L8[7.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-7V/TkrNoiseOccupancy_LayerY2_333001099.tnt";
    L8[15.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-15V/TkrNoiseOccupancy_LayerY2_333001098.tnt";
    L8[40.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-40V/TkrNoiseOccupancy_LayerY2_333001097.tnt";
    L8[60.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-60V/TkrNoiseOccupancy_LayerY2_333001096.tnt";
    L8[90.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-90V/TkrNoiseOccupancy_LayerY2_333001092.tnt";
    L8[120.0] = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30/TkrNoiseOccupancy_LayerY2_333001082.tnt";
    EseguiPlotLayer(L8, "Y2", 8, cMain);

    // LAYER 9: Y3
    std::map<double, std::string> L9;
    L9[3.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-3V/TkrNoiseOccupancy_LayerY3_333001100.tnt";
    L9[7.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-7V/TkrNoiseOccupancy_LayerY3_333001099.tnt";
    L9[15.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-15V/TkrNoiseOccupancy_LayerY3_333001098.tnt";
    L9[40.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-40V/TkrNoiseOccupancy_LayerY3_333001097.tnt";
    L9[60.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-60V/TkrNoiseOccupancy_LayerY3_333001096.tnt";
    L9[90.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-90V/TkrNoiseOccupancy_LayerY3_333001092.tnt";
    L9[120.0] = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30/TkrNoiseOccupancy_LayerY3_333001082.tnt";
    EseguiPlotLayer(L9, "Y3", 9, cMain);

    // LAYER 10: Y4
    std::map<double, std::string> L10;
    L10[3.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-3V/TkrNoiseOccupancy_LayerY4_333001100.tnt";
    L10[7.0]   = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-7V/TkrNoiseOccupancy_LayerY4_333001099.tnt";
    L10[15.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-15V/TkrNoiseOccupancy_LayerY4_333001098.tnt";
    L10[40.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-40V/TkrNoiseOccupancy_LayerY4_333001097.tnt";
    L10[60.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-60V/TkrNoiseOccupancy_LayerY4_333001096.tnt";
    L10[90.0]  = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30-90V/TkrNoiseOccupancy_LayerY4_333001092.tnt";
    L10[120.0] = "/mnt/c/Users/Ludovica/Desktop/uni/MAGISTRALE/Lab int fond/Fermi-Glast/perry_blast/14-04-NO-30/TkrNoiseOccupancy_LayerY4_333001082.tnt";
    EseguiPlotLayer(L10, "Y4", 10, cMain);

    cMain->Update();
    std::cout << "Fatto! Griglia 2x5 generata." << std::endl;
}