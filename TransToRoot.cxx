#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"

// 定义一个结构体来存储每个事件的数据
struct EventData {
    int event_id;
    double b;
    int npart;
    double mult;
    double e2, e3, e4, e5;
    double x2, y2, xy, x2y2;
    double energy;
};

void convertDatToRoot(const std::string& dat_filename, const std::string& root_filename) {
    // 打开.dat文件
    std::ifstream dat_file(dat_filename);
    if (!dat_file.is_open()) {
        std::cerr << "Error: Could not open .dat file: " << dat_filename << std::endl;
        return;
    }

    // 创建ROOT文件
    TFile* root_file = new TFile(root_filename.c_str(), "RECREATE");
    if (!root_file->IsOpen()) {
        std::cerr << "Error: Could not create ROOT file: " << root_filename << std::endl;
        dat_file.close();
        return;
    }

    // 创建TTree
    TTree* tree = new TTree("event_tree", "Collision Event Data");

    // 创建一个EventData对象来存储每一行数据
    EventData event;

    // 将TTree的分支与结构体中的变量关联起来
    tree->Branch("event_id", &event.event_id, "event_id/I");
    tree->Branch("b", &event.b, "b/D");
    tree->Branch("npart", &event.npart, "npart/I");
    tree->Branch("mult", &event.mult, "mult/D");
    tree->Branch("e2", &event.e2, "e2/D");
    tree->Branch("e3", &event.e3, "e3/D");
    tree->Branch("e4", &event.e4, "e4/D");
    tree->Branch("e5", &event.e5, "e5/D");
    tree->Branch("x2", &event.x2, "x2/D");
    tree->Branch("y2", &event.y2, "y2/D");
    tree->Branch("xy", &event.xy, "xy/D");
    tree->Branch("x2y2", &event.x2y2, "x2y2/D");
    tree->Branch("energy", &event.energy, "energy/D");

    // 逐行读取.dat文件并填充TTree
    std::cout << "Reading " << dat_filename << " and writing to " << root_filename << "..." << std::endl;
    while (dat_file >> event.event_id >> event.b >> event.npart >> event.mult >> 
           event.e2 >> event.e3 >> event.e4 >> event.e5 >> 
           event.x2 >> event.y2 >> event.xy >> event.x2y2>> event.energy) {
        tree->Fill();
    }

    // 将TTree写入ROOT文件并关闭文件
    root_file->Write();  

    std::cout << "Conversion complete! " << tree->GetEntries() << " events written." << std::endl;

    root_file->Close();
    dat_file.close();
    delete root_file;

}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_dat_file> <output_root_file>" << std::endl;
        return 1;
    }
    
    std::string dat_file = argv[1];
    std::string root_file = argv[2];
    
    convertDatToRoot(dat_file, root_file);

    return 0;
}