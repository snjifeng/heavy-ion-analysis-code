#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TString.h"
#include "TH1D.h"
#include "TAxis.h"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <limits> // For numeric_limits

// ===========================================
// 全局指针声明
// ===========================================
TProfile *epsilon2Nch = nullptr;
TProfile *epsilon22Nch = nullptr;
TProfile *epsilon24Nch = nullptr;
TProfile *epsilon26Nch = nullptr;
TProfile *epsilon2Cen = nullptr;
TProfile *epsilon22Cen = nullptr;
TProfile *epsilon24Cen = nullptr;
TProfile *epsilon26Cen = nullptr;
TProfile *epsilon2b = nullptr;
TProfile *epsilon22b = nullptr;
TProfile *epsilon24b = nullptr;
TProfile *epsilon26b = nullptr;

TProfile *epsilon3Nch = nullptr;
TProfile *epsilon32Nch = nullptr;
TProfile *epsilon34Nch = nullptr;
TProfile *epsilon3Cen = nullptr;
TProfile *epsilon32Cen = nullptr;
TProfile *epsilon34Cen = nullptr;
TProfile *epsilon3b = nullptr;
TProfile *epsilon32b = nullptr;
TProfile *epsilon34b = nullptr;

TProfile *epsilon4Nch = nullptr;
TProfile *epsilon42Nch = nullptr;
TProfile *epsilon44Nch = nullptr;
TProfile *epsilon4Cen = nullptr;
TProfile *epsilon42Cen = nullptr;
TProfile *epsilon44Cen = nullptr;
TProfile *epsilon4b = nullptr;
TProfile *epsilon42b = nullptr;
TProfile *epsilon44b = nullptr;

TProfile *epsilon2242Nch = nullptr;
TProfile *epsilon2242Cen = nullptr;
TProfile *epsilon2242b = nullptr;

TProfile *epsilon2232Nch = nullptr;
TProfile *epsilon2232Cen = nullptr;
TProfile *epsilon2232b = nullptr;

TProfile *dc1Nch = nullptr; 
TProfile *dc1Cen = nullptr; 
TProfile *dc1b = nullptr;

TProfile *dc2Nch = nullptr; 
TProfile *dc2Cen = nullptr; 
TProfile *dc2b = nullptr;

TProfile *dc3Nch = nullptr; 
TProfile *dc3Cen = nullptr; 
TProfile *dc3b = nullptr;

TProfile *dc4Nch = nullptr; 
TProfile *dc4Cen = nullptr; 
TProfile *dc4b = nullptr;

TProfile *ddc12Nch = nullptr; 
TProfile *ddc12Cen = nullptr; 
TProfile *ddc12b = nullptr;

TProfile *ddc22Nch = nullptr; 
TProfile *ddc22Cen = nullptr; 
TProfile *ddc22b = nullptr;

TProfile *ddc32Nch = nullptr; 
TProfile *ddc32Cen = nullptr; 
TProfile *ddc32b = nullptr;

TProfile *ddc42Nch = nullptr; 
TProfile *ddc42Cen = nullptr; 
TProfile *ddc42b = nullptr;

TProfile *ddc13Nch = nullptr; 
TProfile *ddc13Cen = nullptr; 
TProfile *ddc13b = nullptr;

TProfile *ddc23Nch = nullptr; 
TProfile *ddc23Cen = nullptr; 
TProfile *ddc23b = nullptr;

TProfile *ddc33Nch = nullptr; 
TProfile *ddc33Cen = nullptr; 
TProfile *ddc33b = nullptr;

TProfile *ddc43Nch = nullptr; 
TProfile *ddc43Cen = nullptr; 
TProfile *ddc43b = nullptr;

TProfile *ddc14Nch = nullptr; 
TProfile *ddc14Cen = nullptr; 
TProfile *ddc14b = nullptr;

TProfile *ddc24Nch = nullptr; 
TProfile *ddc24Cen = nullptr; 
TProfile *ddc24b = nullptr;

TProfile *ddc34Nch = nullptr; 
TProfile *ddc34Cen = nullptr; 
TProfile *ddc34b = nullptr;

TProfile *ddc44Nch = nullptr; 
TProfile *ddc44Cen = nullptr; 
TProfile *ddc44b = nullptr;

TProfile *epsilon22ddc1Nch = nullptr;
TProfile *epsilon22ddc2Nch = nullptr;
TProfile *epsilon22ddc3Nch = nullptr;
TProfile *epsilon22ddc4Nch = nullptr;

TProfile *epsilon22ddc1Cen = nullptr;
TProfile *epsilon22ddc2Cen = nullptr;
TProfile *epsilon22ddc3Cen = nullptr;
TProfile *epsilon22ddc4Cen = nullptr;

TProfile *epsilon22ddc1b = nullptr;
TProfile *epsilon22ddc2b = nullptr;
TProfile *epsilon22ddc3b = nullptr;
TProfile *epsilon22ddc4b = nullptr;

TProfile *epsilon22ddc12Nch = nullptr;
TProfile *epsilon22ddc22Nch = nullptr;
TProfile *epsilon22ddc32Nch = nullptr;
TProfile *epsilon22ddc42Nch = nullptr;

TProfile *epsilon22ddc12Cen = nullptr;
TProfile *epsilon22ddc22Cen = nullptr;
TProfile *epsilon22ddc32Cen = nullptr;
TProfile *epsilon22ddc42Cen = nullptr;

TProfile *epsilon22ddc12b = nullptr;
TProfile *epsilon22ddc22b = nullptr;
TProfile *epsilon22ddc32b = nullptr;
TProfile *epsilon22ddc42b = nullptr;

TProfile *epsilon32ddc1Nch = nullptr;
TProfile *epsilon32ddc2Nch = nullptr;
TProfile *epsilon32ddc3Nch = nullptr;
TProfile *epsilon32ddc4Nch = nullptr;

TProfile *epsilon32ddc1Cen = nullptr;
TProfile *epsilon32ddc2Cen = nullptr;
TProfile *epsilon32ddc3Cen = nullptr;
TProfile *epsilon32ddc4Cen = nullptr;

TProfile *epsilon32ddc1b = nullptr;
TProfile *epsilon32ddc2b = nullptr;
TProfile *epsilon32ddc3b = nullptr;
TProfile *epsilon32ddc4b = nullptr;

TProfile *epsilon32ddc12Nch = nullptr;
TProfile *epsilon32ddc22Nch = nullptr;
TProfile *epsilon32ddc32Nch = nullptr;
TProfile *epsilon32ddc42Nch = nullptr;

TProfile *epsilon32ddc12Cen = nullptr;
TProfile *epsilon32ddc22Cen = nullptr;
TProfile *epsilon32ddc32Cen = nullptr;
TProfile *epsilon32ddc42Cen = nullptr;

TProfile *epsilon32ddc12b = nullptr;
TProfile *epsilon32ddc22b = nullptr;
TProfile *epsilon32ddc32b = nullptr;
TProfile *epsilon32ddc42b = nullptr;

// ===========================================
// 函数声明
// ===========================================
void writeHistograms(TFile* FileOutput);
void deleteHistograms();
double get_centrality_percent(double mult);


// ===========================================
// 函数定义
// ===========================================

// 根据 mult 值和阈值返回中心度百分比
double get_centrality_percent(double mult) {
    // === 请将实际计算出的阈值粘贴到这里 ===
    // 键是中心度百分比的上限，值是对应的 mult 阈值
    // 这个map的定义在你的原始代码中是错误的，这里我进行了修正。
    std::map<double, double> mult_thresholds = {
        {0, 999.5},
        {0.2, 246.5},
        {0.5, 230.5},
        {1, 221.5},
        {2, 209.5},
        {3, 200.5},
        {4, 193.5},
        {5, 186.5},
        {6, 180.5},
        {7, 174.5},
        {8, 169.5},
        {9, 163.5},
        {10, 158.5},
        {15, 135.5},
        {20, 116.5},
        {30, 84.5},
        {40, 60.5},
        {50, 42.5},
        {60, 29.5},
        {70, 19.5},
        {80, 12.5},
        {90, 6.5},
    };
    
    // 定义中心度区间的上下限
    std::vector<std::pair<double, double>> centrality_bins = {
        {0.0, 0.2}, 
        {0.2, 0.5}, 
        {0.5, 1.0}, 
        {1.0, 2.0}, 
        {2.0, 3.0}, 
        {3.0, 4.0}, 
        {4.0, 5.0},
        {5.0, 6.0}, 
        {6.0, 7.0}, 
        {7.0, 8.0}, 
        {8.0, 9.0}, 
        {9.0, 10.0},
        {10.0, 15.0},
        {15.0, 20.0}, 
        {20.0, 30.0}, 
        {30.0, 40.0}, 
        {40.0, 50.0},
        {50.0, 60.0}, 
        {60.0, 70.0}, 
        {70.0, 80.0}, 
        {80.0, 90.0}, 
        {90.0, 100.0}

    };

    for (const auto& bin : centrality_bins) {
        double lower_percentile = bin.first;
        double upper_percentile = bin.second;

        // 确保 map 中存在键，以防抛出异常
        if (mult_thresholds.find(lower_percentile) == mult_thresholds.end() || 
            mult_thresholds.find(upper_percentile) == mult_thresholds.end()) {
            continue;
        }

        double lower_mult = mult_thresholds.at(lower_percentile);
        double upper_mult = mult_thresholds.at(upper_percentile);

        // 核心逻辑：判断 mult 落在哪个区间
        if (mult <= lower_mult && mult > upper_mult) {
            return (lower_percentile + upper_percentile) / 2.0;
        }
    }

    // 如果 mult 值超出所有阈值范围，返回-1.0
    return 999.0; 
}


// 主分析函数
void analyze_and_write() {
    // TFile* f = TFile::Open("initial.root");
    // if (!f || f->IsZombie()) {
    //     printf("Error: Could not open file 'initial.root'.\n");
    //     return;
    // }

    // TTree* tree1 = (TTree*)f->Get("event_tree");
    // if (!tree1) {
    //     printf("Error: Could not find TTree 'event_tree' in the file.\n");
    //     f->Close();
    //     return;
    // }


    TFile* f1 = TFile::Open("./initial.root");
    TFile* f2 = TFile::Open("../ws_AV8/initial.root");
    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie()) {
        printf("Error opening files.\n");
        return;
    }

    TTree* tree1 = (TTree*)f1->Get("event_tree");
    TTree* tree2 = (TTree*)f2->Get("event_tree");
    if (!tree1 || !tree2) {
        printf("Error getting trees.\n");
        f1->Close(); f2->Close();
        return;
    }

    // 创建一个新的 TChain 来合并
    TChain chain("event_tree");
    chain.Add(f1->GetName());
    chain.Add(f2->GetName());

    // ===========================================
    // 创建所有 TProfile 直方图
    // ===========================================
    int NBch = 20, NBcen = 1000, NBb = 20*10; // 调整Nch的bin数以覆盖更大的范围
    double max_cen = 100;
    double max_mult = 400; // 假设 mult 最大值在10000左右，请根据你的数据调整
    double max_b = 20;     // 假设 b 最大值在20左右

    epsilon2Nch = new TProfile("epsilon2Nch", "e_{2} vs. Nch", NBch, 0, max_mult);
    epsilon2Nch->Sumw2();
    epsilon22Nch = new TProfile("epsilon22Nch", "e_{2}^{2} vs. Nch", NBch, 0, max_mult);
    epsilon22Nch->Sumw2();
    epsilon24Nch = new TProfile("epsilon24Nch", "e_{2}^{4} vs. Nch", NBch, 0, max_mult);
    epsilon24Nch->Sumw2();
    epsilon26Nch = new TProfile("epsilon26Nch", "e_{2}^{6} vs. Nch", NBch, 0, max_mult);
    epsilon26Nch->Sumw2();

    epsilon2Cen = new TProfile("epsilon2Cen", "e_{2} vs. Centrality", NBcen, 0, max_cen);
    epsilon2Cen->Sumw2();
    epsilon22Cen = new TProfile("epsilon22Cen", "e_{2}^{2} vs. Centrality", NBcen, 0, max_cen);
    epsilon22Cen->Sumw2();
    epsilon24Cen = new TProfile("epsilon24Cen", "e_{2}^{4} vs. Centrality", NBcen, 0, max_cen);
    epsilon24Cen->Sumw2();
    epsilon26Cen = new TProfile("epsilon26Cen", "e_{2}^{6} vs. Centrality", NBcen, 0, max_cen);
    epsilon26Cen->Sumw2();
    epsilon2b = new TProfile("epsilon2b", "e_{2} vs. b", NBb, 0, max_b);
    epsilon2b->Sumw2();
    epsilon22b = new TProfile("epsilon22b", "e_{2}^{2} vs. b", NBb, 0, max_b);
    epsilon22b->Sumw2();
    epsilon24b = new TProfile("epsilon24b", "e_{2}^{4} vs. b", NBb, 0, max_b);
    epsilon24b->Sumw2();
    epsilon26b = new TProfile("epsilon26b", "e_{2}^{6} vs. b", NBb, 0, max_b);
    epsilon26b->Sumw2();

    epsilon3Nch = new TProfile("epsilon3Nch", "e_{3} vs. Nch", NBch, 0, max_mult);
    epsilon3Nch->Sumw2();
    epsilon32Nch = new TProfile("epsilon32Nch", "e_{3}^{2} vs. Nch", NBch, 0, max_mult);
    epsilon32Nch->Sumw2();
    epsilon34Nch = new TProfile("epsilon34Nch", "e_{3}^{4} vs. Nch", NBch, 0, max_mult);
    epsilon34Nch->Sumw2();
    epsilon3Cen = new TProfile("epsilon3Cen", "e_{3} vs. Centrality", NBcen, 0, max_cen);
    epsilon3Cen->Sumw2();
    epsilon32Cen = new TProfile("epsilon32Cen", "e_{3}^{2} vs. Centrality", NBcen, 0, max_cen);
    epsilon32Cen->Sumw2();
    epsilon34Cen = new TProfile("epsilon34Cen", "e_{3}^{4} vs. Centrality", NBcen, 0, max_cen);
    epsilon34Cen->Sumw2();
    epsilon3b = new TProfile("epsilon3b", "e_{3} vs. b", NBb, 0, max_b);
    epsilon3b->Sumw2();
    epsilon32b = new TProfile("epsilon32b", "e_{3}^{2} vs. b", NBb, 0, max_b);
    epsilon32b->Sumw2();
    epsilon34b = new TProfile("epsilon34b", "e_{3}^{4} vs. b", NBb, 0, max_b);
    epsilon34b->Sumw2();

    epsilon4Nch = new TProfile("epsilon4Nch", "e_{4} vs. Nch", NBch, 0, max_mult);
    epsilon4Nch->Sumw2();
    epsilon42Nch = new TProfile("epsilon42Nch", "e_{4}^{2} vs. Nch", NBch, 0, max_mult);
    epsilon42Nch->Sumw2();
    epsilon44Nch = new TProfile("epsilon44Nch", "e_{4}^{4} vs. Nch", NBch, 0, max_mult);
    epsilon44Nch->Sumw2();
    epsilon4Cen = new TProfile("epsilon4Cen", "e_{4} vs. Centrality", NBcen, 0, max_cen);
    epsilon4Cen->Sumw2();
    epsilon42Cen = new TProfile("epsilon42Cen", "e_{4}^{2} vs. Centrality", NBcen, 0, max_cen);
    epsilon42Cen->Sumw2();
    epsilon44Cen = new TProfile("epsilon44Cen", "e_{4}^{4} vs. Centrality", NBcen, 0, max_cen);
    epsilon44Cen->Sumw2();
    epsilon4b = new TProfile("epsilon4b", "e_{4} vs. b", NBb, 0, max_b);
    epsilon4b->Sumw2();
    epsilon42b = new TProfile("epsilon42b", "e_{4}^{2} vs. b", NBb, 0, max_b);
    epsilon42b->Sumw2();
    epsilon44b = new TProfile("epsilon44b", "e_{4}^{4} vs. b", NBb, 0, max_b);
    epsilon44b->Sumw2();

    epsilon2232Nch = new TProfile("epsilon2232Nch", "e_{2}^{2}e_{3}^{2} vs. Nch", NBch, 0, max_mult);
    epsilon2232Nch->Sumw2();
    epsilon2232Cen = new TProfile("epsilon2232Cen", "e_{2}^{2}e_{3}^{2} vs. Centrality", NBcen, 0, max_cen);
    epsilon2232Cen->Sumw2();
    epsilon2232b = new TProfile("epsilon2232b", "e_{2}^{2}e_{3}^{2} vs. b", NBb, 0, max_b);
    epsilon2232b->Sumw2();

    epsilon2242Nch = new TProfile("epsilon2242Nch", "e_{2}^{2}e_{4}^{2} vs. Nch", NBch, 0, max_mult);
    epsilon2242Nch->Sumw2();
    epsilon2242Cen = new TProfile("epsilon2242Cen", "e_{2}^{2}e_{4}^{2} vs. Centrality", NBcen, 0, max_cen);
    epsilon2242Cen->Sumw2();
    epsilon2242b = new TProfile("epsilon2242b", "e_{2}^{2}e_{4}^{2} vs. b", NBb, 0, max_b);
    epsilon2242b->Sumw2();

    dc1Nch = new TProfile("dc1Nch", "dc1 vs. Nch", NBch, 0, max_mult);
    dc1Nch->Sumw2();
    dc1Cen = new TProfile("dc1Cen", "dc1vs. Centrality", NBcen, 0, max_cen);
    dc1Cen->Sumw2();
    dc1b = new TProfile("dc1b", "dc1 vs. b", NBb, 0, max_b);
    dc1b->Sumw2();

    dc2Nch = new TProfile("dc2Nch", "dc2 vs. Nch", NBch, 0, max_mult);
    dc2Nch->Sumw2();
    dc2Cen = new TProfile("dc2Cen", "dc2 vs. Centrality", NBcen, 0, max_cen);
    dc2Cen->Sumw2();
    dc2b = new TProfile("dc2b", "dc2 vs. b", NBb, 0, max_b);
    dc2b->Sumw2();

    dc3Nch = new TProfile("dc3Nch", "dc3 vs. Nch", NBch, 0, max_mult);
    dc3Nch->Sumw2();
    dc3Cen = new TProfile("dc3Cen", "dc3 vs. Centrality", NBcen, 0, max_cen);
    dc3Cen->Sumw2();
    dc3b = new TProfile("dc3b", "dc3 vs. b", NBb, 0, max_b);
    dc3b->Sumw2();

    dc4Nch = new TProfile("dc4Nch", "dc4 vs. Nch", NBch, 0, max_mult);
    dc4Nch->Sumw2();
    dc4Cen = new TProfile("dc4Cen", "dc4 vs. Centrality", NBcen, 0, max_cen);
    dc4Cen->Sumw2();
    dc4b = new TProfile("dc4b", "dc4 vs. b", NBb, 0, max_b);
    dc4b->Sumw2();

    ddc12Nch = new TProfile("ddc12Nch", "ddc12 vs. Nch", NBch, 0, max_mult);
    ddc12Nch->Sumw2();
    ddc12Cen = new TProfile("ddc12Cen", "ddc12 vs. Centrality", NBcen, 0, max_cen);
    ddc12Cen->Sumw2();
    ddc12b = new TProfile("ddc12b", "ddc12 vs. b", NBb, 0, max_b);
    ddc12b->Sumw2();

    ddc22Nch = new TProfile("ddc22Nch", "ddc22 vs. Nch", NBch, 0, max_mult);
    ddc22Nch->Sumw2();
    ddc22Cen = new TProfile("ddc22Cen", "ddc22 vs. Centrality", NBcen, 0, max_cen);
    ddc22Cen->Sumw2();
    ddc22b = new TProfile("ddc22b", "ddc22 vs. b", NBb, 0, max_b);
    ddc22b->Sumw2();

    ddc32Nch = new TProfile("ddc32Nch", "ddc32 vs. Nch", NBch, 0, max_mult);
    ddc32Nch->Sumw2();
    ddc32Cen = new TProfile("ddc32Cen", "ddc32 vs. Centrality", NBcen, 0, max_cen);
    ddc32Cen->Sumw2();
    ddc32b = new TProfile("ddc32b", "ddc32 vs. b", NBb, 0, max_b);
    ddc32b->Sumw2();

    ddc42Nch = new TProfile("ddc42Nch", "ddc42 vs. Nch", NBch, 0, max_mult);
    ddc42Nch->Sumw2();
    ddc42Cen = new TProfile("ddc42Cen", "ddc42 vs. Centrality", NBcen, 0, max_cen);
    ddc42Cen->Sumw2();
    ddc42b = new TProfile("ddc42b", "ddc42 vs. b", NBb, 0, max_b);
    ddc42b->Sumw2();

    ddc13Nch = new TProfile("ddc13Nch", "ddc13 vs. Nch", NBch, 0, max_mult);
    ddc13Nch->Sumw2();
    ddc13Cen = new TProfile("ddc13Cen", "ddc13 vs. Centrality", NBcen, 0, max_cen);
    ddc13Cen->Sumw2();
    ddc13b = new TProfile("ddc13b", "ddc13 vs. b", NBb, 0, max_b);
    ddc13b->Sumw2();

    ddc23Nch = new TProfile("ddc23Nch", "ddc23 vs. Nch", NBch, 0, max_mult);
    ddc23Nch->Sumw2();
    ddc23Cen = new TProfile("ddc23Cen", "ddc23 vs. Centrality", NBcen, 0, max_cen);
    ddc23Cen->Sumw2();
    ddc23b = new TProfile("ddc23b", "ddc23 vs. b", NBb, 0, max_b);
    ddc23b->Sumw2();

    ddc33Nch = new TProfile("ddc33Nch", "ddc33 vs. Nch", NBch, 0, max_mult);
    ddc33Nch->Sumw2();
    ddc33Cen = new TProfile("ddc33Cen", "ddc33 vs. Centrality", NBcen, 0, max_cen);
    ddc33Cen->Sumw2();
    ddc33b = new TProfile("ddc33b", "ddc33 vs. b", NBb, 0, max_b);
    ddc33b->Sumw2();

    ddc43Nch = new TProfile("ddc43Nch", "ddc43 vs. Nch", NBch, 0, max_mult);
    ddc43Nch->Sumw2();
    ddc43Cen = new TProfile("ddc43Cen", "ddc43 vs. Centrality", NBcen, 0, max_cen);
    ddc43Cen->Sumw2();
    ddc43b = new TProfile("ddc43b", "ddc43 vs. b", NBb, 0, max_b);
    ddc43b->Sumw2();

    ddc14Nch = new TProfile("ddc14Nch", "ddc14 vs. Nch", NBch, 0, max_mult);
    ddc14Nch->Sumw2();
    ddc14Cen = new TProfile("ddc14Cen", "ddc14 vs. Centrality", NBcen, 0, max_cen);
    ddc14Cen->Sumw2();
    ddc14b = new TProfile("ddc14b", "ddc14 vs. b", NBb, 0, max_b);
    ddc14b->Sumw2();

    ddc24Nch = new TProfile("ddc24Nch", "ddc24 vs. Nch", NBch, 0, max_mult);
    ddc24Nch->Sumw2();
    ddc24Cen = new TProfile("ddc24Cen", "ddc24 vs. Centrality", NBcen, 0, max_cen);
    ddc24Cen->Sumw2();
    ddc24b = new TProfile("ddc24b", "ddc24 vs. b", NBb, 0, max_b);
    ddc24b->Sumw2();

    ddc34Nch = new TProfile("ddc34Nch", "ddc34 vs. Nch", NBch, 0, max_mult);
    ddc34Nch->Sumw2();
    ddc34Cen = new TProfile("ddc34Cen", "ddc34 vs. Centrality", NBcen, 0, max_cen);
    ddc34Cen->Sumw2();
    ddc34b = new TProfile("ddc34b", "ddc34 vs. b", NBb, 0, max_b);
    ddc34b->Sumw2();

    ddc44Nch = new TProfile("ddc44Nch", "ddc44 vs. Nch", NBch, 0, max_mult);
    ddc44Nch->Sumw2();
    ddc44Cen = new TProfile("ddc44Cen", "ddc44 vs. Centrality", NBcen, 0, max_cen);
    ddc44Cen->Sumw2();
    ddc44b = new TProfile("ddc44b", "ddc44 vs. b", NBb, 0, max_b);
    ddc44b->Sumw2();

    epsilon22ddc1Nch= new TProfile("epsilon22ddc1Nch", "e_{2}^{2}ddc1 vs. Nch", NBch, 0, max_mult);
    epsilon22ddc1Nch->Sumw2();
    epsilon22ddc1Cen= new TProfile("epsilon22ddc1Cen", "e_{2}^{2}ddc1 vs. Centrality", NBcen, 0, max_cen);
    epsilon22ddc1Cen->Sumw2();
    epsilon22ddc1b= new TProfile("epsilon22ddc1b", "e_{2}^{2}ddc1 vs. b", NBb, 0, max_b);
    epsilon22ddc1b->Sumw2();

    epsilon22ddc2Nch= new TProfile("epsilon22ddc2Nch", "e_{2}^{2}ddc2 vs. Nch", NBch, 0, max_mult);
    epsilon22ddc2Nch->Sumw2();
    epsilon22ddc2Cen= new TProfile("epsilon22ddc2Cen", "e_{2}^{2}ddc2 vs. Centrality", NBcen, 0, max_cen);
    epsilon22ddc2Cen->Sumw2();
    epsilon22ddc2b= new TProfile("epsilon22ddc2b", "e_{2}^{2}ddc2 vs. b", NBb, 0, max_b);
    epsilon22ddc2b->Sumw2();

    epsilon22ddc3Nch= new TProfile("epsilon22ddc3Nch", "e_{2}^{2}ddc3 vs. Nch", NBch, 0, max_mult);
    epsilon22ddc3Nch->Sumw2();
    epsilon22ddc3Cen= new TProfile("epsilon22ddc3Cen", "e_{2}^{2}ddc3 vs. Centrality", NBcen, 0, max_cen);
    epsilon22ddc3Cen->Sumw2();
    epsilon22ddc3b= new TProfile("epsilon22ddc3b", "e_{2}^{2}ddc3 vs. b", NBb, 0, max_b);
    epsilon22ddc3b->Sumw2();

    epsilon22ddc4Nch= new TProfile("epsilon22ddc4Nch", "e_{2}^{2}ddc4 vs. Nch", NBch, 0, max_mult);
    epsilon22ddc4Nch->Sumw2();
    epsilon22ddc4Cen= new TProfile("epsilon22ddc4Cen", "e_{2}^{2}ddc4 vs. Centrality", NBcen, 0, max_cen);
    epsilon22ddc4Cen->Sumw2();
    epsilon22ddc4b= new TProfile("epsilon22ddc4b", "e_{2}^{2}ddc4 vs. b", NBb, 0, max_b);
    epsilon22ddc4b->Sumw2();

    epsilon32ddc1Nch= new TProfile("epsilon32ddc1Nch", "e_{3}^{2}ddc1 vs. Nch", NBch, 0, max_mult);
    epsilon32ddc1Nch->Sumw2();
    epsilon32ddc1Cen= new TProfile("epsilon32ddc1Cen", "e_{3}^{2}ddc1 vs. Centrality", NBcen, 0, max_cen);
    epsilon32ddc1Cen->Sumw2();
    epsilon32ddc1b= new TProfile("epsilon32ddc1b", "e_{3}^{2}ddc1 vs. b", NBb, 0, max_b);
    epsilon32ddc1b->Sumw2();

    epsilon32ddc2Nch= new TProfile("epsilon32ddc2Nch", "e_{3}^{2}ddc2 vs. Nch", NBch, 0, max_mult);
    epsilon32ddc2Nch->Sumw2();
    epsilon32ddc2Cen= new TProfile("epsilon32ddc2Cen", "e_{3}^{2}ddc2 vs. Centrality", NBcen, 0, max_cen);
    epsilon32ddc2Cen->Sumw2();
    epsilon32ddc2b= new TProfile("epsilon32ddc2b", "e_{3}^{2}ddc2 vs. b", NBb, 0, max_b);
    epsilon32ddc2b->Sumw2();

    epsilon32ddc3Nch= new TProfile("epsilon32ddc3Nch", "e_{3}^{2}ddc3 vs. Nch", NBch, 0, max_mult);
    epsilon32ddc3Nch->Sumw2();
    epsilon32ddc3Cen= new TProfile("epsilon32ddc3Cen", "e_{3}^{2}ddc3 vs. Centrality", NBcen, 0, max_cen);
    epsilon32ddc3Cen->Sumw2();
    epsilon32ddc3b= new TProfile("epsilon32ddc3b", "e_{3}^{2}ddc3 vs. b", NBb, 0, max_b);
    epsilon32ddc3b->Sumw2();

    epsilon32ddc4Nch= new TProfile("epsilon32ddc4Nch", "e_{3}^{2}ddc4 vs. Nch", NBch, 0, max_mult);
    epsilon32ddc4Nch->Sumw2();
    epsilon32ddc4Cen= new TProfile("epsilon32ddc4Cen", "e_{3}^{2}ddc4 vs. Centrality", NBcen, 0, max_cen);
    epsilon32ddc4Cen->Sumw2();
    epsilon32ddc4b= new TProfile("epsilon32ddc4b", "e_{3}^{2}ddc4 vs. b", NBb, 0, max_b);
    epsilon32ddc4b->Sumw2();

    epsilon22ddc12Nch= new TProfile("epsilon22ddc12Nch", "e_{2}^{2}ddc12 vs. Nch", NBch, 0, max_mult);
    epsilon22ddc12Nch->Sumw2();
    epsilon22ddc12Cen= new TProfile("epsilon22ddc12Cen", "e_{2}^{2}ddc12 vs. Centrality", NBcen, 0, max_cen);
    epsilon22ddc12Cen->Sumw2();
    epsilon22ddc12b= new TProfile("epsilon22ddc12b", "e_{2}^{2}ddc12 vs. b", NBb, 0, max_b);
    epsilon22ddc12b->Sumw2();

    epsilon22ddc22Nch= new TProfile("epsilon22ddc22Nch", "e_{2}^{2}ddc22 vs. Nch", NBch, 0, max_mult);
    epsilon22ddc22Nch->Sumw2();
    epsilon22ddc22Cen= new TProfile("epsilon22ddc22Cen", "e_{2}^{2}ddc22 vs. Centrality", NBcen, 0, max_cen);
    epsilon22ddc22Cen->Sumw2();
    epsilon22ddc22b= new TProfile("epsilon22ddc22b", "e_{2}^{2}ddc22 vs. b", NBb, 0, max_b);
    epsilon22ddc22b->Sumw2();

    epsilon22ddc32Nch= new TProfile("epsilon22ddc32Nch", "e_{2}^{2}ddc32 vs. Nch", NBch, 0, max_mult);
    epsilon22ddc32Nch->Sumw2();
    epsilon22ddc32Cen= new TProfile("epsilon22ddc32Cen", "e_{2}^{2}ddc32 vs. Centrality", NBcen, 0, max_cen);
    epsilon22ddc32Cen->Sumw2();
    epsilon22ddc32b= new TProfile("epsilon22ddc32b", "e_{2}^{2}ddc32 vs. b", NBb, 0, max_b);
    epsilon22ddc32b->Sumw2();

    epsilon22ddc42Nch= new TProfile("epsilon22ddc42Nch", "e_{2}^{2}ddc42 vs. Nch", NBch, 0, max_mult);
    epsilon22ddc42Nch->Sumw2();
    epsilon22ddc42Cen= new TProfile("epsilon22ddc42Cen", "e_{2}^{2}ddc42 vs. Centrality", NBcen, 0, max_cen);
    epsilon22ddc42Cen->Sumw2();
    epsilon22ddc42b= new TProfile("epsilon22ddc42b", "e_{2}^{2}ddc42 vs. b", NBb, 0, max_b);
    epsilon22ddc42b->Sumw2();

    epsilon32ddc12Nch= new TProfile("epsilon32ddc12Nch", "e_{3}^{2}ddc12 vs. Nch", NBch, 0, max_mult);
    epsilon32ddc12Nch->Sumw2();
    epsilon32ddc12Cen= new TProfile("epsilon32ddc12Cen", "e_{3}^{2}ddc12 vs. Centrality", NBcen, 0, max_cen);
    epsilon32ddc12Cen->Sumw2();
    epsilon32ddc12b= new TProfile("epsilon32ddc12b", "e_{3}^{2}ddc12 vs. b", NBb, 0, max_b);
    epsilon32ddc12b->Sumw2();

    epsilon32ddc22Nch= new TProfile("epsilon32ddc22Nch", "e_{3}^{2}ddc22 vs. Nch", NBch, 0, max_mult);
    epsilon32ddc22Nch->Sumw2();
    epsilon32ddc22Cen= new TProfile("epsilon32ddc22Cen", "e_{3}^{2}ddc22 vs. Centrality", NBcen, 0, max_cen);
    epsilon32ddc22Cen->Sumw2();
    epsilon32ddc22b= new TProfile("epsilon32ddc22b", "e_{3}^{2}ddc22 vs. b", NBb, 0, max_b);
    epsilon32ddc22b->Sumw2();

    epsilon32ddc32Nch= new TProfile("epsilon32ddc32Nch", "e_{3}^{2}ddc32 vs. Nch", NBch, 0, max_mult);
    epsilon32ddc32Nch->Sumw2();
    epsilon32ddc32Cen= new TProfile("epsilon32ddc32Cen", "e_{3}^{2}ddc32 vs. Centrality", NBcen, 0, max_cen);
    epsilon32ddc32Cen->Sumw2();
    epsilon32ddc32b= new TProfile("epsilon32ddc32b", "e_{3}^{2}ddc32 vs. b", NBb, 0, max_b);
    epsilon32ddc32b->Sumw2();

    epsilon32ddc42Nch= new TProfile("epsilon32ddc42Nch", "e_{3}^{2}ddc42 vs. Nch", NBch, 0, max_mult);
    epsilon32ddc42Nch->Sumw2();
    epsilon32ddc42Cen= new TProfile("epsilon32ddc42Cen", "e_{3}^{2}ddc42 vs. Centrality", NBcen, 0, max_cen);
    epsilon32ddc42Cen->Sumw2();
    epsilon32ddc42b= new TProfile("epsilon32ddc42b", "e_{3}^{2}ddc42 vs. b", NBb, 0, max_b);
    epsilon32ddc42b->Sumw2();

    double mult, e2, e3, e4, b;
    double x2, y2, xy, x2y2,energy; 
    double dc1,dc2,dc3,dc4;
    double ddc1,ddc2,ddc3,ddc4;    
    chain.SetBranchAddress("mult", &mult);
    chain.SetBranchAddress("e2", &e2);
    chain.SetBranchAddress("e3", &e3);
    chain.SetBranchAddress("e4", &e4);
    chain.SetBranchAddress("b", &b);
    chain.SetBranchAddress("x2", &x2);
    chain.SetBranchAddress("y2", &y2);
    chain.SetBranchAddress("xy", &xy);
    chain.SetBranchAddress("x2y2", &x2y2);
    chain.SetBranchAddress("energy", &energy);


    for (Long64_t i = 0; i < chain.GetEntries(); i++) {
        chain.GetEntry(i);
        
        double centrality = get_centrality_percent(mult);
        
        epsilon2Nch->Fill(mult, e2);
        epsilon22Nch->Fill(mult, e2 * e2);
        epsilon24Nch->Fill(mult, e2 * e2 * e2 * e2);
        epsilon24Nch->Fill(mult, e2 * e2 * e2 * e2 * e2 * e2);
        epsilon3Nch->Fill(mult, e3);
        epsilon32Nch->Fill(mult, e3 * e3);
        epsilon34Nch->Fill(mult, e3 * e3 * e3 * e3);
        epsilon4Nch->Fill(mult, e4);
        epsilon42Nch->Fill(mult, e4 * e4);
        epsilon44Nch->Fill(mult, e4 * e4 * e4 * e4);
        
        epsilon2232Nch->Fill(mult, e2 *e2 * e3 *e3);
        epsilon2242Nch->Fill(mult, e2 *e2 * e4 *e4);

        if (centrality >= 0) {
            epsilon2Cen->Fill(centrality, e2);
            epsilon22Cen->Fill(centrality, e2 * e2);
            epsilon24Cen->Fill(centrality, e2 * e2 * e2 * e2);
            epsilon26Cen->Fill(centrality, e2 * e2 * e2 * e2 * e2 * e2);
            epsilon3Cen->Fill(centrality, e3);
            epsilon32Cen->Fill(centrality, e3 * e3);
            epsilon34Cen->Fill(centrality, e3 * e3 * e3 * e3);
            epsilon4Cen->Fill(centrality, e4);
            epsilon42Cen->Fill(centrality, e4 * e4);
            epsilon44Cen->Fill(centrality, e4 * e4 * e4 * e4);

            epsilon2232Cen->Fill(centrality, e2 *e2 * e3 *e3);
            epsilon2242Cen->Fill(centrality, e2 *e2 * e4 *e4);

        }
        
        epsilon2b->Fill(b, e2);
        epsilon22b->Fill(b, e2 * e2);
        epsilon24b->Fill(b, e2 * e2 * e2 * e2);
        epsilon24b->Fill(b, e2 * e2 * e2 * e2 * e2 * e2);
        epsilon3b->Fill(b, e3);
        epsilon32b->Fill(b, e3 * e3);
        epsilon34b->Fill(b, e3 * e3 * e3 * e3);
        epsilon4b->Fill(b, e4);
        epsilon42b->Fill(b, e4 * e4);
        epsilon44b->Fill(b, e4 * e4 * e4 * e4);

        epsilon2232b->Fill(b, e2 *e2 * e3 *e3);
        epsilon2242b->Fill(b, e2 *e2 * e4 *e4);

        dc1=0;dc2=0;dc3=0,dc4=0;
        if ((x2*y2-xy*xy)>1e-10)  {dc1=1/sqrt(sqrt(x2*y2-xy*xy));}
        if ((x2*y2)>1e-10)        {dc2=1/sqrt(sqrt(x2*y2));}
        if ((x2y2-xy*xy)>1e-10)   {dc3=1/sqrt(sqrt(x2y2-xy*xy));}
        if ((mult)>1e-10)         {dc4=energy/mult;}

        dc1Nch->Fill(mult, dc1);
        dc1b->Fill(b, dc1);
        dc2Nch->Fill(mult, dc2);
        dc2b->Fill(b, dc2);
        dc3Nch->Fill(mult, dc3);
        dc3b->Fill(b, dc3);
        dc4Nch->Fill(mult, dc4);
        dc4b->Fill(b, dc4);

        if (centrality >= 0) {
            dc1Cen->Fill(centrality, dc1);
            dc2Cen->Fill(centrality, dc2);
            dc3Cen->Fill(centrality, dc3);
            dc4Cen->Fill(centrality, dc4);
        }
    }
    
    chain.ResetBranchAddresses(); // 可选，但有助于确保状态干净
    chain.SetBranchAddress("mult", &mult);
    chain.SetBranchAddress("e2", &e2);
    chain.SetBranchAddress("e3", &e3);
    chain.SetBranchAddress("e4", &e4);
    chain.SetBranchAddress("b", &b);
    chain.SetBranchAddress("x2", &x2);
    chain.SetBranchAddress("y2", &y2);
    chain.SetBranchAddress("xy", &xy);
    chain.SetBranchAddress("x2y2", &x2y2);
    chain.SetBranchAddress("energy", &energy);

    for (Long64_t i = 0; i < chain.GetEntries(); i++) {
        chain.GetEntry(i);
        double centrality = get_centrality_percent(mult);

        dc1=0;dc2=0;dc3=0;
        if ((x2*y2-xy*xy)>1e-10)  {dc1=1/sqrt(sqrt(x2*y2-xy*xy));}
        if ((x2*y2)>1e-10)        {dc2=1/sqrt(sqrt(x2*y2));}
        if ((x2y2-xy*xy)>1e-10)   {dc3=1/sqrt(sqrt(x2y2-xy*xy));}
        if ((mult)>1e-10)         {dc4=energy/mult;}

        int binNumber_dc1Nch = dc1Nch->FindBin(mult);
        double value_dc1Nch = dc1Nch->GetBinContent(binNumber_dc1Nch);
        int binNumber_dc2Nch = dc2Nch->FindBin(mult);
        double value_dc2Nch = dc2Nch->GetBinContent(binNumber_dc2Nch);
        int binNumber_dc3Nch = dc3Nch->FindBin(mult);
        double value_dc3Nch = dc3Nch->GetBinContent(binNumber_dc3Nch);
        int binNumber_dc4Nch = dc4Nch->FindBin(mult);
        double value_dc4Nch = dc4Nch->GetBinContent(binNumber_dc4Nch);

        int binNumber_dc1Cen = dc1Cen->FindBin(centrality);
        double value_dc1Cen = dc1Cen->GetBinContent(binNumber_dc1Cen);
        int binNumber_dc2Cen = dc2Cen->FindBin(centrality);
        double value_dc2Cen = dc2Cen->GetBinContent(binNumber_dc2Cen);
        int binNumber_dc3Cen = dc3Cen->FindBin(centrality);
        double value_dc3Cen = dc3Cen->GetBinContent(binNumber_dc3Cen);
        int binNumber_dc4Cen = dc4Cen->FindBin(centrality);
        double value_dc4Cen = dc4Cen->GetBinContent(binNumber_dc4Cen);

        int binNumber_dc1b = dc1b->FindBin(b);
        double value_dc1b = dc1b->GetBinContent(binNumber_dc1b);
        int binNumber_dc2b = dc2b->FindBin(b);
        double value_dc2b = dc2b->GetBinContent(binNumber_dc2b);
        int binNumber_dc3b = dc3b->FindBin(b);
        double value_dc3b = dc3b->GetBinContent(binNumber_dc3b);
        int binNumber_dc4b = dc4b->FindBin(b);
        double value_dc4b = dc4b->GetBinContent(binNumber_dc4b);

        double value_ddc1Nch = (dc1-value_dc1Nch)/value_dc1Nch;
        double value_ddc1Cen = (dc1-value_dc1Cen)/value_dc1Cen;
        double value_ddc1b   = (dc1-value_dc1b)/value_dc1b;

        double value_ddc2Nch = (dc2-value_dc2Nch)/value_dc2Nch;
        double value_ddc2Cen = (dc2-value_dc2Cen)/value_dc2Cen;
        double value_ddc2b   = (dc2-value_dc2b)/value_dc2b;

        double value_ddc3Nch = (dc3-value_dc3Nch)/value_dc3Nch;
        double value_ddc3Cen = (dc3-value_dc3Cen)/value_dc3Cen;
        double value_ddc3b   = (dc3-value_dc3b)/value_dc3b;

        double value_ddc4Nch = (dc4-value_dc4Nch)/value_dc4Nch;
        double value_ddc4Cen = (dc4-value_dc4Cen)/value_dc4Cen;
        double value_ddc4b   = (dc4-value_dc4b)/value_dc4b;

        ddc12Nch->Fill(mult, pow(value_ddc1Nch, 2.0));
        ddc12b->Fill(b, pow(value_ddc1b, 2.0));
        ddc22Nch->Fill(mult, pow(value_ddc2Nch, 2.0));
        ddc22b->Fill(b, pow(value_ddc2b, 2.0));
        ddc32Nch->Fill(mult, pow(value_ddc3Nch, 2.0));
        ddc32b->Fill(b, pow(value_ddc3b, 2.0));
        ddc42Nch->Fill(mult, pow(value_ddc4Nch, 2.0));
        ddc42b->Fill(b, pow(value_ddc4b, 2.0));

        ddc13Nch->Fill(mult, pow(value_ddc1Nch, 3.0));
        ddc13b->Fill(b, pow(value_ddc1b, 3.0));
        ddc23Nch->Fill(mult, pow(value_ddc2Nch, 3.0));
        ddc23b->Fill(b, pow(value_ddc2b, 3.0));
        ddc33Nch->Fill(mult, pow(value_ddc3Nch, 3.0));
        ddc33b->Fill(b, pow(value_ddc3b, 3.0));
        ddc43Nch->Fill(mult, pow(value_ddc4Nch, 3.0));
        ddc43b->Fill(b, pow(value_ddc4b, 3.0));

        ddc14Nch->Fill(mult, pow(value_ddc1Nch, 4.0));
        ddc14b->Fill(b, pow(value_ddc1b, 4.0));
        ddc24Nch->Fill(mult, pow(value_ddc2Nch, 4.0));
        ddc24b->Fill(b, pow(value_ddc2b, 4.0));
        ddc34Nch->Fill(mult, pow(value_ddc3Nch, 4.0));
        ddc34b->Fill(b, pow(value_ddc3b, 4.0));
        ddc44Nch->Fill(mult, pow(value_ddc4Nch, 4.0));
        ddc44b->Fill(b, pow(value_ddc4b, 4.0));

        epsilon22ddc1Nch->Fill(mult, value_ddc1Nch*e2*e2);
        epsilon22ddc1b->Fill(b, value_ddc1Nch*e2*e2);
        epsilon22ddc2Nch->Fill(mult, value_ddc2Nch*e2*e2);
        epsilon22ddc2b->Fill(b, value_ddc2Nch*e2*e2);
        epsilon22ddc3Nch->Fill(mult, value_ddc3Nch*e2*e2);
        epsilon22ddc3b->Fill(b, value_ddc3Nch*e2*e2);
        epsilon22ddc4Nch->Fill(mult, value_ddc4Nch*e2*e2);
        epsilon22ddc4b->Fill(b, value_ddc4Nch*e2*e2);

        epsilon32ddc1Nch->Fill(mult, value_ddc1Nch*e3*e3);
        epsilon32ddc1b->Fill(b, value_ddc1Nch*e3*e3);
        epsilon32ddc2Nch->Fill(mult, value_ddc2Nch*e3*e3);
        epsilon32ddc2b->Fill(b, value_ddc2Nch*e3*e3);
        epsilon32ddc3Nch->Fill(mult, value_ddc3Nch*e3*e3);
        epsilon32ddc3b->Fill(b, value_ddc3Nch*e3*e3);
        epsilon32ddc4Nch->Fill(mult, value_ddc4Nch*e3*e3);
        epsilon32ddc4b->Fill(b, value_ddc4Nch*e3*e3);

        epsilon22ddc12Nch->Fill(mult, value_ddc1Nch*value_ddc1Nch*e2*e2);
        epsilon22ddc12b->Fill(b, value_ddc1b*value_ddc1b*e2*e2);
        epsilon22ddc22Nch->Fill(mult, value_ddc2Nch*value_ddc2Nch*e2*e2);
        epsilon22ddc22b->Fill(b, value_ddc2b*value_ddc2b*e2*e2);
        epsilon22ddc32Nch->Fill(mult, value_ddc3Nch*value_ddc3Nch*e2*e2);
        epsilon22ddc32b->Fill(b, value_ddc3b*value_ddc3b*e2*e2);
        epsilon22ddc42Nch->Fill(mult, value_ddc4Nch*value_ddc4Nch*e2*e2);
        epsilon22ddc42b->Fill(b, value_ddc4b*value_ddc4b*e2*e2);

        epsilon32ddc12Nch->Fill(mult, value_ddc1Nch*value_ddc1Nch*e3*e3);
        epsilon32ddc12b->Fill(b, value_ddc1b*value_ddc1b*e3*e3);
        epsilon32ddc22Nch->Fill(mult, value_ddc2Nch*value_ddc2Nch*e3*e3);
        epsilon32ddc22b->Fill(b, value_ddc2b*value_ddc2b*e3*e3);
        epsilon32ddc32Nch->Fill(mult, value_ddc3Nch*value_ddc3Nch*e3*e3);
        epsilon32ddc32b->Fill(b, value_ddc3b*value_ddc3b*e3*e3);
        epsilon32ddc42Nch->Fill(mult, value_ddc4Nch*value_ddc4Nch*e3*e3);
        epsilon32ddc42b->Fill(b, value_ddc4b*value_ddc4b*e3*e3);

        if (centrality >= 0) {
            ddc12Cen->Fill(centrality, pow(value_ddc1Cen, 2.0));
            ddc22Cen->Fill(centrality, pow(value_ddc2Cen, 2.0));
            ddc32Cen->Fill(centrality, pow(value_ddc3Cen, 2.0));
            ddc42Cen->Fill(centrality, pow(value_ddc4Cen, 2.0));

            ddc13Cen->Fill(centrality, pow(value_ddc1Cen, 3.0));
            ddc23Cen->Fill(centrality, pow(value_ddc2Cen, 3.0));
            ddc33Cen->Fill(centrality, pow(value_ddc3Cen, 3.0));
            ddc43Cen->Fill(centrality, pow(value_ddc4Cen, 3.0));

            ddc14Cen->Fill(centrality, pow(value_ddc1Cen, 4.0));
            ddc24Cen->Fill(centrality, pow(value_ddc2Cen, 4.0));
            ddc34Cen->Fill(centrality, pow(value_ddc3Cen, 4.0));
            ddc44Cen->Fill(centrality, pow(value_ddc4Cen, 4.0));

            epsilon22ddc1Cen->Fill(centrality, value_ddc1Cen*e2*e2);
            epsilon22ddc2Cen->Fill(centrality, value_ddc2Cen*e2*e2);
            epsilon22ddc3Cen->Fill(centrality, value_ddc3Cen*e2*e2);
            epsilon22ddc4Cen->Fill(centrality, value_ddc4Cen*e2*e2);

            epsilon32ddc1Cen->Fill(centrality, value_ddc1Cen*e3*e3);
            epsilon32ddc2Cen->Fill(centrality, value_ddc2Cen*e3*e3);
            epsilon32ddc3Cen->Fill(centrality, value_ddc3Cen*e3*e3);
            epsilon32ddc4Cen->Fill(centrality, value_ddc4Cen*e3*e3);

            epsilon22ddc12Cen->Fill(centrality, value_ddc1Cen*value_ddc1Cen*e2*e2);
            epsilon22ddc22Cen->Fill(centrality, value_ddc2Cen*value_ddc2Cen*e2*e2);
            epsilon22ddc32Cen->Fill(centrality, value_ddc3Cen*value_ddc3Cen*e2*e2);
            epsilon22ddc42Cen->Fill(centrality, value_ddc4Cen*value_ddc4Cen*e2*e2);

            epsilon32ddc12Cen->Fill(centrality, value_ddc1Cen*value_ddc1Cen*e3*e3);
            epsilon32ddc22Cen->Fill(centrality, value_ddc2Cen*value_ddc2Cen*e3*e3);
            epsilon32ddc32Cen->Fill(centrality, value_ddc3Cen*value_ddc3Cen*e3*e3);
            epsilon32ddc42Cen->Fill(centrality, value_ddc4Cen*value_ddc4Cen*e3*e3);
        }

    }

    TFile* FileOutput = new TFile("out.root", "RECREATE");
    writeHistograms(FileOutput);
    
    FileOutput->Close();
    delete FileOutput;

    deleteHistograms();

    f1->Close();
    f2->Close();
    delete f1;
    delete f2;
}

void writeHistograms(TFile* FileOutput) {
    if (FileOutput && FileOutput->IsOpen()) {
        epsilon2Nch->Write();
        epsilon22Nch->Write();
        epsilon24Nch->Write();
        epsilon26Nch->Write();
        epsilon2Cen->Write();
        epsilon22Cen->Write();
        epsilon24Cen->Write();
        epsilon26Cen->Write();
        epsilon2b->Write();
        epsilon22b->Write();
        epsilon24b->Write();
        epsilon26b->Write();

        epsilon3Nch->Write();
        epsilon32Nch->Write();
        epsilon34Nch->Write();
        epsilon3Cen->Write();
        epsilon32Cen->Write();
        epsilon34Cen->Write();
        epsilon3b->Write();
        epsilon32b->Write();
        epsilon34b->Write();

        epsilon4Nch->Write();
        epsilon42Nch->Write();
        epsilon44Nch->Write();
        epsilon4Cen->Write();
        epsilon42Cen->Write();
        epsilon44Cen->Write();
        epsilon4b->Write();
        epsilon42b->Write();
        epsilon44b->Write();
        

        epsilon2242Nch->Write();
        epsilon2242Cen->Write();
        epsilon2242b->Write();

        epsilon2232Nch->Write();
        epsilon2232Cen->Write();
        epsilon2232b->Write();


        dc1Nch->Write();
        dc1Cen->Write();
        dc1b->Write();
        dc2Nch->Write();
        dc2Cen->Write();
        dc2b->Write();
        dc3Nch->Write();
        dc3Cen->Write();
        dc3b->Write();
        dc4Nch->Write();
        dc4Cen->Write();
        dc4b->Write();

        ddc12Nch->Write();
        ddc12Cen->Write();
        ddc12b->Write();
        ddc22Nch->Write();
        ddc22Cen->Write();
        ddc22b->Write();
        ddc32Nch->Write();
        ddc32Cen->Write();
        ddc32b->Write();
        ddc42Nch->Write();
        ddc42Cen->Write();
        ddc42b->Write();

        ddc13Nch->Write();
        ddc13Cen->Write();
        ddc13b->Write();
        ddc23Nch->Write();
        ddc23Cen->Write();
        ddc23b->Write();
        ddc33Nch->Write();
        ddc33Cen->Write();
        ddc33b->Write();
        ddc43Nch->Write();
        ddc43Cen->Write();
        ddc43b->Write();

        ddc14Nch->Write();
        ddc14Cen->Write();
        ddc14b->Write();
        ddc24Nch->Write();
        ddc24Cen->Write();
        ddc24b->Write();
        ddc34Nch->Write();
        ddc34Cen->Write();
        ddc34b->Write();
        ddc44Nch->Write();
        ddc44Cen->Write();
        ddc44b->Write();

        epsilon22ddc1Nch->Write();
        epsilon22ddc1Cen->Write();
        epsilon22ddc1b->Write();

        epsilon22ddc2Nch->Write();
        epsilon22ddc2Cen->Write();
        epsilon22ddc2b->Write();

        epsilon22ddc3Nch->Write();
        epsilon22ddc3Cen->Write();
        epsilon22ddc3b->Write();

        epsilon22ddc4Nch->Write();
        epsilon22ddc4Cen->Write();
        epsilon22ddc4b->Write();

        epsilon32ddc1Nch->Write();
        epsilon32ddc1Cen->Write();
        epsilon32ddc1b->Write();

        epsilon32ddc2Nch->Write();
        epsilon32ddc2Cen->Write();
        epsilon32ddc2b->Write();

        epsilon32ddc3Nch->Write();
        epsilon32ddc3Cen->Write();
        epsilon32ddc3b->Write();

        epsilon32ddc4Nch->Write();
        epsilon32ddc4Cen->Write();
        epsilon32ddc4b->Write();

        epsilon22ddc12Nch->Write();
        epsilon22ddc12Cen->Write();
        epsilon22ddc12b->Write();

        epsilon22ddc22Nch->Write();
        epsilon22ddc22Cen->Write();
        epsilon22ddc22b->Write();

        epsilon22ddc32Nch->Write();
        epsilon22ddc32Cen->Write();
        epsilon22ddc32b->Write();

        epsilon22ddc42Nch->Write();
        epsilon22ddc42Cen->Write();
        epsilon22ddc42b->Write();

        epsilon32ddc12Nch->Write();
        epsilon32ddc12Cen->Write();
        epsilon32ddc12b->Write();

        epsilon32ddc22Nch->Write();
        epsilon32ddc22Cen->Write();
        epsilon32ddc22b->Write();

        epsilon32ddc32Nch->Write();
        epsilon32ddc32Cen->Write();
        epsilon32ddc32b->Write();

        epsilon32ddc42Nch->Write();
        epsilon32ddc42Cen->Write();
        epsilon32ddc42b->Write();

    }
}

void deleteHistograms() {
    delete epsilon2Nch;
    delete epsilon22Nch;
    delete epsilon24Nch;
    delete epsilon26Nch;
    delete epsilon2Cen;
    delete epsilon22Cen;
    delete epsilon24Cen;
    delete epsilon26Cen;
    delete epsilon2b;
    delete epsilon22b;
    delete epsilon24b;
    delete epsilon26b;

    delete epsilon3Nch;
    delete epsilon32Nch;
    delete epsilon34Nch;
    delete epsilon3Cen;
    delete epsilon32Cen;
    delete epsilon34Cen;
    delete epsilon3b;
    delete epsilon32b;
    delete epsilon34b;

    delete epsilon4Nch;
    delete epsilon42Nch;
    delete epsilon44Nch;
    delete epsilon4Cen;
    delete epsilon42Cen;
    delete epsilon44Cen;
    delete epsilon4b;
    delete epsilon42b;
    delete epsilon44b;
    

    delete epsilon2242Nch;
    delete epsilon2242Cen;
    delete epsilon2242b;
    delete epsilon2232Nch;
    delete epsilon2232Cen;
    delete epsilon2232b;

    delete dc1Nch;
    delete dc1Cen;
    delete dc1b;
    delete dc2Nch;
    delete dc2Cen;
    delete dc2b;
    delete dc3Nch;
    delete dc3Cen;
    delete dc3b;
    delete dc4Nch;
    delete dc4Cen;
    delete dc4b;

    delete ddc12Nch;
    delete ddc12Cen;
    delete ddc12b;
    delete ddc22Nch;
    delete ddc22Cen;
    delete ddc22b;
    delete ddc32Nch;
    delete ddc32Cen;
    delete ddc32b;
    delete ddc42Nch;
    delete ddc42Cen;
    delete ddc42b;

    delete ddc13Nch;
    delete ddc13Cen;
    delete ddc13b;
    delete ddc23Nch;
    delete ddc23Cen;
    delete ddc23b;
    delete ddc33Nch;
    delete ddc33Cen;
    delete ddc33b;
    delete ddc43Nch;
    delete ddc43Cen;
    delete ddc43b;

    delete ddc14Nch;
    delete ddc14Cen;
    delete ddc14b;
    delete ddc24Nch;
    delete ddc24Cen;
    delete ddc24b;
    delete ddc34Nch;
    delete ddc34Cen;
    delete ddc34b;
    delete ddc44Nch;
    delete ddc44Cen;
    delete ddc44b;

    delete epsilon22ddc1Nch;
    delete epsilon22ddc1Cen;
    delete epsilon22ddc1b;
    delete epsilon22ddc2Nch;
    delete epsilon22ddc2Cen;
    delete epsilon22ddc2b;
    delete epsilon22ddc3Nch;
    delete epsilon22ddc3Cen;
    delete epsilon22ddc3b;
    delete epsilon22ddc4Nch;
    delete epsilon22ddc4Cen;
    delete epsilon22ddc4b;

    delete epsilon32ddc1Nch;
    delete epsilon32ddc1Cen;
    delete epsilon32ddc1b;
    delete epsilon32ddc2Nch;
    delete epsilon32ddc2Cen;
    delete epsilon32ddc2b;
    delete epsilon32ddc3Nch;
    delete epsilon32ddc3Cen;
    delete epsilon32ddc3b;
    delete epsilon32ddc4Nch;
    delete epsilon32ddc4Cen;
    delete epsilon32ddc4b;

    delete epsilon22ddc12Nch;
    delete epsilon22ddc12Cen;
    delete epsilon22ddc12b;
    delete epsilon22ddc22Nch;
    delete epsilon22ddc22Cen;
    delete epsilon22ddc22b;
    delete epsilon22ddc32Nch;
    delete epsilon22ddc32Cen;
    delete epsilon22ddc32b;
    delete epsilon22ddc42Nch;
    delete epsilon22ddc42Cen;
    delete epsilon22ddc42b;

    delete epsilon32ddc12Nch;
    delete epsilon32ddc12Cen;
    delete epsilon32ddc12b;
    delete epsilon32ddc22Nch;
    delete epsilon32ddc22Cen;
    delete epsilon32ddc22b;
    delete epsilon32ddc32Nch;
    delete epsilon32ddc32Cen;
    delete epsilon32ddc32b;
    delete epsilon32ddc42Nch;
    delete epsilon32ddc42Cen;
    delete epsilon32ddc42b;

}