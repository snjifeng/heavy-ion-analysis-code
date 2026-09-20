#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <utility>

// 这个函数用于从 TTree 中获取数据，计算中心度区间的多重性阈值
void calculate_centrality() {
    // 1. 打开你的 ROOT 文件
    TFile* f = TFile::Open("initial.root");
    if (!f || f->IsZombie()) {
        std::cerr << "错误：无法打开文件 'out.root'。请检查文件是否存在。" << std::endl;
        return;
    }

    // 2. 从文件中获取 TTree
    TTree* tree = (TTree*)f->Get("event_tree");
    if (!tree) {
        std::cerr << "错误：在文件中找不到名为 'event_tree' 的 TTree。" << std::endl;
        f->Close();
        return;
    }

    // 3. 创建一个 TH1D 直方图来存储多重性分布
    // 这里的 binning (1000, 0, 1000) 需要根据你的 mult 数据范围调整
    TH1D* h_mult = new TH1D("RefMulti", "Multiplicity Distribution", 3000, 0, 3000);

    // 4. 遍历 TTree，用 'mult' 分支的数据填充直方图
    // "goff" 选项表示不立即绘图
    tree->Draw("mult>>RefMulti", "", "goff"); 

    if (h_mult->GetEntries() <= 0) {
        std::cerr << "错误：直方图为空，无法计算中心度。" << std::endl;
        f->Close();
        delete h_mult;
        return;
    }

    // 5. 定义中心度百分比的边界
    // 注意：这里的百分比与物理上的习惯（0% 最中心）相符
    std::vector<double> percentile_boundaries = {
        0.0, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
        15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0
    };
    
    // 6. 计算直方图的总积分（总事件数）
    double total_events = h_mult->Integral(1, h_mult->GetNbinsX(), "width");
    
    if (total_events <= 0) {
        std::cerr << "错误：直方图积分小于或等于0，无法计算中心度。" << std::endl;
        f->Close();
        delete h_mult;
        return;
    }

    // 7. 遍历直方图的 bin，计算累积积分并找到阈值
    // 现在从多重性最高的 bin 开始向低多重性遍历 (物理上越中心 mult 越高)
    double cumulative_events = 0.0;
    std::map<double, double> thresholds;
    size_t boundary_index = 0;
    
    for (int i = h_mult->GetNbinsX(); i >= 1; --i) {
        cumulative_events += h_mult->GetBinContent(i);
        double current_percentage = (cumulative_events / total_events) * 100.0;
        
        if (boundary_index < percentile_boundaries.size()) {
            double target_percent = percentile_boundaries[boundary_index];
            
            // 如果当前累积百分比超过或等于目标百分比，记录阈值
            if (current_percentage >= target_percent) {
                thresholds[target_percent] = h_mult->GetBinCenter(i);
                boundary_index++;
            }
        }
    }
    
    // 8. 打印每个中心度区间的阈值
    std::cout << "\n根据多重性直方图计算的中心度阈值（从最中心到最外围）:" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    
    for (size_t i = 0; i < percentile_boundaries.size() - 1; ++i) {
        double lower_percent = percentile_boundaries[i];
        double upper_percent = percentile_boundaries[i+1];
        
        if (thresholds.count(lower_percent) && thresholds.count(upper_percent)) {
            std::cout << "中心度 " << static_cast<int>(lower_percent) << "-" << static_cast<int>(upper_percent) << "%："
                      << " mult 从 " << thresholds[lower_percent] << " 到 " << thresholds[upper_percent] << std::endl;
        } else if (thresholds.count(lower_percent) && upper_percent == 100.0) {
            std::cout << "中心度 " << static_cast<int>(lower_percent) << "-" << static_cast<int>(upper_percent) << "%："
                      << " mult 从 " << thresholds[lower_percent] << " 到 0 (或更低)" << std::endl;
        } else {
             std::cout << "警告: 未找到中心度 " << static_cast<int>(lower_percent) << "-" << static_cast<int>(upper_percent) << "% 的阈值。" << std::endl;
        }
        
    }

    for (size_t i = 0; i < percentile_boundaries.size() - 1; ++i) {
        double lower_percent = percentile_boundaries[i];
        double upper_percent = percentile_boundaries[i+1];
    if (thresholds.count(lower_percent) && thresholds.count(upper_percent)) {
            std::cout <<" {" << static_cast<int>(lower_percent) << ", " << thresholds[lower_percent] <<"},"<< std::endl;
        }
    }
    
    // 清理
    f->Close();

    delete f;
}