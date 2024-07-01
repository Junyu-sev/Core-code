

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
pd.options.mode.chained_assignment = None  # default='warn'
from lifelines import KaplanMeierFitter
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#设置输入路径，需要有每个蛋白的Cutoff存储路径，蛋白浓度路径，疾病定义路径
cutoffpath = '/share/home/zhangjunyu/Project/IHD_MetS_Proteomic_analysis/Fdr-corrected-pval/Result/KM_plot/Cut_off_determination/Non_MetS_cutoff_HR.csv'
propath = '/share/home/zhangjunyu/Project/IHD_MetS_Proteomic_analysis/Fdr-corrected-pval/Data/Proteomics/231008_data_olink_instance_0.csv'
outcomepath = '/share/home/zhangjunyu/Project/IHD_MetS_Proteomic_analysis/Fdr-corrected-pval/Data/Disease_outcomes/IHD/Non_MetS.csv'

#设置输出路径
outpath = '/share/home/zhangjunyu/Project/IHD_MetS_Proteomic_analysis/Fdr-corrected-pval/Result/KM_plot/KMplot_A_Z/Non_MetS/'

#读取cutoff文件，将HR_p_val变成科学计数法
cutoff_df = pd.read_csv(cutoffpath)
cutoff_df = cutoff_df.sort_values(by='Pro_code')

def scientific_notation_to_latex(value):
    if value == 1.0:  # 处理特殊情况，避免输出1*10^0
        return "1.0"
    exponent = int(f"{value:e}".split('e')[1])
    base = f"{value / (10 ** exponent):.1f}"
    return f"${base} \\times 10^{{{exponent}}}$"

cutoff_df['p_val_new_form'] = cutoff_df['HR_p_val'].apply(scientific_notation_to_latex)

pro_f_lst = cutoff_df.Pro_code.tolist()        #蛋白列表
cut_f_lst = cutoff_df.cutoff_result.tolist()   #cutoff列表
riskdir_f_lst = cutoff_df.HR.tolist()          #HR列表
pval_f_lst = cutoff_df.p_val_new_form.tolist() #p_val列表
hrci_f_lst = cutoff_df.HR_ci.tolist()          #HR置信区间列表

#读取蛋白浓度和疾病定义文件
pro_df = pd.read_csv(propath, usecols=['eid'] + pro_f_lst)
target_df = pd.read_csv(outcomepath, usecols = ['eid', 'target_y', 'BL2Target_yrs'])
mydf = pd.merge(target_df, pro_df, how = 'inner', on = 'eid')

#按字母顺序输出成pdf文件
letter_to_elements = {}
for item in pro_f_lst:
    first_letter = item[0]
    if first_letter not in letter_to_elements:
        letter_to_elements[first_letter] = []
    letter_to_elements[first_letter].append(item)

alphabet = list(letter_to_elements.keys())

i = -1

for letter in alphabet:
    pdf_pages = PdfPages(outpath + 'output_' + letter + '.pdf')
    for element in letter_to_elements[letter]:
        i = i +1
        pro_f, f_cut, risk_dir = pro_f_lst[i], np.round(cut_f_lst[i], 2), riskdir_f_lst[i]
        f_hrci, f_pval = hrci_f_lst[i], pval_f_lst[i]
        plotdf = mydf[['eid', 'target_y', 'BL2Target_yrs'] + [pro_f]]
        plotdf.rename(columns={pro_f: 'target_pro'}, inplace=True)
        rm_idx = plotdf.index[plotdf.target_pro.isnull() == True]
        plotdf = plotdf.drop(rm_idx, axis=0)
        plotdf.reset_index(inplace=True)
        high_risk = (plotdf.target_pro > f_cut)
        fig, ax = plt.subplots(figsize=(7.5, 6))
        kmf_l = KaplanMeierFitter()
        kmf_l.fit(durations=plotdf.BL2Target_yrs[~high_risk], event_observed=plotdf.target_y[~high_risk])
        kmf_l.plot_survival_function(ax=ax, color='#2b8cbe', linewidth=3)
        kmf_h = KaplanMeierFitter()
        kmf_h.fit(durations=plotdf.BL2Target_yrs[high_risk], event_observed=plotdf.target_y[high_risk])
        kmf_h.plot_survival_function(ax=ax, color='#f46d43', linewidth=3)
        add_at_risk_counts(kmf_l, kmf_h, labels=['Low protein level', 'High proetin level'], ax=ax, rows_to_show=['At risk', 'Events'],xticks=np.arange(0, 18, 2))
        ax.set_title(pro_f, fontsize=22)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.set_xlim(-0.3, 16.3)
        ax.set_ylim(0.685, 1.015)
        ax.set_xticks(np.arange(0, 18, 2))
        ax.set_yticks([0.7, 0.8, 0.9, 1.0])
        ax.set_xticklabels(np.arange(0, 18, 2), rotation=0)
        ax.set_yticklabels([0.7, 0.8, 0.9, 1], rotation=0)
        ax.set_xlabel('Timeline', fontsize=18)
        ax.set_ylabel('Survival Probability', fontsize=18)
        ax.legend_.remove()
        text1 = "HR = " + str(f_hrci)
        text2 = "P-value = " + str(f_pval)
        plt.text(0.035, 0.26, text1, fontsize=18, color='black', transform=ax.transAxes)
        plt.text(0.035, 0.12, text2, fontsize=18, color='black', transform=ax.transAxes)
        plt.subplots_adjust(left=0.2, bottom=0.2)
        fig.tight_layout()
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['ps.fonttype'] = 42
        pdf_pages.savefig()
        plt.clf()
    pdf_pages.close()
