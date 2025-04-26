import pandas as pd
# Dosya yolunu belirt
file_path = '/Users/lemannur/Downloads/Expression_Public_24Q4_subsetted.csv'

# Dosyayı oku
df = pd.read_csv(file_path)

# İlk birkaç satırı göster (kontrol amaçlı)
print(df.head())
caco2_df = df[df['cell_line_display_name'] == 'CACO2']
other_cells_df = df[df['cell_line_display_name'] != 'CACO2']
# İstediğin genler ve cell_line_display_name
selected_columns = [
    'cell_line_display_name',
    'BIRC3', 'PPP3CB', 'IL1A', 'PIK3R3', 'IRAK2', 'TP53', 'AKT1', 'PRKACB',
    'ENDOD1', 'ATM', 'BCL2L11', 'CAPN2', 'TP53INP1', 'BCL2L1', 'PIK3CD',
    'BCL2', 'CYCS', 'TNFRSF10D', 'PRKX'
]

# Yeni bir DataFrame oluştur
df_selected = df[selected_columns]
# CACO2 olanlar
caco2_df = df_selected[df_selected['cell_line_display_name'] == 'CACO2']

# Diğer tüm cell-line'lar
other_cells_df = df_selected[df_selected['cell_line_display_name'] != 'CACO2']

from scipy.stats import mannwhitneyu
import pandas as pd

# Gen isimleri: 'cell_line_display_name' dışındakiler
gen_columns = df_selected.columns[1:]

# Sonuçları tutmak için boş liste
results = []

for gene in gen_columns:
    group1 = caco2_df[gene].dropna()
    group2 = other_cells_df[gene].dropna()

    stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')

    results.append({
        'Gene': gene,
        'U_statistic': stat,
        'p_value': p_value
    })

# Sonuçları DataFrame'e çevir
results_df = pd.DataFrame(results)

# p-value küçükten büyüğe sırala
results_df = results_df.sort_values('p_value')

# İlk 10 sonucu göster
print(results_df.head(10))
# p-değeri 0.05'ten küçük olanları filtrele
significant_genes_df = results_df[results_df['p_value'] < 0.05]

# Sonuçları göster
print(significant_genes_df)
print(caco2_df.shape)
print(results_df.sort_values('p_value').head(10))

# Seçtiğimiz gen isimleri (cell_line_display_name hariç)
gen_columns = df_selected.columns[1:]

# CACO2 gen değerleri (tek satır olduğu için squeeze ile seri yapalım)
caco2_values = caco2_df[gen_columns].squeeze()

# Diğer cell-line'ların gen değerlerinin ortalaması
others_mean_values = other_cells_df[gen_columns].mean()

import matplotlib.pyplot as plt
import numpy as np

# Bar plot için ayarlamalar
x = np.arange(len(gen_columns))  # x ekseni için gen isimlerinin sıralı dizisi
width = 0.35  # bar genişliği

fig, ax = plt.subplots(figsize=(14, 6))

# İki grup için çubuk grafiği
rects1 = ax.bar(x - width/2, caco2_values, width, label='CACO2')
rects2 = ax.bar(x + width/2, others_mean_values, width, label='Other Cell Lines Mean')

# Label ve başlıklar
ax.set_ylabel('Expression Level')
ax.set_title('Gene Expression: CACO2 vs Other Cell Lines')
ax.set_xticks(x)
ax.set_xticklabels(gen_columns, rotation=45, ha='right')
ax.legend()

# Layout düzelt
fig.tight_layout()

plt.show()
