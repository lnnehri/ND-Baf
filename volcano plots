import pandas as pd

# Dosya yolu
file_path = '/Users/lemannur/Downloads/ND+Baf-vs-NR-all.xlsx'

# Dosyayı oku
df = pd.read_excel(file_path)

# İlk birkaç satırı göster (kontrol amaçlı)
print(df.head())

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# X ekseni: log2FoldChange
x = df['log2FoldChange']

# Y ekseni: -log10(pvalue)
y = -np.log10(df['pvalue'])

# Şekil ve boyut
plt.figure(figsize=(10, 8))

# Tüm genleri gri renkle çiz
plt.scatter(x, y, color='grey', alpha=0.5)

# Önemli genleri seçelim (cutoff: |log2FC| > 1 ve p-adj < 0.05)
significant = (abs(df['log2FoldChange']) > 1) & (df['padj'] < 0.05)

# Önemli genleri kırmızı renkle göster
plt.scatter(df['log2FoldChange'][significant], -np.log10(df['pvalue'][significant]), color='red', alpha=0.7)

# Eşik çizgileri
plt.axhline(-np.log10(0.05), color='blue', linestyle='dashed')  # p-value = 0.05 çizgisi
plt.axvline(1, color='green', linestyle='dashed')               # log2FC = +1
plt.axvline(-1, color='green', linestyle='dashed')              # log2FC = -1

# Eksen başlıkları ve grafik başlığı
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot: ND+Baf vs NR-all')

# Plot'u göster
plt.show()


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# X ekseni: log2FoldChange
x = df['log2FoldChange']

# Y ekseni: -np.log10(pvalue)
y = -np.log10(df['pvalue'])

# Şekil ve boyut
plt.figure(figsize=(12, 10))

# Tüm genleri gri ile çiz
plt.scatter(x, y, color='grey', alpha=0.5)

# Önemli genleri seçelim
significant = (abs(df['log2FoldChange']) > 1) & (df['padj'] < 0.05)

# Önemli genleri kırmızı renkle göster
plt.scatter(df['log2FoldChange'][significant], -np.log10(df['pvalue'][significant]), color='red', alpha=0.7)

# Eşik çizgileri
plt.axhline(-np.log10(0.05), color='blue', linestyle='dashed')
plt.axvline(1, color='green', linestyle='dashed')
plt.axvline(-1, color='green', linestyle='dashed')

# Gen isimlerini ekleyelim
for i in df[significant].index:
    plt.text(
        df.loc[i, 'log2FoldChange'],
        -np.log10(df.loc[i, 'pvalue']),
        df.index[i],  # Burada gen isimleri index olarak görünüyor
        fontsize=8,
        ha='right'
    )

# Eksen ve başlık
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot with Gene Labels: ND+Baf vs NR-all')

# Layout düzelt
plt.tight_layout()

# Plot'u göster
plt.show()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Eğer daha önce yapmadıysan dosyayı oku
# df = pd.read_excel(file_path)

# X ve Y eksenleri
x = df['log2FoldChange']
y = -np.log10(df['pvalue'])

# Şekil
plt.figure(figsize=(12, 10))

# Tüm genler gri
plt.scatter(x, y, color='grey', alpha=0.5)

# Önemli genler
significant = (abs(df['log2FoldChange']) > 1) & (df['padj'] < 0.05)
plt.scatter(df['log2FoldChange'][significant], -np.log10(df['pvalue'][significant]), color='red', alpha=0.7)

# Eşik çizgileri
plt.axhline(-np.log10(0.05), color='blue', linestyle='dashed')
plt.axvline(1, color='green', linestyle='dashed')
plt.axvline(-1, color='green', linestyle='dashed')

# Gen isimlerini doğru yazalım (Unnamed: 0 kolonundan)
for i in df[significant].index:
    plt.text(
        df.loc[i, 'log2FoldChange'],
        -np.log10(df.loc[i, 'pvalue']),
        df.loc[i, 'Unnamed: 0'],  # <<< Burayı değiştirdik
        fontsize=8,
        ha='right'
    )

# Eksenler ve başlık
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot with Gene Labels: ND+Baf vs NR-all')

plt.tight_layout()
plt.show()


import pandas as pd

# Dosya yolu
file_path = '/Users/lemannur/Downloads/NR+Baf-vs-NR-all.xlsx'

# Dosyayı oku
df_nr = pd.read_excel(file_path)

# İlk birkaç satırı kontrol edelim
print(df_nr.head())

import matplotlib.pyplot as plt
import numpy as np

# X ve Y eksenleri
x = df_nr['log2FoldChange']
y = -np.log10(df_nr['pvalue'])

# Şekil
plt.figure(figsize=(12, 10))

# Tüm genler gri
plt.scatter(x, y, color='grey', alpha=0.5)

# Önemli genleri seçelim
significant = (abs(df_nr['log2FoldChange']) > 1) & (df_nr['padj'] < 0.05)

# Önemli genler kırmızı
plt.scatter(df_nr['log2FoldChange'][significant], -np.log10(df_nr['pvalue'][significant]), color='red', alpha=0.7)

# Eşik çizgileri
plt.axhline(-np.log10(0.05), color='blue', linestyle='dashed')
plt.axvline(1, color='green', linestyle='dashed')
plt.axvline(-1, color='green', linestyle='dashed')

# Gen isimlerini çizelim (m kolonu kullanarak)
for i in df_nr[significant].index:
    plt.text(
        df_nr.loc[i, 'log2FoldChange'],
        -np.log10(df_nr.loc[i, 'pvalue']),
        df_nr.loc[i, 'm'],  # <<< Burada gen ismi 'm' kolonundan çekiliyor
        fontsize=8,
        ha='right'
    )

# Eksen ve başlıklar
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot with Gene Labels: NR+Baf vs NR-all')

plt.tight_layout()
plt.show()

import pandas as pd

# Dosya yolu
file_path = '/Users/lemannur/Downloads/ND-vs-NR-all.xlsx'

# Dosyayı oku
df_nd = pd.read_excel(file_path)

# İlk birkaç satırı kontrol edelim
print(df_nd.head())


import matplotlib.pyplot as plt
import numpy as np

# X ve Y eksenleri
x = df_nd['log2FoldChange']
y = -np.log10(df_nd['pvalue'])

# Şekil
plt.figure(figsize=(12, 10))

# Tüm genler gri
plt.scatter(x, y, color='grey', alpha=0.5)

# Önemli genler (cut-off: |log2FC| > 1 ve padj < 0.05)
significant = (abs(df_nd['log2FoldChange']) > 1) & (df_nd['padj'] < 0.05)

# Önemli genleri kırmızı ile göster
plt.scatter(df_nd['log2FoldChange'][significant], -np.log10(df_nd['pvalue'][significant]), color='red', alpha=0.7)

# Eşik çizgileri
plt.axhline(-np.log10(0.05), color='blue', linestyle='dashed')
plt.axvline(1, color='green', linestyle='dashed')
plt.axvline(-1, color='green', linestyle='dashed')

# Gen isimlerini etiketleyelim
for i in df_nd[significant].index:
    plt.text(
        df_nd.loc[i, 'log2FoldChange'],
        -np.log10(df_nd.loc[i, 'pvalue']),
        df_nd.loc[i, 'Unnamed: 0'],  # Gen ismi doğru yerden
        fontsize=8,
        ha='right'
    )

# Eksen ve başlıklar
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot with Gene Labels: ND vs NR-all')

plt.tight_layout()
plt.show()

import pandas as pd

# Dosya yolu
file_path = '/Users/lemannur/Downloads/ND+Baf-vs-ND-all.xlsx'

# Dosyayı oku
df_ndbaf = pd.read_excel(file_path)

# İlk birkaç satırı kontrol edelim
print(df_ndbaf.head())


import matplotlib.pyplot as plt
import numpy as np

# X ve Y eksenleri
x = df_ndbaf['log2FoldChange']
y = -np.log10(df_ndbaf['pvalue'])

# Şekil
plt.figure(figsize=(12, 10))

# Tüm genler gri
plt.scatter(x, y, color='grey', alpha=0.5)

# Önemli genler (|log2FC| > 1 ve padj < 0.05)
significant = (abs(df_ndbaf['log2FoldChange']) > 1) & (df_ndbaf['padj'] < 0.05)

# Önemli genler kırmızı
plt.scatter(df_ndbaf['log2FoldChange'][significant], -np.log10(df_ndbaf['pvalue'][significant]), color='red', alpha=0.7)

# Eşik çizgileri
plt.axhline(-np.log10(0.05), color='blue', linestyle='dashed')
plt.axvline(1, color='green', linestyle='dashed')
plt.axvline(-1, color='green', linestyle='dashed')

# Gen isimlerini etiketleyelim
for i in df_ndbaf[significant].index:
    plt.text(
        df_ndbaf.loc[i, 'log2FoldChange'],
        -np.log10(df_ndbaf.loc[i, 'pvalue']),
        df_ndbaf.loc[i, 'Unnamed: 0'],  # Gen isimleri burada
        fontsize=8,
        ha='right'
    )

# Eksenler ve başlık
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot with Gene Labels: ND+Baf vs ND-all')

plt.tight_layout()
plt.show()
