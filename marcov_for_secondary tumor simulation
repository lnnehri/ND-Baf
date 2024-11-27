###3 asamada, birinde primary tumor daha az uruyor, 2.sinde daha fazla uruyor, sonuncu da daha fazla uruyor ama daha uzun vadede uruyor, bunlarin kodu:


####death consider edilirse: bunu al: 1. simulasyon: burada primary tumor daha az buyuyor, genelde acliktan oluyor
import numpy as np
import matplotlib.pyplot as plt

# Geçiş Matrisi
P = np.array([
    [0.74, 0.15, 0.0, 0.11],  # PT -> PT, M, ST, Death
    [0.1, 0.6, 0.3, 0.0],     # M -> PT, M, ST, Death
    [0.0, 0.22, 0.7, 0.08],   # ST -> PT, M, ST, Death
    [0.0, 0.0, 0.0, 1.0]      # Death -> Death (absorbing state)
])

# Başlangıç Durumu (Hücre Sayısı)
pi_0 = np.array([100, 0, 0, 0])  # Başlangıçta 100 hücre PT'de, ölüm yok

# Çoğalma oranları
growth_rates = np.array([1, 0, 2, 0])  # PT hücreleri 1, ST hücreleri 2 birim çoğalır, ölümde çoğalma yok

# Zaman adımları ve simülasyon
steps = 10
cell_counts = np.zeros((steps, 4))

for step in range(steps):
    # Geçiş matrisiyle bir adım ilerleme
    pi_0 = np.dot(pi_0, P)
    # Çoğalmayı ekle (ölüm hariç)
    pi_0[:3] += growth_rates[:3] * pi_0[:3]
    # Sonuçları kaydet (gerçek hücre sayıları)
    cell_counts[step] = pi_0

# Zaman adımları
time = np.arange(1, steps + 1)

# Oranları hesapla
proportions = cell_counts / cell_counts.sum(axis=1, keepdims=True)
# Görselleştirme: Hücre Popülasyonu Oranları
plt.figure(figsize=(10, 6))
plt.plot(time, proportions[:, 0], label='Primary Tumor (PT)', marker='o')
plt.plot(time, proportions[:, 1], label='Metastasis (M)', marker='s')
plt.plot(time, proportions[:, 2], label='Secondary Tumor (ST)', marker='^')
plt.plot(time, proportions[:, 3], label='Death', marker='x')

plt.title("Metastatic Cycle Dynamics with Death")
plt.xlabel("Time Steps")
plt.ylabel("Proportion of Cells")
plt.legend()
plt.grid()
plt.savefig("metastatic_cycle_proportions.png", dpi=300)  # PNG olarak kaydet
plt.show()

# Görselleştirme: Toplam Hücre Sayısı
plt.figure(figsize=(10, 6))
plt.plot(time, cell_counts.sum(axis=1), label='Total Cells', marker='o', color='purple')
plt.title("Total Cell Population Over Time")
plt.xlabel("Time Steps")
plt.ylabel("Total Cell Count")
plt.legend()
plt.grid()
plt.savefig("total_cell_population.png", dpi=300)  # PNG olarak kaydet
plt.show()

# Görselleştirme: Ölen Hücrelerin Sayısı
plt.figure(figsize=(10, 6))
plt.plot(time, cell_counts[:, 3], label='Dead Cells', marker='x', color='red')
plt.title("Dead Cell Count Over Time")
plt.xlabel("Time Steps")
plt.ylabel("Number of Dead Cells")
plt.legend()
plt.grid()
plt.savefig("dead_cell_count.png", dpi=300)  # PNG olarak kaydet
plt.show()


#####################2.simulasyon:

###son - her durumda primary tumor proportionu azalacak.
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

# Geçiş Matrisi
P = np.array([
    [0.74, 0.15, 0.0, 0.11],  # PT -> PT, M, ST, Death
    [0.1, 0.6, 0.3, 0.0],     # M -> PT, M, ST, Death
    [0.0, 0.22, 0.7, 0.08],   # ST -> PT, M, ST, Death
    [0.0, 0.0, 0.0, 1.0]      # Death -> Death (absorbing state)
])

# Başlangıç Durumu (Hücre Sayısı)
pi_0 = np.array([10000, 0, 0, 0])  # Başlangıçta 10,000 hücre PT'de

# Çoğalma oranları
growth_rates = np.array([10, 0, 12, 0])  # PT hücreleri daha hızlı büyüsün

# Zaman adımları ve simülasyon
steps = 20  # Daha uzun vadede dengeyi görmek için adım sayısını artırdık
cell_counts = np.zeros((steps, 4))

for step in range(steps):
    # Geçiş matrisiyle bir adım ilerleme
    pi_0 = np.dot(pi_0, P)
    # Çoğalmayı ekle (ölüm hariç)
    pi_0[:3] += growth_rates[:3] * pi_0[:3]
    # Sonuçları kaydet (gerçek hücre sayıları)
    cell_counts[step] = pi_0

# Zaman adımları
time = np.arange(1, steps + 1)

# Oranları hesapla
proportions = cell_counts / cell_counts.sum(axis=1, keepdims=True)
# Görselleştirme: Hücre Popülasyonu Oranları
plt.figure(figsize=(10, 6))
plt.plot(time, proportions[:, 0], label='Primary Tumor (PT)', marker='o')
plt.plot(time, proportions[:, 1], label='Metastasis (M)', marker='s')
plt.plot(time, proportions[:, 2], label='Secondary Tumor (ST)', marker='^')
plt.plot(time, proportions[:, 3], label='Death', marker='x')

plt.title("Metastatic Cycle Dynamics with Enhanced Growth")
plt.xlabel("Time Steps")
plt.ylabel("Proportion of Cells")
plt.legend()
plt.grid()
plt.savefig("proportion_of_cells_primary_reduction.png", dpi=300)  # Farklı isimle kaydet
plt.show()

# Görselleştirme: Toplam Hücre Sayısı
plt.figure(figsize=(10, 6))
plt.plot(time, cell_counts.sum(axis=1), label='Total Cells', marker='o', color='purple')
plt.title("Total Cell Population Over Time")
plt.xlabel("Time Steps")
plt.ylabel("Total Cell Count")
plt.legend()
plt.grid()
plt.savefig("total_cell_population_primary_reduction.png", dpi=300)  # Farklı isimle kaydet
plt.show()

# Görselleştirme: Ölen Hücrelerin Sayısı
plt.figure(figsize=(10, 6))
plt.plot(time, cell_counts[:, 3], label='Dead Cells', marker='x', color='red')
plt.title("Dead Cell Count Over Time")
plt.xlabel("Time Steps")
plt.ylabel("Number of Dead Cells")
plt.legend()
plt.grid()
plt.savefig("dead_cell_count_primary_reduction.png", dpi=300)  # Farklı isimle kaydet
plt.show()


#################3.simulasyon:


####daha uzun simulasyon kodu:
import numpy as np
import matplotlib.pyplot as plt

# Geçiş Matrisi
P = np.array([
    [0.74, 0.15, 0.0, 0.11],  # PT -> PT, M, ST, Death
    [0.1, 0.6, 0.3, 0.0],     # M -> PT, M, ST, Death
    [0.0, 0.22, 0.7, 0.08],   # ST -> PT, M, ST, Death
    [0.0, 0.0, 0.0, 1.0]      # Death -> Death (absorbing state)
])

# Başlangıç Durumu (Hücre Sayısı)
pi_0 = np.array([10000, 0, 0, 0])  # Başlangıçta 10,000 hücre PT'de

# Çoğalma oranları
growth_rates = np.array([10, 0, 12, 0])  # PT hücreleri daha hızlı büyüsün

# Zaman adımları ve simülasyon
steps = 100  # Daha uzun vadede dinamikleri görmek için 100 adım
cell_counts = np.zeros((steps, 4))

for step in range(steps):
    # Geçiş matrisiyle bir adım ilerleme
    pi_0 = np.dot(pi_0, P)
    # Çoğalmayı ekle (ölüm hariç)
    pi_0[:3] += growth_rates[:3] * pi_0[:3]
    # Sonuçları kaydet (gerçek hücre sayıları)
    cell_counts[step] = pi_0

# Zaman adımları
time = np.arange(1, steps + 1)

# Oranları hesapla
proportions = cell_counts / cell_counts.sum(axis=1, keepdims=True)
# Görselleştirme: Hücre Popülasyonu Oranları
plt.figure(figsize=(10, 6))
plt.plot(time, proportions[:, 0], label='Primary Tumor (PT)', marker='o')
plt.plot(time, proportions[:, 1], label='Metastasis (M)', marker='s')
plt.plot(time, proportions[:, 2], label='Secondary Tumor (ST)', marker='^')
plt.plot(time, proportions[:, 3], label='Death', marker='x')

plt.title("Metastatic Cycle Dynamics with Enhanced Growth (Long-Term)")
plt.xlabel("Time Steps")
plt.ylabel("Proportion of Cells")
plt.legend()
plt.grid()
plt.savefig("long_term_proportion_of_cells.png", dpi=300)  # Farklı isimle kaydet
plt.show()

# Görselleştirme: Toplam Hücre Sayısı
plt.figure(figsize=(10, 6))
plt.plot(time, cell_counts.sum(axis=1), label='Total Cells', marker='o', color='purple')
plt.title("Total Cell Population Over Time (Long-Term)")
plt.xlabel("Time Steps")
plt.ylabel("Total Cell Count")
plt.legend()
plt.grid()
plt.savefig("long_term_total_cell_population.png", dpi=300)  # Farklı isimle kaydet
plt.show()

# Görselleştirme: Ölen Hücrelerin Sayısı
plt.figure(figsize=(10, 6))
plt.plot(time, cell_counts[:, 3], label='Dead Cells', marker='x', color='red')
plt.title("Dead Cell Count Over Time (Long-Term)")
plt.xlabel("Time Steps")
plt.ylabel("Number of Dead Cells")
plt.legend()
plt.grid()
plt.savefig("long_term_dead_cell_count.png", dpi=300)  # Farklı isimle kaydet
plt.show()

