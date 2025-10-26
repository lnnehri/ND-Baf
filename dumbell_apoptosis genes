import numpy as np, pandas as pd, matplotlib.pyplot as plt

# --- Veriyi oluştur ---
data = {
    "group": [2, 3, 4],
    "BCL2L11": [1.64245599172252, 0.944772245547861, 1.71307243782403],
    "TP53INP1": [1.61370266713904, 1.41775168183201, 1.74967782125782],
    "BIRC3": [-1.29146768268436, 0, -1.07753546657281],
    "CYCS": [-1.54260215746871, 0, 0],
    "TNFRSF10D": [-1.04784716269947, 0, 0],
    "PIK3CD": [0, 1.41383716088679, 0],
    "PIK3R3": [0, 0, 1.91423388178196],
    "BCL2": [0, 0, -1.49750385118025],
    "BCL2L1": [0, 0, -1.61987175055399],
    "CAPN2": [0, 0, -1.19188452024524],
    "PPP3CB": [0, 0, -1.59985351463316],
    "PRKX": [-0.716682839874494, 0, -1.73740842525519],
    "ATM": [0, 0, 1.01540458680992],
    "IL1A": [0, 0, 1.84156605899359],
    "IRAK2": [0, 0, 1.11311074859038],
    "PRKACB": [0, 0, 1.28136556618915],
    "TP53": [0, 0, 1.1814594944741],
    "AKT1": [0, 0, -1.12537252983054],
}

df_raw = pd.DataFrame(data)
group_map = {2: "ND", 3: "Baf", 4: "ND+Baf"}

# Geni satıra, koşulları sütuna çevir
df_tmp = df_raw.set_index("group")
df_tmp.index = df_tmp.index.map(group_map)
df = df_tmp.T  # index: Gene, columns: ["ND","Baf","ND+Baf"]

# Çizim için: ND+Baf’e göre sırala (istersen ND veya Baf'a göre de sıralayabilirsin)
plot_df = df.sort_values("ND+Baf").copy()

plt.figure(figsize=(9, max(4, 0.5*len(plot_df))))
y = np.arange(len(plot_df))

# İki parça çizgi: ND→Baf ve Baf→ND+Baf
plt.hlines(y, plot_df["ND"], plot_df["Baf"], color="lightgray", linewidth=1.5)
plt.hlines(y, plot_df["Baf"], plot_df["ND+Baf"], color="gray", linewidth=1.5)

# Üç ajan noktası
plt.scatter(plot_df["ND"], y, s=60, edgecolor="black", label="ND")
plt.scatter(plot_df["Baf"], y, s=60, marker="s", edgecolor="black", label="Baf")
plt.scatter(plot_df["ND+Baf"], y, s=70, marker="D", edgecolor="black", label="ND+Baf")

plt.axvline(0, color="black", lw=.8)
plt.yticks(y, plot_df.index)
plt.xlabel("log2FC"); plt.ylabel("Gene")
plt.title("Üç ajan karşılaştırması: ND, Baf, ND+Baf", fontsize=12, fontweight="bold")
plt.legend(frameon=False, ncol=3, loc="lower right")
plt.tight_layout()
plt.savefig("apoptosis_three_agent_dumbbell.png", dpi=300)
plt.show()
