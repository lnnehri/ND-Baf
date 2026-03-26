import numpy as np
import matplotlib.pyplot as plt

# --- Yardımcı fonksiyonlar ---
def normalize(arr):
    arr = np.array(arr, dtype=float)
    mn, mx = np.min(arr), np.max(arr)
    return (arr - mn) / (mx - mn + 1e-12)

def normalize_inv(arr):
    return 1.0 - normalize(arr)

# ===========================================================
# 1️⃣ FAO DATA
# ===========================================================
groups = np.array([2, 3, 4])

FAO = {
    "ABCD1": [1.103806896, 0, 0],
    "ABCD3": [0, 0, -1.1922394],
    "ACAD11": [1.58162722, 0, 1.92399438],
    "ACOT8": [0, 0, 1.06839345],
    "ACOX2": [1.514986236, 1.250884391, 1.661705853],
    "AMACR": [0, 0, 1.11508322],
    "CRAT": [1.606008218, 1.215030617, 1.82216515],
    "MCEE": [0, 0, 1.07960814],
    "PCCB": [-0.81918828, 0, -1.6613709],
    "PEX13": [0, 0, 1.06203812],
    "PPARD": [0, 0, 1.09904979],
    "SESN2": [1.93837683, 1.141918406, 2.255953786],
    "SLC25A17": [-1.1810531, -1.036624976, -1.388890654],
}

fao_scores = np.zeros(len(groups))
for gene, vals in FAO.items():
    fao_scores += normalize(vals)

# ===========================================================
# 2️⃣ ROS DATA
# ===========================================================
ROS = {
    "AKT1": [0, 0, -1.12537253],
    "ATP7A": [0.855987085, 0, 1.71298208],
    "BCL2": [0, 0, -1.4975039],
    "BNIP3": [0, 1.81056293, 0],
    "CYBA": [0, 0, -1.9489035],
    "CYP1A1": [0, 1.08441409, 0],
    "CYP1B1": [-1.5417554, 0, 2.46853471],
    "DDIT4": [2.601605169, 2.961595164, 3.475183789],
    "DHFR": [-1.42470485, 0, -1.9538248],
    "DUOX2": [1.985578193, 1.343795499, 0],
    "EDN1": [-1.992393124, 0, -3.7231337],
    "EGFR": [1.07123324, 0, 0],
    "EPHX2": [1.21480231, 0, 0.96828596],
    "GCLM": [-1.1789052, 0, 0],
    "GPX1": [-1.080537207, 0, -1.1144963],
    "GPX2": [0, 1.37895444, 0],
    "GPX3": [-1.150785991, 0, -1.2201409],
    "HMOX1": [0, 0, 1.20650602],
    "HSP90AA1": [-1.945546836, -1.275173899, -2.087067788],
    "MAOB": [0, 0, -2.7808074],
    "NOS2": [-2.2962272, 0, -2.4195555],
    "NOS3": [0, 2.96993846, 3.33696695],
    "NOSTRIN": [-1.364725461, 0, -1.5761996],
    "NOX1": [-2.98504751, -1.921708293, -4.35222649],
    "NOXA1": [1.810805135, 0, 0],
    "PARK7": [-0.924837713, 0, -1.1716096],
    "PRDX1": [-1.469746004, 0, -1.8694652],
    "PRDX3": [-0.7885859, 0, -1.1297303],
    "PRDX4": [-0.977240804, 0, -1.3425915],
    "PREX1": [1.322603536, 1.101901035, 0],
    "PRKG2": [1.738295712, 0, 1.21980613],
    "RFK": [0, 0, -1.4934019],
    "RORA": [1.174270541, 0, 1.21207093],
    "SH3PXD2A": [1.78368154, 1.88584993, 2.432630493],
    "SH3PXD2B": [0, 0, -1.8178487],
    "SLC7A2": [1.4209887, 0, 1.2661619],
    "SOD1": [-1.123875745, 0, -1.497028],
    "TXN": [-1.455434486, 0, -1.8217603],
    "TXNIP": [-3.100828773, 0, -1.5381986],
    "WASL": [0, 0, 1.00528093],
}

inducers_ros = {"ATP7A","BNIP3","CYBA","CYP1A1","CYP1B1","DDIT4","DUOX2","EDN1","EGFR","EPHX2",
                "HMOX1","MAOB","NOS2","NOX1","NOXA1","PREX1","RFK","SH3PXD2A","SH3PXD2B","AKT1","NOS3"}
suppressors_ros = {"WASL","TXNIP","SLC7A2","HSP90AA1","BCL2","DHFR","GCLM","GPX1","GPX2","GPX3",
                   "NOSTRIN","PARK7","PRDX1","PRDX3","PRDX4","PRKG2","RORA","SOD1","TXN"}

ros_scores = np.zeros(len(groups))
for gene, vals in ROS.items():
    arr = np.array(vals, dtype=float)
    if gene in inducers_ros:
        ros_scores += normalize(arr)
    elif gene in suppressors_ros:
        ros_scores += normalize_inv(arr)
    else:
        ros_scores += normalize(arr)

# ===========================================================
# 3️⃣ MITOPHAGY DATA
# ===========================================================
MITO = {
    "PINK1":[1.168312444,1.209997078,1.641245948],
    "SRC":[1.170263989,0.863352332,1.238465989],
    "ULK1":[1.89798269,1.552060398,2.596375092],
    "MTERF3":[-1.037074875,0,-1.2189022],
    "PGAM5":[-1.062996008,0,0],
    "TOMM40":[-1.314918505,0,-1.268570071],
    "TOMM5":[-1.564804791,0,-1.2770355],
    "TOMM6":[-1.205514409,0,-0.8967824],
    "UBE2N":[-1.608133738,0,-1.6661424],
    "VDAC1":[-1.002400677,0,-1.4883174],
    "MAP1LC3B":[0,0,1.16986854],
    "OPTN":[0.763209199,1.1236032,1.146394842],
    "SQSTM1":[0,0.94590932,1.467901722],
    "FUNDC1":[0,0,-1.1838776],
}
suppressors_mito = {"SRC"}
inducers_mito = set(MITO.keys()) - suppressors_mito

mito_scores = np.zeros(len(groups))
for gene, vals in MITO.items():
    arr = np.array(vals, dtype=float)
    if gene in suppressors_mito:
        mito_scores += normalize_inv(arr)
    else:
        mito_scores += normalize(arr)

# ===========================================================
# 4️⃣ Plot: Combine FAO, ROS, and Mitophagy
# ===========================================================
fig, ax = plt.subplots(figsize=(7, 4))
x = np.arange(len(groups))
width = 0.25

ax.bar(x - width, fao_scores, width, label='FAO score')
ax.bar(x, ros_scores, width, label='ROS score')
ax.bar(x + width, mito_scores, width, label='Mitophagy score')

ax.set_xlabel('Group', fontsize=11, fontweight='bold')
ax.set_ylabel('Adjusted Score (normalized sum)', fontsize=11, fontweight='bold')
ax.set_title('FAO, ROS, and Mitophagy composite scores by group', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(groups)
ax.legend(frameon=False)

plt.tight_layout()
plt.savefig('combined_FAO_ROS_Mitophagy_scores.png', dpi=300, bbox_inches='tight')
plt.show()
