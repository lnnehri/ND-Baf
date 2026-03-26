from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import VariableElimination

# Model oluştur
model = BayesianNetwork([
    ('Starvation', 'Ca2+ Sinyalleşmesi'),
    ('Lysosomal Alkalinization', 'Ca2+ Sinyalleşmesi'),
    ('Ca2+ Sinyalleşmesi', 'Hücre Fenotipi')
])

# CPD'leri tanımla
cpd_starvation = TabularCPD(variable='Starvation', variable_card=2, values=[[0.7], [0.3]])
cpd_alkalinization = TabularCPD(variable='Lysosomal Alkalinization', variable_card=2, values=[[0.6], [0.4]])
cpd_ca2 = TabularCPD(
    variable='Ca2+ Sinyalleşmesi', variable_card=2,
    values=[[0.8, 0.6, 0.7, 0.5], [0.2, 0.4, 0.3, 0.5]],
    evidence=['Starvation', 'Lysosomal Alkalinization'],
    evidence_card=[2, 2]
)
cpd_hucre_fenotipi = TabularCPD(
    variable='Hücre Fenotipi', variable_card=3,
    values=[
        [0.6, 0.2],  # Rounded
        [0.3, 0.7],  # Elongated
        [0.1, 0.1]   # Dead
    ],
    evidence=['Ca2+ Sinyalleşmesi'],
    evidence_card=[2]
)


# CPD'leri ekle
model.add_cpds(cpd_starvation, cpd_alkalinization, cpd_ca2, cpd_hucre_fenotipi)

# Modeli kontrol et
assert model.check_model(), "Model geçerli değil!"
print("Model is valid!")

# Çıkarım motoru
inference = VariableElimination(model)

# Çıkarım sorgusu: Starvation ve Lysosomal Alkalinization durumunda Hücre Fenotipi
result = inference.query(variables=['Hücre Fenotipi'], evidence={'Starvation': 1, 'Lysosomal Alkalinization': 1})
print(result)


###bu kompleks model:

from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import VariableElimination

# Modeli oluştur
model = BayesianNetwork([
    ('Starvation', 'ROS'),
    ('Starvation', 'Lipid Oksidasyonu'),
    ('Lysosomal Alkalinization', 'Ca2+ Sinyalleşmesi'),
    ('Ca2+ Sinyalleşmesi', 'Hücre Fenotipi'),
    ('ROS', 'Hücre Fenotipi'),
    ('Lipid Oksidasyonu', 'Hücre Fenotipi')
])

# CPD'leri tanımla
cpd_starvation = TabularCPD(variable='Starvation', variable_card=2, values=[[0.7], [0.3]])
cpd_alkalinization = TabularCPD(variable='Lysosomal Alkalinization', variable_card=2, values=[[0.6], [0.4]])

cpd_ros = TabularCPD(
    variable='ROS', variable_card=2,
    values=[[0.8, 0.5], [0.2, 0.5]],
    evidence=['Starvation'],
    evidence_card=[2]
)

cpd_lipid = TabularCPD(
    variable='Lipid Oksidasyonu', variable_card=2,
    values=[[0.7, 0.4], [0.3, 0.6]],
    evidence=['Starvation'],
    evidence_card=[2]
)

cpd_ca2 = TabularCPD(
    variable='Ca2+ Sinyalleşmesi', variable_card=2,
    values=[[0.6, 0.3], [0.4, 0.7]],
    evidence=['Lysosomal Alkalinization'],
    evidence_card=[2]
)

cpd_hucre_fenotipi = TabularCPD(
    variable='Hücre Fenotipi', variable_card=3,
    values=[
        [0.3, 0.2, 0.3, 0.2, 0.3, 0.2, 0.2, 0.1],  # Rounded
        [0.2, 0.3, 0.3, 0.5, 0.3, 0.4, 0.5, 0.6],  # Elongated
        [0.5, 0.5, 0.4, 0.3, 0.4, 0.4, 0.3, 0.3]   # Dead
    ],
    evidence=['Ca2+ Sinyalleşmesi', 'ROS', 'Lipid Oksidasyonu'],
    evidence_card=[2, 2, 2]
)


# CPD'leri modele ekle
model.add_cpds(cpd_starvation, cpd_alkalinization, cpd_ros, cpd_lipid, cpd_ca2, cpd_hucre_fenotipi)

# Modeli kontrol et
assert model.check_model(), "Model geçerli değil!"
print("Model is valid!")

# Çıkarım motoru
inference = VariableElimination(model)

# Çıkarım sorgusu: Starvation = 1, Lysosomal Alkalinization = 1 durumunda Hücre Fenotipi
result = inference.query(variables=['Hücre Fenotipi'], evidence={'Starvation': 1, 'Lysosomal Alkalinization': 1})
print(result)


####obur sorgu
# Çıkarım motoru
inference = VariableElimination(model)

# Kanıtlarla Hücre Fenotipi olasılıklarını sorgula
result = inference.query(
    variables=['Hücre Fenotipi'],
    evidence={'Starvation': 1, 'Lysosomal Alkalinization': 0}
)
print(result)
# Bayes ağındaki tüm ilişkileri yazdır
print("Bayes Ağı İlişkileri:")
for edge in model.edges():
    print(f"{edge[0]} → {edge[1]}")

# Bayes ağındaki tüm CPD'leri yazdır
print("\nBayes Ağı CPD'leri:")
for cpd in model.get_cpds():
    print(f"\n{cpd.variable} için CPD:")
    print(cpd)
pip install prettytable
from prettytable import PrettyTable

# İlişkiler tablosu
relation_table = PrettyTable(["Parent", "Child"])
for edge in model.edges():
    relation_table.add_row(edge)

# CPD tablosu
cpd_table = PrettyTable(["Variable", "CPD Table"])
for cpd in model.get_cpds():
    cpd_table.add_row([cpd.variable, cpd.values])

# Tabloyu yazdır
print("Bayes Ağı İlişkileri:")
print(relation_table)
print("\nBayes Ağı CPD'leri:")
print(cpd_table)

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Hücre Fenotipi CPD'sini al
cpd = model.get_cpds('Hücre Fenotipi')

# Çok boyutlu CPD'yi düz tabloya çevir
cpd_flat = cpd.values.reshape(cpd.variable_card, -1)

# DataFrame'e dönüştür
columns = [f"Combo_{i}" for i in range(cpd_flat.shape[1])]
cpd_df = pd.DataFrame(cpd_flat, columns=columns, index=['Rounded', 'Elongated', 'Dead'])

# Isı haritası oluştur
plt.figure(figsize=(12, 6))
sns.heatmap(cpd_df, annot=True, cmap="Blues", fmt=".2f")
plt.title("Hücre Fenotipi Olasılıkları")
plt.ylabel("Hücre Fenotipi")
plt.xlabel("Kanıt Kombinasyonları")
plt.show()

#########ing versiyonunu yaziyorum

from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import VariableElimination

# Create the model
model = BayesianNetwork([
    ('Starvation', 'Ca2+ Signaling'),
    ('Lysosomal Alkalinization', 'Ca2+ Signaling'),
    ('Ca2+ Signaling', 'Cell Phenotype')
])

# Define CPDs
cpd_starvation = TabularCPD(variable='Starvation', variable_card=2, values=[[0.7], [0.3]])
cpd_alkalinization = TabularCPD(variable='Lysosomal Alkalinization', variable_card=2, values=[[0.6], [0.4]])
cpd_ca2 = TabularCPD(
    variable='Ca2+ Signaling', variable_card=2,
    values=[[0.8, 0.6, 0.7, 0.5], [0.2, 0.4, 0.3, 0.5]],
    evidence=['Starvation', 'Lysosomal Alkalinization'],
    evidence_card=[2, 2]
)
cpd_cell_phenotype = TabularCPD(
    variable='Cell Phenotype', variable_card=3,
    values=[
        [0.6, 0.2],  # Rounded
        [0.3, 0.7],  # Elongated
        [0.1, 0.1]   # Dead
    ],
    evidence=['Ca2+ Signaling'],
    evidence_card=[2]
)

# Add CPDs to the model
model.add_cpds(cpd_starvation, cpd_alkalinization, cpd_ca2, cpd_cell_phenotype)

# Validate the model
assert model.check_model(), "The model is invalid!"
print("Model is valid!")

# Inference engine
inference = VariableElimination(model)

# Query: Cell Phenotype given Starvation and Lysosomal Alkalinization
result = inference.query(variables=['Cell Phenotype'], evidence={'Starvation': 1, 'Lysosomal Alkalinization': 1})
print(result)

### Complex model:

from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import VariableElimination

# Create the model
model = BayesianNetwork([
    ('Starvation', 'ROS'),
    ('Starvation', 'Lipid Oxidation'),
    ('Lysosomal Alkalinization', 'Ca2+ Signaling'),
    ('Ca2+ Signaling', 'Cell Phenotype'),
    ('ROS', 'Cell Phenotype'),
    ('Lipid Oxidation', 'Cell Phenotype')
])

# Define CPDs
cpd_starvation = TabularCPD(variable='Starvation', variable_card=2, values=[[0.7], [0.3]])
cpd_alkalinization = TabularCPD(variable='Lysosomal Alkalinization', variable_card=2, values=[[0.6], [0.4]])

cpd_ros = TabularCPD(
    variable='ROS', variable_card=2,
    values=[[0.8, 0.5], [0.2, 0.5]],
    evidence=['Starvation'],
    evidence_card=[2]
)

cpd_lipid = TabularCPD(
    variable='Lipid Oxidation', variable_card=2,
    values=[[0.7, 0.4], [0.3, 0.6]],
    evidence=['Starvation'],
    evidence_card=[2]
)

cpd_ca2 = TabularCPD(
    variable='Ca2+ Signaling', variable_card=2,
    values=[[0.6, 0.3], [0.4, 0.7]],
    evidence=['Lysosomal Alkalinization'],
    evidence_card=[2]
)

cpd_cell_phenotype = TabularCPD(
    variable='Cell Phenotype', variable_card=3,
    values=[
        [0.3, 0.2, 0.3, 0.2, 0.3, 0.2, 0.2, 0.1],  # Rounded
        [0.2, 0.3, 0.3, 0.5, 0.3, 0.4, 0.5, 0.6],  # Elongated
        [0.5, 0.5, 0.4, 0.3, 0.4, 0.4, 0.3, 0.3]   # Dead
    ],
    evidence=['Ca2+ Signaling', 'ROS', 'Lipid Oxidation'],
    evidence_card=[2, 2, 2]
)

# Add CPDs to the model
model.add_cpds(cpd_starvation, cpd_alkalinization, cpd_ros, cpd_lipid, cpd_ca2, cpd_cell_phenotype)

# Validate the model
assert model.check_model(), "The model is invalid!"
print("Model is valid!")

# Inference engine
inference = VariableElimination(model)

# Query: Cell Phenotype given Starvation = 1 and Lysosomal Alkalinization = 1
result = inference.query(variables=['Cell Phenotype'], evidence={'Starvation': 1, 'Lysosomal Alkalinization': 1})
print(result)

#### Another query
# Inference engine
inference = VariableElimination(model)

# Query probabilities of Cell Phenotype given evidence
result = inference.query(
    variables=['Cell Phenotype'],
    evidence={'Starvation': 1, 'Lysosomal Alkalinization': 0}
)
print(result)

# Print all relationships in the Bayes Network
print("Bayes Network Relationships:")
for edge in model.edges():
    print(f"{edge[0]} → {edge[1]}")

# Print all CPDs in the Bayes Network
print("\nBayes Network CPDs:")
for cpd in model.get_cpds():
    print(f"\nCPD for {cpd.variable}:")
    print(cpd)

pip install prettytable
from prettytable import PrettyTable

# Relationship table
relation_table = PrettyTable(["Parent", "Child"])
for edge in model.edges():
    relation_table.add_row(edge)

# CPD table
cpd_table = PrettyTable(["Variable", "CPD Table"])
for cpd in model.get_cpds():
    cpd_table.add_row([cpd.variable, cpd.values])

# Print tables
print("Bayes Network Relationships:")
print(relation_table)
print("\nBayes Network CPDs:")
print(cpd_table)

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Get the CPD for Cell Phenotype
cpd = model.get_cpds('Cell Phenotype')

# Flatten multi-dimensional CPD into a table
cpd_flat = cpd.values.reshape(cpd.variable_card, -1)

# Convert to DataFrame
columns = [f"Combo_{i}" for i in range(cpd_flat.shape[1])]
cpd_df = pd.DataFrame(cpd_flat, columns=columns, index=['Rounded', 'Elongated', 'Dead'])

# Create a heatmap
plt.figure(figsize=(12, 6))
sns.heatmap(cpd_df, annot=True, cmap="Blues", fmt=".2f")
plt.title("Cell Phenotype Probabilities")
plt.ylabel("Cell Phenotype")
plt.xlabel("Evidence Combinations")
plt.show()
