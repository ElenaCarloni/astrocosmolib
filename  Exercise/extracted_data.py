import pandas as pd
import os

# Creare una versione semplificata della tabella senza colonne extra
data_simplified = [
    ["BGS", 0.295, "7.942 ± 0.075", "13.588 ± 0.168", "21.863 ± 0.425"],
    ["LRG1", 0.510, "12.720 ± 0.099", "17.351 ± 0.177", "19.455 ± 0.330"],
    ["LRG2", 0.706, "16.050 ± 0.110", "17.351 ± 0.177", "19.455 ± 0.330"],
    ["LRG3+ELG1", 0.934, "19.721 ± 0.091", "21.576 ± 0.152", "17.641 ± 0.193"],
    ["ELG2", 1.321, "24.252 ± 0.174", "27.601 ± 0.318", "14.176 ± 0.221"],
    ["QSO", 1.484, "26.055 ± 0.398", "30.512 ± 0.760", "12.817 ± 0.516"],
    ["Lya", 2.330, "31.267 ± 0.256", "38.988 ± 0.531", "8.632 ± 0.101"],
]

# Definire le nuove intestazioni di colonna
columns_simplified = ["Tracer", "z_eff", "D_V / r_d", "D_M / r_d", "D_H / r_d"]
df_simplified = pd.DataFrame(data_simplified, columns=columns_simplified)

# Eliminare la colonna "Tracer"
df_simplified = df_simplified.drop(columns=["Tracer"])

# Separare le colonne valore ± errore
columns_to_split = ["D_V / r_d", "D_M / r_d", "D_H / r_d"]
for col in columns_to_split:
    df_simplified[[col, col + "_err"]] = df_simplified[col].str.split(" ± ", expand=True)
    df_simplified[col] = df_simplified[col].astype(float)
    df_simplified[col + "_err"] = df_simplified[col + "_err"].astype(float)

# Percorso della cartella di destinazione
directory = "/home/elena/Esercitazione"  # Cambia con il tuo percorso

# Creare la cartella se non esiste
os.makedirs(directory, exist_ok=True)

# Percorso completo del file CSV
csv_filename_simplified = os.path.join(directory, "data.csv")

# Salvare il file CSV
df_simplified.to_csv(csv_filename_simplified, index=False)

print(f"File salvato in: {csv_filename_simplified}")
