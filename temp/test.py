import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv("predictions/1LNIA_1LNIA.csv")
df = df[df["pdb_id"]=="1LNIA"]
x = df["pred_ddg"]
y = df["ddg"]
plt.plot(x,y,'*')
plt.show()