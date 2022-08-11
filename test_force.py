import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from sader_jarvis_function import sjarvis_deconvolution


z = np.linspace(1, 0.25, num=512)
sigma = 0.235
E = 0.371
df_off = 4*E*((sigma/z)**13 - (sigma/z)**7)
df_on = df_off + sign


print(len(z))
df = pd.DataFrame({'deltaF': df_on,'Z': z})

#cu_ON = np.array()
#cd_OFF =np.array()


#force = sjarvis_deconvolution(df_ON=cu_ON,df_OFF=cd_OFF,A=10,f0=25000,k=1000)

plt.plot(df["Z"],df["deltaF"])
plt.ylabel("Frequency shift | df [Hz]")
plt.title(f"df vs Z")
plt.xlabel("Z[m]")
plt.show()