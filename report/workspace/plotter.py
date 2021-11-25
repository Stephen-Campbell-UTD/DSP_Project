#%%
import pandas as pd
import matplotlib.pyplot as plt

# %%
df = pd.read_csv('spectrum_cv.csv', sep='	', names=["freq",'freq-units',"power", 'power-units' ])
freq = df.loc[:,"freq"]
power = df.loc[:,"power"]
# %%
fig, ax = plt.subplots()
ax.plot(freq,power)
ax.set_title('Spectrum of IF Out')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('Power [dBm]')
plt.savefig('../images/spectrum_cv.png')
# %%
