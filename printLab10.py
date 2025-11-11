import pandas as pd
import matplotlib.pyplot as plt
import sys

data = pd.read_csv("./simulation_data_20_300.000000_0_0.300000.csv")
L = 400

plt.figure(figsize=(10, 4))
plt.plot(data['time'], data['co_ad_count'], label='co_ad_count')
plt.ylim(0, L)
plt.xlabel('Time')
plt.ylabel('Count')
plt.title('co_ad_count vs Time')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 4))
plt.plot(data['time'], data['o_ad_count'], label='o_ad_count')
plt.ylim(0, L)
plt.xlabel('Time')
plt.ylabel('Count')
plt.title('o_ad_count vs Time')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 4))
plt.plot(data['time'], data['o_v_count'], label='o_v_count')
plt.ylim(0, L)
plt.xlabel('Time')
plt.ylabel('Count')
plt.title('o_v_count vs Time')
plt.legend()
plt.grid(True)
plt.show()