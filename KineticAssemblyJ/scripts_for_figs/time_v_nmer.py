import matplotlib.pyplot as plt

# Use a style to enhance the plot aesthetics
plt.style.use('seaborn-white')

# Data
n_sizes = [3, 4, 5, 6, 7, 8]
time_taken = [24.572/1000,105.021/1000,494.365/1000,2.321,10.951,53.637]
plt.figure(figsize=(8, 6))

# Plot with markers
plt.plot(n_sizes, time_taken, marker='o', linestyle='-',lw=5,markersize=15)

# Labels and title
plt.xlabel('N size',size=20)
plt.ylabel('Time (S)',size=20)
plt.title('Time per Iteration vs Nmer Size',size=20)
plt.xticks(size=20)
plt.yticks(size=20)


plt.tight_layout()
plt.savefig('time_v_n.png')