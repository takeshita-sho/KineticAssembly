import matplotlib.pyplot as plt

# Use a style to enhance the plot aesthetics
plt.style.use('seaborn-white')

# Data
n_sizes = [3, 4, 5, 6, 7, 8]
peak_memory_usage = [2012.78125, 2012.78125, 2012.78125, 2037.1640625, 2544.71484375, 4266.52734375]

plt.figure(figsize=(10, 6))

# Plot with markers
plt.plot(n_sizes, peak_memory_usage, marker='o', linestyle='-', color='b')

# Labels and title
plt.xlabel('N size')
plt.ylabel('Peak Memory Usage (MB)')
plt.title('Memory Usage vs Nmer Size')


plt.tight_layout()
plt.savefig('memory_v_n.png')