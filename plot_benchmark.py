import csv
import matplotlib.pyplot as plt
from collections import defaultdict

data = defaultdict(lambda: defaultdict(list))

with open("benchmark_results_7sizes.csv", "r") as file:
    reader = csv.DictReader(file)

    for row in reader:
        p = int(row["p"])
        version = row["version"]
        time = float(row["total_seconds"])
        data[p][version].append(time)

sizes = sorted(data.keys())

sequential = []
parallel = []

for p in sizes:
    sequential.append(sum(data[p]["sequential"]) / len(data[p]["sequential"]))
    parallel.append(sum(data[p]["parallel"]) / len(data[p]["parallel"]))

plt.figure(figsize=(10, 6))
plt.plot(sizes, sequential, marker="o", label="Sequential")
plt.plot(sizes, parallel, marker="o", label="Parallel (4 threads)")

plt.xscale("log")
plt.yscale("log")

plt.xlabel("Modulus p")
plt.ylabel("Average execution time (seconds)")
plt.title("Index-calculus: Sequential vs Parallel")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()

plt.savefig("benchmark_plot.png", dpi=300)
plt.show()