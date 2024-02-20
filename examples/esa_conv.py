import matplotlib.pyplot as plt
import numpy as np
import csv

# n = [36, 38, 41, 44, 46, 48, 51]
n = [31]
files = []
h = []
for n_v in n:
    files.append(
        f"/media/evgen/SecondLinuxDisk/4_level/"
        f"Electromagnetic-Waves-Scattering/vtk_files/examples/plane/sigmas/sigma_{n_v}.csv")
    if n_v != n[-1]:
        h.append(1/(n_v-1))
errors = []
sigmas = []

for i in range(len(files)):
    sigmas.append([])
angle = []

with open(files[0]) as csvfile:
    r = csv.reader(csvfile, delimiter=',')
    next(r)
    for line in r:
        angle.append(float(line[1]))

for i in range(len(files)):
    with open(files[i]) as csvfile:
        r = csv.reader(csvfile, delimiter=',')
        next(r)
        for line in r:
            sigmas[i].append(float(line[0]))

fig, ax = plt.subplots()
for i in range(len(files)):
    ax.plot(angle, sigmas[i], label=f"sigma_{i}")

for i in range(len(files) - 1):
    d = np.array(sigmas[i]) - np.array(sigmas[-1])
    errors.append(np.max(np.abs(d)))
# log_h = np.log(h)
# log_e = np.log(errors)
# k, b = np.polyfit(log_h, log_e, 1)
# print(k)
# ax.plot(log_h, log_e)
# ax.plot(log_h, np.polyval((k, b), log_h))
# ax.scatter(log_h, log_e)
# ax.plot(h, errors)
ax.grid()
plt.show()
