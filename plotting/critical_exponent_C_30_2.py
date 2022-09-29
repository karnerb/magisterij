import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
from scipy import interpolate
from matplotlib import rc

rc("text", usetex=True)


def load_data(filename):
    data = pd.read_csv(filename, delimiter=" ", header=0, skiprows=1)
    data['T'] = 1 / data['beta']
    f = open(filename)
    n = f.readline()
    n = int(n[2:])

    return n, data


def average_over_cycles(data):
    data['Energy_squared'] = data['Energy'] * data['Energy']
    avg = data.groupby(by=['beta'], as_index=False).mean()
    avg['T'] = 1 / avg['beta']
    betas = 1 / avg.index
    return avg


def std_over_cycles(data):
    std = data.groupby(by=['beta'], as_index=False).std()
    std['T'] = 1 / std['beta']
    return std


def plot(data):
    # n, data = load_data(filename)
    # calcP2s(data)

    avg = average_over_cycles(data)

    # plt.xlabel('T')
    # plt.ylabel('P2, P2m')
    # plt.plot(avg['T'], avg['P2'])
    # plt.plot(avg['T'], avg['P2m'])
    # plt.show()
    plt.subplot(3, 1, 1)
    plt.xlabel('T')
    plt.ylabel('E')
    plt.plot(avg['T'], avg['Energy'])
    plt.scatter(data['T'], data['Energy'], s=2)
    plt.plot(avg['T'], avg['Energy'])
    plt.subplot(3, 1, 2)
    plt.xlabel('T')
    plt.ylabel('Cv')
    plt.plot(avg['T'], n ** 3 * (-avg['Energy'] ** 2 + avg['Energy_squared']) * avg['beta'] ** 2)
    plt.subplot(3, 1, 3)
    plt.xlabel('T')
    plt.ylabel('Cv')
    plt.plot(avg['T'], n ** 3 * (-avg['Energy'] ** 2 + avg['Energy_squared']) * avg['beta'] ** 2)
    # plt.plot(avg['T'], np.gradient(avg['Energy'], avg['T']))


def calcP2(data):
    split = [y for x, y in data.groupby(['beta'], as_index=False)]

    betas = []
    P2 = np.zeros((3, 3))
    P2s = []
    P2ms = []
    Nsample = 500
    for group in split:

        betas.append(group['beta'][group.index.min()])
        randints = random.sample(set(group.index.values), Nsample)
        p2temp = []
        p2mtemp = []
        for i in range(Nsample):
            P2[0, 0] = group['P2_00'][randints[i]]
            P2[0, 1] = group['P2_01'][randints[i]]
            P2[0, 2] = group['P2_02'][randints[i]]

            P2[1, 0] = group['P2_10'][randints[i]]
            P2[1, 1] = group['P2_11'][randints[i]]
            P2[1, 2] = group['P2_12'][randints[i]]

            P2[2, 0] = group['P2_20'][randints[i]]
            P2[2, 1] = group['P2_21'][randints[i]]
            P2[2, 2] = group['P2_22'][randints[i]]

            eigvals, eigvecs = np.linalg.eig(P2)
            eigvals = np.sort(eigvals)
            p2temp.append(eigvals[-1])
            p2mtemp.append(eigvals[1])
        P2s.append(np.average(p2temp))
        P2ms.append(np.average(p2mtemp))
    # plt.plot(1/np.array(betas), P2s)
    # plt.plot(1/np.array(betas), P2ms)
    # plt.show()
    return 1 / np.array(betas), P2s, P2ms


files = [  # "test.txt"
    # "cooldown_zoom_cluster.txt",
    # "heatup_zoom_cluster.txt",
    "cluster30cooldown.txt",
    "cluster30heatup.txt",
    # "BW_cooldown.txt",
    # "BW_HEATUP.txt",
    # "heatup_zoom.txt"
    # "00heatup.txt",
    # "01heatup.txt",
    # "02heatup.txt",
    # "03heatup.txt",
    # "04heatup.txt",
    # "05heatup.txt",
    # "06heatup.txt",
    # "07heatup.txt",
    # "08heatup.txt",
    # "09heatup.txt",
    # "00cooldown.txt",
    # "01cooldown.txt",
    # "02cooldown.txt",
    # "03cooldown.txt",
    # "04cooldown.txt",
    # "05cooldown.txt",
    # "06cooldown.txt",
    # "07cooldown.txt",
    # "08cooldown.txt",
    # "09cooldown.txt"
]

folder = "../output/"

dataset = []
dataset_joined = []

labels = [ "ohlajanje cluster",
           "segrevanje cluster",
           "ohlajanje Barker-Watts",
           "segrevanje Barker-Watts",
          "$r=0.0$ heatup",
          "$r=0.1$ heatup",
          "$r=0.2$ heatup",
          "$r=0.3$ heatup",
          "$r=0.4$ heatup",
          "$r=0.5$ heatup",
          "$r=0.6$ heatup",
          "$r=0.7$ heatup",
          "$r=0.8$ heatup",
          "$r=0.9$ heatup",
          "$r=0.0$ cooldown",
          "$r=0.1$ cooldown",
          "$r=0.2$ cooldown",
          "$r=0.3$ cooldown",
          "$r=0.4$ cooldown"
          "$r=0.0$ cooldown",
          "$r=0.1$ cooldown",
          "$r=0.2$ cooldown",
          "$r=0.3$ cooldown",
          "$r=0.4$ cooldown"
          ]
# join two heatup and cooldown data

for file in files:
    n, data = load_data(folder + file)
    data = data[data['beta'].between(1/1.1235, 1/1.11)]
    dataset.append(data)




avgs = []
stds = []
P2s = []
P2ms = []
Ts = []

avgs_joined = []
stds_joined = []
P2s_joined = []
P2ms_joined = []
Ts_joined = []

# calculate averages, standard deviations and P2, P2m
for data in dataset:
    avgs.append(average_over_cycles(data))
    stds.append(std_over_cycles(data))
    T, P2, P2m = calcP2(data)
    P2s.append(P2)
    P2ms.append(P2m)
    Ts.append(T)





print(len(avgs))
#exit()
    # plt.errorbar(avgs[i]['T'], avgs[i]['polar_order'], stds[i]['polar_order'], label=labels[i])

# plot separated data

# plt.scatter(avgs[0]['T'], avgs[0]['Energy'], label=labels[0], s=12, c = 'blue')
# plt.scatter(avgs[1]['T'], avgs[1]['Energy'], label=labels[1], s=12, c='red')
# plt.scatter(avgs[2]['T'], avgs[2]['Energy'], label=labels[2], s=12, c = 'green')
# plt.scatter(avgs[3]['T'], avgs[3]['Energy'], label=labels[3], s=12, c='black')

    # plt.errorbar(avgs[i]['T'], avgs[i]['polar_order'], stds[i]['polar_order'], label=labels[i])


# plt.title('$N=30$')
# plt.xlabel('$T$')
# plt.ylabel('$E$')
# plt.legend()
# plt.show()

Cs_joined = n ** 3 * (-avgs[0]['Energy'] ** 2 + avgs[0]['Energy_squared']) * avgs[0]['beta'] ** 2
for i in range(1, 2):
    Cs_joined = np.concatenate((Cs_joined, n ** 3 * (-avgs[i]['Energy'] ** 2 + avgs[i]['Energy_squared']) * avgs[i]['beta'] ** 2))

#plt.scatter(np.concatenate((avgs[0]['T'],avgs[1]['T'],avgs[2]['T'],avgs[3]['T'])), Cs_joined)
#plt.show()

# plt.scatter(avgs[0]['T'], n ** 3 * (-avgs[0]['Energy'] ** 2 + avgs[0]['Energy_squared']) * avgs[0]['beta'] ** 2 )
# plt.show()


plt.scatter(np.log(-avgs[0]['T']+1.1235), np.log(n ** 3 * (-avgs[0]['Energy'] ** 2 + avgs[0]['Energy_squared']) * avgs[0]['beta'] ** 2), c='blue', s=12)
plt.scatter(np.log(-avgs[1]['T']+1.1235), np.log(n ** 3 * (-avgs[1]['Energy'] ** 2 + avgs[1]['Energy_squared']) * avgs[1]['beta'] ** 2), c='red', s=12)
# plt.scatter(np.log(avgs[2]['T']-1.1235), np.log(n ** 3 * (-avgs[2]['Energy'] ** 2 + avgs[2]['Energy_squared']) * avgs[2]['beta'] ** 2))
# plt.scatter(np.log(avgs[3]['T']-1.1235), np.log(n ** 3 * (-avgs[3]['Energy'] ** 2 + avgs[3]['Energy_squared']) * avgs[3]['beta'] ** 2))


m, b = np.polyfit(np.log(-np.concatenate((avgs[0]['T'],avgs[1]['T']))+1.1235), np.log(Cs_joined), 1)
print(m,b)
plt.plot(np.log(-avgs[0]['T']+1.1235), m*np.log(-avgs[0]['T']+1.1235)+b, label = r'$\tilde{\alpha}' + ' = {0:.2f}$'.format(-m), c='black')
plt.xlabel('$\log(T-T_{NI})$')
plt.ylabel('$\log(C_V)$')
plt.title(r'$C_V\propto (T_{NI}-T)^{-\tilde{\alpha}}$')
plt.legend()
plt.show()

for i in range(len(avgs)-1):
    plt.scatter(avgs[i]['T'], n ** 3 * (-avgs[i]['Energy'] ** 2 + avgs[i]['Energy_squared']) * avgs[i]['beta'] ** 2,
                label=labels[i], s=12, c = 'blue')
    plt.scatter(avgs[i+1]['T'], n ** 3 * (-avgs[i+1]['Energy'] ** 2 + avgs[i+1]['Energy_squared']) * avgs[i+1]['beta'] ** 2,
                label=labels[i+1], s=12, c='red')


plt.title('$N=30$, cluster algorithm')
plt.xlabel('$T$')
plt.ylabel('$C_v$')
plt.legend()
plt.show()


for i in range(len(P2s_joined)):
    plt.plot(Ts_joined[i], P2s_joined[i], c='black')
    plt.plot(Ts_joined[i], P2ms_joined[i], c='black')

for i in range(len(P2s)-1):
    plt.scatter(np.array(Ts[i]), P2s[i], s=12, c='blue')
    plt.scatter(np.array(Ts[i]), P2ms[i], label=labels[i], s=12, c='blue')
    plt.scatter(np.array(Ts[i+1]), P2s[i+1], s=12, c='red')
    plt.scatter(np.array(Ts[i+1]), P2ms[i+1], label=labels[i+1], s=12, c='red')
    # plt.plot(np.array(Ts[i]), P2ms[i], label=labels[i])

plt.title('$N=30$, cluster algorithm')
plt.xlabel('$T$')
plt.ylabel('$P_2$, $P_{2m}$')
plt.legend()
plt.show()






