import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
from scipy import interpolate
from matplotlib import rc

rc("text", usetex=True)
#
#
# def load_data(filename):
#     data = pd.read_csv(filename, delimiter=" ", header=0, skiprows=1)
#     data['T'] = 1 / data['beta']
#     f = open(filename)
#     n = f.readline()
#     n = int(n[2:])
#
#     return n, data
#
#
# def average_over_cycles(data):
#     data['Energy_squared'] = data['Energy'] * data['Energy']
#     avg = data.groupby(by=['beta'], as_index=False).mean()
#     avg['T'] = 1 / avg['beta']
#     betas = 1 / avg.index
#     return avg
#
#
# def std_over_cycles(data):
#     std = data.groupby(by=['beta'], as_index=False).std()
#     std['T'] = 1 / std['beta']
#     return std
#
#
# def plot(data):
#     # n, data = load_data(filename)
#     # calcP2s(data)
#
#     avg = average_over_cycles(data)
#
#     # plt.xlabel('T')
#     # plt.ylabel('P2, P2m')
#     # plt.plot(avg['T'], avg['P2'])
#     # plt.plot(avg['T'], avg['P2m'])
#     # plt.show()
#     plt.subplot(3, 1, 1)
#     plt.xlabel('T')
#     plt.ylabel('E')
#     plt.plot(avg['T'], avg['Energy'])
#     plt.scatter(data['T'], data['Energy'], s=2)
#     plt.plot(avg['T'], avg['Energy'])
#     plt.subplot(3, 1, 2)
#     plt.xlabel('T')
#     plt.ylabel('Cv')
#     plt.plot(avg['T'], n ** 3 * (-avg['Energy'] ** 2 + avg['Energy_squared']) * avg['beta'] ** 2)
#     plt.subplot(3, 1, 3)
#     plt.xlabel('T')
#     plt.ylabel('Cv')
#     plt.plot(avg['T'], n ** 3 * (-avg['Energy'] ** 2 + avg['Energy_squared']) * avg['beta'] ** 2)
#     # plt.plot(avg['T'], np.gradient(avg['Energy'], avg['T']))
#
#
# def calcP2(data):
#     split = [y for x, y in data.groupby(['beta'], as_index=False)]
#
#     betas = []
#     P2 = np.zeros((3, 3))
#     P2s = []
#     P2ms = []
#     Nsample = 500
#     for group in split:
#
#         betas.append(group['beta'][group.index.min()])
#         randints = random.sample(set(group.index.values), Nsample)
#         p2temp = []
#         p2mtemp = []
#         for i in range(Nsample):
#             P2[0, 0] = group['P2_00'][randints[i]]
#             P2[0, 1] = group['P2_01'][randints[i]]
#             P2[0, 2] = group['P2_02'][randints[i]]
#
#             P2[1, 0] = group['P2_10'][randints[i]]
#             P2[1, 1] = group['P2_11'][randints[i]]
#             P2[1, 2] = group['P2_12'][randints[i]]
#
#             P2[2, 0] = group['P2_20'][randints[i]]
#             P2[2, 1] = group['P2_21'][randints[i]]
#             P2[2, 2] = group['P2_22'][randints[i]]
#
#             eigvals, eigvecs = np.linalg.eig(P2)
#             eigvals = np.sort(eigvals)
#             p2temp.append(eigvals[-1])
#             p2mtemp.append(eigvals[1])
#         P2s.append(np.average(p2temp))
#         P2ms.append(np.average(p2mtemp))
#     # plt.plot(1/np.array(betas), P2s)
#     # plt.plot(1/np.array(betas), P2ms)
#     # plt.show()
#     return 1 / np.array(betas), P2s, P2ms
#
# file_pairs = [["00heatup.txt",
#                "00cooldown.txt"],
#               ["01heatup.txt",
#                "01cooldown.txt"],
#               ["02heatup.txt",
#                "02cooldown.txt"],
#               ["03heatup.txt",
#                "03cooldown.txt"],
#               ["04heatup.txt",
#                "04cooldown.txt"],
#               ["05heatup.txt",
#                "05cooldown.txt"],
#               ["06heatup.txt",
#                "06cooldown.txt"],
#               ["07heatup.txt",
#                "07cooldown.txt"],
#               ["08heatup.txt",
#                "08cooldown.txt"],
#               ["09heatup.txt",
#                "09cooldown.txt"]
#               ]
#
# folder = "../output/"
#
# dataset = []
# dataset_joined = []
#
# labels_joined=["$r=0.0$",
#                "$r=0.1$",
#                "$r=0.2$",
#                "$r=0.3$",
#                "$r=0.4$",
#                "$r=0.5$",
#                "$r=0.6$",
#                "$r=0.7$",
#                "$r=0.8$",
#                "$r=0.9$",
#
# ]
#
# labels = [ "$r=0.0$ ohlajanje",
#            "$r=0.0$ segrevanje",
#            "$r=0.1$ ohlajanje",
#            "$r=0.1$ segrevanje",
#           "$r=0.3$ heatup",
#           "$r=0.4$ heatup",
#           "$r=0.5$ heatup",
#           "$r=0.6$ heatup",
#           "$r=0.7$ heatup",
#           "$r=0.8$ heatup",
#           "$r=0.9$ heatup",
#           "$r=0.0$ cooldown",
#           "$r=0.1$ cooldown",
#           "$r=0.2$ cooldown",
#           "$r=0.3$ cooldown",
#           "$r=0.4$ cooldown"
#           "$r=0.0$ cooldown",
#           "$r=0.1$ cooldown",
#           "$r=0.2$ cooldown",
#           "$r=0.3$ cooldown",
#           "$r=0.4$ cooldown"
#           ]
#
# prop_cycle = plt.rcParams['axes.prop_cycle']
# colors = prop_cycle.by_key()['color']
#
# j=0
# TsofMaxC=[]
# for files in file_pairs:
#     dataset = []
#     dataset_joined = []
#     for file in files:
#         n, data = load_data(folder + file)
#         dataset.append(data)
#
#     for i in range(len(dataset)//2):
#         dataset_joined.append(pd.concat([dataset[2*i], dataset[2*i+1]], ignore_index=True, sort=False))
#
#     avgs = []
#     stds = []
#     P2s = []
#     P2ms = []
#     Ts = []
#
#     avgs_joined = []
#     stds_joined = []
#     P2s_joined = []
#     P2ms_joined = []
#     Ts_joined = []
#
#     # calculate averages, standard deviations and P2, P2m
#     for data in dataset:
#         avgs.append(average_over_cycles(data))
#         stds.append(std_over_cycles(data))
#         T, P2, P2m = calcP2(data)
#         P2s.append(P2)
#         P2ms.append(P2m)
#         Ts.append(T)
#
#     for data in dataset_joined:
#         avgs_joined.append(average_over_cycles(data))
#         stds_joined.append(std_over_cycles(data))
#         T, P2, P2m = calcP2(data)
#         P2s_joined.append(P2)
#         P2ms_joined.append(P2m)
#         Ts_joined.append(T)
#
#
#
#     for i in range(len(avgs_joined)):
#         plt.figure(1)
#         plt.plot(avgs_joined[i]['T'], avgs_joined[i]['Energy'], label = labels_joined[j], color = colors[j])
#
#
#
#     # plt.errorbar(avgs[i]['T'], avgs[i]['polar_order'], stds[i]['polar_order'], label=labels[i])
#
# # plot separated data
# # for i in range(len(avgs)-1):
# #     plt.scatter(avgs[i]['T'], avgs[i]['Energy'], label=labels[i], s=12, c = 'blue')
# #     plt.scatter(avgs[i+1]['T'], avgs[i+1]['Energy'], label=labels[i+1], s=12, c='red')
# #     # plt.errorbar(avgs[i]['T'], avgs[i]['polar_order'], stds[i]['polar_order'], label=labels[i])
#
#
#     plt.title('$N=10$')
#     plt.xlabel('$T$')
#     plt.ylabel('$E$')
#     plt.legend()
#     # plt.show()
#
# # for i in range(len(avgs)-1):
# #     plt.scatter(avgs[i]['T'], n ** 3 * (-avgs[i]['Energy'] ** 2 + avgs[i]['Energy_squared']) * avgs[i]['beta'] ** 2,
# #                 label=labels[i], s=12, c = 'blue')
# #     plt.scatter(avgs[i+1]['T'], n ** 3 * (-avgs[i+1]['Energy'] ** 2 + avgs[i+1]['Energy_squared']) * avgs[i+1]['beta'] ** 2,
# #                 label=labels[i+1], s=12, c='red')
#
#     for i in range(len(avgs_joined)):
#         plt.figure(2)
#         plt.scatter(avgs_joined[i]['T'], n ** 3 * (-avgs_joined[i]['Energy'] ** 2 + avgs_joined[i]['Energy_squared']) * avgs_joined[i]['beta'] ** 2, label=labels_joined[j], color = colors[j], s = 12)
#         TsofMaxC.append(avgs_joined[i]['T'][np.argmax(n ** 3 * (-avgs_joined[i]['Energy'] ** 2 + avgs_joined[i]['Energy_squared']) * avgs_joined[i]['beta'] ** 2)])
#     plt.title('$N=10$')
#     plt.xlabel('$T$')
#     plt.ylabel('$C_v$')
#     plt.legend()
#     # plt.show()
#
#
#
#
#     for i in range(len(P2s_joined)):
#         plt.figure(3)
#         plt.plot(Ts_joined[i], P2s_joined[i], color = colors[j], label=labels_joined[j])
#         plt.plot(Ts_joined[i], P2ms_joined[i], color = colors[j])
#
# # for i in range(len(P2s)-1):
# #     plt.scatter(np.array(Ts[i]), P2s[i], s=12, c='blue')
# #     plt.scatter(np.array(Ts[i]), P2ms[i], label=labels[i], s=12, c='blue')
# #     plt.scatter(np.array(Ts[i+1]), P2s[i+1], s=12, c='red')
# #     plt.scatter(np.array(Ts[i+1]), P2ms[i+1], label=labels[i+1], s=12, c='red')
# #     # plt.plot(np.array(Ts[i]), P2ms[i], label=labels[i])
#
#     plt.title('$N=10$')
#     plt.xlabel('$T$')
#     plt.ylabel('$P_2$, $P_{2m}$')
#     plt.legend()
#     # plt.show()
#     j+=1
#     print(j)
#
# plt.show()


TsofMaxC = [1.1320002897920742, 1.1080000930720078, 1.0840002644960645, 1.060000254400061, 1.0440004426561875, 1.015999967488001, 1.0, 0.9760001561600251, 0.9599969280098303, 0.9279967334514982]
rs=np.array([0.0,0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9])

m, b = np.polyfit(rs, TsofMaxC, 1)
plt.scatter(rs, TsofMaxC, s=12, color='red')
plt.plot(rs, m*rs+b, c='k', label = '${0:.3f}r+{1:.3f}$'.format(m,b))
plt.title('$N=10$')
plt.xlabel('$r$')
plt.ylabel('$T_{NI}$')
plt.legend()
plt.show()



