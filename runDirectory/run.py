from cProfile import label
import os
from pathlib import Path
import shutil
from matplotlib import pyplot as plt
import numpy as np

LEVELS = 8

rez = Path(Path.cwd(), 'rez')
shutil.rmtree(rez, ignore_errors=True)
os.mkdir(rez)

for rightPart in range(2):
    for ref_levels in range(LEVELS):
        path1 = Path(rez, f'{rightPart}_grid_{ref_levels}.txt')
        path2 = Path(rez, f'{rightPart}_time_{ref_levels}.txt')
        os.system(f'(time ./tutorialRun -rightPart {rightPart} -refine_levels {ref_levels}) > {path1} 2> {path2}')


cells = []
time = []

for i in range(LEVELS):
    path1 = Path(rez, f'1_grid_{i}.txt')
    path2 = Path(rez, f'1_time_{i}.txt')
    with open(path1, 'r') as f:
        line = f.readline()
        tmp = line.rfind(':')
        cells.append(int(line[tmp + 1:]))
    with open(path2, 'r') as f:
        line = (f.readline()).split()[:3]
        elapsed = line[2][:-7]
        q = elapsed.rfind('.')
        time.append([line[0][:-4], line[1][:-6], elapsed[q - 1:]])
time = np.array(time)

timeNames = ['user', 'system', 'elapsed']
for i in range(3):
    tmp = list(map(float, time[:, i]))
    p = np.polyfit(cells, tmp, 1)
    y = np.polyval(p, cells)
    print(cells, tmp, p)
    plt.scatter(cells, tmp, label=timeNames[i])
    a = '%.2e' % p[0]
    b = '%.2e' % p[1]
    plt.plot(cells, y, label = f'{timeNames[i]}: y = ' + a + r'$\cdot x + $' + b)
plt.title('New right part')
plt.xlabel('cells')
plt.ylabel(f'time')
plt.legend()
plt.grid(True)
plt.savefig('rez_new.png', dpi=500)

plt.clf()
    
cells = []
time = []

for i in range(LEVELS):
    path1 = Path(rez, f'0_grid_{i}.txt')
    path2 = Path(rez, f'0_time_{i}.txt')
    with open(path1, 'r') as f:
        line = f.readline()
        tmp = line.rfind(':')
        cells.append(int(line[tmp + 1:]))
    with open(path2, 'r') as f:
        line = (f.readline()).split()[:3]
        elapsed = line[2][:-7]
        q = elapsed.rfind('.')
        time.append([line[0][:-4], line[1][:-6], elapsed[q - 1:]])
time = np.array(time)

timeNames = ['user', 'system', 'elapsed']
for i in range(3):
    tmp = list(map(float, time[:, i]))
    p = np.polyfit(cells, tmp, 1)
    y = np.polyval(p, cells)
    print(cells, tmp, p)
    plt.scatter(cells, tmp, label=timeNames[i])
    a = '%.2e' % p[0]
    b = '%.2e' % p[1]
    plt.plot(cells, y, label = f'{timeNames[i]}: y = ' + a + r'$\cdot x + $' + b)
plt.title('Const right part')
plt.xlabel('cells')
plt.ylabel(f'time')
plt.legend()
plt.grid(True)
plt.savefig('rez_const.png', dpi=500)