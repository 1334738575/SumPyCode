from sympy import *
import sympy
from sympy.utilities.codegen import codegen

a,b,r,x,y,z = symbols('a b r x y z')
center = symbols('center:3')
M0 = symbols('M:3(:4)')
rM0 = symbols('rM:3(:4)')
u, v = symbols('u v')
planeG0 = symbols("planeG:4")
ca = cos(a)
sa = sin(a)
cb = cos(b)
sb = sin(b)
cr = cos(r)
sr = sin(r)
def printMatrix(M, name):
    print(name + ":")
    for v in range(M.rows):
        for j in range(M.cols):
            end = ", "
            if j == M.cols - 1:
                end = '\n'
            print(M.row(v)[j], end = end)
    print("")
    return
# 创建一个 2x2 矩阵
rotate = Matrix([
    [ca * cb,   ca * sb * sr - sa * cr,     ca * sb * cr + sa * sr], 
    [sa * cb,   sa * sb * sr + ca * cr,     sa * sb * cr - ca * sr],
    [-sb,       cb * sr,                    cb * cr]])
T0 = Matrix([[x], [y], [z]])
C = Matrix(center)
T = T0 + C - rotate * C
M = Matrix([M0[0:4], M0[4:8], M0[8:12], [0, 0, 0, 1]])
rM = Matrix([rM0[0:4], rM0[4:8], rM0[8:12], [0, 0, 0, 1]])
planeG = Matrix(planeG0)
def gendR(R, T):
    return Matrix([
        [R[0, 0], R[0, 1], R[0, 2], T[0]],
        [R[1, 0], R[1, 1], R[1, 2], T[1]],
        [R[2, 0], R[2, 1], R[2, 2], T[2]],
        [0, 0, 0, 1]
    ])
dR = gendR(rotate, T)
plane0 = M.transpose() * dR.transpose() * planeG
def cal_d(plane, u, v): 
    return -plane[3] / (plane[0] * u + plane[1] * v + plane[2])
d = cal_d(plane0, u, v)
pg = dR * M * Matrix([u * d, v * d, d, 1])
initValue = {a: 0, b: 0, r: 0, x:0, y:0, z:0}
if False:
    pi4 = symbols("pi:4")
    pi40 = rM.transpose() * Matrix(pi4)
    printMatrix(pi40, "pi40")
    ddd = cal_d(pi40, u, v)
    print(ddd)
    ppp = M * Matrix([u, v, ddd, 1])
    print(ppp)
    exit(0)

with open("_r_difftest.txt", 'w') as f:
    f.write("//predefined\n")
    f.write("typedef double ValueType;\n")
    for i in range(4):
        for j in range(4):
            f.write(f"const auto& M{i}{j} = M({i}, {j});")
    f.write("\n")
    for i in range(4):
        for j in range(4):
            f.write(f"const auto& rM{i}{j} = rM({i}, {j});")
    f.write("\n")
    for i in range(4):
        f.write(f"const auto& planeG{i} = planeG[{i}];")
    f.write("\n")
    for i in range(3):
        f.write(f"const auto& center{i} = center[{i}];")
    f.write("\n")

    f.write("//----------------with  a = 0, b = 0, c = 0, dx = 0, dy = 0, dz = 0\n")
    f.write("//ValueType d = " + str(d.subs(initValue)) + ";\n")
    f.write("ValueType px = " + str(pg[0].subs(initValue)) + ";\n")
    f.write("ValueType py = " + str(pg[1].subs(initValue)) + ";\n")
    f.write("ValueType pz = " + str(pg[2].subs(initValue)) + ";\n")
    f.write("\n")
    def test(idx, nnn):
        f.write(f"ValueType dp{nnn}0 = " + str(diff(pg[idx], a).subs(initValue)) + ";\n")
        f.write(f"ValueType dp{nnn}1 = " + str(diff(pg[idx], b).subs(initValue)) + ";\n")
        f.write(f"ValueType dp{nnn}2 = " + str(diff(pg[idx], r).subs(initValue)) + ";\n")
        f.write(f"ValueType dp{nnn}3 = " + str(diff(pg[idx], x).subs(initValue)) + ";\n")
        f.write(f"ValueType dp{nnn}4 = " + str(diff(pg[idx], y).subs(initValue)) + ";\n")
        f.write(f"ValueType dp{nnn}5 = " + str(diff(pg[idx], z).subs(initValue)) + ";\n")
        f.write("\n")
    test(0, "x")
    test(1, "y")
    test(2, "z")
    if True:
        f.write("//------------------------------full define\n")
        f.write("full_px = " + str(pg[0]) + "\n")
        f.write("full_py = " + str(pg[1]) + "\n")
        f.write("full_pz = " + str(pg[2]) + "\n")
        f.write("\n")
        def test(idx, nnn):
            f.write(f"full_dp{nnn}0 = " + str(diff(pg[idx], a)) + "\n")
            f.write(f"full_dp{nnn}1 = " + str(diff(pg[idx], b)) + "\n")
            f.write(f"full_dp{nnn}2 = " + str(diff(pg[idx], r)) + "\n")
            f.write(f"full_dp{nnn}3 = " + str(diff(pg[idx], x)) + "\n")
            f.write(f"full_dp{nnn}4 = " + str(diff(pg[idx], y)) + "\n")
            f.write(f"full_dp{nnn}5 = " + str(diff(pg[idx], z)) + "\n")
            f.write("\n")
        test(0, "x")
        test(1, "y")
        test(2, "z")
    f.close()
    
def find_all_substring_indices(string, substring):
    start = 0
    indices = []
    while True:
        start = string.find(substring, start)
        if start == -1:
            break
        indices.append(start)
        start += 1  # 移动到下一个字符位置
    return indices

with open("_r_difftest.txt", 'r') as f:
    lines = f.readlines()
    a1 = []
    f.close()
    for line in lines:
        line = line.replace("1.0", "ValueType(1.0f)")
        if "**" not in line:
            a1.append(line)
            continue
        b = ""
        pp = find_all_substring_indices(line, "**2")
        rp = []
        first = 0
        for p in pp:
            fp = p - 1
            current = 0
            while True:
                if line[fp] == ")":
                    current += 1
                    fp -= 1
                    continue
                if line[fp] == '(':
                    current -= 1
                    if current == 0:
                        break
                    fp -= 1
                    continue
                fp -= 1
                continue
            rp.append(line[first:fp])
            rp.append("pow")
            rp.append(line[fp:p - 1])
            rp.append(",2)")
            first = p + 3
        rp.append(line[first:])
        for v in rp:
            b += v
        a1.append(b)
with open("_r_difftest.txt", 'w') as f:
    for v in a1:
        f.write(v)
exit(0)