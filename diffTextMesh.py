from sympy import *
import sympy
from sympy.utilities.codegen import codegen


r1, p1, w1, x1, y1, z1, r2, p2, w2, x2, y2, z2 = symbols('rA pA wA xA yA zA rB pB wB xB yB zB')
tw1 = symbols('twA:3')
Rw1 = symbols('RwA:3(:3)')
tw2 = symbols('twB:3')
Rw2 = symbols('RwB:3(:3)')
u1,v1 = symbols('uA vA')
u2,v2 = symbols('uB vB')
Planew = symbols('Planew:4')
# invK = symbols('invK:3(:3)')
n1 = symbols('nA:3')
n2 = symbols('nB:3')
T1invK = symbols('TAinvK:3(:4)')

n1Matrix = Matrix([[n1[0]], [n1[1]], [n1[2]]])
n2Matrix = Matrix([[n2[0]], [n2[1]], [n2[2]]])
def Rt2T(R, t):
    return Matrix([
        [R[0, 0], R[0, 1], R[0, 2], t[0]],
        [R[1, 0], R[1, 1], R[1, 2], t[1]],
        [R[2, 0], R[2, 1], R[2, 2], t[2]],
        [0, 0, 0, 1]
    ])
Rw1Matrix = Matrix([Rw1[0:3], Rw1[3:6], Rw1[6:9]])
Rw2Matrix = Matrix([Rw2[0:3], Rw2[3:6], Rw2[6:9]])
tw1Matrix = Matrix([[tw1[0]],[tw1[1]],[tw1[2]]])
tw2Matrix = Matrix([[tw2[0]],[tw2[1]],[tw2[2]]])
# Tw1 = Rt2T(Rw1Matrix, tw1Matrix)
# Tw2 = Rt2T(Rw2Matrix, tw2Matrix)

def EulerToRot(r, p, w):
    cr = cos(r)
    sr = sin(r)
    cp = cos(p)
    sp = sin(p)
    cw = cos(w)
    sw = sin(w)
    rotate = Matrix([
    [cr * cp,   cr * sp * sw - sr * cw,     cr * sp * cw + sr * sw], 
    [sr * cp,   sr * sp * sw + cr * cw,     sr * sp * cw - cr * sw],
    [-sp,       cp * sw,                    cp * cw]])
    return rotate
dR1 = EulerToRot(r1,p1,w1)
# dt1Tmp = Matrix([[x1], [y1], [z1]])
# dt1 = tw1Matrix + dt1Tmp - dR1 * tw1Matrix
dt1 = Matrix([[x1], [y1], [z1]])

dR2 = EulerToRot(r2,p2,w2)
# dt2Tmp = Matrix([[x2], [y2], [z2]])
# dt2 = tw2Matrix + dt2Tmp - dR2 * tw2Matrix
dt2 = Matrix([[x2], [y2], [z2]])

dT1 = Rt2T(dR1, dt1)
# dT2 = Rt2T(dR2, dt2)


planeMatrix = Matrix(Planew)
# invKplus = Matrix([
#         [invK[0], invK[1], invK[2], 0],
#         [invK[3], invK[4], invK[5], 0],
#         [invK[6], invK[6], invK[8], 0],
#         [0, 0, 0, 1]
#     ])
T1invKMatrix = Matrix([T1invK[0:4], T1invK[4:8], T1invK[8:12], [0, 0, 0, 1]])
plane0 = T1invKMatrix.transpose() * dT1.transpose() * planeMatrix
def cal_d(plane, u, v): 
    return -plane[3] / (plane[0] * u + plane[1] * v + plane[2])
d = cal_d(plane0, u1, v1)
pg = dT1 * T1invKMatrix * Matrix([u1 * d, v1 * d, d, 1])
PwPlane = Matrix([[pg[0]], [pg[1]], [pg[2]]])

# invKMatrix = Matrix([invK[0:3], invK[3:6], invK[6:9]])
# def PixToNormal(invK, u, v):
#     n = Matrix([[u], [v], [1]])
#     return invK * n
# nc1 = PixToNormal(invKMatrix, u1, v1)
# nc2 = PixToNormal(invKMatrix, u2, v2)
nw1 = dR1 * Rw1Matrix * n1Matrix
nw2 = dR2 * Rw2Matrix * n2Matrix

def NormalCross(n1, n2):
    skew_symmetrix = Matrix([
                             [0, -n1[2,0], n1[1,0]],
                             [n1[2,0], 0, -n1[0,0]],
                             [-n1[1,0], n1[0,0], 0]
                             ])
    return skew_symmetrix * n2
nwx = NormalCross(nw1, nw2)
nw22 = NormalCross(nwx, nw1)
nw12 = NormalCross(nwx, nw2)

# initValue = {r1:0, p1:0, w1:0, x1:0, y1:0, z1:0, r2:0, p2:0, w2:0, x2:0, y2:0, z2:0}
# with open("_r_diffTextMesh.txt", 'w') as f:
#     f.write("ValueType px = " + str(nwx[0].subs(initValue)) + ";\n")
#     f.write("ValueType px = " + str(nwx[1].subs(initValue)) + ";\n")
#     f.write("ValueType px = " + str(nwx[2].subs(initValue)) + ";\n")

APw = Matrix([
    [nw22[0,0],nw22[1,0],nw22[2,0]],
    [nwx[0,0],nwx[1,0],nwx[2,0]],
    [nw12[0,0],nw12[1,0],nw12[2,0]]
])
bPw = Matrix([
    [nw22.transpose() * (tw1Matrix + dt1)],
    [nwx.transpose() * (tw1Matrix + dt1)],
    [nw12.transpose() * (tw2Matrix + dt2)]
])
# Pw = APw.inv() * bPw

normalDist = (tw2Matrix + dt2 - tw1Matrix - dt1).transpose() * nwx
# planeDist = planeMatrix.transpose() * Matrix([ [Pw[0,0]], [Pw[1,0]], [Pw[2,0]], [1] ])
PwDist = APw * PwPlane - bPw

initValue = {r1:0, p1:0, w1:0, x1:0, y1:0, z1:0, r2:0, p2:0, w2:0, x2:0, y2:0, z2:0}
with open("_r_diffTextMesh.txt", 'w') as f:
    f.write("ValueType errTriDist = " + str(normalDist[0].subs(initValue)) + ";\n")
    f.write("ValueType errPDist0 = " + str(PwDist[0].subs(initValue)) + ";\n")
    f.write("ValueType errPDist1 = " + str(PwDist[1].subs(initValue)) + ";\n")
    f.write("ValueType errPDist2 = " + str(PwDist[2].subs(initValue)) + ";\n")
    f.write("\n")
    def test2(idx, nnn):
        f.write(f"ValueType jacA{nnn}0 = " + str(diff(normalDist[idx], r1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacA{nnn}1 = " + str(diff(normalDist[idx], p1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacA{nnn}2 = " + str(diff(normalDist[idx], w1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacA{nnn}3 = " + str(diff(normalDist[idx], x1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacA{nnn}4 = " + str(diff(normalDist[idx], y1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacA{nnn}5 = " + str(diff(normalDist[idx], z1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}0 = " + str(diff(normalDist[idx], r2).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}1 = " + str(diff(normalDist[idx], p2).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}2 = " + str(diff(normalDist[idx], w2).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}3 = " + str(diff(normalDist[idx], x2).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}4 = " + str(diff(normalDist[idx], y2).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}5 = " + str(diff(normalDist[idx], z2).subs(initValue)) + ";\n")
        f.write("\n")
    test2(0, "0")
    def test(idx, nnn):
        f.write(f"ValueType jacA{nnn}0 = " + str(diff(PwDist[idx], r1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacA{nnn}1 = " + str(diff(PwDist[idx], p1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacA{nnn}2 = " + str(diff(PwDist[idx], w1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacA{nnn}3 = " + str(diff(PwDist[idx], x1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacA{nnn}4 = " + str(diff(PwDist[idx], y1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacA{nnn}5 = " + str(diff(PwDist[idx], z1).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}0 = " + str(diff(PwDist[idx], r2).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}1 = " + str(diff(PwDist[idx], p2).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}2 = " + str(diff(PwDist[idx], w2).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}3 = " + str(diff(PwDist[idx], x2).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}4 = " + str(diff(PwDist[idx], y2).subs(initValue)) + ";\n")
        f.write(f"ValueType jacB{nnn}5 = " + str(diff(PwDist[idx], z2).subs(initValue)) + ";\n")
        f.write("\n")
    test(0, "1")
    test(1, "2")
    test(2, "3")


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

with open("_r_diffTextMesh.txt", 'r') as f:
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
with open("_r_diffTextMesh_code.txt", 'w') as f:
    for v in a1:
        f.write(v)
exit(0)
