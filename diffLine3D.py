from sympy import *

outputName = "_r_diffLine3D"

r1, p1, w1, x1, y1, z1, r2, p2, w2, x2, y2, z2 = symbols('rA pA wA xA yA zA rB pB wB xB yB zB')
tw1 = symbols('twA:3')
Rw1 = symbols('RwA:3(:3)')
tw2 = symbols('twB:3')
Rw2 = symbols('RwB:3(:3)')
Pc1 = symbols("PcA:3")
Pc2 = symbols("PcB:3")
nc1 = symbols("ncA:3")
nc2 = symbols("ncB:3")
dc1 = symbols("dcA:3")
dc2 = symbols("dcB:3")
center1 = symbols("centerA:3")
center2 = symbols("centerB:3")

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
Tw1 = Rt2T(Rw1Matrix, tw1Matrix)
Tw2 = Rt2T(Rw2Matrix, tw2Matrix)

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
C1 = Matrix(center1)
dt1Tmp = Matrix([[x1], [y1], [z1]])
dt1 = C1 + dt1Tmp - dR1 * C1
dT1 = Rt2T(dR1, dt1)
dR2 = EulerToRot(r2,p2,w2)
C2 = Matrix(center2)
dt2Tmp = Matrix([[x2], [y2], [z2]])
dt2 = C2 + dt2Tmp - dR2 * C2
dT2 = Rt2T(dR2, dt2)

Pc1Matrix = Matrix([Pc1[0], Pc1[1], Pc1[2], 1])
Pc2Matrix = Matrix([Pc2[0], Pc2[1], Pc2[2], 1])
Pc1nc1 = Matrix([nc1[0], nc1[1], nc1[2], 0])
Pc2nc2 = Matrix([nc2[0], nc2[1], nc2[2], 0])
Pc1dc1 = Matrix([dc1[0], dc1[1], dc1[2], 0])
Pc2dc2 = Matrix([dc2[0], dc2[1], dc2[2], 0])

finalTw1 = dT1 * Tw1
finalTw2 = dT2 * Tw2
finalPw1 = finalTw1 * Pc1Matrix
finalPw2 = finalTw2 * Pc2Matrix
# finalRw1 = Matrix([finalTw1[0:3], finalTw1[4:7], finalTw1[8:11]])
# finalRw2 = Matrix([finalTw2[0:3], finalTw2[4:7], finalTw2[8:11]])
finalnw1 = finalTw1 * Pc1nc1
finalnw2 = finalTw2 * Pc2nc2
finaldw1 = finalTw1 * Pc1dc1
finaldw2 = finalTw2 * Pc2dc2

detPw = finalPw1 - finalPw2
errn = detPw.transpose() * finalnw2
errd = detPw.transpose() * finaldw2

initValue = {r1:0, p1:0, w1:0, x1:0, y1:0, z1:0, r2:0, p2:0, w2:0, x2:0, y2:0, z2:0}
with open(outputName+".txt", 'w') as f:
    f.write("ValueType errn = " + str(errn[0].subs(initValue)) + ";\n")
    f.write("ValueType errd = " + str(errd[0].subs(initValue)) + ";\n")
    f.write("\n")
    def test2():
        f.write(f"ValueType jac_n_T1_0 = " + str(diff(errn[0], r1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T1_1 = " + str(diff(errn[0], p1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T1_2 = " + str(diff(errn[0], w1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T1_3 = " + str(diff(errn[0], x1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T1_4 = " + str(diff(errn[0], y1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T1_5 = " + str(diff(errn[0], z1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T2_0 = " + str(diff(errn[0], r2).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T2_1 = " + str(diff(errn[0], p2).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T2_2 = " + str(diff(errn[0], w2).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T2_3 = " + str(diff(errn[0], x2).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T2_4 = " + str(diff(errn[0], y2).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_n_T2_5 = " + str(diff(errn[0], z2).subs(initValue)) + ";\n")

        f.write(f"ValueType jac_d_T1_0 = " + str(diff(errd[0], r1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T1_1 = " + str(diff(errd[0], p1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T1_2 = " + str(diff(errd[0], w1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T1_3 = " + str(diff(errd[0], x1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T1_4 = " + str(diff(errd[0], y1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T1_5 = " + str(diff(errd[0], z1).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T2_0 = " + str(diff(errd[0], r2).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T2_1 = " + str(diff(errd[0], p2).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T2_2 = " + str(diff(errd[0], w2).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T2_3 = " + str(diff(errd[0], x2).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T2_4 = " + str(diff(errd[0], y2).subs(initValue)) + ";\n")
        f.write(f"ValueType jac_d_T2_5 = " + str(diff(errd[0], z2).subs(initValue)) + ";\n")
        f.write("\n")
    test2()


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

with open(outputName+".txt", 'r') as f:
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
with open(outputName+"_code.txt", 'w') as f:
    for v in a1:
        f.write(v)
exit(0)
