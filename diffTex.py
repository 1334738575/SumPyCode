from sympy import *
import sympy
from sympy.utilities.codegen import codegen


r, p, w, x, y, z, Px, Py, Pz = symbols('r p w x y z Px Py Pz')
Twc = symbols('Twc:3(:4)')
u,v = symbols('u v')
fx, fy, cx, cy = symbols('fx fy cx cy')
Pw = symbols('Pw:3')
center = symbols('center:3')
# detPw = symbols('detPw:3')
planew = symbols('planew:4')


def Rt2T(R, t):
    return Matrix([
        [R[0, 0], R[0, 1], R[0, 2], t[0]],
        [R[1, 0], R[1, 1], R[1, 2], t[1]],
        [R[2, 0], R[2, 1], R[2, 2], t[2]],
        [0, 0, 0, 1]
    ])
TwcMatrix = Matrix([Twc[0:4], Twc[4:8], Twc[8:12], [0, 0, 0, 1]])
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
dR = EulerToRot(r,p,w)
C = Matrix(center)
dtTmp = Matrix([[x], [y], [z]])
dt = C + dtTmp - dR * C
dT = Rt2T(dR, dt)

KMatrix = Matrix([
        [fx, 0, cx, 0],
        [0, fy, cy, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
PwMatrix = Matrix([Pw[0], Pw[1], Pw[2], 1])
# detPwMatrix = Matrix([detPw[0], detPw[1], detPw[2], 0])
dPw = Matrix([Px, Py, Pz, 0])

finalTwc = dT * TwcMatrix
# finalPw = PwMatrix - detPwMatrix + dPw
finalPw = PwMatrix + dPw
finalRwc = Matrix([finalTwc[0:3], finalTwc[4:7], finalTwc[8:11]])
finaltwc = Matrix([finalTwc[3], finalTwc[7], finalTwc[11]])
finalRcw = finalRwc.transpose()
finaltcw = -1 * finalRcw * finaltwc
finalTcw = Rt2T(finalRcw, finaltcw)

Pc = finalTcw * finalPw
PcNorm = Pc / Pc[2]
uv = KMatrix * PcNorm
erru = u - uv[0]
errv = v - uv[1]

planewMatrix = Matrix([planew[0], planew[1], planew[2], planew[3]])
errPlane = (planewMatrix.transpose() * finalPw)[0]

initValue = {r:0, p:0, w:0, x:0, y:0, z:0, Px:0, Py:0, Pz:0}
with open("_r_diffText.txt", 'w') as f:
    f.write("ValueType erru = " + str(erru.subs(initValue)) + ";\n")
    f.write("ValueType errv = " + str(errv.subs(initValue)) + ";\n")
    f.write("ValueType errPlane = " + str(errPlane.subs(initValue)) + ";\n")
    f.write("\n")
    def test2():
        f.write(f"ValueType jacuT0 = " + str(diff(erru, r).subs(initValue)) + ";\n")
        f.write(f"ValueType jacuT1 = " + str(diff(erru, p).subs(initValue)) + ";\n")
        f.write(f"ValueType jacuT2 = " + str(diff(erru, w).subs(initValue)) + ";\n")
        f.write(f"ValueType jacuT3 = " + str(diff(erru, x).subs(initValue)) + ";\n")
        f.write(f"ValueType jacuT4 = " + str(diff(erru, y).subs(initValue)) + ";\n")
        f.write(f"ValueType jacuT5 = " + str(diff(erru, z).subs(initValue)) + ";\n")
        f.write(f"ValueType jacuP0 = " + str(diff(erru, Px).subs(initValue)) + ";\n")
        f.write(f"ValueType jacuP1 = " + str(diff(erru, Py).subs(initValue)) + ";\n")
        f.write(f"ValueType jacuP2 = " + str(diff(erru, Pz).subs(initValue)) + ";\n")

        f.write(f"ValueType jacvT0 = " + str(diff(errv, r).subs(initValue)) + ";\n")
        f.write(f"ValueType jacvT1 = " + str(diff(errv, p).subs(initValue)) + ";\n")
        f.write(f"ValueType jacvT2 = " + str(diff(errv, w).subs(initValue)) + ";\n")
        f.write(f"ValueType jacvT3 = " + str(diff(errv, x).subs(initValue)) + ";\n")
        f.write(f"ValueType jacvT4 = " + str(diff(errv, y).subs(initValue)) + ";\n")
        f.write(f"ValueType jacvT5 = " + str(diff(errv, z).subs(initValue)) + ";\n")
        f.write(f"ValueType jacvP0 = " + str(diff(errv, Px).subs(initValue)) + ";\n")
        f.write(f"ValueType jacvP1 = " + str(diff(errv, Py).subs(initValue)) + ";\n")
        f.write(f"ValueType jacvP2 = " + str(diff(errv, Pz).subs(initValue)) + ";\n")

        f.write(f"ValueType jacdT0 = " + str(diff(errPlane, r).subs(initValue)) + ";\n")
        f.write(f"ValueType jacdT1 = " + str(diff(errPlane, p).subs(initValue)) + ";\n")
        f.write(f"ValueType jacdT2 = " + str(diff(errPlane, w).subs(initValue)) + ";\n")
        f.write(f"ValueType jacdT3 = " + str(diff(errPlane, x).subs(initValue)) + ";\n")
        f.write(f"ValueType jacdT4 = " + str(diff(errPlane, y).subs(initValue)) + ";\n")
        f.write(f"ValueType jacdT5 = " + str(diff(errPlane, z).subs(initValue)) + ";\n")
        f.write(f"ValueType jacdP0 = " + str(diff(errPlane, Px).subs(initValue)) + ";\n")
        f.write(f"ValueType jacdP1 = " + str(diff(errPlane, Py).subs(initValue)) + ";\n")
        f.write(f"ValueType jacdP2 = " + str(diff(errPlane, Pz).subs(initValue)) + ";\n")
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

with open("_r_diffText.txt", 'r') as f:
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
with open("_r_diffText_code.txt", 'w') as f:
    for v in a1:
        f.write(v)
exit(0)
