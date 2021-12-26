import sympy as sym
import numpy as np
from sympy.matrices import Matrix

def VolumenControl(nodos, elem, nodosVC, k):
    r, s = sym.Symbol("r"), sym.Symbol("s")
    He = [1 - r - s, r, s]    # Funciones de forma
    dHe = sym.Matrix([[-1, 1, 0], [-1, 0, 1]])  # Derivada de Funciones de forma
    KG = np.zeros([len(nodos.xnod), len(nodos.xnod)])
    RG = np.zeros(len(nodos.xnod))
   # BordeNeuman = nodos.puntos[10][:]

    for i in range(len(nodos.xnod)):
        for j in range(len(nodosVC.ElemNode[i][:])):
            Rp1, Rp2 = PuntosIntegración(nodosVC.posElem[i][j])
            invJe, detJe, xi, yi = Jacobiano(He, r, s, nodos.xnod, i, elem)
            dx, dy, surf, nx, ny = deltas(r, s, xi, yi, Rp1, Rp2)

            # Be * n* ds
            Be = invJe * dHe
            # Be = sym.Transpose(sym.Matrix(invJe) * sym.Matrix([He.diff(r), He.diff(s)]))
            integralRp1 = sym.Transpose(Be) * sym.Matrix([[nx[0]], [ny[0]]]) * surf[0]
            integralRp2 = sym.Transpose(Be) * sym.Matrix([[nx[1]], [ny[1]]]) * surf[1]

            integral_ao = integralRp1.subs([(r, Rp1[1][0]), (s, Rp1[1][1])]) -\
                          integralRp1.subs([(r, Rp1[0][0]), (s, Rp1[0][1])])

            integral_oc = integralRp2.subs([(r, Rp2[1][0]), (s, Rp2[1][1])]) -\
                          integralRp2.subs([(r, Rp2[0][0]), (s, Rp2[0][1])])
            Ke = np.array(k * (integral_ao + integral_oc))

            # Ensamble Matriz Global
            elementoActual = nodosVC.ElemNode[i][j]
            nodosElementos = elem[elementoActual, 1:4]
            # nodosActuales = elem[elementoActual, 1:4]

            pos = np.where(nodosElementos == i + 1)
            if pos[0] == 0:
                nodosActuales = nodosElementos
            elif pos[0] == 1:
                nodosActuales = [nodosElementos[1], nodosElementos[2], nodosElementos[0]]
            else:
                nodosActuales = [nodosElementos[2], nodosElementos[0], nodosElementos[1]]

            #KG[0, 0] += Ke[0]
            KG[nodosVC.ID[i], nodosActuales[0] - 1] += Ke[0]
            KG[nodosVC.ID[i], nodosActuales[1] - 1] += Ke[1]
            KG[nodosVC.ID[i], nodosActuales[2] - 1] += Ke[2]
    return KG, RG

# -----------------------------------------------------------------
# ------------------ Puntos de Integración ------------------------
# -----------------------------------------------------------------
def PuntosIntegración(UbicacionNodo):
    # posxy1 = | xa , ya |          posxy2 = | xo , yo |
    #          | xo , yo |                   | xc , yc |
    if UbicacionNodo == 0:
        posxy1, posxy2 = [[0.5, 0.0], [1/3, 1/3]], [[1/3, 1/3], [0.0, 0.5]]
    elif UbicacionNodo == 1:
        posxy1, posxy2 = [[0.5, 0.5], [1/3, 1/3]], [[1/3, 1/3], [0.5, 0.0]]
    else:
        posxy1, posxy2 = [[0.0, 0.5], [1/3, 1/3]], [[1/3, 1/3], [0.5, 0.5]]

    return posxy1, posxy2

# -----------------------------------------------------------------
# ----------------------- JACOBIANO -------------------------------
# -----------------------------------------------------------------
def Jacobiano(H, r, s, POSnodos, IDnodo, elem):
    # Elemento triangular
    xy = [elem[IDnodo, 1], elem[IDnodo, 2], elem[IDnodo, 3]]
    x = H[0] * POSnodos[xy[0] - 1][0] + H[1] * POSnodos[xy[1] - 1][0] + H[2] * POSnodos[xy[2] - 1][0]
    y = H[0] * POSnodos[xy[0] - 1][1] + H[1] * POSnodos[xy[1] - 1][1] + H[2] * POSnodos[xy[2] - 1][1]

    # Jacobiano elemental
    dx_dr, dx_ds = x.diff(r), x.diff(s)
    dy_dr, dy_ds = y.diff(r), y.diff(s)
    Je = sym.Matrix([[dx_dr, dy_dr], [dx_ds, dy_ds]])

    # Inversa Jacobiano elemental
    invJacobiano = Je.inv()
    # Determinante Jacobiano elemental
    detJacobiano = Je.det()

    return invJacobiano, detJacobiano, x, y

# -----------------------------------------------------------------
# --------------------- Función Deltas ----------------------------
# -----------------------------------------------------------------
def deltas(r, s, x, y, Rp1, Rp2):
    ddx1 = - x.subs([(r, Rp1[0][0]), (s, Rp1[0][1])]) + x.subs([(r, Rp1[1][0]), (s, Rp1[1][1])])
    ddx2 = - x.subs([(r, Rp2[0][0]), (s, Rp2[0][1])]) + x.subs([(r, Rp2[1][0]), (s, Rp2[1][1])])
    dx = [ddx1, ddx2]
    # dx = [sym.Abs(ddx1), sym.Abs(ddx2)]

    ddy1 = - y.subs([(r, Rp1[0][0]), (s, Rp1[0][1])]) + y.subs([(r, Rp1[1][0]), (s, Rp1[1][1])])
    ddy2 = - y.subs([(r, Rp2[0][0]), (s, Rp2[0][1])]) + y.subs([(r, Rp2[1][0]), (s, Rp2[1][1])])
    dy = [ddy1, ddy2]
    # dy = [sym.Abs(ddy1), sym.Abs(ddy2)]

    surf = [(dx[0]**2 + dy[0]**2)**0.5, (dx[1]**2 + dy[1]**2)**0.5]

    nx = [dy[0]/surf[0], dy[1]/surf[1]]
    ny = [dx[0] / surf[0], dx[1] / surf[1]]

    return dx, dy, surf, nx, ny
