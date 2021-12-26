import numpy as np

def VolumenControl(nodos, elem, nodosVC, k, tipo):
    KG = np.zeros([len(nodos.xnod), len(nodos.xnod)])
    RG = np.zeros(len(nodos.xnod))

    for i in range(len(nodos.xnod)):
        for j in range(len(nodosVC.ElemNode[i][:])):

            # ---------------- TRIANGULAR ---------------
            if tipo == 1:
                elementoActual = nodosVC.ElemNode[i][j]
                nodosElementos = elem[elementoActual, 1:4]

                # ---- posición del Nodo " i ", dentro del elemento ----
                pos = np.where(nodosElementos == i + 1)
                if pos[0] == 0:
                    nodosActuales = nodosElementos
                elif pos[0] == 1:
                    nodosActuales = [nodosElementos[1], nodosElementos[2], nodosElementos[0]]
                else:
                    nodosActuales = [nodosElementos[2], nodosElementos[0], nodosElementos[1]]

                # x1, x2, x3   ---   y1, y2, y3
                x = [nodos.xnod[nodosActuales[0] - 1, 0], nodos.xnod[nodosActuales[1] - 1, 0], nodos.xnod[nodosActuales[2] - 1, 0]]
                y = [nodos.xnod[nodosActuales[0] - 1, 1], nodos.xnod[nodosActuales[1] - 1, 1], nodos.xnod[nodosActuales[2] - 1, 1]]

                # Volumen Elemental, siendo espesor Unitario
                Ve = 0.5 * ((x[1]*y[2] - x[2]*y[1]) - x[0]*(y[2] - y[1]) + y[0]*(x[2] - x[1]))

                # Derivadas de Funciones de Forma
                N1x = (y[1] - y[2]) * 0.5 / Ve
                N1y = (x[2] - x[1]) * 0.5 / Ve
                N2x = (y[2] - y[0]) * 0.5 / Ve
                N2y = (x[0] - x[2]) * 0.5 / Ve
                N3x = (y[0] - y[1]) * 0.5 / Ve
                N3y = (x[1] - x[0]) * 0.5 / Ve

                # deltas
                dxf1 = x[2] / 3 - x[1] / 6 - x[0] / 6
                dxf2 = - x[1] / 3 + x[2] / 6 + x[0] / 6
                dyf1 = y[2] / 3 - y[1] / 6 - y[0] / 6
                dyf2 = - y[1] / 3 + y[2] / 6 + y[0] / 6

                # Areas
                af1 = (dxf1 ** 2 + dyf1 ** 2) ** 0.5
                af2 = (dxf2 ** 2 + dyf2 ** 2) ** 0.5

                # matriz elemental a global
                a1k = - N1x * dyf1 + N1y * dxf1 - N1x * dyf2 + N1y * dxf2
                a2k = N2x * dyf1 - N2y * dxf1 + N2x * dyf2 - N2y * dxf2
                a3k = N3x * dyf1 - N3y * dxf1 + N3x * dyf2 - N3y * dxf2

                #KG[0, 0] += Ke[0]
                KG[nodosVC.ID[i], nodosActuales[0] - 1] -= k * a1k
                KG[nodosVC.ID[i], nodosActuales[1] - 1] += k * a2k
                KG[nodosVC.ID[i], nodosActuales[2] - 1] += k * a3k

            # ------------------ CUAD -------------------
            else:
                elementoActual = nodosVC.ElemNode[i][j]
                nodosElementos = elem[elementoActual, 1:5]

                # ---- posición del Nodo " i ", dentro del elemento ----
                pos = np.where(nodosElementos == i + 1)
                if pos[0] == 0:
                    nodosActuales = nodosElementos
                    numeracion = [0, 1, 2, 3]
                    #                      x  ,  y      x  ,  y
                    limite_ao = np.array([[0.0, 1.0], [0.0, 0.0]])
                    limite_ob = np.array([[0.0, 0.0], [1.0, 0.0]])
                    df1, df2, x, y = deltas(limite_ao, limite_ob, nodosElementos, nodos)
                elif pos[0] == 1:
                    nodosActuales = [nodosElementos[1], nodosElementos[2], nodosElementos[3], nodosElementos[0]]
                    numeracion = [1, 2, 3, 0]
                    limite_ao = np.array([[-1.0, 0.0], [0.0, 0.0]])
                    limite_ob = np.array([[0.0, 0.0], [0.0, 1.0]])
                    df1, df2, x, y = deltas(limite_ao, limite_ob, nodosElementos, nodos)
                elif pos[0] == 2:
                    nodosActuales = [nodosElementos[2], nodosElementos[3], nodosElementos[0], nodosElementos[1]]
                    numeracion = [2, 3, 0, 1]
                    limite_ao = np.array([[0.0, -1.0], [0.0, 0.0]])
                    limite_ob = np.array([[0.0, 0.0], [-1.0, 0.0]])
                    df1, df2, x, y = deltas(limite_ao, limite_ob, nodosElementos, nodos)
                else:
                    nodosActuales = [nodosElementos[3], nodosElementos[0], nodosElementos[1], nodosElementos[2]]
                    numeracion = [3, 0, 1, 2]
                    limite_ao = np.array([[1.0, 0.0], [0.0, 0.0]])
                    limite_ob = np.array([[0.0, 0.0], [0.0, -1.0]])
                    df1, df2, x, y = deltas(limite_ao, limite_ob, nodosElementos, nodos)

                # ----- Derivadas de funciones de forma -----
                # en r = 0   y   s = 0
                dH = np.array([[0.25, -0.25, -0.25, 0.25], [0.25, 0.25, -0.25, -0.25]])

                # ----- Jacobiano -----
                invJe, detJe = JacobianoElemental(x, y, dH)

                # ----- Gradiente de Funciones de Forma -----
                Be = np.dot(invJe, dH)

                # matriz elemental a global
                ak = k * (np.dot(np.transpose(Be), df1) + np.dot(np.transpose(Be), df2))

                # KG[0, 0] += Ke[0]
                KG[nodosVC.ID[i], nodosActuales[0] - 1] += ak[numeracion[0]]
                KG[nodosVC.ID[i], nodosActuales[1] - 1] += ak[numeracion[1]]
                KG[nodosVC.ID[i], nodosActuales[2] - 1] += ak[numeracion[2]]
                KG[nodosVC.ID[i], nodosActuales[3] - 1] += ak[numeracion[3]]
    return KG, RG


# ---------------------------------------------------------
# -----------------  Jacobiano  ---------------------------
# ---------------------------------------------------------
def JacobianoElemental(X, Y, dH):

    dxdr = X[0] * dH[0, 0] + X[1] * dH[0, 1] + X[2] * dH[0, 2] + X[3] * dH[0, 3]
    dxds = X[0] * dH[1, 0] + X[1] * dH[1, 1] + X[2] * dH[1, 2] + X[3] * dH[1, 3]
    dydr = Y[0] * dH[0, 0] + Y[1] * dH[0, 1] + Y[2] * dH[0, 2] + Y[3] * dH[0, 3]
    dyds = Y[0] * dH[1, 0] + Y[1] * dH[1, 1] + Y[2] * dH[1, 2] + Y[3] * dH[1, 3]

    Je = np.array([[dxdr, dydr], [dxds, dyds]])

    invJe = np.linalg.inv(Je)
    detJe = np.linalg.det(Je)
    return invJe, detJe

# ---------------------------------------------------------
# -------------  Superficie de Integración  ---------------
# ---------------------------------------------------------
def deltas(lim_ao,lim_ob, nodosActuales, nodos):

    # x1, x2, x3, x4   ---   y1, y2, y3, y4
    x = np.array([nodos.xnod[nodosActuales[0] - 1, 0], nodos.xnod[nodosActuales[1] - 1, 0],
                  nodos.xnod[nodosActuales[2] - 1, 0], nodos.xnod[nodosActuales[3] - 1, 0]])
    y = np.array([nodos.xnod[nodosActuales[0] - 1, 1], nodos.xnod[nodosActuales[1] - 1, 1],
                  nodos.xnod[nodosActuales[2] - 1, 1], nodos.xnod[nodosActuales[3] - 1, 1]])

    # ----- deltas -----
    h1_lim_ao = np.array([[0.25 * (1 + lim_ao[0, 0]) * (1 + lim_ao[0, 1])],
                          [0.25 * (1 + lim_ao[1, 0]) * (1 + lim_ao[1, 1])]])
    h2_lim_ao = np.array([[0.25 * (1 - lim_ao[0, 0]) * (1 + lim_ao[0, 1])],
                          [0.25 * (1 - lim_ao[1, 0]) * (1 + lim_ao[1, 1])]])
    h3_lim_ao = np.array([[0.25 * (1 - lim_ao[0, 0]) * (1 - lim_ao[0, 1])],
                          [0.25 * (1 - lim_ao[1, 0]) * (1 - lim_ao[1, 1])]])
    h4_lim_ao = np.array([[0.25 * (1 + lim_ao[0, 0]) * (1 - lim_ao[0, 1])],
                          [0.25 * (1 + lim_ao[1, 0]) * (1 - lim_ao[1, 1])]])

    h1_lim_ob = np.array([[0.25 * (1 + lim_ob[0, 0]) * (1 + lim_ob[0, 1])],
                          [0.25 * (1 + lim_ob[1, 0]) * (1 + lim_ob[1, 1])]])
    h2_lim_ob = np.array([[0.25 * (1 - lim_ob[0, 0]) * (1 + lim_ob[0, 1])],
                          [0.25 * (1 - lim_ob[1, 0]) * (1 + lim_ob[1, 1])]])
    h3_lim_ob = np.array([[0.25 * (1 - lim_ob[0, 0]) * (1 - lim_ob[0, 1])],
                          [0.25 * (1 - lim_ob[1, 0]) * (1 - lim_ob[1, 1])]])
    h4_lim_ob = np.array([[0.25 * (1 + lim_ob[0, 0]) * (1 - lim_ob[0, 1])],
                          [0.25 * (1 + lim_ob[1, 0]) * (1 - lim_ob[1, 1])]])

    xa = x[0] * h1_lim_ao[0] + x[1] * h2_lim_ao[0] + x[2] * h3_lim_ao[0] + x[3] * h4_lim_ao[0]
    xo = x[0] * h1_lim_ao[1] + x[1] * h2_lim_ao[1] + x[2] * h3_lim_ao[1] + x[3] * h4_lim_ao[1]
    xb = x[0] * h1_lim_ob[1] + x[1] * h2_lim_ob[1] + x[2] * h3_lim_ob[1] + x[3] * h4_lim_ob[1]

    ya = y[0] * h1_lim_ao[0] + y[1] * h2_lim_ao[0] + y[2] * h3_lim_ao[0] + y[3] * h4_lim_ao[0]
    yo = y[0] * h1_lim_ao[1] + y[1] * h2_lim_ao[1] + y[2] * h3_lim_ao[1] + y[3] * h4_lim_ao[1]
    yb = y[0] * h1_lim_ob[1] + y[1] * h2_lim_ob[1] + y[2] * h3_lim_ob[1] + y[3] * h4_lim_ob[1]

    dxf1 = xo - xa
    dyf1 = yo - ya
    dxf2 = xb - xo
    dyf2 = yb - yo

    df1 = np.array([dyf1[0], -dxf1[0]])
    df2 = np.array([dyf2[0], -dxf2[0]])

    # Areas
    af1 = (dxf1[0] ** 2 + dyf1[0] ** 2) ** 0.5
    af2 = (dxf2[0] ** 2 + dyf2[0] ** 2) ** 0.5

    return df1, df2, x, y
