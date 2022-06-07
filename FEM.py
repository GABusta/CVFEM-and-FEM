import numpy as np
# otro cmabio
def FiniteElement(nodos, elem, k, tipo):
    # puntos de integración
    rg, sg, w, H, dH = Gauss_dif_Forma(tipo)
    Ntotales = len(nodos.xnod[:, 0])
    KG = np.zeros((Ntotales, Ntotales))
    RG = np.zeros(Ntotales)

    for i in range(len(elem[:, 0])):

        # ---------------- TRIANGULAR ---------------
        if tipo == 1:
            Elemento = elem[i, 1:4] - 1
            nodosX = np.array([nodos.xnod[Elemento[0], 0], nodos.xnod[Elemento[1], 0],
                               nodos.xnod[Elemento[2], 0]])
            nodosY = np.array([nodos.xnod[Elemento[0], 1], nodos.xnod[Elemento[1], 1],
                               nodos.xnod[Elemento[2], 1]])

            invJe, detJe = JacobianoElemental(nodosX, nodosY, dH, tipo)
            # ----- Matriz Elemental -----
            ke11 = np.zeros((3, 3))
            for ig in range(len(rg)):
                # Matriz "Be" Elemental
                Be = np.dot(invJe, dH)
                # k11 =  BeT * k * Be * |Je|
                k11 = k * detJe * np.dot(np.transpose(Be), Be)
                ke11 += w[ig] * k11
                ke = ke11 * 0.5
                # ----- Matriz Global -----
                for j in range(len(rg)):
                    KG[Elemento[0], Elemento[j]] += ke[0, j]
                    KG[Elemento[1], Elemento[j]] += ke[1, j]
                    KG[Elemento[2], Elemento[j]] += ke[2, j]


        # ------------------ CUAD -------------------
        else:
            Elemento = elem[i, 1:5] - 1
            nodosX = np.array([nodos.xnod[Elemento[0], 0], nodos.xnod[Elemento[1], 0],
                               nodos.xnod[Elemento[2], 0], nodos.xnod[Elemento[3], 0]])
            nodosY = np.array([nodos.xnod[Elemento[0], 1], nodos.xnod[Elemento[1], 1],
                               nodos.xnod[Elemento[2], 1], nodos.xnod[Elemento[3], 1]])

            # ----- Matriz Elemental -----
            ke11 = np.zeros((4, 4))
            for ig in range(len(rg)):
                invJe, detJe = JacobianoElemental(nodosX, nodosY, dH, tipo)
                # Matriz "Be" Elemental
                Be = np.dot(invJe, dH)
                # k11 =  BeT * k * Be * |Je|
                k11 = k * detJe * np.dot(np.transpose(Be), Be)
                ke11 += w[ig] * k11
                ke = ke11

            # ----- Matriz Global -----
            for j in range(len(rg)):
                KG[Elemento[0], Elemento[j]] += ke[0, j]
                KG[Elemento[1], Elemento[j]] += ke[1, j]
                KG[Elemento[2], Elemento[j]] += ke[2, j]
                KG[Elemento[3], Elemento[j]] += ke[3, j]

    return KG, RG


# -------------------------------------------------------
# --------------- Puntos de Gauss -----------------------
# -------------------------------------------------------
def Gauss_dif_Forma(tipo):
    # --------------------------------
    # ------ TRIANGULAR ELEMENT ------
    # --------------------------------
    if tipo == 1:
        rg = np.array([1 / 6, 2 / 3, 1 / 6])
        sg = np.array([rg[0], rg[0], rg[1]])
        w = np.array([1 / 3, 1 / 3, 1 / 3])

        # H( 3 , 2 , 3)
        # ----- funciones de forma -----
        H = np.array([[[1.0 - rg[0] - sg[0], rg[0], sg[0]],
                      [1.0 - rg[0] - sg[0], rg[0], sg[0]]],

                     [[1.0 - rg[1] - sg[1], rg[1], sg[1]],
                      [1.0 - rg[1] - sg[1], rg[1], sg[1]]],

                     [[1.0 - rg[2] - sg[2], rg[2], sg[2]],
                      [1.0 - rg[2] - sg[2], rg[2], sg[2]]]])

        # ----- Derivadas de funciones de forma -----
        dH = np.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]])

    # ---------------------------
    # ------ CUAD. ELEMENT ------
    # ---------------------------
    else:
        pg = 0.577350269189626
        rg = np.array([pg, -pg, -pg, pg])
        sg = np.array([pg, pg, -pg, -pg])
        w = np.array([1.0, 1.0, 1.0, 1.0])

        # H( 4 , 2 , 3)
        H = np.zeros((4, 2, 4))
        for i in range(4):
            # ----- funciones de forma -----
            h0 = np.array([[0.25 * (1 + rg[i])*(1 + sg[i]), 0.25 * (1 - rg[i])*(1 + sg[i]),
                           0.25 * (1 - rg[i])*(1 - sg[i]), 0.25 * (1 + rg[i])*(1 - sg[i])],

                          [0.25 * (1 + rg[i]) * (1 + sg[i]), 0.25 * (1 - rg[i]) * (1 + sg[i]),
                           0.25 * (1 - rg[i]) * (1 - sg[i]), 0.25 * (1 + rg[i]) * (1 - sg[i])]])

            H[i, :, :] = h0[:, :]

            # ----- Derivadas de funciones de forma -----
            # en r = 0   y   s = 0
            dH = np.array([[0.25, -0.25, -0.25, 0.25], [0.25, 0.25, -0.25, -0.25]])

    return rg, sg, w, H, dH

# -------------------------------------------------------
# ------------- Jacobiano Elemental ---------------------
# -------------------------------------------------------
def JacobianoElemental(X, Y, dH, tipo):
    # --------------------------------
    # ------ TRIANGULAR ELEMENT ------
    # --------------------------------
    if tipo == 1:
        dxdr = X[0] * dH[0, 0] + X[1] * dH[0, 1] + X[2] * dH[0, 2]
        dxds = X[0] * dH[1, 0] + X[1] * dH[1, 1] + X[2] * dH[1, 2]
        dydr = Y[0] * dH[0, 0] + Y[1] * dH[0, 1] + Y[2] * dH[0, 2]
        dyds = Y[0] * dH[1, 0] + Y[1] * dH[1, 1] + Y[2] * dH[1, 2]
    # ---------------------------
    # ------ CUAD. ELEMENT ------
    # ---------------------------
    else:
        dxdr = X[0] * dH[0, 0] + X[1] * dH[0, 1] + X[2] * dH[0, 2] + X[3] * dH[0, 3]
        dxds = X[0] * dH[1, 0] + X[1] * dH[1, 1] + X[2] * dH[1, 2] + X[3] * dH[1, 3]
        dydr = Y[0] * dH[0, 0] + Y[1] * dH[0, 1] + Y[2] * dH[0, 2] + Y[3] * dH[0, 3]
        dyds = Y[0] * dH[1, 0] + Y[1] * dH[1, 1] + Y[2] * dH[1, 2] + Y[3] * dH[1, 3]

    Je = np.array([[dxdr, dydr], [dxds, dyds]])

    invJe = np.linalg.inv(Je)
    detJe = np.linalg.det(Je)
    return invJe, detJe
