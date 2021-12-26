import numpy as np
import os

def matrices(nodos, elem, tipo, Ktotal, RutaTexto, tiempo):

    determinante = np.zeros(len(elem[:, 0]))
    for i in range(len(elem[:, 0])):
        # --------------------------------
        # ------ TRIANGULAR ELEMENT ------
        # --------------------------------
        if tipo == 1:
            dH = np.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]])
            X = np.array([nodos.xnod[elem[i, 1] - 1, 0], nodos.xnod[elem[i, 2] - 1, 0], nodos.xnod[elem[i, 3] - 1, 0]])
            Y = np.array([nodos.xnod[elem[i, 1] - 1, 1], nodos.xnod[elem[i, 2] - 1, 1], nodos.xnod[elem[i, 3] - 1, 1]])
            dxdr = X[0] * dH[0, 0] + X[1] * dH[0, 1] + X[2] * dH[0, 2]
            dxds = X[0] * dH[1, 0] + X[1] * dH[1, 1] + X[2] * dH[1, 2]
            dydr = Y[0] * dH[0, 0] + Y[1] * dH[0, 1] + Y[2] * dH[0, 2]
            dyds = Y[0] * dH[1, 0] + Y[1] * dH[1, 1] + Y[2] * dH[1, 2]
        # ---------------------------
        # ------ CUAD. ELEMENT ------
        # ---------------------------
        else:
            # en r = 0   y   s = 0
            dH = np.array([[0.25, -0.25, -0.25, 0.25], [0.25, 0.25, -0.25, -0.25]])
            X = np.array([nodos.xnod[elem[i, 1] - 1, 0], nodos.xnod[elem[i, 2] - 1, 0],
                          nodos.xnod[elem[i, 3] - 1, 0], nodos.xnod[elem[i, 4] - 1, 0]])
            Y = np.array([nodos.xnod[elem[i, 1] - 1, 1], nodos.xnod[elem[i, 2] - 1, 1],
                          nodos.xnod[elem[i, 3] - 1, 1], nodos.xnod[elem[i, 4] - 1, 1]])
            dxdr = X[0] * dH[0, 0] + X[1] * dH[0, 1] + X[2] * dH[0, 2] + X[3] * dH[0, 3]
            dxds = X[0] * dH[1, 0] + X[1] * dH[1, 1] + X[2] * dH[1, 2] + X[3] * dH[1, 3]
            dydr = Y[0] * dH[0, 0] + Y[1] * dH[0, 1] + Y[2] * dH[0, 2] + Y[3] * dH[0, 3]
            dyds = Y[0] * dH[1, 0] + Y[1] * dH[1, 1] + Y[2] * dH[1, 2] + Y[3] * dH[1, 3]


        Je = np.array([[dxdr, dydr], [dxds, dyds]])

        invJe = np.linalg.inv(Je)
        detJe = np.linalg.det(Je)
        determinante[i] = detJe

    # --------- Condicion de la matriz ---------
    condicion = np.linalg.cond(Ktotal)

    # --------- Salida de Texto ---------

    file = open(RutaTexto, "w")
    file.write("Determinante = ")
    file.write('%d' % np.min(determinante))
    file.write(os.linesep)
    file.write("Condición Matriz Global = ")
    file.write('%d' % condicion)
    file.write(os.linesep)
    file.write("Nº de Nodos = ")
    file.write('%d' % len(nodos.xnod[:, 0]))
    file.write(os.linesep)
    file.write("Nº de Elementos = ")
    file.write('%d' % len(elem[:, 0]))
    file.write(os.linesep)
    file.write("Tiempo de Simulación [seg.] = ")
    file.write('%f' % tiempo)
    file.write(os.linesep)
    file.close()

    return determinante