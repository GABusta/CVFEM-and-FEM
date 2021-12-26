import numpy as np

def BoundaryConditions(KG, RG, nodos, temp):
    # ---- Condiciones en Puntos ----
    l1 = np.array([[1, 2], [temp[1], temp[1]]])
    l3 = np.array([[3, 5, 7, 9, 11, 13, 15, 17],
                  [temp[2], temp[2], temp[2], temp[2], temp[2], temp[2], temp[2], temp[2]]])
    l2 = np.array([[4, 6, 8, 10, 12, 14, 16, 18],
                  [temp[0], temp[0], temp[0], temp[0], temp[0], temp[0], temp[0], temp[0]]])
    condiciones = np.append(l1, l2, axis=1)
    condiciones = np.append(condiciones, l3, axis=1)

    # ---- Condiciones en Bordes ----
    T = np.zeros(20)
    T[0], T[1:10], T[11:20] = temp[1], temp[0], temp[2]

    for i in range(20):
        nodosEncontrados = np.zeros((2, len(nodos.bordes[i][:])))
        if i != 10:
            nodosEncontrados[0, :] = np.array(nodos.bordes[i][:])
            nodosEncontrados[1, :] = T[i]
            condiciones = np.append(condiciones, nodosEncontrados, axis=1)

    arreglo = np.lexsort((condiciones[1, :], condiciones[0, :]))

    # ---- Aplicación de las Condiciones de Contorno ----
    for j in range(len(condiciones[0, :])):
        NodoCondicion = arreglo[j]
        posNodo = int(condiciones[0, NodoCondicion]) - 1
        for k in range(len(nodos.xnod[:, 0])):
            RG[k] -= KG[k, posNodo] * condiciones[1, NodoCondicion]
        KG[posNodo, :] = 0.0
        KG[:, posNodo] = 0.0
        KG[posNodo, posNodo] = 1.0
        RG[posNodo] = condiciones[1, NodoCondicion]

    return KG, RG
