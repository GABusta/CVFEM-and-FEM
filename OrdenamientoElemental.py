import numpy as np
# ---------------------------------------------------
# --------- Corrección de la numeración -------------
# ---------------------------------------------------
def elemOrdenado(Ne, coord):
    # ---------------------------------
    # ------ Elemento Triangular ------
    # ---------------------------------
    normal = np.zeros((len(Ne[:, 0]), 3))
    for i in range(len(Ne[:, 0])):
        # Cálculo de las normales elementales
        p1 = np.array([coord[Ne[i, 1] - 1, 0], coord[Ne[i, 1] - 1, 1]])
        p2 = np.array([coord[Ne[i, 2] - 1, 0], coord[Ne[i, 2] - 1, 1]])
        p3 = np.array([coord[Ne[i, 3] - 1, 0], coord[Ne[i, 3] - 1, 1]])

        v1 = p2 - p1
        v2 = p3 - p2
        v3 = p1 - p3

        normal[i, :] = np.array([np.cross(v1, v2), np.cross(v2, v3), np.cross(v3, v1)])

        nn = Ne[i, 1:4]
        if normal[i, 0] < 0.0 or normal[i, 1] < 0.0 or normal[i, 2] < 0.0:
            Ne[i, 1:4] = nn[::-1]

        # -----------------------------------
        # ------ Elemento Cuadrangular ------
        # -----------------------------------
        normal = np.zeros((len(Ne[:, 0]), 4))
        for i in range(len(Ne[:, 0])):
            # Cálculo de las normales elementales
            p1 = np.array([coord[Ne[i, 1] - 1, 0], coord[Ne[i, 1] - 1, 1]])
            p2 = np.array([coord[Ne[i, 2] - 1, 0], coord[Ne[i, 2] - 1, 1]])
            p3 = np.array([coord[Ne[i, 3] - 1, 0], coord[Ne[i, 3] - 1, 1]])
            p3 = np.array([coord[Ne[i, 4] - 1, 0], coord[Ne[i, 4] - 1, 1]])

            v1 = p2 - p1
            v2 = p3 - p2
            v3 = p1 - p3

            normal[i, :] = np.array([np.cross(v1, v2), np.cross(v2, v3), np.cross(v3, v1)])

            nn = Ne[i, 1:5]
            if normal[i, 0] < 0.0 or normal[i, 1] < 0.0 or normal[i, 2] < 0.0:
                Ne[i, 1:5] = nn[::-1]
    return Ne, normal
