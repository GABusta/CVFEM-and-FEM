
import numpy as np

# Nodos = ObjetoNodo(puntos, puntos_xnod, bordes, bordes_xnod, sup, sup_xnod, xnod)
class ObjetoNodo:
    def __init__(self, puntos, puntos_xnod, bordes, bordes_xnod, sup, sup_xnod, xnod):
        self.puntos = puntos
        self.puntos_xnod = puntos_xnod
        self.bordes = bordes
        self.bordes_xnod = bordes_xnod
        self.sup = sup
        self.sup_xnod = sup_xnod
        self.xnod = xnod

# Nodes to Control Volume
class Nodo_VC:
    def __init__(self, ID, ElemNode, posElem):
        self.ID = ID
        self.ElemNode = ElemNode
        self.posElem = posElem


def mallado(archivo, tipo, dim=2):
    # ----------------------------------------------------------------------------
    # ---------------------  Lectura de Archivo  ---------------------------------
    # ----------------------------------------------------------------------------
    m = open(archivo)
    malla_ini = m.readlines()
    malla_ini = [linea.rstrip('\n') for linea in malla_ini]

    # --- Se corta el archivo ---
    inicio_nodos, fin_nodos = malla_ini.index('$Nodes'), malla_ini.index('$EndNodes')
    inicio_elem, fin_elem = malla_ini.index('$Elements'), malla_ini.index('$EndElements')

    malla = malla_ini[inicio_nodos + 1:fin_nodos]  # NODOS
    malla2 = malla_ini[inicio_elem + 1:fin_elem]   # ELEMENTOS

    # ----------------------------------------------------------------------------
    # ----------------------  Lectura de Nodos  ----------------------------------
    # ----------------------------------------------------------------------------
    nblocks, nnodos = [int(n) for n in malla[0].split()[0:2]]

    #  -----     puntos    -----------    bordes      -----------       superfice
    nodos_puntos, xnod_puntos, nodos_bordes, xnod_bordes, nodos_superf, xnod_superf= [], [], [], [], [], []

    for j in range(1, len(malla)):
        # --- Lectura de Línea de Nodos ---
        line = malla[j]
        if len(line.split()) == 4:  # Reporte de cada bloque
            tipo_ent = int(line.split()[0])                 # Punto, borde o superf.
            nno_bloque = int(line.split()[-1])              # No. de nodos del bloque
            nodos_bloque = malla[j + 1:j + 1 + nno_bloque]  # Lista de nodos del bloque
            nodos_bloque = [int(n) for n in nodos_bloque]   # Se convierten a enteros

            xnod_b = malla[j + 1 + nno_bloque:j + 1 + 2 * nno_bloque]

            # --- Se guarda la Linea como "matriz" dinamica ---
            xnod_bloque = []
            for l in xnod_b:
                coord = [float(n) for n in l.split()]
                xnod_bloque.append(coord)

            # --- Se guarda información (punto, línea o superficie) ---
            if tipo_ent == 0:
                nodos_puntos.append(nodos_bloque)
                xnod_puntos.append(xnod_bloque)
            elif tipo_ent == 1:
                nodos_bordes.append(nodos_bloque)
                xnod_bordes.append(xnod_bloque)
            elif tipo_ent == 2:
                nodos_superf.append(nodos_bloque)
                xnod_superf.append(xnod_bloque)

    # --- Ensamble Xnod ---
    xnod = np.empty((nnodos, 3))

    # Primero se adicionan los nodos correspondientes a los puntos:
    for i in range(len(nodos_puntos)):
        idx = np.array(nodos_puntos[i]) - 1
        xnod[idx, :] = np.array(xnod_puntos[i])

    # Luego se agregan los nodos correspondientes a los bordes:
    for i in range(len(nodos_bordes)):
        idx = np.array(nodos_bordes[i]) - 1
        xnod[idx, :] = np.array(xnod_bordes[i])

    # Finalmente se agregan los nodos correspondientes a la superficie:
    for i in range(len(nodos_superf)):
        idx = np.array(nodos_superf[i]) - 1
        xnod[idx, :] = np.array(xnod_superf[i])

    Nodos = ObjetoNodo(nodos_puntos, xnod_puntos, nodos_bordes, xnod_bordes, nodos_superf, xnod_superf, xnod)

    # ----------------------------------------------------------------------------
    # -------------------  Lectura de Elementos  ---------------------------------
    # ----------------------------------------------------------------------------
    nblocks, nelem = [int(n) for n in malla2[0].split()[0:2]]
    matrices = []   # Guardará las matrices LaG asociadas a cada entidad
    tags = []       # Tag correspondiente a cada entidad
    fila = 1        # Contador de filas recorridas
    nmat = 1
    for i in range(nblocks):
        linea = malla2[fila]

        dim, tag, tipo_ef, nef = [int(n) for n in linea.split()]
        if dim == 2:  # solo se toman nodos pertenecientes a superficies
            tags.extend([nmat] * nef)
            nmat += 1
            matrices.append(malla2[fila + 1:fila + 1 + nef])
        fila += (1 + nef)

    Elem = []  # Se inicializa la matriz LaG

    for matriz in matrices:
        for i in range(len(matriz)):
            linea = [int(n) for n in matriz[i].split()]
            Elem.append(linea)
    Elem = np.array(Elem)

    # Se reporta la superficie a la que pertenece cada elemento finito
    Elem[:, 0] = np.array(tags)

    # ----------------------------------------------------------------------------
    # ------------------  Nodos to Control Volume  -------------------------------
    # ----------------------------------------------------------------------------
    nodoID, ElemNodo, NodoPosElem = [], [], []
    for i in range(len(Nodos.xnod[:, 0])):
        elementos = np.array(Elem[:, 1:4]) if tipo == 1 else np.array(Elem[:, 1:5])
        posE_N = np.argwhere(elementos == i+1)

        #ElemNodo, NodoPosElem = np.append(posE_N[:, 0]), posE_N[:, 1]
        nodoID.append(i)
        ElemNodo.append(posE_N[:, 0])
        NodoPosElem.append(posE_N[:, 1])
    Elem_To_nodo  = Nodo_VC(nodoID, ElemNodo, NodoPosElem)

    return Nodos, Elem, Elem_To_nodo
