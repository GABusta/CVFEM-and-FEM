from leer_GMSH import mallado
from CVFEM22 import VolumenControl
import time
from CondicionContorno import BoundaryConditions
import numpy as np
from FEM import FiniteElement
from OrdenamientoElemental import elemOrdenado
from Graficas import ploteos
from CaracteristicasMatriz import matrices

ruta1 = ['/Users/Gustavo/Downloads/DOCTORADO/CURSADAS/ANAV/TP.N2/codigos/MalladoTriangular_struc',
         '/Users/Gustavo/Downloads/DOCTORADO/CURSADAS/ANAV/TP.N2/codigos/MalladoTriangular_unstruc',
         '/Users/Gustavo/Downloads/DOCTORADO/CURSADAS/ANAV/TP.N2/codigos/MalladoCuad_struc',
         '/Users/Gustavo/Downloads/DOCTORADO/CURSADAS/ANAV/TP.N2/codigos/MalladoCuad_unstruc']

mallasCreadas = [['/trian_10.msh', '/trian_05.msh', '/trian_07.msh', '/trian_03.msh'],
                 ['/trian_10.msh', '/trian_05.msh', '/trian_07.msh', '/trian_03.msh'],
                 ['/cuad_10.msh', '/cuad_05.msh', '/cuad_07.msh', '/cuad_03.msh'],
                 ['/cuad_10.msh', '/cuad_05.msh', '/cuad_07.msh', '/cuad_03.msh']]

dim = 2
k = 50.0                             # constante difusividad
Temp = np.array([18.0, 20.0, 16.0])  # Condiciones de contorno
for i in range(4):      # ----- Rutas
    for j in range(4):  # ----- mallas
        arr = mallasCreadas[i][j]
        archivo = ruta1[i] + arr
        directorio = arr[1:len(arr)-4]
        tipo = 1 if i <= 1 else 2    # tipo 1 = triang   ,    tipo 2 = cuad.

        nodos, elem, nodosVC = mallado(archivo, tipo, dim)

        # ------ Corrección Normal Elemental ------
        # elem, normales = elemOrdenado(elem, nodos.xnod)

        # ----------------------------
        # ------ CVFEM solution ------
        # ----------------------------
        TiempoInicio = time.time()
        estudio = '/CVFEM'
        Ktotal, Rtotal = VolumenControl(nodos, elem, nodosVC, k, tipo)
        Ktotal, Rtotal = BoundaryConditions(Ktotal, Rtotal, nodos, Temp)
        solucion = np.dot(np.linalg.inv(Ktotal), Rtotal)
        TiempoFin = time.time()
        TiempoSimulacion = TiempoFin - TiempoInicio
        RutaTexto = ruta1[i] + estudio + '-' + directorio + '.txt'
        det = matrices(nodos, elem, tipo, Ktotal, RutaTexto, TiempoSimulacion)
        # ------ Graficas ------
        b = ploteos(nodos, solucion, elem, tipo, ruta1[i], directorio, estudio)

        # ----------------------------
        # ------- FEM solution -------
        # ----------------------------
        TiempoInicio = time.time()
        estudio = '/FEM'
        Ktotal, Rtotal = FiniteElement(nodos, elem, k, tipo)
        Ktotal, Rtotal = BoundaryConditions(Ktotal, Rtotal, nodos, Temp)
        solucion = np.dot(np.linalg.inv(Ktotal), Rtotal)
        TiempoFin = time.time()
        TiempoSimulacion = TiempoFin - TiempoInicio
        RutaTexto = ruta1[i] + estudio + '-' + directorio + '.txt'
        det = matrices(nodos, elem, tipo, Ktotal, RutaTexto, TiempoSimulacion)
        # ------ Graficas ------
        b = ploteos(nodos, solucion, elem, tipo, ruta1[i], directorio, estudio)
