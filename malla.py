
import numpy as np
import gmsh

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)  # GMSH imprime mensajes en la consola
gmsh.model.add("modelo_2")                    # Nombre al modelo
filename = 'trian_02.msh'
tm = 2.0                   # Tamaño de malla a utilizar en todos los puntos
type = 1                # malla triangular = 1, cuadratica = 2

# -----------------------------------------------------------
# ---------------------      Puntos     ---------------------
# -----------------------------------------------------------
xcm, ycm = 107.5, 300.0  # Centro del perfil
p = np.zeros((28, 2))
p[[0, 2, 4, 6, 8], :] = [[0.0, 0.0], [0.0, 11.4], [16.4, 32.4], [86.8, 50.0], [96.7, 62.6]]

for i in range(1, 6, 1):
    c = i * 2 - 1
    p[c, :] = p[c-1, :]
    p[c, 0] = xcm * 2 - p[c-1, 0]

c1 = 9
for i in range(1, 6, 1):
     c = i * 2 - 1
     p[c + c1, :] = p[c1 - c, :]
     p[c + c1, 1] = ycm * 2 - p[c1 - c, 1]

     p[c + c1 + 1, :] = p[c1 - c + 1, :]
     p[c + c1 + 1, 1] = ycm * 2 - p[c1 - c + 1, 1]

r1, r2 = 21.6, 13.0  # radio del círculo
p[[20, 21, 22, 23], :] = [[21.6, 11.4], [xcm * 2 - 21.6, 11.4], [21.6, ycm * 2 - 11.4], [xcm * 2 - 21.6, ycm * 2 - 11.4]]
p[[24, 25, 26, 27], :] = [[83.7, 62.6], [xcm * 2 - 83.7, 62.6], [83.7, ycm * 2 - 62.6], [xcm * 2 - 83.7, ycm * 2 - 62.6]]

# Sintaxis: gmsh.model.geo.addPoint(x, y, z, tm, tag)
for i in range(0, 28, 1):
    gmsh.model.geo.addPoint(p[i, 0],  p[i, 1],  0.0, tm, i + 1)

# -----------------------------------------------------------
# ---------------------      Lineas    ----------------------
# -----------------------------------------------------------
# Sintaxis: gmsh.model.geo.addLine(punto inicial, punto final, tag)
# Sintaxis: gmsh.model.geo.addCircle(p. inicial, p. centro, p. final, tag)
gmsh.model.geo.addLine(1, 2, 1)  # Éste sería el borde inferior
gmsh.model.geo.addLine(2, 4, 2)  # ...
gmsh.model.geo.addCircleArc(4, 22, 6, 3)  # Arco r2
gmsh.model.geo.addLine(6, 8, 4)
gmsh.model.geo.addCircleArc(8, 26, 10, 5)  # Arco r6
gmsh.model.geo.addLine(10, 12, 6)
gmsh.model.geo.addCircleArc(12, 28, 14, 7)  # Arco r8
gmsh.model.geo.addLine(14, 16, 8)
gmsh.model.geo.addCircleArc(16, 24, 18, 9)  # Arco r4
gmsh.model.geo.addLine(18, 20, 10)
gmsh.model.geo.addLine(20, 19, 11)
gmsh.model.geo.addLine(19, 17, 12)
gmsh.model.geo.addCircleArc(17, 23, 15, 13)  # Arco r3
gmsh.model.geo.addLine(15, 13, 14)
gmsh.model.geo.addCircleArc(13, 27, 11, 15)  # Arco r7
gmsh.model.geo.addLine(11, 9, 16)
gmsh.model.geo.addCircleArc(9, 25, 7, 17)  # Arco r5
gmsh.model.geo.addLine(7, 5, 18)
gmsh.model.geo.addCircleArc(5, 21, 3, 19)  # Arco r1
gmsh.model.geo.addLine(3, 1, 20)

# -----------------------------------------------------------
# -------------------      Superficie    --------------------
# -----------------------------------------------------------
# Para definir superficies, se deben definir 'Curve Loops' que las limiten.
# En este caso un Curve Loop será el borde exterior
curvas = np.linspace(1, 20, 20, dtype=int)
gmsh.model.geo.addCurveLoop(curvas[:], 21)

# Sintaxis: gmsh.model.geo.addPlaneSurface([Lista de Curve Loops], tag)
gmsh.model.geo.addPlaneSurface([21], 22)

# Por defecto, si hay grupos físicos definidos, GMSH solo reporta elementos fini-
# tos que pertenezcan a algún grupo físico. En este caso se crearán dos:
#    - Una superficie física que contenga nuestra superficie creada
#    - Una curva física que contenga el borde inferior

# Sintaxis: gmsh.model.addPhysicalGroup(dimensión, lista de entidades, tag)
gmsh.model.addPhysicalGroup(1, [1], 101)  # T2
gmsh.model.setPhysicalName(1, 101, "T2")

gmsh.model.addPhysicalGroup(1, [2, 3, 4, 5, 6, 7, 8, 9, 10], 102)  # T1
gmsh.model.setPhysicalName(1, 102, "T1")

gmsh.model.addPhysicalGroup(1, [11], 103)  # G1
gmsh.model.setPhysicalName(1, 103, "G1")

gmsh.model.addPhysicalGroup(1, [12, 13, 14, 15, 16, 17, 18, 19, 20], 104)  # T3
gmsh.model.setPhysicalName(1, 104, "T3")

s = gmsh.model.addPhysicalGroup(2, [22])  # En este caso no se especifica tag
gmsh.model.setPhysicalName(2, s, "superficie")  # Se puede definir un nombre

# -----------------------------------------------------------
# --------------------      Mallado    ----------------------
# -----------------------------------------------------------
gmsh.model.geo.synchronize()   # emparejar Modelo con malla

if (type == 1):
    gmsh.model.mesh.generate(2)  # Unstructured triangle

elif (type == 2):
    gmsh.model.geo.mesh.setTransfiniteSurface(22)
    gmsh.model.mesh.setRecombine(2, 22)
    gmsh.model.mesh.generate(2)  # Unstructured Quad

gmsh.option.setNumber('Mesh.SurfaceFaces', 1)   # Ver las "caras" de los elementos finitos 2D

# Y finalmente guardar la malla
gmsh.write(filename)

# -----------------------------------------------------------
# ---------------      Mallado visualizar   -----------------
# -----------------------------------------------------------
# Podemos visualizar el resultado en la interfaz gráfica de GMSH
gmsh.fltk.run()
gmsh.finalize()
