import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

def ploteos(nodos, solucion, elem, tipo, ruta, directorio, estudio):
    # --------------------------------
    # ------ TRIANGULAR ELEMENT ------
    # --------------------------------
    if tipo == 1:
        # --------- Gráfica del Mallado ---------
        grid = plt.GridSpec(2, 2, wspace=0.008, hspace=0.2)
        plot1 = plt.subplot(grid[0:, 0])
        plot2 = plt.subplot(grid[0, 1])
        plot3 = plt.subplot(grid[1, 1])
        for i in range(len(elem)):
            plot1.fill(nodos.xnod[elem[i, 1:4] - 1, 0], nodos.xnod[elem[i, 1:4] - 1, 1],
                       edgecolor='black', fill=False, linewidth=0.2)
            plot2.fill(nodos.xnod[elem[i, 1:4] - 1, 0], nodos.xnod[elem[i, 1:4] - 1, 1],
                       edgecolor='black', fill=False, linewidth=0.2)
            plot3.fill(nodos.xnod[elem[i, 1:4] - 1, 0], nodos.xnod[elem[i, 1:4] - 1, 1],
                       edgecolor='black', fill=False, linewidth=0.2)

        plot1.set_aspect('equal', adjustable='box')
        plot1.set_ylabel('Distancia[mm]')
        plot1.set_ylim(0, 600)

        plot2.set_aspect('equal', adjustable='box')
        plot2.set_ylim(450, 600)
        plot3.set_aspect('equal', adjustable='box')
        plot3.set_ylim(0, 150)

        name = ruta + estudio + '- mallado - ' + directorio + '.png'
        plt.savefig(name, dpi=600)

        # --------- Distribución de Temperaturas ---------
        grid1 = plt.GridSpec(2, 2, wspace=0.008, hspace=0.2)
        plot4 = plt.subplot(grid1[0:, 1])
        plot5 = plt.subplot(grid1[0, 0])
        plot6 = plt.subplot(grid1[1, 0])

        X = np.array(nodos.xnod[:, 0])
        Y = np.array(nodos.xnod[:, 1])
        triangulos = elem[:, 1:4] - 1

        trian = mtri.Triangulation(X, Y, triangulos)
        Z = np.asarray(solucion)

        plot5.tricontourf(trian, Z)
        plot5.set_aspect('equal', adjustable='box')
        plot5.set_ylim(450, 600)
        plot6.tricontourf(trian, Z)
        plot6.set_aspect('equal', adjustable='box')
        plot6.set_ylim(0, 150)

        a1 = plot4.tricontourf(trian, Z)
        plot4.set_aspect('equal', adjustable='box')
        plot4.set_ylabel('Distancia[mm]')
        plt.colorbar(ax=plot4, mappable = a1)

        name = ruta + estudio + '- Temperaturas - ' + directorio + '.png'
        plt.savefig(name, dpi=600)

    # ---------------------------
    # ------ CUAD. ELEMENT ------
    # ---------------------------
    else:
        # --------- Gráfica del Mallado ---------
        grid = plt.GridSpec(2, 2, wspace=0.008, hspace=0.2)
        plot1 = plt.subplot(grid[0:, 0])
        plot2 = plt.subplot(grid[0, 1])
        plot3 = plt.subplot(grid[1, 1])
        for i in range(len(elem[:, 0])):
            # plot1.plot(nodos.xnod[elem[i, 1:5] - 1, 0], nodos.xnod[elem[i, 1:5] - 1, 1],
            #           'k', linewidth=0.2)
            plot1.fill(nodos.xnod[elem[i, 1:5] - 1, 0], nodos.xnod[elem[i, 1:5] - 1, 1],
                       edgecolor='black', fill=False, linewidth=0.2)
            plot2.fill(nodos.xnod[elem[i, 1:5] - 1, 0], nodos.xnod[elem[i, 1:5] - 1, 1],
                       edgecolor='black', fill=False, linewidth=0.2)
            plot3.fill(nodos.xnod[elem[i, 1:5] - 1, 0], nodos.xnod[elem[i, 1:5] - 1, 1],
                       edgecolor='black', fill=False, linewidth=0.2)

        plot1.set_aspect('equal', adjustable='box')
        plot1.set_ylabel('Distancia[mm]')
        plot1.set_ylim(0, 600)

        plot2.set_aspect('equal', adjustable='box')
        plot2.set_ylim(450, 600)
        plot3.set_aspect('equal', adjustable='box')
        plot3.set_ylim(0, 150)

        name = ruta + estudio + '- mallado - ' + directorio + '.png'
        plt.savefig(name, dpi=600)

        # --------- Distribución de Temperaturas ---------
        grid1 = plt.GridSpec(2, 2, wspace=0.008, hspace=0.2)
        plot4 = plt.subplot(grid1[0:, 1])
        plot5 = plt.subplot(grid1[0, 0])
        plot6 = plt.subplot(grid1[1, 0])

        Elementos_triangulos = quads_to_tri(elem)

        X = np.array(nodos.xnod[:, 0])
        Y = np.array(nodos.xnod[:, 1])
        trian = mtri.Triangulation(X, Y, Elementos_triangulos)

        # plot_mesh(X, Y, elem)
        Z = np.asarray(solucion)

        plot5.tricontourf(trian, Z)
        plot5.set_aspect('equal', adjustable='box')
        plot5.set_ylim(450, 600)
        plot6.tricontourf(trian, Z)
        plot6.set_aspect('equal', adjustable='box')
        plot6.set_ylim(0, 150)

        a1 = plot4.tricontourf(trian, Z)
        plot4.set_aspect('equal', adjustable='box')
        plot4.set_ylabel('Distancia[mm]')
        plt.colorbar(ax=plot4, mappable=a1)

        name = ruta + estudio + '- Temperaturas - ' + directorio + '.png'
        plt.savefig(name, dpi=600)

    a = 1
    return a


# ---------------------------------------------------
# ----------- Conversión quad a triangulo -----------
# ---------------------------------------------------
def quads_to_tri(elem):
    tris = [[None for j in range(3)] for i in range(2*len(elem[:, 0]))]
    for i in range(len(elem[:, 0])):
        j = 2 * i
        n0 = elem[i, 1]
        n1 = elem[i, 2]
        n2 = elem[i, 3]
        n3 = elem[i, 4]
        tris[j][0] = n0 - 1
        tris[j][1] = n1 - 1
        tris[j][2] = n2 - 1
        tris[j + 1][0] = n2 - 1
        tris[j + 1][1] = n3 - 1
        tris[j + 1][2] = n0 - 1
    return tris
