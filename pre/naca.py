import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import Geom.PyTree as D
import Converter.Internal as I

cellSizes = [1, 10, 20, 40, 80, 160]
for cellsize in cellSizes:
    print(f'Generating mesh with cell size = {cellsize} ...')
    size = 1 / cellsize
    mesh11 = G.cartr2((0,0,0), (size,size,0), (1.,1.,0.0), (5.,5.,0.))
    mesh12 = G.cartr2((5.,0,0), (size,size,0), (2.,1.0,0.0), (150.,5.,0.))
    mesh21 = G.cartr2((0,5.,0), (size,size,0), (1.0,2.,0.0), (5.,150.,0.))
    mesh22 = G.cartr2((5.,5.,0), (size,size,0), (2.,2.,0.0), (150.,150.,0.))

    mesh11_12 = T.join(mesh11, mesh12)
    mesh21_22 = T.join(mesh21, mesh22)
    mesh1 = T.join(mesh11_12, mesh21_22)
    mesh2 = T.rotate(mesh1, (0,0,0), (0,0,1), 90)
    mesh3 = T.rotate(mesh1, (0,0,0), (0,0,1), 180)
    mesh4 = T.rotate(mesh1, (0,0,0), (0,0,1), 270)
    mesh5 = T.join(mesh1, mesh2)
    mesh6 = T.join(mesh3, mesh4)
    mesh = T.join(mesh5, mesh6)
    mesh = T.reorder(mesh, (-2,-1,3))
    mesh = T.translate(mesh, (0.5,0.,0.))

    # I.printTree(mesh)

    zones = I.getZones(mesh)
    print("Zones before renaming:")
    for z in zones:
        print("  ", z[0])

    for i, z in enumerate(zones):
        mesh = I.renameNode(mesh, z[0], f'Block_{i+1}')
    C.convertPyTree2File(mesh, f'cylinder/cartesian_mesh_{cellsize}.cgns')


# cellSizes = {0.1: '0-1', 0.05: '0-05', 0.01: '0-01', 0.005: '0-005', 0.0005: '0-0005'}
# for cellsize in cellSizes.keys():
#     mesh1 = G.cartr2((0,0,0), (cellsize,cellsize,0), (1.015,1.015,0.0), (150.,150.,0.))
#     mesh2 = G.cartr2((0,0,0), (cellsize,cellsize,0), (1.015,1.015,0.0), (150.,150.,0.))
#     mesh2 = T.rotate(mesh2, (0,0,0), (0,0,1), 90)
#     mesh3 = T.join(mesh1, mesh2)
#     mesh4 = T.rotate(mesh3, (0,0,0), (0,0,1), 180)
#     mesh = T.join(mesh3, mesh4)
#     mesh = T.reorder(mesh, (-2,-1,3))
#     zones = I.getZones(mesh)
#     for i, z in enumerate(zones):
#         mesh = I.renameNode(mesh, z[0], f'Block_{i+1}')
#     C.convertPyTree2File(mesh, f'cartesian_mesh_{cellSizes[cellsize]}.cgns')

# naca = D.naca('0012', N=10001)
# zones = I.getZones(naca)
# naca = I.renameNode(naca, zones[0][0], 'Geometry')
# C.convertPyTree2File(naca, 'naca0012_geometry.cgns')