import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import Geom.PyTree as D
import Converter.Internal as I

circle = D.circle((0.5,0,0), 0.5, N=4096)
zones = I.getZones(circle)
circle = I.renameNode(circle, zones[0][0], 'Geometry')
C.convertPyTree2File(circle, 'cylinder_geometry.cgns')