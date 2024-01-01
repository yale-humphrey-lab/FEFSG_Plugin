import numpy as np
import pyvista as pv
import vtk

def create_quadratic_hexahedron_mesh(nodes, connectivity):
    # Create points
    points = vtk.vtkPoints()
    for node in nodes:
        points.InsertNextPoint(node)

    # Create the hexahedron cells
    hexahedra = vtk.vtkCellArray()
    for hexa in connectivity:
        hexahedra.InsertNextCell(20, hexa)  # 20 is the number of points in a quadratic hexahedron

    # Create a VTK unstructured grid
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.SetCells(vtk.VTK_QUADRATIC_HEXAHEDRON, hexahedra)

    return grid

def save_vtk_mesh(grid, filename):
    # Write the mesh to a VTK file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(grid)
    writer.Write()


print("Initializing cylindrical vessel...")

numCirc = 20*2 #Must be divisible by 4!
numLen = 1*2
numRad = 4*2
radius = 6.47e-01
thickness = 4.10e-02
length = 6.47e-02

points = np.empty([(numCirc+1)*(numLen+1)*(numRad+1),3])
cells = []
point_ids = []
fix_x = []
fix_y = []
fix_z = []
inner_surf = []
axis_a = []
axis_d = []

num = 1

for i in range(numLen+1):
    for j in range(numCirc+1):
        for k in range(numRad+1):

            if j % 2 == 1:
                xPt1 = (radius + thickness*k/numRad)*np.cos(2.0*(j-1)*np.pi/numCirc)
                yPt1 = (radius + thickness*k/numRad)*np.sin(2.0*(j-1)*np.pi/numCirc)

                xPt2 = (radius + thickness*k/numRad)*np.cos(2.0*(j+1)*np.pi/numCirc)
                yPt2 = (radius + thickness*k/numRad)*np.sin(2.0*(j+1)*np.pi/numCirc)

                xPt = (xPt1 + xPt2)/2.0
                yPt = (yPt1 + yPt2)/2.0

            else:
                xPt = (radius + thickness*k/numRad)*np.cos(2.0*j*np.pi/numCirc)
                yPt = (radius + thickness*k/numRad)*np.sin(2.0*j*np.pi/numCirc)

            zPt = length*i/numLen - length/2.0

            points[i*(numCirc+1)*(numRad+1) + j*(numRad+1) + k,:] = [xPt, yPt, zPt]
            point_ids.append(num)
            num+=1

coords=[[0, 1, 0],
        [0, 1, 1],
        [1, 1, 1],
        [1, 1, 0],
        [0, 0, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 0, 0]]

hex_coords = [
        [0, 0, 0],
        [2, 0, 0],
        [2, 2, 0],
        [0, 2, 0],
        [0, 0, 2],
        [2, 0, 2],
        [2, 2, 2],
        [0, 2, 2],
        [1, 0, 0],
        [2, 1, 0],
        [1, 2, 0],
        [0, 1, 0],
        [1, 0, 2],
        [2, 1, 2],
        [1, 2, 2],
        [0, 1, 2],
        [0, 0, 1],
        [2, 0, 1],
        [2, 2, 1],
        [0, 2, 1]]

quad_coords =  [
        [0, 0],
        [0, 2],
        [2, 2],
        [2, 0],
        [0, 1],
        [1, 2],
        [2, 1],
        [1, 0]]


for i in range(0, numLen, 2):
    for j in range(0, numCirc, 2):
        for k in range(0, numRad, 2):

            cellPts = []

            for coord in hex_coords:
                cellPts.append(point_ids[(i+coord[1])*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k+coord[2])])

            cells.append(cellPts)

            cellCenter = points[(i+1)*(numCirc+1)*(numRad+1) + (j+1)*(numRad+1) + (k+1),:]

            if  j == 0 or j == 2*numCirc/4:
                quadPts = []
                for coord in quad_coords:
                    quadPts.append(point_ids[(i+coord[0])*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k+coord[1])])
                fix_y.append(quadPts)

            if j == 1*numCirc/4 or j == 3*numCirc/4:
                quadPts = []
                for coord in quad_coords:
                    quadPts.append(point_ids[(i+coord[0])*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k+coord[1])])
                fix_x.append(quadPts)


            if i == 0:
                quadPts = []
                for coord in quad_coords:
                    quadPts.append(point_ids[(i)*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k+coord[1])])
                fix_z.append(quadPts)

            if i == numLen-2:
                quadPts = []
                for coord in quad_coords:
                    quadPts.append(point_ids[(i + 2)*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k+coord[1])])
                fix_z.append(quadPts)

            if k == 0:
                quadPts = []
                for coord in quad_coords:
                    quadPts.append(point_ids[(i+coord[1])*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k)])
                inner_surf.append(quadPts)


            theta = 2.0*(j+1)*np.pi/numCirc

            xPt = np.cos(theta)
            yPt = np.sin(theta)

            axis_a.append([xPt,yPt,0])
            axis_d.append([-yPt,xPt,0])


with open("output_file.xml", "w") as f:
    f.write('<Geometry>\n')

    # Write Nodes
    f.write('\t<Nodes name="Object01">\n')
    for i, node in enumerate(points, start=1):
        str_val = ', '.join(map(str, node))
        f.write(f'\t\t<node id="{i}">{str_val}</node>\n')
    f.write('\t</Nodes>\n')

    # Write Elements
    f.write('\t<Elements type="hex20" mat="1" name="Part1">\n')
    for i, element in enumerate(cells, start=1):
        str_val = ', '.join(map(str, element))
        f.write(f'\t\t<elem id="{i}">{str_val}</elem>\n')
    f.write('\t</Elements>\n')

    # Write Surfaces
    f.write('\t<Surface name="FixXs">\n')
    for i, surface in enumerate(fix_x, start=1):
        str_val = ', '.join(map(str, surface))
        f.write(f'\t\t<quad8 id="{i}">{str_val}</quad8>\n')
    f.write('\t</Surface>\n')

    f.write('\t<Surface name="FixYs">\n')
    for i, surface in enumerate(fix_y, start=1):
        str_val = ', '.join(map(str, surface))
        f.write(f'\t\t<quad8 id="{i}">{str_val}</quad8>\n')
    f.write('\t</Surface>\n')

    f.write('\t<Surface name="FixZs">\n')
    for i, surface in enumerate(fix_z, start=1):
        str_val = ', '.join(map(str, surface))
        f.write(f'\t\t<quad8 id="{i}">{str_val}</quad8>\n')
    f.write('\t</Surface>\n')

    f.write('\t<Surface name="PressureLoad1">\n')
    for i, surface in enumerate(inner_surf, start=1):
        str_val = ', '.join(map(str, surface))
        f.write(f'\t\t<quad8 id="{i}">{str_val}</quad8>\n')
    f.write('\t</Surface>\n')
    f.write('</Geometry>\n')
    f.write('\n')


    f.write('\t<ElementData var="mat_axis" elem_set="Part1">\n')

    for i, surface in enumerate(axis_a):
        str_val1 = ', '.join(map(str, axis_a[i]))
        str_val2 = ', '.join(map(str, axis_d[i]))
        f.write(f'\t\t<elem lid="{i+1}">\n')
        f.write(f'\t\t\t<a>{str_val1}</a>\n')
        f.write(f'\t\t\t<d>{str_val2}</d>\n')
        f.write(f'\t\t</elem>\n')


    f.write('\t</ElementData>\n')


# Create the quadratic hexahedron mesh
quadratic_hex_mesh = create_quadratic_hexahedron_mesh(points, np.array(cells)-1)

# Save the mesh to a VTK file
save_vtk_mesh(quadratic_hex_mesh, 'quadratic_hex_mesh.vtu')