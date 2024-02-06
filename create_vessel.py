import numpy as np
import pyvista as pv
import vtk
import xml.etree.ElementTree as ET
import argparse

def update_geometry(xml_file_path, xml_items):
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    # Find the Geometry tag
    geometry_tag = root.find(".//Mesh")

    # Remove existing content under Geometry tag
    geometry_tag.clear()

    for item in xml_items:
        geometry_tag.append(item)

    ET.indent(tree, '  ')
    tree.write(xml_file_path, encoding='utf-8', xml_declaration=True) 

def create_quadratic_hexahedron_mesh(nodes, connectivity):
    # Create points
    points = vtk.vtkPoints()
    for node in nodes:
        points.InsertNextPoint(node)

    # Create the hexahedron cells
    hexahedra = vtk.vtkCellArray()
    for hexa in connectivity:
        hexahedra.InsertNextCell(27, hexa)  # 27 is the number of points in a quadratic hexahedron

    # Create a VTK unstructured grid
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.SetCells(vtk.VTK_TRIQUADRATIC_HEXAHEDRON, hexahedra)

    return grid

def save_vtk_mesh(grid, filename):
    # Write the mesh to a VTK file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(grid)
    writer.Write()

def getAneurysmValue(z, theta):
    """                                                                                                                                                                                                      
    Get the value of vessel behavior based on point location and radius                                                                                                                                      
    """

    zod = 0.3
    zapex = 0.0

    thetaod = np.pi
    thetaapex = np.pi

    vend = 0.1
    vapex = 1.0
    vz = 2.0
    vtheta = 2.0

    vesselValue = (
        vend + (vapex - vend) * np.exp(-np.abs((z - zapex) / zod) ** vz)
        * np.exp(-np.abs((theta - thetaapex) / thetaod) ** vtheta)
    )

    #zPt = point[2]
    #vesselValue = 0.65*np.exp(-abs(zPt/(radius*4.0))**2)
    return vesselValue

def getGeometry():
    print("Initializing cylindrical vessel...")

    numCirc = 20 #Must be divisible by 4!
    numLen = 10
    numRad = 1
    radius = 6.468e-01
    thickness = 4.02e-02
    length = 7.5

    half = False
    quarter = True
    linear = False
    len_half = True

    if linear == False:
        numCirc = numCirc*2 #Must be divisible by 4!
        numLen = numLen*2
        numRad = numRad*2

    points = [] #np.empty([(numCirc+1)*(numLen+1)*(numRad+1),3])
    point_ids = {} #[]
    cells = []
    fix_x = []
    fix_y = []
    fix_z = []
    inner_surf = []
    axis_a = []
    axis_d = []
    aneurysm_val = []

    num = 1


    maxCirc = numCirc
    if half:
        maxCirc = numCirc//2 + 1
    elif quarter:
        maxCirc = numCirc//4 + 1

    maxLen = numLen + 1
    if len_half:
        maxLen = numLen//2 + 1

    for i in range(maxLen):
        for j in range(maxCirc):
            for k in range(numRad+1):

                xPt = (radius + thickness*k/numRad)*np.cos(2.0*j*np.pi/numCirc)
                yPt = (radius + thickness*k/numRad)*np.sin(2.0*j*np.pi/numCirc)

                zPt = length*i/numLen - length/2.0

                points.append([xPt, yPt, zPt])
                point_ids[i*(numCirc+1)*(numRad+1) + j*(numRad+1) + k] = num
                num+=1

                if half:
                    if  (j == 0 or j == 2*numCirc/4):
                        fix_y.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                    if (j == 1*numCirc/4 or j == 3*numCirc/4) and (k == 0) and (i == 0 or i == numLen):
                        fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                elif quarter:
                    if  (j == 0 or j == 2*numCirc/4):
                        fix_y.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                    if (j == 1*numCirc/4 or j == 3*numCirc/4):
                        fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                else:
                    if  (j == 0 or j == 2*numCirc/4) and (k == 0) and (i == 0 or i == numLen):
                        fix_y.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                    if (j == 1*numCirc/4 or j == 3*numCirc/4) and (k == 0) and (i == 0 or i == numLen):
                        fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                

    linear_coords=[
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1]]

    linear_quad_coords =  [
            [0, 0],
            [0, 1],
            [1, 1],
            [1, 0]]

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
            [0, 2, 1],
            [1, 0, 1],
            [2, 1, 1],
            [1, 2, 1],
            [0, 1, 1],
            [1, 1, 0],
            [1, 1, 2],
            [1, 1, 1]]

    hex_quad_coords =  [
            [0, 0],
            [0, 2],
            [2, 2],
            [2, 0],
            [0, 1],
            [1, 2],
            [2, 1],
            [1, 0],
            [1, 1]]


    maxCirc = numCirc
    numJump = 2
    elem_coords = hex_coords
    quad_coords = hex_quad_coords

    if linear == True:
        numJump = 1
        elem_coords = linear_coords
        quad_coords = linear_quad_coords

    if half:
        maxCirc = numCirc//2
    elif quarter:
        maxCirc = numCirc//4


    maxLen = numLen
    if len_half:
        maxLen = numLen//2

    for i in range(0, maxLen, numJump):
        for j in range(0, maxCirc, numJump):
            for k in range(0, numRad, numJump):

                cellPts = []

                for coord in elem_coords:
                    cellPts.append(point_ids[(i+coord[1])*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k+coord[2])])

                cells.append(cellPts)

                if i == 0:
                    quadPts = []
                    for coord in quad_coords:
                        quadPts.append(point_ids[(i)*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k+coord[1])])
                    fix_z.append(quadPts)

                if i == maxLen-numJump:
                    quadPts = []
                    for coord in quad_coords:
                        quadPts.append(point_ids[(i + numJump)*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k+coord[1])])
                    fix_z.append(quadPts)

                if k == 0:
                    quadPts = []
                    for coord in quad_coords:
                        quadPts.append(point_ids[(i+coord[1])*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k)])
                    inner_surf.append(quadPts)

                theta = 2.0*(j+numJump/2.)*np.pi/numCirc

                xPt = np.cos(theta)
                yPt = np.sin(theta)

                zPt = length*(i+numJump/2.)/numLen - length/2.0
                aneurysm_val.append(getAneurysmValue(zPt,theta))


    xml_content = []

    # Write Nodes
    xml_object = ''
    xml_object += '\t<Nodes name="Object01">\n'
    for i, node in enumerate(points, start=1):
        str_val = ', '.join(map(str, node))
        xml_object += f'\t\t<node id="{i}">{str_val}</node>\n'
    xml_object += '\t</Nodes>\n'
    xml_content.append(ET.fromstring(xml_object))

    # Write Elements
    xml_object = ''
    if linear:
        xml_object += '\t<Elements type="hex8" mat="1" name="Part1">\n'
    else:
        xml_object += '\t<Elements type="hex27" mat="1" name="Part1">\n'
    for i, element in enumerate(cells, start=1):
        str_val = ', '.join(map(str, element))
        xml_object += f'\t\t<elem id="{i}">{str_val}</elem>\n'
    xml_object += '\t</Elements>\n'
    xml_content.append(ET.fromstring(xml_object))

    # Write Surfaces
    xml_object = ''
    xml_object += '\t<NodeSet name="FixXs">\n'
    result_string = ', '.join(str(x) for x in fix_x)
    xml_object += f'\t\t{result_string}\n'
    xml_object += '\t</NodeSet>\n'
    xml_content.append(ET.fromstring(xml_object))


    xml_object = ''
    xml_object += '\t<NodeSet name="FixYs">\n'
    result_string = ', '.join(str(x) for x in fix_y)
    xml_object += f'\t\t{result_string}\n'
    xml_object += '\t</NodeSet>\n'
    xml_content.append(ET.fromstring(xml_object))


    xml_object = ''
    xml_object += '\t<Surface name="FixZs">\n'
    for i, surface in enumerate(fix_z, start=1):
        str_val = ', '.join(map(str, surface))
        if linear:
            xml_object += f'\t\t<quad4 id="{i}">{str_val}</quad4>\n'
        else:
            xml_object += f'\t\t<quad9 id="{i}">{str_val}</quad9>\n'
    xml_object += '\t</Surface>\n'
    xml_content.append(ET.fromstring(xml_object))


    xml_object = ''
    xml_object += '\t<Surface name="PressureLoad1">\n'
    for i, surface in enumerate(inner_surf, start=1):
        str_val = ', '.join(map(str, surface))
        if linear:
            xml_object += f'\t\t<quad4 id="{i}">{str_val}</quad4>\n'
        else:
            xml_object += f'\t\t<quad9 id="{i}">{str_val}</quad9>\n'
    xml_object += '\t</Surface>\n'
    xml_content.append(ET.fromstring(xml_object))

    return xml_content

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Update Geometry tag in XML file.')
    parser.add_argument('xml_file', help='Path to the XML file')
    args = parser.parse_args()

    update_geometry(args.xml_file, getGeometry())
