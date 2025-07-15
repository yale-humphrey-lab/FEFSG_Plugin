import numpy as np
import pyvista as pv
import vtk
import xml.etree.ElementTree as ET
import argparse
from scipy.special import erfinv
from scipy.spatial.distance import cdist, squareform
from numpy.linalg import pinv
from scipy.stats import multivariate_normal, norm
import matplotlib.pyplot as plt

def update_geometry(xml_file_path, xml_items1, xml_items2, xml_items3):
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    # Find the Geometry tag
    geometry_tag = root.find(".//Mesh")
    # Remove existing content under Geometry tag
    geometry_tag.clear()
    for item in xml_items1:
        geometry_tag.append(item)

    # Find the Boundary tag
    boundary_tag = root.find(".//Boundary")
    # Remove existing content under Geometry tag
    boundary_tag.clear()
    for item in xml_items2:
        boundary_tag.append(item)

    # Find the Boundary tag
    boundary_tag = root.find(".//MeshData")
    # Remove existing content under Geometry tag
    boundary_tag.clear()
    for item in xml_items3:
        boundary_tag.append(item)


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


# Ensure semi-positive definiteness (SPD)
def nearestSPD(A):
    """Find the nearest symmetric positive definite matrix."""
    B = (A + A.T) / 2
    _, s, V = np.linalg.svd(B)
    H = np.dot(V.T, np.dot(np.diag(s), V))
    A2 = (B + H) / 2
    A3 = (A2 + A2.T) / 2
    if isSPD(A3):
        return A3
    spacing = np.spacing(np.linalg.norm(A))
    I = np.eye(A.shape[0])
    k = 1
    while not isSPD(A3):
        mineig = np.min(np.real(np.linalg.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1
    return A3


def isSPD(A):
    """Check if a matrix is symmetric positive definite."""
    try:
        np.linalg.cholesky(A)
        return True
    except np.linalg.LinAlgError:
        return False


def getRandomInjuryField(T, Z, r):

    overall_perturbed = 0.15
    slope = 0.6
    L_t = 2.0
    L_z = 1.5
    np.random.seed(125)

    mu = 1/2 - 1/slope/np.sqrt(np.pi) * np.exp(-erfinv(1 - 2 * overall_perturbed)**2) * erfinv(1 - 2 * overall_perturbed)
    sigma = 1/slope/np.sqrt(2 * np.pi) * np.exp(-erfinv(1 - 2 * overall_perturbed)**2)

    # Pairwise distances in theta [mm], taking into account the vessel circumference
    D_t = cdist(T.flatten()[:, np.newaxis] + np.pi,T.flatten()[:, np.newaxis] + np.pi)*r
    # To enforce periodicity in the covariance function (using the periodic kernel)
    D_t = 2 * np.sin(np.pi * D_t / (2 * np.pi * r))
    # Pairwise distances in z [mm]
    D_z = cdist(Z.flatten()[:, np.newaxis],Z.flatten()[:, np.newaxis])

    # Overall standardized squared distance metric
    D_squared = (D_t * r / L_t) ** 2 + (D_z / L_z) ** 2

    # Covariance matrix of mechanosensing (M*) values at all surface points
    SIGMA = sigma**2 * np.exp(-0.5 * D_squared)

    # Correlation matrix
    diag_SIGMA = np.sqrt(np.diag(SIGMA))
    CORR = SIGMA / np.outer(diag_SIGMA, diag_SIGMA)

    # Number of standard deviations away to place M* at the boundaries
    n_std = 2

    # Non-boundary indices
    idx1 = np.where(~((Z.flatten() == np.min(Z)) | (Z.flatten() == np.max(Z))))[0]

    # Boundary indices
    idx2 = np.where((Z.flatten() == np.min(Z)) | (Z.flatten() == np.max(Z)))[0]

    # Matrix of "regression coefficients"
    regmat = SIGMA[np.ix_(idx1, idx2)] @ pinv(SIGMA[np.ix_(idx2, idx2)])

    # Conditional mean
    MU_cond = np.full(T.size, np.nan)  # Initialize with NaN values
    MU_cond[idx2] = mu - n_std * sigma
    MU_cond[idx1] = mu + regmat @ (mu - n_std * sigma - mu * np.ones(len(idx2)))
    MU_cond = MU_cond.T  # Ensure it's a row vector

    # Conditional covariance matrix
    SIGMA_cond = np.zeros_like(SIGMA)
    SIGMA_cond[np.ix_(idx1, idx1)] = SIGMA[np.ix_(idx1, idx1)] - regmat @ SIGMA[np.ix_(idx2, idx1)]

    # Apply nearestSPD to ensure positive semi-definiteness
    SIGMA_cond_SPD = nearestSPD(SIGMA_cond)

    # Generate random field (latent Gaussian process)
    M_star = multivariate_normal.rvs(mean=MU_cond, cov=SIGMA_cond_SPD).reshape(T.shape)

    # Transform to match Normal(mu, sigma) exactly using inverse CDF (norminv equivalent)
    # In Python, we use np.histogram for ksdensity approximation and norm.ppf for inverse CDF
    sorted_M_star = np.sort(M_star.flatten())
    
    # Approximate CDF using empirical distribution (like ksdensity)
    cdf_vals = np.linspace(0, 1, len(sorted_M_star), endpoint=False)

    # Apply the inverse CDF (percent-point function)
    M_star_transformed = norm.ppf(cdf_vals, loc=mu, scale=sigma)
    
    # Reshape to the original size
    M_star = np.interp(M_star.flatten(), sorted_M_star, M_star_transformed).reshape(M_star.shape)

    # Mechanosensing field in [0, 1]
    M = np.clip(M_star, 0, 1)

    return M

def getGeometry():
    print("Initializing cylindrical vessel...")

    vesselType = "cylinder" # Can be "torus" or "cylinder"
    torusFraction = 0.25
    torusRadius = 10.0
    numCirc = 20 # Number of circumfrential elements # Must be divisible by 4
    numLen = 20 # Number of axial elements
    numRad = 1 # Number of radial elements
    radius = 6.468e-01 #0.809 #6.468e-01
    thickness = 4.02e-02 #0.041 #4.02e-02
    length = 15. #0.041 #15.0


    # Add aneurysm conditions
    random_field = False

    # Symmetry conditions
    half_circumfrence = False
    quarter_circumfrence = True
    half_length = False

    # Element types (default quadratic hexehdron)
    hex_8 = False
    hex_20 = False
    

    if hex_8 == False:
        numCirc = numCirc*2 #Must be divisible by 4!
        numLen = numLen*2
        numRad = numRad*2

    points = [] #np.empty([(numCirc+1)*(numLen+1)*(numRad+1),3])
    point_ids = {} #[]
    cells = []
    fix_x = []
    fix_y = []
    fix_z = []
    fix_x_quad = []
    fix_y_quad = []
    fix_z_quad = []
    inner_surf = []
    injury_val = []

    num = 1

    maxCirc = numCirc
    if half_circumfrence:
        maxCirc = numCirc//2 + 1
    elif quarter_circumfrence:
        maxCirc = numCirc//4 + 1

    maxLen = numLen + 1

    if torusFraction == 1.0:
        maxLen = numLen

    if half_length:
        maxLen = numLen//2 + 1

    for i in range(maxLen):
        theta = torusFraction * 2 * np.pi * i / numLen  # Azimuthal angle (circumference)
        for j in range(maxCirc):
            phi = -2 * np.pi * j / numCirc
            for k in range(numRad+1):

                if vesselType == "cylinder":
                    xPt = (radius + thickness*k/numRad)*np.cos(2.0*j*np.pi/numCirc)
                    yPt = (radius + thickness*k/numRad)*np.sin(2.0*j*np.pi/numCirc)
                    zPt = length*i/numLen - length/2.0

                if vesselType == "torus":
                    xPt = radius * np.sin(phi) + (np.sin(phi)*(thickness*k/numRad))
                    yPt = (torusRadius + radius * np.cos(phi)) * np.cos(theta) + ((np.cos(phi) * np.cos(theta))*(thickness*k/numRad))
                    zPt = (torusRadius + radius * np.cos(phi)) * np.sin(theta) + ((np.cos(phi) * np.sin(theta))*(thickness*k/numRad))

                points.append([xPt, yPt, zPt])
                point_ids[i*(numCirc+1)*(numRad+1) + j*(numRad+1) + k] = num
                num+=1

                if vesselType == "cylinder":
                    if half_circumfrence:
                        if  (j == 0 or j == 2*numCirc/4):
                            fix_y.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        if (j == 1*numCirc/4 or j == 3*numCirc/4) and (k == 0) and (i == 0 or i == numLen):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                    elif quarter_circumfrence:
                        if  (j == 0 or j == 2*numCirc/4):
                            fix_y.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        if (j == 1*numCirc/4 or j == 3*numCirc/4):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                    else:
                        if  (j == 0 or j == 2*numCirc/4) and (k == 0) and (i == 0 or i == numLen):
                            fix_y.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        if (j == 1*numCirc/4 or j == 3*numCirc/4) and (k == 0) and (i == 0 or i == numLen):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])

                elif vesselType == "torus":

                    if (theta == 0):
                        if (phi == 0.0):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        if (phi == -np.pi):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        fix_z.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])

                    if  (theta == np.pi/2.0):
                        if (phi == 0.0):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        if (phi == -np.pi):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        fix_y.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])


                    if  (theta == np.pi):
                        if (phi == 0.0):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        if (phi == -np.pi):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        fix_z.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])

                    if  (theta == 3.0*np.pi/2.0):
                        if (phi == 0.0):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        if (phi == -np.pi):
                            fix_x.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                        fix_y.append(point_ids[(i)*(numCirc+1)*(numRad+1) + (j)*(numRad+1) + (k)])
                

    hex_8_coords=[
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1]]

    hex_8_quad_coords =  [
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

    hex_20_coords = [
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

    hex_20_quad_coords =  [
            [0, 0],
            [0, 2],
            [2, 2],
            [2, 0],
            [0, 1],
            [1, 2],
            [2, 1],
            [1, 0]]


    maxCirc = numCirc
    numJump = 2
    elem_coords = hex_coords
    quad_coords = hex_quad_coords

    if hex_8 == True:
        numJump = 1
        elem_coords = hex_8_coords
        quad_coords = hex_8_quad_coords

    if hex_20 == True:
        elem_coords = hex_20_coords
        quad_coords = hex_20_quad_coords

    if half_circumfrence:
        maxCirc = numCirc//2
    elif quarter_circumfrence:
        maxCirc = numCirc//4


    maxLen = numLen
    if half_length:
        maxLen = numLen//2


    T, Z = np.meshgrid(np.linspace(0, 2.0*np.pi, maxCirc), np.linspace(-length/2.0, length/2.0, maxLen))
    
    if random_field:
        M = getRandomInjuryField(T, Z, radius)

    for i in range(0, maxLen, numJump):
        for j in range(0, maxCirc, numJump):
            for k in range(0, numRad, numJump):

                cellPts = []

                if torusFraction == 1.0:
                    for coord in elem_coords:
                        cellPts.append(point_ids[((i+coord[1])%(numLen))*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k+coord[2])])
                else:
                    for coord in elem_coords:
                        cellPts.append(point_ids[(i+coord[1])*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k+coord[2])])

                cells.append(cellPts)


                if k == 0:
                    quadPts = []
                    for coord in quad_coords:
                        if torusFraction == 1.0:
                            quadPts.append(point_ids[((i+coord[1])%(numLen))*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k)])
                        else:
                            quadPts.append(point_ids[(i+coord[1])*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k)])
                    inner_surf.append(quadPts)


                if vesselType == "cylinder":
                    if i == 0:
                        quadPts = []
                        for coord in quad_coords:
                            quadPts.append(point_ids[(i)*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k+coord[1])])
                        fix_z_quad.append(quadPts)

                    if i == maxLen-numJump:
                        quadPts = []
                        for coord in quad_coords:
                            quadPts.append(point_ids[(i + numJump)*(numCirc+1)*(numRad+1) + ((j+coord[0])%(numCirc))*(numRad+1) + (k+coord[1])])
                        fix_z_quad.append(quadPts)


                theta = 2.0*(j+numJump/2.)*np.pi/numCirc
                xPt = np.cos(theta)
                yPt = np.sin(theta)
                zPt = length*(i+numJump/2.)/numLen - length/2.0

                if random_field:
                    injury_val.append(M[i,j])



    xml_content1 = []
    xml_content2 = []
    xml_content3 = []

    # Write Nodes
    xml_object = ''
    xml_object += '\t<Nodes name="Object01">\n'
    for i, node in enumerate(points, start=1):
        str_val = ', '.join(map(str, node))
        xml_object += f'\t\t<node id="{i}">{str_val}</node>\n'
    xml_object += '\t</Nodes>\n'
    xml_content1.append(ET.fromstring(xml_object))

    # Write Elements
    xml_object = ''
    if hex_8:
        xml_object += '\t<Elements type="hex8" mat="1" name="Part1">\n'
    elif hex_20:
        xml_object +=  '\t<Elements type="hex20" mat="1" name="Part1">\n'
    else:
        xml_object += '\t<Elements type="hex27" mat="1" name="Part1">\n'
    for i, element in enumerate(cells, start=1):
        str_val = ', '.join(map(str, element))
        xml_object += f'\t\t<elem id="{i}">{str_val}</elem>\n'
    xml_object += '\t</Elements>\n'
    xml_content1.append(ET.fromstring(xml_object))

    if fix_x:
        # Write Surfaces
        xml_object = ''
        xml_object += '\t<NodeSet name="FixXs">\n'
        result_string = ', '.join(str(x) for x in fix_x)
        xml_object += f'\t\t{result_string}\n'
        xml_object += '\t</NodeSet>\n'
        xml_content1.append(ET.fromstring(xml_object))

    if fix_y:
        xml_object = ''
        xml_object += '\t<NodeSet name="FixYs">\n'
        result_string = ', '.join(str(x) for x in fix_y)
        xml_object += f'\t\t{result_string}\n'
        xml_object += '\t</NodeSet>\n'
        xml_content1.append(ET.fromstring(xml_object))

    if fix_z:
        xml_object = ''
        xml_object += '\t<NodeSet name="FixZs">\n'
        result_string = ', '.join(str(x) for x in fix_z)
        xml_object += f'\t\t{result_string}\n'
        xml_object += '\t</NodeSet>\n'
        xml_content1.append(ET.fromstring(xml_object))


    if fix_x_quad:
        xml_object = ''
        xml_object += '\t<Surface name="FixXs">\n'
        for i, surface in enumerate(fix_x_quad, start=1):
            str_val = ', '.join(map(str, surface))
            if hex_8:
                xml_object += f'\t\t<quad4 id="{i}">{str_val}</quad4>\n'
            elif hex_20:
                xml_object += f'\t\t<quad8 id="{i}">{str_val}</quad8>\n'
            else:
                xml_object += f'\t\t<quad9 id="{i}">{str_val}</quad9>\n'
        xml_object += '\t</Surface>\n'
        xml_content1.append(ET.fromstring(xml_object))


    if fix_y_quad:
        xml_object = ''
        xml_object += '\t<Surface name="FixYs">\n'
        for i, surface in enumerate(fix_y_quad, start=1):
            str_val = ', '.join(map(str, surface))
            if hex_8:
                xml_object += f'\t\t<quad4 id="{i}">{str_val}</quad4>\n'
            elif hex_20:
                xml_object += f'\t\t<quad8 id="{i}">{str_val}</quad8>\n'
            else:
                xml_object += f'\t\t<quad9 id="{i}">{str_val}</quad9>\n'
        xml_object += '\t</Surface>\n'
        xml_content1.append(ET.fromstring(xml_object))


    if fix_z_quad:
        xml_object = ''
        xml_object += '\t<Surface name="FixZs">\n'
        for i, surface in enumerate(fix_z_quad, start=1):
            str_val = ', '.join(map(str, surface))
            if hex_8:
                xml_object += f'\t\t<quad4 id="{i}">{str_val}</quad4>\n'
            elif hex_20:
                xml_object += f'\t\t<quad8 id="{i}">{str_val}</quad8>\n'
            else:
                xml_object += f'\t\t<quad9 id="{i}">{str_val}</quad9>\n'
        xml_object += '\t</Surface>\n'
        xml_content1.append(ET.fromstring(xml_object))

    xml_object = ''
    xml_object += '\t<Surface name="PressureLoad1">\n'
    for i, surface in enumerate(inner_surf, start=1):
        str_val = ', '.join(map(str, surface))
        if hex_8:
            xml_object += f'\t\t<quad4 id="{i}">{str_val}</quad4>\n'
        elif hex_20:
            xml_object += f'\t\t<quad8 id="{i}">{str_val}</quad8>\n'
        else:
            xml_object += f'\t\t<quad9 id="{i}">{str_val}</quad9>\n'
    xml_object += '\t</Surface>\n'
    xml_content1.append(ET.fromstring(xml_object))


    if fix_x:
        xml_object = ''
        xml_object += '\t<bc name="FixXs" type="zero displacement" node_set="FixXs">\n'
        xml_object += f'\t\t<x_dof>1</x_dof>\n'
        xml_object += f'\t\t<y_dof>0</y_dof>\n'
        xml_object += f'\t\t<z_dof>0</z_dof>\n'
        xml_object += '\t</bc>\n'
        xml_content2.append(ET.fromstring(xml_object))
    if fix_y:
        xml_object = ''
        xml_object += '\t<bc name="FixYs" type="zero displacement" node_set="FixYs">\n'
        xml_object += f'\t\t<x_dof>0</x_dof>\n'
        xml_object += f'\t\t<y_dof>1</y_dof>\n'
        xml_object += f'\t\t<z_dof>0</z_dof>\n'
        xml_object += '\t</bc>\n'
        xml_content2.append(ET.fromstring(xml_object))
    if fix_z:
        xml_object = ''
        xml_object += '\t<bc name="FixZs" type="zero displacement" node_set="FixZs">\n'
        xml_object += f'\t\t<x_dof>0</x_dof>\n'
        xml_object += f'\t\t<y_dof>0</y_dof>\n'
        xml_object += f'\t\t<z_dof>1</z_dof>\n'
        xml_object += '\t</bc>\n'
        xml_content2.append(ET.fromstring(xml_object))

    if fix_x_quad:
        xml_object = ''
        xml_object += '\t<bc name="FixXs" type="zero displacement" node_set="@surface:FixXs">\n'
        xml_object += f'\t\t<x_dof>1</x_dof>\n'
        xml_object += f'\t\t<y_dof>0</y_dof>\n'
        xml_object += f'\t\t<z_dof>0</z_dof>\n'
        xml_object += '\t</bc>\n'
        xml_content2.append(ET.fromstring(xml_object))
    if fix_y_quad:
        xml_object = ''
        xml_object += '\t<bc name="FixYs" type="zero displacement" node_set="@surface:FixYs">\n'
        xml_object += f'\t\t<x_dof>0</x_dof>\n'
        xml_object += f'\t\t<y_dof>1</y_dof>\n'
        xml_object += f'\t\t<z_dof>0</z_dof>\n'
        xml_object += '\t</bc>\n'
        xml_content2.append(ET.fromstring(xml_object))
    if fix_z_quad:
        xml_object = ''
        xml_object += '\t<bc name="FixZs" type="zero displacement" node_set="@surface:FixZs">\n'
        xml_object += f'\t\t<x_dof>0</x_dof>\n'
        xml_object += f'\t\t<y_dof>0</y_dof>\n'
        xml_object += f'\t\t<z_dof>1</z_dof>\n'
        xml_object += '\t</bc>\n'
        xml_content2.append(ET.fromstring(xml_object))


    """
    if injury_val:
        xml_object = ''
        xml_object += '\t<NodeData name="injury_map" node_set="Object01">\n'
        for i, val in enumerate(injury_val, start=1):
            str_val = str(val)
            xml_object += f'\t\t<node lid="{i}">{str_val}</node>\n'
        xml_object += '\t</NodeData>\n'
        xml_content3.append(ET.fromstring(xml_object))
    """


    if injury_val:
        xml_object = ''
        xml_object += '\t<ElementData name="injury_map" elem_set="Part1">\n'
        for i, val in enumerate(injury_val, start=1):
            str_val = str(val)
            xml_object += f'\t\t<elem lid="{i}">{str_val}</elem>\n'
        xml_object += '\t</ElementData>\n'
        xml_content3.append(ET.fromstring(xml_object))



    return xml_content1, xml_content2, xml_content3

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Update Geometry tag in XML file.')
    parser.add_argument('xml_file', help='Path to the XML file')
    args = parser.parse_args()

    xml_content1, xml_content2, xml_content3 = getGeometry()

    update_geometry(args.xml_file, xml_content1, xml_content2, xml_content3)
