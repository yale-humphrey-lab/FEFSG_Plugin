import math

def generate_circle_data(num_elements):
    data = '<ElementData var="mat_axis" elem_set="Part1">\n'
    
    for i in range(1, num_elements + 1):
        theta = 2 * math.pi * (i + 0.5) / num_elements

        a_x = round(math.cos(theta), 6)
        a_y = round(math.sin(theta), 6)
        d_x = round(-a_y, 6)
        d_y = round(a_x, 6)
        
        data += f'\t<elem lid="{i}">\n'
        data += f'\t\t<a>{a_x}, {a_y}, 0</a>\n'
        data += f'\t\t<d>{d_x}, {d_y}, 0</d>\n'
        data += '\t</elem>\n'
    
    data += '</ElementData>'
    
    return data

# Set the number of elements in the circle
num_elements = 800

# Generate circle data
circle_data = generate_circle_data(num_elements)

# Print the generated XML data
print(circle_data)