import os

from abaqus import *
from abaqusConstants import *
import random
import mesh
import regionToolset
import numpy as np
from odbAccess import *


def AssignPropertiesRun(ABD_list,Lx,Ly,element_size):
    #for modelName in mdb.models.keys():
    #    del mdb.models[modelName]
    model_name = 'Model-1'
    model = mdb.Model(name=model_name)

    # Create a 2D part (plate)
    plate_width = Lx
    plate_height = Ly
    sketch = model.ConstrainedSketch(name='PlateSketch', sheetSize=2*Lx)
    sketch.rectangle(point1=(0.0, 0.0), point2=(plate_width, plate_height))
    plate_part = model.Part(name='Part-1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    plate_part.BaseShell(sketch=sketch)

    # Create a random material
    random_material = model.Material(name='RandomMaterial')
    elastic_modulus = random.uniform(100e3, 210e3)  # MPa
    poisson_ratio = random.uniform(0.2, 0.35)
    random_material.Elastic(table=((elastic_modulus, poisson_ratio),))

    # Mesh the part
    plate_part.seedPart(size=element_size, deviationFactor=0.1, minSizeFactor=0.1)
    plate_part.generateMesh()
    sectionBaseName = 'ShellSection'

    
    for i, element in enumerate(plate_part.elements):
        #nodes = element.connectivity # Get the nodes of the element
        #coords = [plate_part.nodes[node_id].coordinates for node_id in nodes]
        #centroid = np.mean(coords, axis=0) # Calculate the centroid as the mean of the nodal coordinates # FOR CHECKING IF COORD AGREE WITH STOCHASTIC GRID
        ABD_tuple = tuple(ABD_list[i,:])
        sectionName = sectionBaseName + str(i)
        if sectionName not in model.sections.keys():
            shell_section = model.GeneralStiffnessSection(name=sectionName, referenceTemperature=None, stiffnessMatrix=ABD_tuple,applyThermalStress=0, poissonDefinition=DEFAULT, useDensity=OFF)
            # Assign the section to the element
            region = regionToolset.Region(elements=plate_part.elements[i:i+1])
            plate_part.SectionAssignment(region=region, sectionName=sectionName)

    # Create an assembly and instance
    assembly = model.rootAssembly
    instance = assembly.Instance(name='PlateInstance', part=plate_part, dependent=ON)

    # Apply boundary conditions 
    fixed_edge = assembly.instances['PlateInstance'].edges.findAt(((0.0, 25.0, 0.0),))  # left edge
    assembly.Set(edges=fixed_edge, name='FixedEdge') # Create a set for the fixed edge
    model.DisplacementBC(name='FixEdge', createStepName='Initial', region=assembly.sets['FixedEdge'], u1=0.0, u2=0.0, u3 = 0.0, ur1= 0.0, ur2 = 0.0, ur3=0.0) 

    stepName = 'LoadStep'
    model.StaticStep(name=stepName, previous='Initial', timePeriod=1.0, initialInc=0.1, minInc=1e-5, maxInc=0.1)
    
    # Define the region for the distributed load (on a specific face)
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['PlateInstance'].faces
    side1Faces1 = s1.findAt(((33.333333, 33.333333, 0.0), ))
    region = regionToolset.Region(side1Faces=side1Faces1)

    # Create a pressure load
    model.Pressure(name='DistributedLoad', createStepName='LoadStep', region=region, magnitude=0.001)  

    # Job Parameters
    job_name = "my_simulation"  # Change to your job name
    my_job = mdb.Job(name=job_name, model=model_name, description='Simulation Job')
    try:
        my_job.submit(consistencyChecking=OFF)
        print("Job submitted successfully.")
    except Exception as e:
        print("Error submitting job")
    my_job.waitForCompletion()



def ReadSaveOutput(ID,Lx,Ly,Lw):
    odb_path = 'my_simulation.odb'  
    odb = openOdb(path=odb_path)
    #for instance_name in odb.rootAssembly.instances.keys():
    #    print instance_name
    #    #print(f"Instance found: {instance_name}")
    step_name = 'LoadStep'      # Access the first step and its last frame
    step = odb.steps[step_name]
    last_frame = step.frames[-1]  # Access the last frame (final state)

    displacement_field = last_frame.fieldOutputs['U'] # Access the displacement field output (U)
    instance_name = 'PLATEINSTANCE'  # Replace with the correct instance name
    instance = odb.rootAssembly.instances[instance_name]
    node_displacements = []
    for value in displacement_field.values:
        node_label = value.nodeLabel  # Node label
        u3 = value.data[2]  # U3 displacement (Z-direction)
        node_coords = instance.nodes[node_label - 1].coordinates  # Adjust for zero-based indexing # THIS SEEMS CORRECT, BUT NOT IN ABAQUS  INPUT
        node_displacements.append((node_label, node_coords, u3))
    folder_path = 'Output_L' + str(Lx) + 'x' + str(Ly) + '/Cantilever/Lw' + str(Lw)
    if not os.path.isdir(folder_path):
        os.makedirs(folder_path)
    with open('Output_L' + str(Lx) + 'x' + str(Ly) + '/Cantilever/Lw' + str(Lw) + '/node_displacements_with_coords_' + str(ID) + '.txt', 'w') as file:     
        for node, coords, u3 in node_displacements:
            #file.write(f'{node}, {coords[0]}, {coords[1]}, {coords[2]}, {u3}\n')
            file.write('%d, %f, %f, %f, %f\n' % (node, coords[0], coords[1], coords[2], u3))
    odb.close()

def is_pos_def(x): # for 6x6 abd matrix
    return np.all(np.linalg.eigvals(x) > 0)

def GenerateRealization(eigenvalue_vector,eigenvector_matrix,Params,n):    

    m=len(eigenvalue_vector)
    nw = int(m/8)
    nonpos_def_count =1
    num_attempts = 0
    while  nonpos_def_count and  num_attempts<10:
        nonpos_def_count =0
        # Generate multivariate random field realization
        DC = np.diag(eigenvalue_vector) 
        UC= np.random.standard_normal(m)
        XC = np.dot(eigenvector_matrix, np.dot(np.sqrt(DC), UC))
        #
        ABD = np.zeros((nw, n))
        for i in range(1, n + 1):  
            ABD[:, i - 1] = XC[(i - 1) * nw : i * nw]
        ABD_list = np.zeros((nw,21))
        for i in range(0,nw):
            A11= Params[0,0] + Params[0,1]*ABD[i,0]
            A12= Params[1,0] + Params[1,1]*ABD[i,1]
            A22= Params[2,0] + Params[2,1]*ABD[i,2]
            A33= Params[3,0] + Params[3,1]*ABD[i,3]
            D11= Params[4,0] + Params[4,1]*ABD[i,4]
            D12= Params[5,0] + Params[5,1]*ABD[i,5]
            D22= Params[6,0] + Params[6,1]*ABD[i,6]
            D33= Params[7,0] + Params[7,1]*ABD[i,7]
            ABD_list[i,:] = [A11,A12,A22,0,0,A33,0,0,0,D11,0,0,0,D12,D22,0,0,0,0,0,D33]
            #print([A11,A12,A22,0,0,A33,0,0,0,D11,0,0,0,D12,D22,0,0,0,0,0,D33])
            symm_matrix = np.array([
        [A11,         A12,        0,         0,         0,                 0],
        [A12,         A22,        0,         0,         0,                 0],
        [    0,        0,        A33,        0,         0,                 0],
        [    0,        0,         0,        D11,       D12,                0],
        [    0,        0,         0,        D12,       D22,                0],
        [    0,        0,         0,         0,          0,            D33]])
            if is_pos_def(symm_matrix) == False: # check if pos def, otherwise ABD matrix invalid
                print("NOT pos def" + str(i))
                nonpos_def_count = nonpos_def_count +1
        num_attempts = num_attempts +1
    return ABD_list



