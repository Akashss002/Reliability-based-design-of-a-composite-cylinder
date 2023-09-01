# from abaqus import *
# from abaqusConstants import *
# import section
# import regionToolset
# import displayGroupMdbToolset as dgm
# import part
# import material
# import assembly
# import step
# import interaction
# import load
# import mesh
# import optimization
# import job
# import sketch
# import visualization
# import xyPlot
# import displayGroupOdbToolset as dgo
# import connectorBehavior

# import odbAccess
# from odbAccess import *
# import pandas as pd
# import glob
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import stats


def model():
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=6.0)
    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
        point2=(0.3, 0.1))
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-1'].BaseShell(sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']
    mdb.models['Model-1'].Material(name='Composite1')
    mdb.models['Model-1'].materials['Composite1'].Density(table=((1600, ), ))
    mdb.models['Model-1'].materials['Composite1'].Elastic(table=((259.7, 6.9, 
        0.253, 5.2, 2.5, 2.5), ), type=LAMINA)
    mdb.models['Model-1'].parts['Part-1'].Surface(name='Surf-1', side1Faces=
        mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#1 ]', 
        ), ))
    mdb.models['Model-1'].parts['Part-1'].Set(edges=
        mdb.models['Model-1'].parts['Part-1'].edges.getSequenceFromMask(('[#f ]', 
        ), ), name='Set-1')
    mdb.models['Model-1'].parts['Part-1'].CompositeLayup(description='', 
        elementType=SHELL, name='CompositeLayup-1', offsetType=MIDDLE_SURFACE, 
        symmetric=False, thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1'].compositeLayups['CompositeLayup-1'].Section(
        integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
        temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
    mdb.models['Model-1'].parts['Part-1'].compositeLayups['CompositeLayup-1'].CompositePly(
        additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
        , axis=AXIS_3, material='Composite1', numIntPoints=3, orientationType=
        SPECIFY_ORIENT, orientationValue=0.0, plyName='Ply-1', region=Region(
        faces=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(
        mask=('[#1 ]', ), )), suppressed=False, thickness=0.001, thicknessType=
        SPECIFY_THICKNESS)
    mdb.models['Model-1'].parts['Part-1'].compositeLayups['CompositeLayup-1'].CompositePly(
        additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
        , axis=AXIS_3, material='Composite1', numIntPoints=3, orientationType=
        SPECIFY_ORIENT, orientationValue=45.0, plyName='Ply-2', region=Region(
        faces=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(
        mask=('[#1 ]', ), )), suppressed=False, thickness=0.001, thicknessType=
        SPECIFY_THICKNESS)
    mdb.models['Model-1'].parts['Part-1'].compositeLayups['CompositeLayup-1'].CompositePly(
        additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
        , axis=AXIS_3, material='Composite1', numIntPoints=3, orientationType=
        SPECIFY_ORIENT, orientationValue=90.0, plyName='Ply-3', region=Region(
        faces=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(
        mask=('[#1 ]', ), )), suppressed=False, thickness=0.001, thicknessType=
        SPECIFY_THICKNESS)
    mdb.models['Model-1'].parts['Part-1'].compositeLayups['CompositeLayup-1'].ReferenceOrientation(
        additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
        , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
        localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
        normalAxisRegion=mdb.models['Model-1'].parts['Part-1'].surfaces['Surf-1'], 
        orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
        AXIS_1, primaryAxisRegion=
        mdb.models['Model-1'].parts['Part-1'].sets['Set-1'], stackDirection=
        STACK_3)
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
        part=mdb.models['Model-1'].parts['Part-1'])
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
    mdb.models['Model-1'].FieldOutputRequest(createStepName='Step-1', 
        layupLocationMethod=SPECIFIED, layupNames=('Part-1-1.CompositeLayup-1', ), 
        name='F-Output-2', outputAtPlyBottom=False, outputAtPlyMid=True, 
        outputAtPlyTop=False, rebar=EXCLUDE, variables=('S', 'MISES'))
    mdb.models['Model-1'].rootAssembly.Set(edges=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.getSequenceFromMask(
        ('[#8 ]', ), ), name='Set-1')
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-1', region=mdb.models['Model-1'].rootAssembly.sets['Set-1'], u1=0.0, 
        u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)
    mdb.models['Model-1'].rootAssembly.Surface(name='Surf-1', side1Edges=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.getSequenceFromMask(
        ('[#2 ]', ), ))
    mdb.models['Model-1'].ShellEdgeLoad(createStepName='Step-1', distributionType=
        UNIFORM, field='', localCsys=None, magnitude=-100.0, name='Load-1', region=
        mdb.models['Model-1'].rootAssembly.surfaces['Surf-1'], traction=TRANSVERSE)
    mdb.models['Model-1'].parts['Part-1'].seedPart(deviationFactor=0.1, 
        minSizeFactor=0.1, size=0.005)
    mdb.models['Model-1'].parts['Part-1'].generateMesh()
    mdb.models['Model-1'].rootAssembly.regenerate()

#ply_buckle
def job_creationE1(job_name, E1,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['IM7/ MTM45 facesheet'].Elastic(table=((E1[i], 8.62e+9, 
            0.36, 5.31e+9, 5.31e+9, 2.65e+9), ), type=LAMINA)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationCE1(job_name, CE1,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['Core'].Elastic(table=((CE1[i], 262000.0, 
        413690000.0, 0.45, 0.0001, 0.0001, 120000.0, 203400000.0, 82740000.0), ), 
        type=ENGINEERING_CONSTANTS)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationE2(job_name, E2,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['IM7/ MTM45 facesheet'].Elastic(table=((137.89e+9, E2[i], 
            0.36, 5.31e+9, 5.31e+9, 2.65e+9), ), type=LAMINA)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationCE2(job_name, CE2,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['Core'].Elastic(table=((344000, CE2[i], 
        413690000.0, 0.45, 0.0001, 0.0001, 120000.0, 203400000.0, 82740000.0), ), 
        type=ENGINEERING_CONSTANTS)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationCE3(job_name, CE3,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['Core'].Elastic(table=((344000, 262000, 
        CE3[i], 0.45, 0.0001, 0.0001, 120000.0, 203400000.0, 82740000.0), ), 
        type=ENGINEERING_CONSTANTS)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationCv(job_name, Cv,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['Core'].Elastic(table=((344000,262000.0, 
        413690000.0, Cv[i], 0.0001, 0.0001, 120000.0, 203400000.0, 82740000.0), ), 
        type=ENGINEERING_CONSTANTS)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationCG1(job_name, CG1,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['Core'].Elastic(table=((344000,262000.0, 
        413690000.0, 0.45, 0.0001, 0.0001, CG1[i], 203400000.0, 82740000.0), ), 
        type=ENGINEERING_CONSTANTS)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationCG2(job_name, CG2,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['Core'].Elastic(table=((344000,262000.0, 
        413690000.0, 0.45, 0.0001, 0.0001, 120000.0, CG2[i], 82740000.0), ), 
        type=ENGINEERING_CONSTANTS)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationCG3(job_name, CG3,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['Core'].Elastic(table=((344000,262000.0, 
        413690000.0, Cv[i], 0.0001, 0.0001, 120000.0, 203400000.0, CG3[i]), ), 
        type=ENGINEERING_CONSTANTS)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationv(job_name, v,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['IM7/ MTM45 facesheet'].Elastic(table=((137.89e+9, 8.62e+9, 
            v[i], 5.31e+9, 5.31e+9, 2.65e+9), ), type=LAMINA)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationG1(job_name, G1,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['IM7/ MTM45 facesheet'].Elastic(table=((137.89e+9, 8.62e+9, 
            0.36, G1[i], 5.31e+9, 2.65e+9), ), type=LAMINA)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationG2(job_name, G2,model):
    for i in range(0,10):
        mdb.models[str(model)].materials['IM7/ MTM45 facesheet'].Elastic(table=((137.89e+9, 8.62e+9, 
            0.36, 5.31e+9, G2[i], 2.65e+9), ), type=LAMINA)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(model), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def job_creationG3(job_name, G3,m):
    for i in range(0,10):
        mdb.models[str(m)].materials['IM7/ MTM45 facesheet'].Elastic(table=((137.89e+9, 8.62e+9, 
            0.36, 5.31e+9, 5.31e+9, G3[i]), ), type=LAMINA)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(m), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def post_processing(variable_name, limit):
    S11_margin = []
    S12_margin = []
    S22_margin = []
    #S11
    for i in range(1,11):
        o3 = session.openOdb(name='C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/job'+str(variable_name)+str(i)+'.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
        session.viewports['Viewport: 1'].makeCurrent()
        odb = session.odbs['C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/job'+str(variable_name)+str(i)+'.odb']
        S11_max = []
        for j in range(1,20):
            if j==10:
                session.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
                variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S11'), ), ('CORE', 
                MIDDLE)), ), elementSets=("PART-1-1.COMPOSITELAYUP-1-1", ))
                xy = []
                for k in range(1,2521):
                    x = session.xyDataObjects['S:S11 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(k)+' IP: 1']
                    xy.append(x)
                S11 = np.maximum(abs(np.max(xy)),abs(np.min(xy)))
                S11_max.append(S11)
                for l in range(1,2521):
                    del session.xyDataObjects['S:S11 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(l)+' IP: 1']
            else:
                session.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
                    variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S11'), ), ('PLY-'+str(j), 
                    MIDDLE)), ), elementSets=("PART-1-1.COMPOSITELAYUP-1-1", ))
                xy = []
                for k in range(1,2521):
                    x = session.xyDataObjects['S:S11 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(k)+' IP: 1']
                    xy.append(x)
                S11 = np.maximum(abs(np.max(xy)),abs(np.min(xy)))
                S11_max.append(S11)
                for l in range(1,2521):
                    del session.xyDataObjects['S:S11 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(l)+' IP: 1']
        margin = (limit/np.max(S11_max)) - 1
        del(S11_max)
        S11_margin.append(margin)
    csv_file_path = "C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/output/S"+str(variable_name)+"S11.csv"
    with open(csv_file_path, mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(S11_margin)
    del(S11_margin)
  
    #S22    
    for i in range(1,11):
        o3 = session.openOdb(name='C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/job'+str(variable_name)+str(i)+'.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
        session.viewports['Viewport: 1'].makeCurrent()
        odb = session.odbs['C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/job'+str(variable_name)+str(i)+'.odb']
        S22_max = []
        for j in range(1,20):
            if j==10:
                session.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
                variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S22'), ), ('CORE', 
                MIDDLE)), ), elementSets=("PART-1-1.COMPOSITELAYUP-1-1", ))
                xy = []
                for k in range(1,2521):
                    x = session.xyDataObjects['S:S22 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(k)+' IP: 1']
                    xy.append(x)
                S22 = np.maximum(abs(np.max(xy)),abs(np.min(xy)))
                S22_max.append(S22)
                for l in range(1,2521):
                    del session.xyDataObjects['S:S22 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(l)+' IP: 1']
            else:
                session.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
                    variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S22'), ), ('PLY-'+str(j), 
                    MIDDLE)), ), elementSets=("PART-1-1.COMPOSITELAYUP-1-1", ))
                xy = []
                for k in range(1,2521):
                    x = session.xyDataObjects['S:S22 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(k)+' IP: 1']
                    xy.append(x)
                S22 = np.maximum(abs(np.max(xy)),abs(np.min(xy)))
                S22_max.append(S22)
                for l in range(1,2521):
                    del session.xyDataObjects['S:S22 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(l)+' IP: 1']
        margin = (limit/np.max(S22_max)) - 1
        del(S22_max)
        S22_margin.append(margin)
    csv_file_path2 = "C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/output/S"+str(variable_name)+"S22.csv"
    with open(csv_file_path2, mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(S22_margin)
    del(S22_margin)
    #S12
    for i in range(1,11):
        o3 = session.openOdb(name='C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/job'+str(variable_name)+str(i)+'.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
        session.viewports['Viewport: 1'].makeCurrent()
        odb = session.odbs['C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/job'+str(variable_name)+str(i)+'.odb']
        S12_max = []
        for j in range(1,20):
            if j==10:
                session.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
                variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S12'), ), ('CORE', 
                MIDDLE)), ), elementSets=("PART-1-1.COMPOSITELAYUP-1-1", ))
                xy = []
                for k in range(1,2521):
                    x = session.xyDataObjects['S:S12 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(k)+' IP: 1']
                    xy.append(x)
                S12 = np.maximum(abs(np.max(xy)),abs(np.min(xy)))
                S12_max.append(S12)
                for l in range(1,2521):
                    del session.xyDataObjects['S:S12 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(l)+' IP: 1']
            else:
                session.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
                    variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S12'), ), ('PLY-'+str(j), 
                    MIDDLE)), ), elementSets=("PART-1-1.COMPOSITELAYUP-1-1", ))
                xy = []
                for k in range(1,2521):
                    x = session.xyDataObjects['S:S12 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(k)+' IP: 1']
                    xy.append(x)
                S12 = np.maximum(abs(np.max(xy)),abs(np.min(xy)))
                S12_max.append(S12)
                for l in range(1,2521):
                    del session.xyDataObjects['S:S12 SP:'+str(3*(j-1)+2)+' PI: PART-1-1 E: '+str(l)+' IP: 1']
        margin = (limit/np.max(S12_max)) - 1
        del(S12_max)
        S12_margin.append(margin)
    
    csv_file_path3 = "C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/output/S"+str(variable_name)+"S12.csv"
    with open(csv_file_path3, mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(S12_margin)
    del(S12_margin)

def normal_fit(file_name, xlabel):
    data = []
    csv_file_path='C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/output/'+str(file_name)+'.csv'
    with open(csv_file_path, mode='r') as file:
        reader = csv.reader(file)
        row = next(reader)
        data = [float(cell) for cell in row]
    mu, sigma = stats.norm.fit(data)

    # Generate the probability density function (PDF) based on the fitted parameters
    x = np.linspace(min(data), max(data), 100)
    pdf = stats.norm.pdf(x, mu, sigma)

    # Generate the cumulative distribution function (CDF) based on the fitted parameters
    cdf = stats.norm.cdf(x, mu, sigma)

    # Plot the histogram of the data and the fitted PDF
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.hist(data, bins='auto', density=True, alpha=0.7, label='Data')
    plt.plot(x, pdf, 'r', linewidth=2, label='Fitted PDF')
    plt.xlabel(xlabel)
    plt.ylabel('Probability')
    plt.title('PDF for Varying '+str(xlabel))
    plt.legend()

    # Plot the CDF of the data
    plt.subplot(1, 2, 2)
    plt.plot(x, cdf, 'b', linewidth=2, label='Fitted CDF')
    plt.xlabel(xlabel)
    plt.ylabel('Cumulative Probability')
    plt.title('CDF for Varying '+str(xlabel))        
    plt.legend()

    plt.tight_layout()
    plt.savefig("C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/Output/"+str(file_name)+"normal"+".pdf", format="pdf", bbox_inches="tight")

def normal_fit_data(k, xlabel):
    mu, sigma = stats.norm.fit(k)

    # Generate the probability density function (PDF) based on the fitted parameters
    x = np.linspace(min(k), max(k), 100)
    pdf = stats.norm.pdf(x, mu, sigma)

    # Generate the cumulative distribution function (CDF) based on the fitted parameters
    cdf = stats.norm.cdf(x, mu, sigma)

    # Plot the histogram of the data and the fitted PDF
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.hist(k, bins='auto', density=True, alpha=0.7, label='Data')
    plt.plot(x, pdf, 'r', linewidth=2, label='Fitted PDF')
    plt.xlabel(xlabel)
    plt.ylabel('Probability')
    plt.title('PDF for Varying '+str(xlabel))
    plt.legend()

    # Plot the CDF of the data
    plt.subplot(1, 2, 2)
    plt.plot(x, cdf, 'b', linewidth=2, label='Fitted CDF')
    plt.xlabel(xlabel)
    plt.ylabel('Cumulative Probability')
    plt.title('CDF for Varying '+str(xlabel))        
    plt.legend()

    plt.tight_layout()
    plt.savefig("C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/Output/"+str(xlabel)+"normal"+".pdf", format="pdf", bbox_inches="tight")

def log_normal_fit(variable_name, variable_array):
    for j in range(1,10):
        # Read data from CSV file
        data = pd.read_csv('C:/Users/akash/Documents/Akash/CMSE Internship/Trial 4/output'+str(variable_name)+'.csv').iloc[j-1, 1:].values
        
        # Remove non-numeric values
        data = pd.to_numeric(data, errors='coerce')
        data = data[np.isfinite(data)]
        # Fit data to a log-normal distribution
        shape, loc, scale = stats.lognorm.fit(data)

        # Generate the probability density function (PDF) based on the fitted parameters
        x = np.linspace(min(data), max(data), 100)
        pdf = stats.lognorm.pdf(x, shape, loc=loc, scale=scale)

        # Generate the cumulative distribution function (CDF) based on the fitted parameters
        cdf = stats.lognorm.cdf(x, shape, loc=loc, scale=scale)

        # Plot the histogram of the data and the fitted PDF
        plt.figure(figsize=(10, 4))
        plt.subplot(1, 3, 1)
        plt.plot(x, pdf, 'r', linewidth=2, label='Fitted PDF')
        if j-1 == 0:
            plt.xlabel('S11 in ply 1')
        if j-1 == 1:
            plt.xlabel('S11 in ply 2')
        if j-1 == 2:
            plt.xlabel('S11 in ply 3')
        if j-1 == 3:
            plt.xlabel('S22 in ply 1')
        if j-1 == 4:
            plt.xlabel('S22 in ply 2')
        if j-1 == 5:
            plt.xlabel('S22 in ply 3')
        if j-1 == 6:
            plt.xlabel('S12 in ply 1')
        if j-1 == 7:
            plt.xlabel('S12 in ply 2')
        if j-1 == 8:
            plt.xlabel('S12 in ply 3')
        plt.ylabel('Probability')
        plt.title('PDF for Varying '+str(variable_name))
        plt.legend()

        plt.subplot(1, 3, 2)
        plt.hist(data, bins='auto', density=True, alpha=0.7, label='Data')
        if j-1 == 0:
            plt.xlabel('S11 in ply 1')
        if j-1 == 1:
            plt.xlabel('S11 in ply 2')
        if j-1 == 2:
            plt.xlabel('S11 in ply 3')
        if j-1 == 3:
            plt.xlabel('S22 in ply 1')
        if j-1 == 4:
            plt.xlabel('S22 in ply 2')
        if j-1 == 5:
            plt.xlabel('S22 in ply 3')
        if j-1 == 6:
            plt.xlabel('S12 in ply 1')
        if j-1 == 7:
            plt.xlabel('S12 in ply 2')
        if j-1 == 8:
            plt.xlabel('S12 in ply 3')
        plt.ylabel('Probability')
        plt.title('Data for Varying '+str(variable_name))
        plt.legend()
        # Plot the CDF of the data
        plt.subplot(1, 3, 3)
        plt.plot(x, cdf, 'b', linewidth=2, label='Fitted CDF')
        if j-1 == 0:
            plt.xlabel('S11 in ply 1')
        if j-1 == 1:
            plt.xlabel('S11 in ply 2')
        if j-1 == 3:
            plt.xlabel('S11 in ply 3')
        if j-1 == 4:
            plt.xlabel('S22 in ply 1')
        if j-1 == 5:
            plt.xlabel('S22 in ply 2')
        if j-1 == 6:
            plt.xlabel('S22 in ply 3')
        if j-1 == 7:
            plt.xlabel('S12 in ply 1')
        if j-1 == 8:
            plt.xlabel('S12 in ply 2')
        if j-1 == 9:
            plt.xlabel('S12 in ply 3')
        plt.ylabel('Cumulative Probability')
        plt.title('CDF for Varying '+str(variable_name))        
        plt.legend()
        plt.tight_layout()
        plt.savefig("C:/Users/akash/Documents/Akash/CMSE Internship/Trial 4/output/"+str(variable_name)+"lognormal"+str(j)+".pdf", format="pdf", bbox_inches="tight")

def gamma_fit(variable_name, variable_array):
    for j in range(1,10):
        # Read data from CSV file
        data = pd.read_csv('C:/Users/akash/Documents/Akash/CMSE Internship/Trial 4/output'+str(variable_name)+'.csv').iloc[j-1, 1:].values
        
        # Remove non-numeric values
        data = pd.to_numeric(data, errors='coerce')
        data = data[np.isfinite(data)]
        # Fit data to a gamma distribution
        shape, loc, scale = stats.gamma.fit(data)

        # Generate the probability density function (PDF) based on the fitted parameters
        x = np.linspace(min(data), max(data), 100)
        pdf = stats.gamma.pdf(x, shape, loc=loc, scale=scale)

        # Generate the cumulative distribution function (CDF) based on the fitted parameters
        cdf = stats.gamma.cdf(x, shape, loc=loc, scale=scale)

        # Plot the histogram of the data and the fitted PDF
        plt.figure(figsize=(10, 4))
        plt.subplot(1, 2, 1)
        plt.hist(data, bins='auto', density=True, alpha=0.7, label='Data')
        plt.plot(x, pdf, 'r', linewidth=2, label='Fitted PDF')
        if j-1 == 0:
            plt.xlabel('S11 in ply 1')
        if j-1 == 1:
            plt.xlabel('S11 in ply 2')
        if j-1 == 2:
            plt.xlabel('S11 in ply 3')
        if j-1 == 3:
            plt.xlabel('S22 in ply 1')
        if j-1 == 4:
            plt.xlabel('S22 in ply 2')
        if j-1 == 5:
            plt.xlabel('S22 in ply 3')
        if j-1 == 6:
            plt.xlabel('S12 in ply 1')
        if j-1 == 7:
            plt.xlabel('S12 in ply 2')
        if j-1 == 8:
            plt.xlabel('S12 in ply 3')
        plt.ylabel('Probability')
        plt.title('PDF for Varying '+str(variable_name))
        plt.legend()

        # Plot the CDF of the data
        plt.subplot(1, 2, 2)
        plt.plot(x, cdf, 'b', linewidth=2, label='Fitted CDF')
        if j-1 == 0:
            plt.xlabel('S11 in ply 1')
        if j-1 == 1:
            plt.xlabel('S11 in ply 2')
        if j-1 == 3:
            plt.xlabel('S11 in ply 3')
        if j-1 == 4:
            plt.xlabel('S22 in ply 1')
        if j-1 == 5:
            plt.xlabel('S22 in ply 2')
        if j-1 == 6:
            plt.xlabel('S22 in ply 3')
        if j-1 == 7:
            plt.xlabel('S12 in ply 1')
        if j-1 == 8:
            plt.xlabel('S12 in ply 2')
        if j-1 == 9:
            plt.xlabel('S12 in ply 3')
        plt.ylabel('Cumulative Probability')
        plt.title('CDF for Varying '+str(variable_name))
        plt.legend()

        plt.tight_layout()
        plt.savefig("C:/Users/akash/Documents/Akash/CMSE Internship/Trial 4/output/"+str(variable_name)+"gamma"+str(j)+".pdf", format="pdf", bbox_inches="tight")

def weibull_fit(variable_name, variable_array):
    for j in range(1,10):
        # Read data from CSV file
        data = pd.read_csv('C:/Users/akash/Documents/Akash/CMSE Internship/Trial 4/output'+str(variable_name)+'.csv').iloc[j-1, 1:].values
        
        # Remove non-numeric values
        data = pd.to_numeric(data, errors='coerce')
        data = data[np.isfinite(data)]
        # Fit data to a Weibull distribution
        shape, loc, scale = stats.weibull_min.fit(data)

        # Generate the probability density function (PDF) based on the fitted parameters
        x = np.linspace(min(data), max(data), 100)
        pdf = stats.weibull_min.pdf(x, shape, loc=loc, scale=scale)

        # Generate the cumulative distribution function (CDF) based on the fitted parameters
        cdf = stats.weibull_min.cdf(x, shape, loc=loc, scale=scale)

        # Plot the histogram of the data and the fitted PDF
        plt.figure(figsize=(10, 4))
        plt.subplot(1, 2, 1)
        plt.hist(data, bins='auto', density=True, alpha=0.7, label='Data')
        plt.plot(x, pdf, 'r', linewidth=2, label='Fitted PDF')
        if j-1 == 0:
            plt.xlabel('S11 in ply 1')
        if j-1 == 1:
            plt.xlabel('S11 in ply 2')
        if j-1 == 2:
            plt.xlabel('S11 in ply 3')
        if j-1 == 3:
            plt.xlabel('S22 in ply 1')
        if j-1 == 4:
            plt.xlabel('S22 in ply 2')
        if j-1 == 5:
            plt.xlabel('S22 in ply 3')
        if j-1 == 6:
            plt.xlabel('S12 in ply 1')
        if j-1 == 7:
            plt.xlabel('S12 in ply 2')
        if j-1 == 8:
            plt.xlabel('S12 in ply 3')
        plt.ylabel('Probability')
        plt.title('PDF for Varying '+str(variable_name))
        plt.legend()

        # Plot the CDF of the data
        plt.subplot(1, 2, 2)
        plt.plot(x, cdf, 'b', linewidth=2, label='Fitted CDF')
        if j-1 == 0:
            plt.xlabel('S11 in ply 1')
        if j-1 == 1:
            plt.xlabel('S11 in ply 2')
        if j-1 == 3:
            plt.xlabel('S11 in ply 3')
        if j-1 == 4:
            plt.xlabel('S22 in ply 1')
        if j-1 == 5:
            plt.xlabel('S22 in ply 2')
        if j-1 == 6:
            plt.xlabel('S22 in ply 3')
        if j-1 == 7:
            plt.xlabel('S12 in ply 1')
        if j-1 == 8:
            plt.xlabel('S12 in ply 2')
        if j-1 == 9:
            plt.xlabel('S12 in ply 3')
        plt.ylabel('Cumulative Probability')
        plt.title('CDF for Varying '+str(variable_name)) 
        plt.legend()

        plt.tight_layout()
        plt.savefig("C:/Users/akash/Documents/Akash/CMSE Internship/Trial 4/output/"+str(variable_name)+"weibull"+str(j)+".pdf", format="pdf", bbox_inches="tight")

def uniform_fit(variable_name, variable_array):
    for j in range(1,10):
        # Read data from CSV file
        data = pd.read_csv('C:/Users/akash/Documents/Akash/CMSE Internship/Trial 4/output'+str(variable_name)+'.csv').iloc[j-1, 1:].values
        
        # Remove non-numeric values
        data = pd.to_numeric(data, errors='coerce')
        data = data[np.isfinite(data)]
        # Fit data to a uniform distribution
        loc = min(data)
        scale = max(data) - min(data)

        # Generate the probability density function (PDF) based on the fitted parameters
        x = np.linspace(min(data), max(data), 100)
        pdf = stats.uniform.pdf(x, loc=loc, scale=scale)

        # Generate the cumulative distribution function (CDF) based on the fitted parameters
        cdf = stats.uniform.cdf(x, loc=loc, scale=scale)

        # Plot the histogram of the data and the fitted PDF
        plt.figure(figsize=(10, 4))
        plt.subplot(1, 2, 1)
        plt.hist(data, bins='auto', density=True, alpha=0.7, label='Data')
        plt.plot(x, pdf, 'r', linewidth=2, label='Fitted PDF')
        if j-1 == 0:
            plt.xlabel('S11 in ply 1')
        if j-1 == 1:
            plt.xlabel('S11 in ply 2')
        if j-1 == 2:
            plt.xlabel('S11 in ply 3')
        if j-1 == 3:
            plt.xlabel('S22 in ply 1')
        if j-1 == 4:
            plt.xlabel('S22 in ply 2')
        if j-1 == 5:
            plt.xlabel('S22 in ply 3')
        if j-1 == 6:
            plt.xlabel('S12 in ply 1')
        if j-1 == 7:
            plt.xlabel('S12 in ply 2')
        if j-1 == 8:
            plt.xlabel('S12 in ply 3')
        plt.ylabel('Probability')
        plt.title('PDF for Varying '+str(variable_name))
        plt.legend()

        # Plot the CDF of the data
        plt.subplot(1, 2, 2)
        plt.plot(x, cdf, 'b', linewidth=2, label='Fitted CDF')
        if j-1 == 0:
            plt.xlabel('S11 in ply 1')
        if j-1 == 1:
            plt.xlabel('S11 in ply 2')
        if j-1 == 3:
            plt.xlabel('S11 in ply 3')
        if j-1 == 4:
            plt.xlabel('S22 in ply 1')
        if j-1 == 5:
            plt.xlabel('S22 in ply 2')
        if j-1 == 6:
            plt.xlabel('S22 in ply 3')
        if j-1 == 7:
            plt.xlabel('S12 in ply 1')
        if j-1 == 8:
            plt.xlabel('S12 in ply 2')
        if j-1 == 9:
            plt.xlabel('S12 in ply 3')
        plt.ylabel('Cumulative Probability')
        plt.title('CDF for Varying '+str(variable_name)) 
        plt.legend()

        plt.tight_layout()
        plt.savefig("C:/Users/akash/Documents/Akash/CMSE Internship/Trial 4/output/"+str(variable_name)+"uniform"+str(j)+".pdf", format="pdf", bbox_inches="tight")

def fit_dist(variable_name, variable_array):
    normal_fit(variable_name, variable_array)
    log_normal_fit(variable_name, variable_array)
    gamma_fit(variable_name, variable_array)
    weibull_fit(variable_name, variable_array)
    uniform_fit(variable_name, variable_array)

def post_processing_buckling(variable_name):
    CBF = []
    for i in range(1,11):
        odb_path = 'C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/jobB'+str(variable_name)+str(i)+'.odb'
        a1 = odbAccess.openOdb(odb_path)
        a1.steps['Step-1'].frames[1].mode

        eigen_mode = a1.steps['Step-1'].frames[1].description
        eiger_value_1 = float(eigen_mode.split('=')[-1].strip())
        e1 = eiger_value_1 * 0.65
        CBF.append(e1)
    csv_file_path = "C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/output/B"+str(variable_name)+".csv"
    with open(csv_file_path, mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(CBF)

def plot_var(file_name, xlabel, ylabel, variable_array):
    data = []
    csv_file_path='C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/output/'+str(file_name)+'.csv'
    with open(csv_file_path, mode='r') as file:
        reader = csv.reader(file)
        row = next(reader)
        data = [float(cell) for cell in row]
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    plt.plot(variable_array,data,marker ='.')
    plt.savefig("C:/Users/akash/Documents/Akash/CMSE Internship/Cylinder2/output/"+str(file_name)+".pdf", format="pdf", bbox_inches="tight")
    plt.clf()

def jobcreationt_fs(job_name,t,model):
    for i in range(0,10):
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region1=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region2=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region3=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region4=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region5=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region6=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region7=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region8=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region9=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region10=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region11=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region12=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region13=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region14=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region15=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region16=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region17=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region18=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region19=regionToolset.Region(faces=faces)
        compositeLayup = mdb.models[str(model)].parts['Part-1'].compositeLayups['CompositeLayup-1']
        compositeLayup.deletePlies()
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=t[i], orientationType=SPECIFY_ORIENT, 
            orientationValue=45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region2, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=t[i], orientationType=SPECIFY_ORIENT, 
            orientationValue=-45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-3', region=region3, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=t[i], orientationType=SPECIFY_ORIENT, 
            orientationValue=45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-4', region=region4, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=t[i], orientationType=SPECIFY_ORIENT, 
            orientationValue=-45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-5', region=region5, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=t[i], orientationType=SPECIFY_ORIENT, 
            orientationValue=0.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-6', region=region6, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=t[i], orientationType=SPECIFY_ORIENT, 
            orientationValue=0.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-7', region=region7, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=t[i], orientationType=SPECIFY_ORIENT, 
            orientationValue=0.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-8', region=region8, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=t[i], orientationType=SPECIFY_ORIENT, 
            orientationValue=45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-9', region=region9, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=t[i], orientationType=SPECIFY_ORIENT, 
            orientationValue=-45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Core', region=region10, 
            material='Core', thicknessType=SPECIFY_THICKNESS, thickness=0.00635, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-11', 
            region=region11, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=t[i], 
            orientationType=SPECIFY_ORIENT, orientationValue=-45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-12', 
            region=region12, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=t[i], 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-13', 
            region=region13, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=t[i], 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-14', 
            region=region14, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=t[i], 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-15', 
            region=region15, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=t[i], 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-16', 
            region=region16, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=t[i], 
            orientationType=SPECIFY_ORIENT, orientationValue=-45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-17', 
            region=region17, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=t[i], 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-18', 
            region=region18, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=t[i], 
            orientationType=SPECIFY_ORIENT, orientationValue=-45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-19', 
            region=region19, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=t[i], 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(model), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

def jobCreationt_c(job_name,t,model):
    for i in range(0,10):
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region1=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region2=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region3=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region4=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region5=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region6=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region7=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region8=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region9=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region10=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region11=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region12=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region13=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region14=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region15=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region16=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region17=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region18=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region19=regionToolset.Region(faces=faces)
        compositeLayup = mdb.models[str(model)].parts['Part-1'].compositeLayups['CompositeLayup-1']
        compositeLayup.deletePlies()
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region2, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=-45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-3', region=region3, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-4', region=region4, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=-45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-5', region=region5, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=0.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-6', region=region6, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=0.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-7', region=region7, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=0.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-8', region=region8, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-9', region=region9, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=-45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Core', region=region10, 
            material='Core', thicknessType=SPECIFY_THICKNESS, thickness=t[i], 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-11', 
            region=region11, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=-45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-12', 
            region=region12, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-13', 
            region=region13, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-14', 
            region=region14, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-15', 
            region=region15, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-16', 
            region=region16, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=-45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-17', 
            region=region17, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-18', 
            region=region18, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=-45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-19', 
            region=region19, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(model), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()
            
def jobCreationang(job_name,t,model):
    for i in range(0,10):
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region1=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region2=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region3=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region4=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region5=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region6=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region7=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region8=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region9=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region10=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region11=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region12=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region13=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region14=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region15=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region16=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region17=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region18=regionToolset.Region(faces=faces)
        p = mdb.models[str(model)].parts['Part-1']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region19=regionToolset.Region(faces=faces)
        compositeLayup = mdb.models[str(model)].parts['Part-1'].compositeLayups['CompositeLayup-1']
        compositeLayup.deletePlies()
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region2, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=-45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-3', region=region3, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-4', region=region4, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=-45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-5', region=region5, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=t[i], additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-6', region=region6, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=t[i], additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-7', region=region7, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=t[i], additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-8', region=region8, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-9', region=region9, 
            material='IM7/ MTM45 facesheet', thicknessType=SPECIFY_THICKNESS, 
            thickness=0.0001321, orientationType=SPECIFY_ORIENT, 
            orientationValue=-45.0, additionalRotationType=ROTATION_NONE, 
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Core', region=region10, 
            material='Core', thicknessType=SPECIFY_THICKNESS, thickness=0.00635, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-11', 
            region=region11, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=-45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-12', 
            region=region12, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-13', 
            region=region13, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=t[i], 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-14', 
            region=region14, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=t[i], 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-15', 
            region=region15, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=t[i], 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-16', 
            region=region16, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=-45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-17', 
            region=region17, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-18', 
            region=region18, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=-45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-19', 
            region=region19, material='IM7/ MTM45 facesheet', 
            thicknessType=SPECIFY_THICKNESS, thickness=0.0001321, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(model), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=str(job_name)+str(i+1), nodalOutputPrecision=SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        mdb.jobs[str(job_name)+str(i+1)].submit(consistencyChecking=OFF)
        mdb.jobs[str(job_name)+str(i+1)].waitForCompletion()

E1 = [129.85e+9,132.55e+9,132.57e+9,134.25e+9,135.92e+9,136.98e+9,137.62e+9,138.57e+9,138.73e+9,144.63e+9]
E2 = [8.33e+9,8.46e+9,8.47e+9,8.54e+9,8.57e+9,8.65e+9,8.66e+9,8.76e+9,8.79e+9,8.8e+9]
v = [0.34,0.351,0.355,0.358,0.361,0.362,0.364,0.366,0.369,0.37]
G1 = [5.25e+9,5.28e+9,5.29e+9,5.295e+9,5.33e+9,5.36e+9,5.39e+9,5.46e+9,5.48e+9,5.54e+9]
G2 = [5.25e+9,5.28e+9,5.29e+9,5.295e+9,5.33e+9,5.36e+9,5.39e+9,5.46e+9,5.48e+9,5.54e+9]
G3 = [2.6e+9,2.63e+9,2.635e+9,2.655e+9,2.66e+9,2.67e+9,2.69e+9,2.7e+9,2.72e+9,3.77e+9]

CE1 = [0.339e+6,0.3394e+6,0.342e+6,0.347e+6,0.348e+6,0.349e+6,0.352e+6,0.3526e+6,0.353e+6,0.355e+6]
CE2 = [0.260e+6,0.2603e+6,0.262e+6,0.263e+6,0.2632e+6,0.266e+6,0.267e+6,0.269e+6,0.272e+6,0.275e+6]
CE3 = [394.13e+6,394.184e+6,403.27e+6,407.05e+6,408.37e+6,412.37e+6,415.78e+6,417.01e+6,417.55e+6,426.41e+6]
Cv = [0.43,0.432,0.441,0.444,0.445,0.45,0.46,0.461,0.47,0.472]
CG1 = [0.113e+6,0.117e+6,0.1171e+6,0.118e+6,0.1181e+6,0.1183e+6,0.119e+6,0.121e+6,0.1211e+6,0.122e+6]
CG2 = [193.96e+6,200.9e+6,203.83e+6,204.57e+6,204.62e+6,205.21e+6,205.28e+6,207.43e+6,207.49e+6,208.23e+6]
CG3 = [79.62e+6,81.16e+6,82.23e+6,82.34e+6,82.62e+6,83.2e+6,84.05e+6,84.2e+6,85.08e+6,85.081e+6]

t_fs = [0.000129,0.0001317,0.0001320,0.0001321,0.00013215,0.00013217,0.0001322,0.0001325,0.000134, 0.0001382]
t_core = [0.00618,0.00619,0.0062,0.00622,0.00624,0.00632,0.00635,0.00641,0.00643,0.00647]

theta = [1,2,3,4,5,6,7,8,9,10]

# model()

#buckling results:

# job_creationE1('jobBfsE1', E1)
# post_processing_buckling('fsE1')
# plot_var('BfsE1','E11 of face sheet','Critical Buckling factor',E1)

# job_creationE2('jobBfsE2', E2)
# post_processing_buckling('fsE2')
#plot_var('BfsE2','E22 of face sheet','Critical Buckling factor',E2)

# job_creationv('jobBfsv', v)
# post_processing_buckling('fsv')
# plot_var('Bfsv','v of face sheet','Critical Buckling factor',v)

# job_creationG1('jobBfsG1', G1)
# post_processing_buckling('fsG1')
# plot_var('BfsG1','G12 of face sheet','Critical Buckling factor',G1)

# job_creationG2('jobBfsG2', G2,'Model-1')
# post_processing_buckling('fsG2')
# plot_var('BfsG2','G13 of face sheet','Critical Buckling factor',G2)

# job_creationG3('jobBfsG3', G3)
# post_processing_buckling('fsG3')
# plot_var('BfsG3','G23 of face sheet','Critical Buckling factor',G3)

# jobcreationt_fs('jobBtfs',t_fs,'Model-1')
# post_processing_buckling('tfs')
# plot_var('Btfs','Thickness of face sheet','Critical Buckling factor',t_fs)
# normal_fit('Btfs','Critical Buckling factor')

# jobCreationang('jobBang',theta,'Model-1')
# post_processing_buckling('ang')
# plot_var('Bang','Fiber orientation angle','Critical Buckling factor',theta)
# normal_fit('Bang','Critical Buckling factor')

# job_creationCE1('jobBcE1',CE1,'Model-1')
# post_processing_buckling('cE1')
# plot_var('BcE1','E1 of core','Critical Buckling factor',CE1)

# job_creationCE2('jobBcE2',CE2)
# post_processing_buckling('cE2')
# plot_var('BcE2','E2 of core','Critical Buckling factor',CE2)

# job_creationCE3('jobBcE3',CE3)
# post_processing_buckling('cE3')
# plot_var('BcE3','E3 of core','Critical Buckling factor',CE3)

# job_creationCv('jobBcv',Cv)
# post_processing_buckling('cv')
# plot_var('Bcv','v of core','Critical Buckling factor',Cv)

# job_creationCG1('jobBcG1',CG1)
# post_processing_buckling('cG1')
# plot_var('BcG1','G1 of core','Critical Buckling factor',CG1)

# job_creationCG2('jobBcG2',CG2)
# post_processing_buckling('cG2')
# plot_var('BcG2','G2 of core','Critical Buckling factor',CG2)

# job_creationCG3('jobBcG3',CG3)
# post_processing_buckling('cG3')
# plot_var('BcG3','G3 of core','Critical Buckling factor',CG3)

# jobCreationt_c('jobBtc',t_core,'Model-1')
# post_processing_buckling('tc')
# plot_var('Btc','Thickness of core', 'Critical Buckling factor', t_core)
# normal_fit('Btc','Critical Buckling factor')

#Margin results:

#job_creationE1('jobSfsE1',E1)
#post_processing('SfsE1',506500000)
# plot_var('SSfsE1S11','E1 of face sheet','Margin in S11',E1)
# plot_var('SSfsE1S22','E1 of face sheet','Margin in S22',E1)
# plot_var('SSfsE1S12','E1 of face sheet','Margin in S12',E1)

#job_creationE2('jobSfsE2',E2)
#post_processing('SfsE2',506500000)
# plot_var('SSfsE2S11','E2 of face sheet','Margin in S11',E2)
# plot_var('SSfsE2S22','E2 of face sheet','Margin in S22',E2)
# plot_var('SSfsE2S12','E2 of face sheet','Margin in S12',E2)

# job_creationv('jobSfsv',v)
# post_processing('Sfsv',506500000)
# plot_var('SSfsvS11','v of face sheet','Margin in S11',v)
# plot_var('SSfsvS22','v of face sheet','Margin in S22',v)
# plot_var('SSfsvS12','v of face sheet','Margin in S12',v)

# job_creationG1('jobSfsG1',G1)
# post_processing('SfsG1',506500000)
# plot_var('SSfsG1S11','G12 of face sheet','Margin in S11',G1)
# plot_var('SSfsG1S22','G12 of face sheet','Margin in S22',G1)
# plot_var('SSfsG1S12','G12 of face sheet','Margin in S12',G1)

# job_creationG2('jobSfsG2',G2)
# post_processing('SfsG2',506500000)
# plot_var('SSfsG2S11','G13 of face sheet','Margin in S11',G2)
# plot_var('SSfsG2S22','G13 of face sheet','Margin in S22',G2)
# plot_var('SSfsG2S12','G13 of face sheet','Margin in S12',G2)

#job_creationG3('jobSfsG3',G3)
#post_processing('SfsG3',506500000)
# plot_var('SSfsG3S11','G23 of face sheet','Margin in S11',G3)
# plot_var('SSfsG3S22','G23 of face sheet','Margin in S22',G3)
# plot_var('SSfsG3S12','G23 of face sheet','Margin in S12',G3)

# job_creationCE1('jobScE1',CE1,'Model-1-Static')
# post_processing('ScE1',506500000)
# plot_var('SScE1S11','E11 of core','Margin in S11',CE1)
# plot_var('SScE1S22','E11 of core','Margin in S22',CE1)
# plot_var('SScE1S12','E11 of core','Margin in S12',CE1)
# normal_fit('SScE1S11','Margin in S11')
# normal_fit('SScE1S22','Margin in S22')
# normal_fit('SScE1S12','Margin in S12')

# job_creationCE2('jobScE2',CE2,'Model-1-Static')
# post_processing('ScE2',506500000)
# plot_var('SScE2S11','E22 of core','Margin in S11',CE2)
# plot_var('SScE2S22','E22 of core','Margin in S22',CE2)
# plot_var('SScE2S12','E22 of core','Margin in S12',CE2)
# normal_fit('SScE2S11','Margin in S11')
# normal_fit('SScE2S22','Margin in S22')
# normal_fit('SScE2S12','Margin in S12')

# job_creationCE3('jobScE3',CE3,'Model-1-Static')
# post_processing('ScE3',506500000)
# plot_var('SScE3S11','E33 of core','Margin in S11',CE3)
# plot_var('SScE3S22','E33 of core','Margin in S22',CE3)
# plot_var('SScE3S12','E33 of core','Margin in S12',CE3)
# normal_fit('SScE3S11','Margin in S11')
# normal_fit('SScE3S22','Margin in S22')
# normal_fit('SScE3S12','Margin in S12')

# job_creationCv('jobScv',Cv,'Model-1-Static')
# post_processing('Scv',506500000)
# plot_var('SScvS11','v of core','Margin in S11',Cv)
# plot_var('SScvS22','v of core','Margin in S22',Cv)
# plot_var('SScvS12','v of core','Margin in S12',Cv)
# normal_fit('SScvS11','Margin in S11')
# normal_fit('SScvS22','Margin in S22')
# normal_fit('SScvS12','Margin in S12')

# job_creationCG1('jobScG1',CG1,'Model-1-Static')
# post_processing('ScG1',506500000)
# plot_var('SScG1S11','G12 of core','Margin in S11',CG1)
# plot_var('SScG1S22','G12 of core','Margin in S22',CG1)
# plot_var('SScG1S12','G12 of core','Margin in S12',CG1)
# normal_fit('SScG1S11','Margin in S11')
# normal_fit('SScG1S22','Margin in S22')
# normal_fit('SScG1S12','Margin in S12')

# job_creationCG2('jobScG2',CG2,'Model-1-Static')
# post_processing('ScG2',506500000)
# plot_var('SScG2S11','G13 of core','Margin in S11',CG2)
# plot_var('SScG2S22','G13 of core','Margin in S22',CG2)
# plot_var('SScG2S12','G13 of core','Margin in S12',CG2)
# normal_fit('SScG2S11','Margin in S11')
# normal_fit('SScG2S22','Margin in S22')
# normal_fit('SScG2S12','Margin in S12')

# job_creationCG3('jobScG3',CG3,'Model-1-Static')
# post_processing('ScG3',506500000)
# plot_var('SScG3S11','G23 of core','Margin in S11',CG3)
# plot_var('SScG3S22','G23 of core','Margin in S22',CG3)
# plot_var('SScG3S12','G23 of core','Margin in S12',CG3)
# normal_fit('SScG3S11','Margin in S11')
# normal_fit('SScG3S22','Margin in S22')
# normal_fit('SScG3S12','Margin in S12')

# jobcreationt_fs('jobStfs',t_fs,'Model-1-Static')
# # post_processing('Stfs',506500000)
# plot_var('SStfsS11','Thickness of face sheet','Margin in S11',t_fs)
# plot_var('SStfsS22','Thickness of face sheet','Margin in S22',t_fs)
# plot_var('SStfsS12','Thickness of face sheet','Margin in S12',t_fs)
# normal_fit('SStfsS11','Margin in S11')
# normal_fit('SStfsS22','Margin in S22')
# normal_fit('SStfsS12','Margin in S12')

# jobCreationt_c('jobStc',t_core,'Model-1-Static')
# post_processing('Stc',506500000)
# plot_var('SStcS11','Thickness of core','Margin in S11',t_core)
# plot_var('SStcS22','Thickness of core','Margin in S22',t_core)
# plot_var('SStcS12','Thickness of core','Margin in S12',t_core)
# normal_fit('SStcS11','Margin in S11')
# normal_fit('SStcS22','Margin in S22')
# normal_fit('SStcS12','Margin in S12')

# jobCreationt_c('jobSang',theta,'Model-1-Static')
# post_processing('Sang',506500000)
# plot_var('SSangS11','Fiber orientation angle','Margin in S11',theta)
# plot_var('SSangS22','Fiber orientation angle','Margin in S22',theta)
# plot_var('SSangS12','Fiber orientation angle','Margin in S12',theta)
# normal_fit('SSangS11','Margin in S11')
# normal_fit('SSangS22','Margin in S22')
# normal_fit('SSangS12','Margin in S12')

# normal_fit_data(t_fs,'Thickness of face sheet')