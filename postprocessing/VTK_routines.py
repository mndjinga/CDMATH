#!/usr/bin/env python
# -*-coding:utf-8 -*-

import os
import numpy as np
import vtk
from vtk.util import numpy_support as npvtk 
# do I need to kill the pipeline?

def Extract_VTK_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution):

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    probeLine = vtk.vtkLineSource()
    probeLine.SetPoint1(point1)
    probeLine.SetPoint2(point2)
    probeLine.SetResolution(resolution)
    
    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(probeLine.GetOutputPort())
    probe.SetSourceData(reader.GetOutput())
    probe.Update()

    vtkarray = probe.GetOutput().GetPointData().GetArray(0) # or Slice1.GetCellData() # or Clip1.GetCellData()
    numpy_array = npvtk.vtk_to_numpy(vtkarray)

    return numpy_array
    
def Extract_VTK_data_over_line_to_txt_file(inputFileName, outputFileName, point1, point2, resolution):

    numpy_array = Extract_VTK_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution)

    np.savetxt(outputFileName, numpy_array, delimiter=" ")
   
def Extract_field_data_over_line_to_numpyArray(field, point1, point2, resolution):

    inputFileName = field.getName()#os.getcwd()+field.get_name()
    field.writeVTK(inputFileName)

    numpy_array = Extract_VTK_data_over_line_to_numpyArray(inputFileName+"_0.vtu", point1, point2, resolution)

    os.remove(inputFileName+"_0.vtu")
    return numpy_array

def Extract_field_data_over_line_to_txt_file(field, point1, point2, resolution, outputFileName):

    numpy_array = Extract_field_data_over_line_to_numpyArray(field, point1, point2, resolution)

    np.savetxt(outputFileName, numpy_array, delimiter=" ")

def Slice_VTK_data_to_numpyArray(inputFileName,
                                 point, normal,
                                 resolution
                                           ):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    plane = vtk.vtkPlane()
    plane.SetOrigin(point)
    plane.SetNormal(normal)

    cutter = vtk.vtkFiltersCorePython.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(reader.GetOutputPort())
    cutter.Update()

    vtkarray = cutter.GetOutput().GetPointData().GetArray(0)
    numpy_array = npvtk.vtk_to_numpy(vtkarray)
    
    return numpy_array

    
def Slice_VTK_data_to_txt_file(inputFileName, outputFileName,
                                           point, normal,
                                           resolution
                                           ):
    numpy_array =   Slice_VTK_data_to_numpyArray(inputFileName, point, normal, resolution )  
    
    np.savetxt(outputFileName, numpy_array, delimiter=" ")
    
     
def Slice_field_data_to_numpyArray(field,
                                   point, normal,
                                   resolution
                                   ):
    inputFileName = field.getName()
    field.writeVTK(inputFileName)
 
    numpy_array = Slice_VTK_data_to_numpyArray(inputFileName+"_0.vtu", point, normal, resolution)

    os.remove(inputFileName+"_0.vtu")
    return numpy_array

def Slice_field_data_to_txt_file(field, outputFileName,
                                        point, normal,
                                        resolution):
    numpy_array = Slice_field_data_to_numpyArray(field, point, normal, resolution)

    np.savetxt(outputFileName, numpy_array, delimiter=" ")

def Slice_VTK_data_to_VTK(inputFileName,
                             outputFileName,
                                 point, normal,
                                 resolution
                                           ):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    plane = vtk.vtkPlane()
    plane.SetOrigin(point)
    plane.SetNormal(normal)

    cutter = vtk.vtkFiltersCorePython.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(reader.GetOutputPort())
    cutter.Update()

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(cutter.GetOutput())
    writer.SetFileName(outputFileName)
    writer.Write()

def Clip_VTK_data_to_VTK(inputFileName,
                             outputFileName,
                                 point, normal,
                                 resolution
                                           ):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    plane = vtk.vtkPlane()
    plane.SetOrigin(point)
    plane.SetNormal(normal)

    clipper = vtk.vtkClipDataSet()
    clipper.SetClipFunction(plane)
    clipper.SetInputConnection(reader.GetOutputPort())
    clipper.Update()

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(clipper.GetOutput())
    writer.SetFileName(outputFileName)
    writer.Write()

def Save_VTK_data_to_picture_file(inputFileName,
                             outputFileName
                             ):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()

#-------------------------------------------------------------------------------
    # create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    
    # create source
    source = vtk.vtkSphereSource()
    source.SetCenter(0,0,0)
    source.SetRadius(5.0)
 
    # mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(source.GetOutput())
 
    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
 
    # color the actor
    actor.GetProperty().SetColor(1,0,0) # (R,G,B)
 
    # assign actor to the renderer
    ren.AddActor(actor)
    
    renWin.Render()
    
    wif = vtk.vtkWindowToImageFilter()
    wif.SetInput(renWin)
    wif.Update()

    writer = vtk.vtkPNGWriter()
    #writer.SetInputConnection(reader.GetOutputPort())
    #writer.SetInputData(reader.GetOutput().GetPointData())
    writer.SetInputConnection(wif.GetOutputPort())
    writer.SetFileName(outputFileName+"bis.png")
    writer.Write()

#-------------------------------------------------------------------------------------------------------
    renwin = vtk.vtkRenderWindow()
    renderer = vtk.vtkRenderer() 
    renwin.AddRenderer(renderer) 
    renwin.Render() 
    #mapper = vtk.vtkUnstructuredGridVolumeRayCastMapper() 
    mapper = vtk.vtkDataSetMapper() 
    mapper.SetInputConnection(reader.GetOutputPort())
    #mapper.SetInput(id)
    
    actor = vtk.vtkActor() 
    actor.SetMapper(mapper) 
    renderer.AddViewProp(actor) 
    renwin.Render() 

    #mapper2D = vtk.vtkImageMapper()
    #actor2D = vtk.vtkActor2D()
    #actor2D.SetMapper(mapper2D)
    #renwin.AddActor2D(actor2D)
    
    wif = vtk.vtkWindowToImageFilter()
    wif.SetInput(renwin)
    wif.Update()
    

#--------------------------------------------------------------------------------
    
    
    writer = vtk.vtkPNGWriter()
    #writer.SetInputConnection(reader.GetOutputPort())
    #writer.SetInputData(reader.GetOutput().GetPointData())
    writer.SetInputConnection(wif.GetOutputPort())
    writer.SetFileName(outputFileName+".png")
    writer.Write()
