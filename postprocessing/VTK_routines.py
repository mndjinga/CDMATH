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

    vtkarray = probe.GetOutput().GetCellData().GetArray(0) # or Slice1.GetCellData() # or Clip1.GetCellData()
    numpy_array = npvtk.vtk_to_numpy(vtkarray)

    return numpy_array
    
def Extract_VTK_data_over_line_to_csv_file(inputFileName, outputFileName, point1, point2, resolution):

    numpy_array = Extract_VTK_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution)

    np.savetxt(outputFileName, numpy_array, delimiter=" ")
   
def Extract_field_data_over_line_to_numpyArray(field, point1, point2, resolution):
    field.writeVTK(field.get_name())
    inputFileName = os.getcwd()+field.get_name()+"_0.vtu"

    numpy_array = Extract_VTK_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution)

    os.remove(inputFileName)
    return numpy_array

def Extract_field_data_over_line_to_csv_file(field, point1, point2, resolution, outputFileName):
    field.writeVTK(field.get_name())
    inputFileName = os.getcwd()+field.get_name()+"_0.vtu"

    numpy_array = Extract_VTK_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution)

    os.remove(inputFileName)
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

    vtkarray = cutter.GetOutput().GetCellData().GetArray(0)
    numpy_array = npvtk.vtk_to_numpy(vtkarray)
    
    return numpy_array

    
def Slice_VTK_data_to_csv_file(inputFileName, outputFileName,
                                           point, normal,
                                           resolution
                                           ):
    numpy_array =   Slice_VTK_data_to_numpyArray(inputFileName, point, normal, resolution )  
    
    np.savetxt(outputFileName, numpy_array, delimiter=" ")
    
     
def Slice_field_data_over_line_to_numpyArray(field,
                                           point, normal,
                                           resolution
                                           ):
    field.writeVTK(field.get_name())
    inputFileName = os.getcwd()+field.get_name()+"_0.vtu"
 
    numpy_array = Slice_VTK_to_numpyArray(inputFileName, point, normal, resolution)

    os.remove(inputFileName)
    return numpy_array

def Slice_field_data_over_line_to_csv_file(field, outputFileName,
                                        point, normal,
                                        resolution):
    field.writeVTK(field.get_name())
    inputFileName = os.getcwd()+field.get_name()+"_0.vtu"
 
    numpy_array = Slice_VTK_to_numpyArray(inputFileName, point, normal, resolution)

    os.remove(inputFileName)
    np.savetxt(outputFileName, numpy_array, delimiter=" ")

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
    #Update or Write ?
    writer.Write()

def Save_VTK_data_to_picture_file(inputFileName,
                             outputFileName
                             ):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    writer = vtk.vtkPNGWriter()
    writer.SetInputConnection(cdata_vtu.GetOutput())
    writer.SetFileName(outputFileName+".png")
    writer.Write()
