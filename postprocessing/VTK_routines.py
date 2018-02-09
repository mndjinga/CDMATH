import os
import time
import numpy as np
import vtk
from vtk.util import numpy_support as npvtk 
# do I need to kill the pipeline?

def Extract_VTK_data_over_line_to_csv_file(inputFileName, outputFileName,
                                           point, normal,
                                           resolution
                                           ):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    plane = vtk.vtkPlane()
    plane.SetOrigin(point)
    plane.SetNormal(normal)

    cutter = vtk.vtkFiltersCorePython.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(reader.GetOutputPort())
    cutter.Update()

    SaveData('outputFileName', proxy=cutter)
    
def Extract_VTK_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution):

    data_vtu = vtk.vtkUnstructuredGridReader(FileName=[inputFileName])

    PlotOverLine1 = pvs.PlotOverLine(data_vtu, Source="High Resolution Line Source" )
    PlotOverLine1.Source.Point1 = point1
    PlotOverLine1.Source.Point2 = point2
    PlotOverLine1.Source.Resolution = resolution

    vtkarray = PlotOverLine1.GetCellData() # or Slice1.GetCellData() # or Clip1.GetCellData()
    numpy_array = npvtk.vtk_to_numpy(vtkarray)

    return numpy_array
    
def Extract_field_data_over_line_to_numpyArray(field, point1, point2, resolution):
    field.writeVTK(field.get_name())
    inputFileName = os.getcwd()+field.get_name()+".vtu"

    result = Extract_VTK_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution)

    os.remove(inputFileName)
    return result

def Slice_VTK_to_numpyArray(inputFileName, point, normal):

    data_vtu = vtk.vtkUnstructuredGridReader(FileName=[inputFileName])
    Slice1 = pvs.Slice(SliceType="Plane")
    Slice1.SliceOffsetValues = [0.0]
    Slice1.SliceType.Origin = point
    Slice1.SliceType.Normal = normal
    CellCenters1 = pvs.CellCenters()    

    vtkarray = Slice1.GetCellData() # or PlotOverLine1.GetCellData() # or Clip1.GetCellData()
    numpy_array = npvtk.vtk_to_numpy(vtkarray)

    return numpy_array

def Slice_MED_data_over_line_to_csv_file(inputFileName,
                                        outputFileName,
                                        point,
                                        normal):

    data_vtu = vtk.vtkUnstructuredGridReader(FileName=[inputFileName])
    Slice1 = pvs.Slice(SliceType="Plane")
    Slice1.SliceOffsetValues = [0.0]
    Slice1.SliceType.Origin = point
    Slice1.SliceType.Normal = normal
    CellCenters1 = pvs.CellCenters()    
    writer = pvs.CreateWriter(outputFileName, CellCenters1)
    writer.Precision=30
    writer.FieldAssociation = "Cells"  # or "Points" (default in cdmath is finite volumes)
    writer.UpdatePipeline()
 
def Slice_field_data_over_line_to_numpyArray(field,
                                        outputFileName,
                                        point,
                                        normal):
    field.writeVTK(field.get_name())
    inputFileName = os.getcwd()+field.get_name()+".vtu"
 
    result = Slice_VTK_to_numpyArray(inputFileName, point1, point2, resolution)

    os.remove(inputFileName)
    return result

def Save_VTK_to_picture_file(inputFileName,
                             outputFileName,
                             ):
    data_vtu = vtk.vtkUnstructuredGridReader(FileName=[inputFileName])
    Show(data_vtu)
    Render()
    SaveScreenshot(outputFileName) 
