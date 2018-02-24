import paraview.simple as pvs
import os
import time
import numpy as np
import vtk.util.numpy_support as vn
# do I need to kill the pipeline?

def Extract_VTU_data_over_line_to_csv_file(inputFileName, outputFileName,
                                           point1, point2,
                                           resolution
                                           ):
    pvs._DisableFirstRenderCameraReset()
    data_vtu = pvs.XMLUnstructuredGridReader(FileName=[inputFileName])
    PlotOverLine1 = pvs.PlotOverLine(Source="High Resolution Line Source"
                                     )
    PlotOverLine1.Source.Point1 = point1
    PlotOverLine1.Source.Point2 = point2
    PlotOverLine1.Source.Resolution = resolution
    writer = pvs.CreateWriter(outputFileName, PlotOverLine1)
    writer.FieldAssociation = "Points"  # or "Cells"
    writer.UpdatePipeline()
    
def Extract_VTU_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution):
    dirName = os.path.dirname(inputFileName)
    outputFileName = os.path.join(dirName, "tmp." + str(os.getpid()) + str(time.clock()) + ".csv")
    Extract_VTU_data_over_line_to_csv_file(inputFileName, outputFileName, point1, point2, resolution)
    var, x1, x2, x3 = np.loadtxt(outputFileName, delimiter=',',
                                skiprows=1, usecols=(0, 3, 4, 5), unpack=True)
    os.remove(outputFileName)
    return x1, x2, x3, var

# TODO : this routine 
# should be rewritten in ordre to skip the tmp file
# and directly convert the Slide vtk field elt to a numpy array
# 
def Slice_VTU_data_over_line_to_numpyArray(inputFileName, point, normal):
    dirName = os.path.dirname(inputFileName)
    outputFileName = os.path.join(dirName, "tmp." + str(os.getpid()) + str(time.clock()) + ".csv")
    Slice_VTU_data_over_line_to_csv_file(inputFileName, outputFileName, point, normal)
    var, x1, x2, x3 = np.loadtxt(outputFileName, delimiter=',',
                                skiprows=1, usecols=(0, 1, 2, 3), unpack=True)
    os.remove(outputFileName)
    return x1, x2, x3, var

def Slice_VTU_data_over_line_to_csv_file(inputFileName,
                                        outputFileName,
                                        point,
                                        normal):
    pvs._DisableFirstRenderCameraReset()
    data_vtu = pvs.XMLUnstructuredGridReader(FileName=[inputFileName])
    Slice1 = pvs.Slice(SliceType="Plane")
    Slice1.SliceOffsetValues = [0.0]
    Slice1.SliceType.Origin = point
    Slice1.SliceType.Normal = normal
    CellCenters1 = pvs.CellCenters()    
    writer = pvs.CreateWriter(outputFileName, CellCenters1)
    writer.Precision=30
    writer.FieldAssociation = "Points"  # or "Cells"
    writer.UpdatePipeline()
 

