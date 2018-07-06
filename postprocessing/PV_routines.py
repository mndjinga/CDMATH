import paraview.simple as pvs
import os
import time
import numpy as np
import vtk.util.numpy_support as vn
# do I need to kill the pipeline?

def Extract_PV_data_over_line_to_txt_file(inputFileName, outputFileName,
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
    
def Extract_PV_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution):
    dirName = os.path.dirname(inputFileName)
    outputFileName = os.path.join(dirName, "tmp." + str(os.getpid()) + str(time.clock()) + ".csv")
    Extract_PV_data_over_line_to_txt_file(inputFileName, outputFileName, point1, point2, resolution)
    var, x1, x2, x3 = np.loadtxt(outputFileName, delimiter=',',
                                skiprows=1, usecols=(0, 3, 4, 5), unpack=True)
    os.remove(outputFileName)
    return x1, x2, x3, var

# TODO : this routine 
# should be rewritten in order to skip the tmp file
# and directly convert the Slide vtk field elt to a numpy array
# 
def Slice_PV_data_to_numpyArray(inputFileName, point, normal, resolution):
    dirName = os.path.dirname(inputFileName)
    outputFileName = os.path.join(dirName, "tmp." + str(os.getpid()) + str(time.clock()) + ".csv")
    Slice_PV_data_to_txt_file(inputFileName, outputFileName, point, normal, resolution)
    var = np.loadtxt(outputFileName, delimiter=',',
                                skiprows=1, usecols=(0), unpack=True)
    os.remove(outputFileName)
    return var

def Slice_PV_data_to_txt_file(inputFileName,
                                        outputFileName,
                                        point,
                                        normal,resolution):
    pvs._DisableFirstRenderCameraReset()
    data_vtu = pvs.XMLUnstructuredGridReader(FileName=[inputFileName])
    Slice1 = pvs.Slice(SliceType="Plane")
    Slice1.SliceOffsetValues = [0.0]
    Slice1.SliceType.Origin = point
    Slice1.SliceType.Normal = normal
    CellCenters1 = pvs.CellCenters()    
    writer = pvs.CreateWriter(outputFileName, CellCenters1)
    writer.Precision=resolution
    writer.FieldAssociation = "Points"  # or "Cells"
    writer.UpdatePipeline()
 
def Slice_PV_field_data_to_numpyArray(field,
                                   point, normal,
                                   resolution
                                   ):
    inputFileName = field.getName()
    field.writeVTK(inputFileName)
 
    numpy_array = Slice_PV_data_to_numpyArray(inputFileName+"_0.vtu", point, normal, resolution)

    os.remove(inputFileName+"_0.vtu")
    return numpy_array

def Slice_PV_field_data_to_txt_file(field, outputFileName,
                                        point, normal,
                                        resolution):
    numpy_array = Slice_PV_field_data_to_numpyArray(field, point, normal, resolution)

    np.savetxt(outputFileName, numpy_array, delimiter=" ")



def Save_PV_data_to_picture_file(inputFileName, field_name,
                             node_or_cell,outputFileName
                             ):
    pvs._DisableFirstRenderCameraReset()

    # create a new 'XML Unstructured Grid Reader'
    reader = pvs.XMLUnstructuredGridReader(FileName=[inputFileName])
    if node_or_cell== 'CELLS':
        reader.CellArrayStatus = [field_name]
    elif node_or_cell== 'NODES':
        reader.PointArrayStatus = [field_name]
    else:
        raise ValueError("unknown type : should be CELLS or NODES")
        
    # get active view
    renderView1 = pvs.GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1057, 499]

    # show data in view
    display = pvs.Show(reader, renderView1);
    # trace defaults for the display properties.
    display.ColorArrayName = [None, '']
    display.GlyphType = 'Arrow'
    display.ScalarOpacityUnitDistance = 0.02234159571242408

    # reset view to fit data
    renderView1.ResetCamera()

    # set scalar coloring
    if node_or_cell== 'CELLS':
        pvs.ColorBy(display, ('CELLS', field_name))
    elif node_or_cell== 'NODES':
        pvs.ColorBy(display, ('POINTS', field_name))
    else:
        raise ValueError("unknown type : should be CELLS or NODES")

    # rescale color and/or opacity maps used to include current data range
    display.RescaleTransferFunctionToDataRange(True)

    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)

    pvs.SaveScreenshot(outputFileName+".png", magnification=1, quality=100, view=renderView1)
    display.SetScalarBarVisibility(renderView1, False)
    pvs.Delete()
    

