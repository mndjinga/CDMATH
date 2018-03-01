import paraview.simple as pvs
import os
import time
import numpy as np
import vtk.util.numpy_support as vn
# do I need to kill the pipeline?

def Extract_PV_data_over_line_to_csv_file(inputFileName, outputFileName,
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
    Extract_PV_data_over_line_to_csv_file(inputFileName, outputFileName, point1, point2, resolution)
    var, x1, x2, x3 = np.loadtxt(outputFileName, delimiter=',',
                                skiprows=1, usecols=(0, 3, 4, 5), unpack=True)
    os.remove(outputFileName)
    return x1, x2, x3, var

# TODO : this routine 
# should be rewritten in ordre to skip the tmp file
# and directly convert the Slide vtk field elt to a numpy array
# 
def Slice_PV_data_over_line_to_numpyArray(inputFileName, point, normal):
    dirName = os.path.dirname(inputFileName)
    outputFileName = os.path.join(dirName, "tmp." + str(os.getpid()) + str(time.clock()) + ".csv")
    Slice_PV_data_over_line_to_csv_file(inputFileName, outputFileName, point, normal)
    var, x1, x2, x3 = np.loadtxt(outputFileName, delimiter=',',
                                skiprows=1, usecols=(0, 1, 2, 3), unpack=True)
    os.remove(outputFileName)
    return x1, x2, x3, var

def Slice_PV_data_over_line_to_csv_file(inputFileName,
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

def Clip_VTK_data_to_VTK(inputFileName,field_name,POINTS_or_CELLS,
                             outputFileName,
                                 point, normal,
                                 resolution
                                           ):
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'XML Unstructured Grid Reader'
    reader = XMLUnstructuredGridReader(FileName=[inputFileName])
    reader.PointArrayStatus = [field_name]

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1057, 499]

    # show data in view
    display = Show(finiteElementsResultField_3D0vtu, renderView1)
    # trace defaults for the display properties.
    display.ColorArrayName = [None, '']
    display.GlyphType = 'Arrow'
    display.ScalarOpacityUnitDistance = 0.02234159571242408

    # reset view to fit data
    renderView1.ResetCamera()

    # set scalar coloring
    ColorBy(display, (POINTS_or_CELLS, field_name))

    # rescale color and/or opacity maps used to include current data range
    finiteElementsResultField_3D0vtuDisplay.RescaleTransferFunctionToDataRange(True)

    # show color bar/color legend
    finiteElementsResultField_3D0vtuDisplay.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'Resultfield'
    resultfieldLUT = GetColorTransferFunction(field_name)

    # get opacity transfer function/opacity map for 'Resultfield'
    resultfieldPWF = GetOpacityTransferFunction(field_name)

    # create a new 'Clip'
    clip1 = Clip(Input=reader)
    clip1.ClipType = 'Plane'
    clip1.Scalars = [POINTS_or_CELLS, field_name]
    clip1.Value = 0.5000234246253967

    # init the 'Plane' selected for 'ClipType'
    clip1.ClipType.Origin = [0.5, 0.5, 0.5]

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=clip1)

    # show data in view
    clip1Display = Show(clip1, renderView1)
    # trace defaults for the display properties.
    clip1Display.ColorArrayName = [POINTS_or_CELLS, field_name]
    clip1Display.LookupTable = resultfieldLUT
    clip1Display.GlyphType = 'Arrow'
    clip1Display.ScalarOpacityUnitDistance = 0.02440176703101865

    # hide data in view
    Hide(reader, renderView1)

    # show color bar/color legend
    clip1Display.SetScalarBarVisibility(renderView1, True)

    # save animation geometry from a view
    WriteAnimationGeometry(outputFileName, view=renderView1)

    #### saving camera placements for all active views

    # current camera placement for renderView1
    renderView1.CameraPosition = [-0.14244381190851818, 1.6398475304062, 3.579637312846104]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
    renderView1.CameraViewUp = [0.4538877555358159, 0.862319473450665, -0.22447946695060406]
    renderView1.CameraParallelScale = 0.8660254037844386

    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).

def Save_VTK_data_to_picture_file(inputFileName, field_name,POINTS_or_CELLS,
                             outputFileName
                             ):
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'XML Unstructured Grid Reader'
    finiteElementsResultField_3D0vtu = XMLUnstructuredGridReader(FileName=[inputFileName])
    finiteElementsResultField_3D0vtu.PointArrayStatus = [field_name]

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1057, 499]

    # show data in view
    display = Show(finiteElementsResultField_3D0vtu, renderView1)
    # trace defaults for the display properties.
    display.ColorArrayName = [None, '']
    display.GlyphType = 'Arrow'
    display.ScalarOpacityUnitDistance = 0.02234159571242408

    # reset view to fit data
    renderView1.ResetCamera()

    # set scalar coloring
    ColorBy(display, (POINTS_or_CELLS, field_name))

    # rescale color and/or opacity maps used to include current data range
    display.RescaleTransferFunctionToDataRange(True)

    # show color bar/color legend
    finiteElementsResultField_3D0vtuDisplay.SetScalarBarVisibility(renderView1, True)

    resultfieldLUT = GetColorTransferFunction(field_name)

    resultfieldPWF = GetOpacityTransferFunction(field_name)

    renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
    renderView1.CameraParallelScale = 0.8660254037844386

    SaveScreenshot(outputFileName, magnification=1, quality=100, view=renderView1)
    

