__author__ = 'ymelniko'

from vtk import *
import numpy as np



def read_exodus_mesh(filename):
    """
    Reads unstructured grid from an Exodus file (.e, .exo)
    :param filename:
    :return:
    """

    reader =vtk.vtkExodusIIReader()
    reader.SetFileName(filename)
    reader.Update() # Needed because of GetScalarRange
    exodus_structure = reader.GetOutput()

    geom =vtk.vtkCompositeDataGeometryFilter()
    geom.SetInput(exodus_structure)
    geom.Update()
    grid= geom.GetOutput()

    return grid



def read_unstructured_grid(filepath):
    """
    Reads unstructured grid from a VTU file
    :param filepath: Path to *.VTU file, containing unstructured grid
    :return: instanse of vtkXMLUnstructuredGridReader, containing pre-processed grid
    """
    reader =vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filepath)
    reader.Update()
    grid = reader.GetOutput()
    return grid




def apply_quality(grid,metric):
    """
    Applies vtkMeshQuality filter to an input grid, for a more detailed description
    visit http://www.vtk.org/doc/nightly/html/classvtkMeshQuality.html#details
    :param grid: input grid, instance of vtkUnstructuredGrid,
           metric: metric name,string
    :return: instance of vtkMeshQuality
    """

    qualityFilter=vtk.vtkMeshQuality()
    qualityFilter.SetInput(grid)

    metric_fun_name = 'SetTetQualityMeasureTo'+metric

    getattr(qualityFilter, metric_fun_name)()

    qualityFilter.Update()


    return qualityFilter


def filter_quality(grid, qmin=0.0, qmax=float("inf"), array="Quality"):
    """
    Finds elements which quality measure is qmin<q<qmax
    :param grid: input grid, instance of vtkUnstructuredGrid
    :param qmin: quality lower bound
    :param qmax: quality upper bound
    :param array: name of mesh quality cellular scalar field. By default, it is the
           same, as vtkMeshQuality produce
    :return:
    """
    threshold = vtk.vtkThreshold()
    threshold.SetInput(grid)
    threshold.ThresholdBetween(qmin, qmax)
    threshold.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, array)
    threshold.Update()
    return threshold.GetOutput()



def yield_background_grid_actor(grid):


    mapper = vtk.vtkDataSetMapper()
    mapper.SetInput(grid)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.0,1.0,0.0)
    actor.GetProperty().SetOpacity(1.0)
    actor.GetProperty().LightingOff()
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetEdgeColor(0.36, 0.36, 0.36)
    actor.GetProperty().SetLineWidth(0.1)
    return actor


# see http://www.kennethmoreland.com/color-maps/
def yield_lookup_table():
    lookup_table = vtk.vtkLookupTable()
    num_colors = 256
    lookup_table.SetNumberOfTableValues(num_colors)

    transfer_func = vtk.vtkColorTransferFunction()
    transfer_func.SetColorSpaceToDiverging()
    transfer_func.AddRGBPoint(0, 0.630, 0.299,  0.954)
    transfer_func.AddRGBPoint(1, 0.206, 0.700, 0.150)

    for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
        cc = transfer_func.GetColor(ss)
        lookup_table.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)

    return lookup_table


def yield_filtered_grid_actor(grid, quality_min, quality_max):
    lookup_table = yield_lookup_table()

    mapper = vtk.vtkDataSetMapper()
    mapper.SetInput(grid)
    mapper.SetScalarRange(quality_min, quality_max)
    mapper.SetLookupTable(lookup_table)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().EdgeVisibilityOn()

    legend_actor = vtk.vtkScalarBarActor()
    legend_actor.SetLookupTable(lookup_table)

    return actor, legend_actor


def yield_corner_annotation(metric, qmin, qmax, num_cells_showed):
    metric_info = "Quality metric is: %s" % metric
    filter_info = "Filtering criteria is: %s <= quality <= %s. Number of cells found is %s" %\
        (qmin, qmax, num_cells_showed)

    corner_annotation = vtk.vtkCornerAnnotation()
    corner_annotation.SetLinearFontScaleFactor(2)
    corner_annotation.SetNonlinearFontScaleFactor(1)
    corner_annotation.SetMaximumFontSize(20)
    corner_annotation.SetText(0, metric_info)
    corner_annotation.SetText(2, filter_info)
    #color = vtk.vtkNamedColors().GetColor3d('warm_grey')
    color = [0.5,0.5,0.5]
    corner_annotation.GetTextProperty().SetColor(color)
    return corner_annotation


def visualise_quality(grid, filtered_grid,qmin,qmax, metric):

    filtered_grid_actor, legend_actor = yield_filtered_grid_actor(filtered_grid, qmin, qmax)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(yield_background_grid_actor(grid))
    renderer.AddActor(filtered_grid_actor)
    renderer.AddActor(legend_actor)
    renderer.SetBackground(1, 1, 1)
    renderer.AddViewProp(yield_corner_annotation(metric, qmin, qmax, filtered_grid.GetNumberOfCells()))

    renderer_window = vtk.vtkRenderWindow()
    renderer_window.AddRenderer(renderer)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)

    interactor_style = vtk.vtkInteractorStyleTrackballCamera()  # ParaView-like interaction
    interactor.SetInteractorStyle(interactor_style)

    interactor.Initialize()
    interactor.Start()

    return



def get_mesh_quality(filename,metric,qmin=None,qmax=None, n_int=None):

    # Get vtk grid
    filename_ext = filename.split('.')[-1]
    if filename_ext=='vtu':
        grid = read_unstructured_grid(filename)
    else:
        if filename_ext=='e' or filename_ext=='exo':

           grid = read_exodus_mesh(filename)
        else:
            print 'Wrong mesh format'
            raise Exception


    quality = apply_quality(grid, metric)

    quality_min = quality.GetOutput().GetFieldData()\
        .GetArray('Mesh Tetrahedron Quality').GetComponent(0, 0)
    quality_max = quality.GetOutput().GetFieldData()\
        .GetArray('Mesh Tetrahedron Quality').GetComponent(0, 2)

    if qmin is None:
        qmin = quality_min
    if qmax is None:
        qmax = quality_max
    if qmin > qmax:
        print 'Lower bound qmin is bigger than upper bound qmax. Check your input.'
        raise Exception
    if n_int is None:
       n_int = 5

    print qmin, qmax, quality_min, quality_max



    # Find statistics for the metric

    values_range=[qmin,qmax]
    divide_range = np.linspace(qmin,qmax,n_int+1)
    intervals=[ (a,b) for a,b in zip (divide_range[:-1], divide_range[1:])]

    grid_size =grid.GetCellData().GetNumberOfTuples()

    filtered_grid =  filter_quality(quality.GetOutput(), qmin, qmax)
    qualityArray=    filtered_grid.GetCellData().GetArray("Quality")
    num_of_elements_filter = len([qualityArray.GetValue(i) for i in range(qualityArray.GetNumberOfTuples())])

    f = open('mesh_quality'+filename+'.txt','w')
    f.write('Mesh from %s \n' %filename)

    f.write('Metric name: %s  \t  Total number of elements: %s\n' %(metric, grid_size))


    f.write('Range  \t \t  Number of elements \t Percentage of elements \n')

    f.write('[%10.2f %10.2f]  \t %s  \t %10.2f \n' %(qmin,qmax ,num_of_elements_filter, num_of_elements_filter/float(grid_size)*100.0))

    f.write('Split range \n')
    # Interval statistics
    for interval in intervals:
        qmin_int = interval[0]
        qmax_int = interval[1]

        filtered_grid_int = filter_quality(quality.GetOutput(), qmin_int, qmax_int)
        qualityArray=filtered_grid_int.GetCellData().GetArray("Quality")
        quality_values = [qualityArray.GetValue(i) for i in range(qualityArray.GetNumberOfTuples())]
        num_of_elements = len(quality_values)

        print interval, len(quality_values)
        f.write('[%10.2f %10.2f]  \t %s  \t %10.2f \n' %(qmin_int,qmax_int,num_of_elements, num_of_elements/float(grid_size)*100.0))


    f.close()

    # Visualise mesh quality

    visualise_quality(grid, filtered_grid, qmin, qmax, metric)

    return







