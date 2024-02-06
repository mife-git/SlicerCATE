import os
from typing import Annotated, Optional, Union
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import csv
from csv import writer

import vtk
import qt
import traceback

import slicer
from slicer.i18n import tr as _
from slicer.i18n import translate
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
from slicer.parameterNodeWrapper import (
    parameterNodeWrapper,
    WithinRange,
)
from slicer import (
    vtkMRMLModelNode,
    vtkMRMLSegmentationNode,
    vtkMRMLMarkupsFiducialNode,
    vtkMRMLTableNode
)


#
# TortuosityEvaluator
#


class TortuosityEvaluator(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("Tortuosity Evaluator")
        self.parent.categories = [translate("qSlicerAbstractCoreModule", "SlicerCATE")]
        self.parent.dependencies = [""]  # TODO: add here list of module names that this module requires
        self.parent.contributors = ["Michela Ferrari (University of Pavia & Fondazione IRCCS Policlino San Matteo)"]
        # TODO: update with short description of the module and a link to online module documentation
        # _() function marks text as translatable to other languages
        self.parent.helpText = _("""
            This is an example of scripted loadable module bundled in an extension.
            See more information in <a href="https://github.com/organization/projectname#TortuosityEvaluator">module documentation</a>.
        """)
        # TODO: replace with organization, grant and thanks
        self.parent.acknowledgementText = _("""
            This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
                and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
            """)

        # Additional initialization step after application startup is complete
        slicer.app.connect("startupCompleted()", registerSampleData)

#
# Register sample data sets in Sample Data module
#


def registerSampleData():
    """Add data sets to Sample Data module."""
    # It is always recommended to provide sample data for users to make it easy to try the module,
    # but if no sample data is available then this method (and associated startupCompeted signal connection) can be removed.

    import SampleData

    iconsPath = os.path.join(os.path.dirname(__file__), "Resources/Icons")

    # To ensure that the source code repository remains small (can be downloaded and installed quickly)
    # it is recommended to store data sets that are larger than a few MB in a Github release.

    # TortuosityEvaluator1
    SampleData.SampleDataLogic.registerCustomSampleDataSource(
        # Category and sample name displayed in Sample Data module
        category="TortuosityEvaluator",
        sampleName="TortuosityEvaluator1",
        # Thumbnail should have size of approximately 260x280 pixels and stored in Resources/Icons folder.
        # It can be created by Screen Capture module, "Capture all views" option enabled, "Number of images" set to "Single".
        thumbnailFileName=os.path.join(iconsPath, "TortuosityEvaluator1.png"),
        # Download URL and target file name
        uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
        fileNames="TortuosityEvaluator1.nrrd",
        # Checksum to ensure file integrity. Can be computed by this command:
        #  import hashlib; print(hashlib.sha256(open(filename, "rb").read()).hexdigest())
        checksums="SHA256:998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
        # This node name will be used when the data set is loaded
        nodeNames="TortuosityEvaluator1",
    )

    # TortuosityEvaluator2
    SampleData.SampleDataLogic.registerCustomSampleDataSource(
        # Category and sample name displayed in Sample Data module
        category="TortuosityEvaluator",
        sampleName="TortuosityEvaluator2",
        thumbnailFileName=os.path.join(iconsPath, "TortuosityEvaluator2.png"),
        # Download URL and target file name
        uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
        fileNames="TortuosityEvaluator2.nrrd",
        checksums="SHA256:1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
        # This node name will be used when the data set is loaded
        nodeNames="TortuosityEvaluator2",
    )


#
# TortuosityEvaluatorParameterNode
#


@parameterNodeWrapper
class TortuosityEvaluatorParameterNode:
    """
    The parameters needed by module.

    inputModel - Surface model for centerline computation.
    centerlineSeeds - Markups for centerline computation (at least two seeds representing start and end of the vessel
    are required).
    outputModel - The output model that will contain the extracted centerline.
    curveSamplingDistance - Sampling distance between the centerline points (millimeters).
    targetPointsCount - Number of points in the preprocessed surface model.
    decimationAggressiveness - Decimation factor to reduce triangle count of the surface model while maintaining
    the visual appearance.
    subdivideInputSurface - Subdivide the input surface to make centerline computation more accurate (avoid problems
    because of slim models).
    smoothingFactor - Centerline smoothing factor.
    smoothingIterations - Number of centerline smoothing iterations.
    centerlineModel - Input centerline for tortuosity evaluation.
    outputTable - Output table that will contain the tortuosity indexes.
    """

    # Centerline extraction
    inputModel: Union[vtkMRMLModelNode, None, vtkMRMLSegmentationNode, None]
    centerlineSeeds: vtkMRMLMarkupsFiducialNode
    outputModel: vtkMRMLModelNode
    curveSamplingDistance: Annotated[float, WithinRange(0.00, 5.00)] = 0.50  # millimeters
    targetPointsCount: Annotated[float, WithinRange(0.1, 300.0)] = 5.0
    decimationAggressiveness: Annotated[float, WithinRange(0.00, 15.00)] = 4.00
    subdivideInputSurface: bool = False
    smoothingIterations: Annotated[int, WithinRange(0, 500)] = 100
    smoothingFactor: Annotated[float, WithinRange(0.00, 1.00)] = 0.1

    # Tortuosity metrics calculation
    centerlineModel: vtkMRMLModelNode
    outputTable: vtkMRMLTableNode


#
# TortuosityEvaluatorWidget
#


class TortuosityEvaluatorWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent=None) -> None:
        """Called when the user opens the module the first time and the widget is initialized."""
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)  # needed for parameter node observation
        self.logic = None
        self._parameterNode = None
        self._parameterNodeGuiTag = None

        # Define global variables
        self.reportFolderPath = os.path.join(os.path.expanduser("~"), "SlicerCATEReportFolder")
        self.trajPath = os.path.join(self.reportFolderPath, "Trajectories.png")
        self.colormapsPath = os.path.join(self.reportFolderPath, "Colormaps.png")

        # Create report folder if it does not exist
        if not os.path.exists(self.reportFolderPath):
            os.makedirs(self.reportFolderPath)
        self.imageWidget = qt.QLabel()

    def setup(self) -> None:
        """Called when the user opens the module the first time and the widget is initialized."""
        ScriptedLoadableModuleWidget.setup(self)

        # Load widget from .ui file (created by Qt Designer).
        # Additional widgets can be instantiated manually and added to self.layout.
        uiWidget = slicer.util.loadUI(self.resourcePath("UI/TortuosityEvaluator.ui"))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)

        # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
        # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
        # "setMRMLScene(vtkMRMLScene*)" slot.
        uiWidget.setMRMLScene(slicer.mrmlScene)

        # Create logic class. Logic implements all computations that should be possible to run
        # in batch mode, without a graphical user interface.
        self.logic = TortuosityEvaluatorLogic()

        # The analysis should start from centerline extraction, so other steps are collapsed
        self.ui.advancedCollapsibleButton.collapsed = True
        self.ui.tortuosityCollapsibleButton.collapsed = True

        # Add all the options for artery selection
        arteries = ["Select the artery name", "Left anterior descending artery", "Left circumflex artery",
                    "Right coronary artery"]
        for idx, artery in enumerate(arteries):
            self.ui.arterySelector.addItem(artery, idx)

        # Connections

        # These connections ensure that we update parameter node when scene is closed
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

        # Buttons
        self.ui.centerlineApplyButton.connect("clicked(bool)", self.onCenterlineApplyButton)
        self.ui.calculateTortuosityButton.connect('clicked(bool)', self.onCalculateTortuosityButton)
        self.ui.saveDataButton.connect('clicked(bool)', self.onSaveDataButton)
        self.ui.changeModuleButton.connect('clicked(bool)', self.onChangeModuleButton)

        # Make sure parameter node is initialized (needed for module reload)
        self.initializeParameterNode()

        # Make sure module dependencies are installed
        self.importRequiredModules()

    def cleanup(self) -> None:
        """Called when the application closes and the module widget is destroyed."""
        # Remove observers
        self.removeObservers()
        # Remove temp files
        self.removeTempFiles()

    def enter(self) -> None:
        """Called each time the user opens this module."""
        # Make sure parameter node exists and observed
        self.initializeParameterNode()

    def exit(self) -> None:
        """Called each time the user opens a different module."""
        # Do not react to parameter node changes (GUI will be updated when the user enters into the module)
        if self._parameterNode:
            self._parameterNode.disconnectGui(self._parameterNodeGuiTag)
            self._parameterNodeGuiTag = None
            self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self._checkCanApply)

    def onSceneStartClose(self, caller, event) -> None:
        """Called just before the scene is closed."""
        # Parameter node will be reset, do not use it anymore
        self.setParameterNode(None)
        # Remove temp files
        self.removeTempFiles()

    def onSceneEndClose(self, caller, event) -> None:
        """Called just after the scene is closed."""
        # If this module is shown while the scene is closed then recreate a new parameter node immediately
        if self.parent.isEntered:
            self.initializeParameterNode()

    def initializeParameterNode(self) -> None:
        """Ensure parameter node exists and observed."""
        # Parameter node stores all user choices in parameter values, node selections, etc.
        # so that when the scene is saved and reloaded, these settings are restored.
        self.setParameterNode(self.logic.getParameterNode())

        # Select default input nodes if nothing is selected yet to save a few clicks for the user
        if not self._parameterNode.inputModel:
            firstModelNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLModelNode")
            if firstModelNode:
                self._parameterNode.inputModel = firstModelNode

    def setParameterNode(self, inputParameterNode: Optional[TortuosityEvaluatorParameterNode]) -> None:
        """
        Set and observe parameter node.
        Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
        """

        if self._parameterNode and self.hasObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self._checkCanApply):
            self._parameterNode.disconnectGui(self._parameterNodeGuiTag)
            self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self._checkCanApply)
        self._parameterNode = inputParameterNode
        if self._parameterNode:
            # Note: in the .ui file, a Qt dynamic property called "SlicerParameterName" is set on each
            # ui element that needs connection.
            self._parameterNodeGuiTag = self._parameterNode.connectGui(self.ui)
            self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self._checkCanApply)
            self._checkCanApply()

    def _checkCanApply(self, caller=None, event=None) -> None:
        """
        Check if all the needed inputs are available before enabling "apply" buttons.
        """

        # The input for centerline computation can be of two types: model node or segmentation node
        centerlineInput = self._parameterNode.inputModel
        if centerlineInput and centerlineInput.IsA("vtkMRMLSegmentationNode"):
            self.ui.inputSegmentSelectorWidget.setVisible(True)  # allow the user to select the segment
        else:
            self.ui.inputSegmentSelectorWidget.setVisible(False)

        # Check the presence of the needed inputs before enabling the "Extract centerline" button
        if self._parameterNode and self._parameterNode.inputModel and self._parameterNode.centerlineSeeds:
            self.ui.centerlineApplyButton.toolTip = _("Compute output centerline.")
            self.ui.centerlineApplyButton.enabled = True
        else:
            self.ui.centerlineApplyButton.toolTip = _("Select input model and place centerline seeds.")
            self.ui.centerlineApplyButton.enabled = False

        # Check the presence of the needed inputs before enabling the "Calculate tortuosity metrics" button
        if self._parameterNode and self._parameterNode.centerlineModel:
            self.ui.calculateTortuosityButton.toolTip = _("Compute tortuosity metrics.")
            self.ui.calculateTortuosityButton.enabled = True
        else:
            self.ui.calculateTortuosityButton.toolTip = _("Select input centerline.")
            self.ui.calculateTortuosityButton.enabled = False

        # Check the presence of the needed inputs before enabling the "Save data" button
        if self._parameterNode and self._parameterNode.outputTable and self._parameterNode.centerlineModel:
            self.ui.saveDataButton.toolTip = _("Save data in csv and pdf report.")
            self.ui.saveDataButton.enabled = True
        else:
            self.ui.saveDataButton.toolTip = _("Select centerline and report table you want to save.")
            self.ui.saveDataButton.enabled = False

        # Check if checkboxes are checked
        if self.ui.subdivideInputSurfaceModelCheckBox.checked:
            self._parameterNode.subdivideInputSurface = True

    def onCenterlineApplyButton(self):
        """
        Run processing when user clicks "Extract centerline" button.
        Code adapted from:
        https://github.com/vmtk/SlicerExtension-VMTK/blob/master/ExtractCenterline/ExtractCenterline.py
        """
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            # Import ExtractCenterline logic from "SlicerVMTK"
            import ExtractCenterline
            ecl = ExtractCenterline.ExtractCenterlineLogic()

            # Get required inputs
            inputSurfaceModelNode = self._parameterNode.inputModel
            endPointsMarkupsNode = self._parameterNode.centerlineSeeds
            segmentID = self.ui.inputSegmentSelectorWidget.currentSegmentID()
            targetNumberOfPoints = float(self._parameterNode.targetPointsCount)*1000
            decimationAggressiveness = float(self._parameterNode.decimationAggressiveness)
            curveSamplingDistance = float(self._parameterNode.curveSamplingDistance)
            subdivideInputSurface = self._parameterNode.subdivideInputSurface
            centerlineModelNode = self._parameterNode.outputModel
            if not centerlineModelNode:  # create new node to store centerline data
                newCenterlineModelNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLModelNode")
                newCenterlineModelNode.UnRegister(None)
                newCenterlineModelNode.SetName(slicer.mrmlScene.GetUniqueNameByString("CenterlineModel"))
                centerlineModelNode = slicer.mrmlScene.AddNode(newCenterlineModelNode)
                self.ui.outputCenterlineModelSelector.setCurrentNode(centerlineModelNode)

            # Perform pre-processing
            inputSurfacePolyData = ecl.polyDataFromNode(inputSurfaceModelNode, segmentID)
            preprocessedPolyData = ecl.preprocess(inputSurfacePolyData, targetNumberOfPoints, decimationAggressiveness,
                                                  subdivideInputSurface)
            # Remove non-manifold edges
            finalModel = self.logic.removeNonManifoldEdges(preprocessedPolyData)

            # Extract centerline
            # The "curveSamplingDistance" parameter is one of the inputs, but the function sets the
            # "SetCenterlineResampling" parameter of vtkvmtkComputationalGeometry.vtkvmtkPolyDataCenterlines() to False
            # so resampling is not performed
            centerlinePolyData, voronoiDiagramPolyData = ecl.extractCenterline(finalModel,
                                                                               endPointsMarkupsNode,
                                                                               curveSamplingDistance)

            # Perform post-processing (merging, smoothing) on extracted centerline, here resampling is performed
            finalCenterlinePolyData = self.logic.postprocessCenterline(centerlinePolyData, curveSamplingDistance)

            if centerlineModelNode:
                centerlineModelNode.SetAndObservePolyData(finalCenterlinePolyData)
                if not centerlineModelNode.GetDisplayNode():
                    centerlineModelNode.CreateDefaultDisplayNodes()
                    centerlineModelNode.GetDisplayNode().SetColor(0.0, 1.0, 0.0)
                    centerlineModelNode.GetDisplayNode().SetLineWidth(3)
                    inputSurfaceModelNode.GetDisplayNode().SetOpacity(0.5)

        except Exception as e:
            slicer.util.errorDisplay("Failed to compute results: " + str(e))
            traceback.print_exc()

        qt.QApplication.restoreOverrideCursor()
        slicer.util.showStatusMessage("Centerline extraction complete.", 3000)

    def onCalculateTortuosityButton(self):
        """
        Run processing when user clicks "Calculate tortuosity metrics" button.
        """
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            # Get required inputs
            centerlineModelNode = self._parameterNode.centerlineModel
            centerlinePropertiesTableNode = self._parameterNode.outputTable
            smoothingIterations = self._parameterNode.smoothingIterations
            smoothingFactor = self._parameterNode.smoothingFactor
            if not centerlinePropertiesTableNode:
                newTableNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLTableNode")
                newTableNode.UnRegister(None)
                newTableNode.SetName(slicer.mrmlScene.GetUniqueNameByString("CenterlineQuantification"))
                centerlinePropertiesTableNode = slicer.mrmlScene.AddNode(newTableNode)
                self.ui.outputCenterlinePropertiesTableSelector.setCurrentNode(centerlinePropertiesTableNode)

            # Extract tortuosity indexes and create table to store them
            centerline = centerlineModelNode.GetPolyData()
            centerlineArray = self.logic.createPropertiesTable(centerline, centerlinePropertiesTableNode,
                                                               smoothingFactor, smoothingIterations)

            # Create 2D plots
            self.logic.createPlot2D(centerlineArray - centerlineArray[0, :], self.reportFolderPath)
            # Create static image view for plots
            pm = qt.QPixmap(self.trajPath)
            self.imageWidget.setPixmap(pm)
            self.imageWidget.setScaledContents(True)
            self.imageWidget.show()

        except Exception as e:
            slicer.util.errorDisplay("Failed to compute results: "+str(e))
            traceback.print_exc()

        qt.QApplication.restoreOverrideCursor()
        slicer.util.showStatusMessage("Tortuosity analysis complete.", 3000)

    def onSaveDataButton(self):
        """
        Save patient data as pdf (visual report) and csv (indexes).
        """
        try:
            # Get exam ID defined by the user
            examID = self.ui.IDTextEdit.toPlainText()

            # Get centerline
            centerlineModel = self._parameterNode.centerlineModel

            # Get report table and save it as .csv file in the report folder
            reportTable = self._parameterNode.outputTable
            reportTablePath = os.path.join(self.reportFolderPath, examID + "_reportTable.csv")
            slicer.util.saveNode(reportTable, reportTablePath)

            # Get table rows and columns
            nRows = reportTable.GetNumberOfRows()
            nCols = reportTable.GetNumberOfColumns()
            # Create an array with all properties stored in the report table (used to create a unique row for all
            # indexes)
            propertiesArray = [examID]
            for i in range(0, nRows):
                for j in range(1, nCols):  # skip first column because it has segments names
                    text = reportTable.GetCellText(i, j)
                    propertiesArray.append(text)

            self.logic.saveData(propertiesArray, nRows, self.reportFolderPath)
            self.logic.createColorPlot(float(propertiesArray[7]), self.colormapsPath)
            self.logic.capture3DWiews(self.reportFolderPath, centerlineModel)
            self.logic.createPdf(self.reportFolderPath, propertiesArray[0])

        except Exception as e:
            slicer.util.errorDisplay("Failed to compute results: " + str(e))
            traceback.print_exc()

        qt.QApplication.restoreOverrideCursor()
        slicer.util.showStatusMessage("Saving complete.", 3000)

    def importRequiredModules(self):
        try:
            import ExtractCenterline
        except ModuleNotFoundError:
            if slicer.util.confirmOkCancelDisplay(
                    "This module requires 'SlicerVMTK' Extension. Click OK to install it now. "
                    "Note that installation requires application restart; otherwise you can click Cancel and install "
                    "it manually using the Extension Manager."):
                extensionName = 'SlicerVMTK'
                em = slicer.app.extensionsManagerModel()
                em.interactive = False  # prevent display of popups
                restart = True
                if not em.installExtensionFromServer(extensionName, restart):
                    raise ValueError(f"Failed to install {extensionName} extension")

    def removeTempFiles(self):
        """
        Remove temporary files.
        """
        if os.path.exists(self.trajPath):
            os.remove(self.trajPath)
        if os.path.exists(self.colormapsPath):
            os.remove(self.colormapsPath)
        for axis in range(0, 6):
            if os.path.exists(os.path.join(self.reportFolderPath, "view3D_axis_" + str(axis) + ".png")):
                os.remove(os.path.join(self.reportFolderPath, "view3D_axis_" + str(axis) + ".png"))

    def onChangeModuleButton(self):
        self.logic.switchModule('ArterySegmentation')


#
# TortuosityEvaluatorLogic
#


class TortuosityEvaluatorLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self) -> None:
        """Called when the logic class is instantiated. Can be used for initializing member variables."""
        ScriptedLoadableModuleLogic.__init__(self)

    def getParameterNode(self):
        return TortuosityEvaluatorParameterNode(super().getParameterNode())

    def removeNonManifoldEdges(self, preprocessedPolyData):

        """
        Removes non-manifold edges from the preprocessed centerline model.
        Code source: https://github.com/SlicerIGT/SlicerBoneReconstructionPlanner/issues/35#issuecomment-1193403647
        """

        # Create filter to extract non-manifold edges (used by three or more polygons)
        featureEdges = vtk.vtkFeatureEdges()
        featureEdges.SetInputData(preprocessedPolyData)
        featureEdges.BoundaryEdgesOff()
        featureEdges.FeatureEdgesOff()
        featureEdges.ManifoldEdgesOff()
        featureEdges.NonManifoldEdgesOn()
        featureEdges.Update()
        nonManifoldEdgesPolyData = featureEdges.GetOutput()

        ids = vtk.vtkIdTypeArray()
        ids.SetNumberOfComponents(1)

        for i in range(nonManifoldEdgesPolyData.GetNumberOfPoints()):
            pointOfNonManifoldEdges = nonManifoldEdgesPolyData.GetPoints().GetPoint(i)
            pointIdInOriginal = preprocessedPolyData.FindPoint(pointOfNonManifoldEdges)
            ids.InsertNextValue(pointIdInOriginal)

        selectionNode = vtk.vtkSelectionNode()
        selectionNode.SetFieldType(vtk.vtkSelectionNode.POINT)
        selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
        selectionNode.SetSelectionList(ids)
        selectionNode.GetProperties().Set(vtk.vtkSelectionNode.CONTAINING_CELLS(), 1)
        selectionNode.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)

        selection = vtk.vtkSelection()
        selection.AddNode(selectionNode)

        # Create filter to extract non-manifold edges points from centerline polydata
        extractSelection = vtk.vtkExtractSelection()
        extractSelection.SetInputData(0, preprocessedPolyData)
        extractSelection.SetInputData(1, selection)
        extractSelection.Update()

        # Create filter to remove those points
        geometryFilter = vtk.vtkGeometryFilter()
        geometryFilter.SetInputData(extractSelection.GetOutput())
        geometryFilter.Update()

        # Save the model without non-manifold edges
        return geometryFilter.GetOutput()

    def postprocessCenterline(self, centerlinePolyData, curveSamplingDistance):
        """
        Perform branch splitting and centerline merging (with resampling).
        Why branch splitting? http://www.vmtk.org/tutorials/BranchSplitting.html
        Why merging? To combine multiple centerlines which lie within the same branch of a split and grouped centerline.
        Code adapted from:
        https://github.com/vmtk/SlicerExtension-VMTK/blob/master/ExtractCenterline/ExtractCenterline.py
        """

        # Import vmtk functions
        try:
            import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
        except ImportError:
            raise ImportError("VMTK library is not found")

        branchExtractor = vtkvmtkComputationalGeometry.vtkvmtkCenterlineBranchExtractor()
        branchExtractor.SetInputData(centerlinePolyData)
        branchExtractor.SetBlankingArrayName('Blanking')
        branchExtractor.SetRadiusArrayName('Radius')
        branchExtractor.SetGroupIdsArrayName('GroupIds')
        branchExtractor.SetCenterlineIdsArrayName('CenterlineIds')
        branchExtractor.SetTractIdsArrayName('TractIds')
        branchExtractor.Update()
        centerlines = branchExtractor.GetOutput()

        mergeCenterlines = vtkvmtkComputationalGeometry.vtkvmtkMergeCenterlines()
        mergeCenterlines.SetInputData(centerlines)
        mergeCenterlines.SetBlankingArrayName('Blanking')
        mergeCenterlines.SetRadiusArrayName('Radius')
        mergeCenterlines.SetGroupIdsArrayName('GroupIds')
        mergeCenterlines.SetCenterlineIdsArrayName('CenterlineIds')
        mergeCenterlines.SetTractIdsArrayName('TractIds')
        mergeCenterlines.SetResamplingStepLength(curveSamplingDistance)
        mergeCenterlines.SetMergeBlanked(True)
        mergeCenterlines.Update()
        mergedCenterlines = mergeCenterlines.GetOutput()

        return mergedCenterlines

    def createPropertiesTable(self, finalCenterlines, centerlinePropertiesTableNode,
                              smoothingFactor, smoothingIterations):
        """
        This function creates a table node containing the tortuosity indexes of the complete centerline
        and for every segment.
        """

        # Create a table with all indexes
        centerlinePropertiesTableNode.RemoveAllColumns()
        indicesNames = ['Centerline', 'Length', 'TI', 'ICM', 'SOAM', 'avgAbsCurv', 'RMSCurv', 'SqrDerCurv',
                        'avgAbsTors', '2DTX', '2DTY', '2DTZ', 'CorT45', 'CorT90']
        for name in indicesNames:
            col = centerlinePropertiesTableNode.AddColumn()
            col.SetName(name)

        # Get centerline metrics
        centerlinePoints, Curvature, Torsion, Tangent, Normal, Binormal, centerlineSegments, segmentsCurvature, \
            segmentsTorsion, TangentArray, NormalArray, BinormalArray = self.getCenterlineMetrics(finalCenterlines,
                                                                                                  smoothingFactor,
                                                                                                  smoothingIterations)

        # Get indexes
        indicesArray = self.calculateIndices(centerlinePoints, Curvature, Torsion, Tangent, Normal, Binormal)
        centerlinePropertiesTableNode.AddEmptyRow()
        centerlinePropertiesTableNode.SetCellText(0, 0, 'Complete')

        for i, index in enumerate(indicesArray):
            centerlinePropertiesTableNode.SetCellText(0, i + 1, str(round(index, 3)))

        # Segments' indexes work if there is only one centerline with two or more segments specified by user seeds
        for i in range(len(centerlineSegments)):
            indicesArray = self.calculateIndices(centerlineSegments[i], segmentsCurvature[i], segmentsTorsion[i],
                                                 TangentArray[i], NormalArray[i], BinormalArray[i])
            centerlinePropertiesTableNode.AddEmptyRow()
            for j, index in enumerate(indicesArray):
                centerlinePropertiesTableNode.SetCellText(i+1, 0, 'Segment {}'.format(i + 1))
                centerlinePropertiesTableNode.SetCellText(i+1, j+1, str(round(index, 3)))

        centerlinePropertiesTableNode.GetTable().Modified()
        return centerlinePoints

    def getCenterlineMetrics(self, mergedCenterlines, smoothingFactor, smoothingIterations):
        """
        This function calculates the following metrics from the input centerline (and for every segment speciefied by
        the input seeds):
            1. Curvature
            2. Torsion
            3. Frenet Frame (Normal, Binormal and Tangent array)
        """

        # Import vmtk functions
        try:
            import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
        except ImportError:
            raise ImportError("VMTK library is not found")

        centerlineSmoothing = vtkvmtkComputationalGeometry.vtkvmtkCenterlineSmoothing()
        centerlineGeometry = vtkvmtkComputationalGeometry.vtkvmtkCenterlineGeometry()

        # Get number of centerline segments (it depends on the number of input seeds)
        numberOfCells = mergedCenterlines.GetNumberOfCells()
        if numberOfCells < 1 or numberOfCells > 100:
            raise ValueError("cannot generate table, check the input")

        # Now get curvature and torsion for each segment
        segments, segmentsCurvature, segmentsTorsion = [], [], []
        segmentsTangent, segmentsNormal, segmentsBinormal = [], [], []
        for i in range(numberOfCells):
            # Get single segment
            centerlineSegment = mergedCenterlines.GetCell(i).GetPoints()
            # Apply smoothing
            smoothedSegment = vtk.vtkPoints()
            centerlineSmoothing.SmoothLine(centerlineSegment, smoothedSegment,
                                           smoothingIterations, smoothingFactor)
            # Compute Frenet system arrays
            TangentArray = vtk.vtkDoubleArray()
            NormalArray = vtk.vtkDoubleArray()
            BinormalArray = vtk.vtkDoubleArray()
            centerlineGeometry.ComputeLineFrenetReferenceSystem(smoothedSegment, TangentArray,
                                                                NormalArray, BinormalArray)
            # Compute curvature
            curvatureArray = vtk.vtkDoubleArray()
            centerlineGeometry.ComputeLineCurvature(smoothedSegment, curvatureArray)
            # Compute torsion
            torsionArray = vtk.vtkDoubleArray()
            centerlineGeometry.ComputeLineTorsion(smoothedSegment, torsionArray)

            # Convert to numpy and remove duplicated points before concatenating arrays
            segments.append(vtk_to_numpy(smoothedSegment.GetData()))
            if i == (numberOfCells - 1):
                segmentsCurvature.append(vtk_to_numpy(curvatureArray))
                segmentsTorsion.append(vtk_to_numpy(torsionArray))
                segmentsTangent.append(vtk_to_numpy(TangentArray))
                segmentsNormal.append(vtk_to_numpy(NormalArray))
                segmentsBinormal.append(vtk_to_numpy(BinormalArray))
            else:
                segmentsCurvature.append(vtk_to_numpy(curvatureArray)[:-1])
                segmentsTorsion.append(vtk_to_numpy(torsionArray)[:-1])
                segmentsTangent.append(vtk_to_numpy(TangentArray)[:-1])
                segmentsNormal.append(vtk_to_numpy(NormalArray)[:-1])
                segmentsBinormal.append(vtk_to_numpy(BinormalArray)[:-1])

        # Concatenate segments arrays to create the total array
        tmp = np.vstack(segments)
        _, idx, count = np.unique(tmp, axis=0, return_index=True, return_counts=True)
        total = tmp[np.sort(idx)]
        totalCurvature = np.concatenate(segmentsCurvature)
        totalTorsion = np.concatenate(segmentsTorsion)
        totalTangent = np.concatenate(segmentsTangent)
        totalNormal = np.concatenate(segmentsNormal)
        totalBinormal = np.concatenate(segmentsBinormal)

        # Add curvature and torsion arrays to the complete centerline
        totalCurvature_vtk = numpy_to_vtk(totalCurvature)
        totalCurvature_vtk.SetName("Curvature")
        mergedCenterlines.GetPointData().AddArray(totalCurvature_vtk)
        totalTorsion_vtk = numpy_to_vtk(totalTorsion)
        totalTorsion_vtk.SetName("Torsion")
        mergedCenterlines.GetPointData().AddArray(totalTorsion_vtk)

        return total, totalCurvature, totalTorsion, totalTangent, \
            totalNormal, totalBinormal, segments, segmentsCurvature, segmentsTorsion, \
            segmentsTangent, segmentsNormal, segmentsBinormal

    def calculateIndices(self, curve, Curvature, Torsion, Tangent, Normal, Binormal):
        """
        This function calculates the all the tortuosity indexes.
        """

        # Tortuosity Index
        L = self.calculateCurveLength(curve)
        D = np.sqrt(sum((curve[-1, :] - curve[0, :]) ** 2))
        TI = L/D

        # Inflection Count Metric
        deltaN = Normal[1:-1, :] - Normal[:-2, :]
        dotN = np.sum(deltaN * deltaN, axis=1)
        ind = np.where(dotN > 1)
        coord = curve[ind, :]
        n = len(coord[0, :, :])  # number of inflection points
        ICM = n*TI

        # Sum Of Angles Metric
        T1, T2, T3 = [], [], []
        for i in range(1, len(curve)-2):
            T1.append(curve[i, :] - curve[i-1, :])
            T2.append(curve[i+1, :] - curve[i, :])
            T3.append(curve[i+2, :] - curve[i+1, :])

        T1_norm = np.array(T1)/np.linalg.norm(T1, axis=1, keepdims=True)
        T2_norm = np.array(T2)/np.linalg.norm(T2, axis=1, keepdims=True)
        cos_sum1 = np.sum(T1_norm * T2_norm, axis=1)
        err = np.where(cos_sum1 > 1)  # remove errors due to noise (when there is an angle with a cosine > 1)
        cos_sum1[err] = 1
        T1T2_cross = np.cross(T1, T2)
        T1T2_cross_norm = np.linalg.norm(T1T2_cross, axis=1, keepdims=True)
        T2T3_cross = np.cross(T2, T3)
        T2T3_cross_norm = np.linalg.norm(T2T3_cross, axis=1, keepdims=True)
        cos_sum2 = np.sum(T1T2_cross_norm * T2T3_cross_norm, axis=1)
        err = np.where(cos_sum2 > 1)
        cos_sum2[err] = 1

        IP = np.arccos(cos_sum1)  # In Plane angles
        TP = np.arccos(cos_sum2)  # Torsional angles
        CP = np.sqrt((IP * IP) + (TP * TP))
        SOAM = sum(CP)/L

        # Average absolute curvature
        avgAbsCurvature = np.sum(Curvature)/L  # no need to do abs because it's already done in vmtk

        # Average absolute torsion
        avgAbsTorsion = np.sum(np.abs(Torsion))/L

        # Root mean square curvature
        RMSCurvature = np.sqrt(np.sum(Curvature**2)/L)

        # Average squared-derivative-curvature
        firstDer = np.diff(Curvature)
        avgSquareDerCurvature = np.sum(firstDer**2)/L * 1000  # multiply by 1000 because it's a small value

        # CorT and local TI (LT)
        # Find the first point after 0.5 cm (5 mm) from the ostium and last point 0.5 cm before the end
        s = self.calculateSpaceCoordinates(curve)
        dist = 5  # mm
        start = np.argmax(s >= dist)  # argmax stops at first index that evaluates condition to True
        stop = np.argmax(s > s[-1] - dist)
        angles = []
        LT = []
        for i in range(start, stop):
            prev = np.argmax(s >= s[i] - dist)
            nex = np.argmax(s >= s[i] + dist)
            up = self.leastSquareFit3D(curve[prev:i, :])
            down = self.leastSquareFit3D(curve[i:nex, :])

            if up @ (curve[i] - curve[prev]) < 0:
                up = -up
            if down @ (curve[nex - 1] - curve[i]) < 0:
                down = -down

            angle = np.degrees(np.arccos(up @ down))
            angles.append(angle)

        angles = np.array(angles)
        intervals_45 = []
        intervals_90 = []
        start_idx_45 = None
        start_idx_90 = None
        for i, angle in enumerate(angles):
            if angle >= 45 and start_idx_45 is None:
                start_idx_45 = i
            elif angle >= 90 and start_idx_90 is None:
                start_idx_90 = i
            elif angle < 45 and start_idx_45 is not None:
                intervals_45.append((start_idx_45, i - 1))
                start_idx_45 = None
            elif angle < 90 and start_idx_90 is not None:
                intervals_90.append((start_idx_90, i - 1))
                start_idx_90 = None
        if start_idx_45 is not None:  # if last interval continues until the end of the list
            intervals_45.append((start_idx_45, len(angles) - 1))
        if start_idx_90 is not None:  # if last interval continues until the end of the list
            intervals_90.append((start_idx_90, len(angles) - 1))

        corT45 = len(intervals_45) - len(intervals_90)
        corT90 = len(intervals_90)

        # 2D Trajectories
        critx, crit_ptsx, crity, crit_ptsy, critz, crit_ptsz = self.calculateMaxMin2D(curve)
        inflx, infl_ptsx, infly, infl_ptsy, inflz, infl_ptsz = self.calculateInflectionPoints2D(curve)
        anglesX = self.calculateAngles(curve[:, 0], s, inflx, critx)
        anglesY = self.calculateAngles(curve[:, 1], s, infly, crity)
        anglesZ = self.calculateAngles(curve[:, 2], s, inflz, critz)
        LX = self.calculateCurveLength((np.vstack((s, curve[:, 0]))).T)
        LY = self.calculateCurveLength((np.vstack((s, curve[:, 1]))).T)
        LZ = self.calculateCurveLength((np.vstack((s, curve[:, 2]))).T)
        _2DTX = sum(anglesX) / LX
        _2DTY = sum(anglesY) / LY
        _2DTZ = sum(anglesZ) / LZ

        return [L, TI, ICM, SOAM, avgAbsCurvature, RMSCurvature, avgSquareDerCurvature, avgAbsTorsion,
                _2DTX, _2DTY, _2DTZ, corT45, corT90]

    def calculateCurveLength(self, curve):
        """
        This function calculates the length of a curve.
        """
        diffs = curve[1:, :] - curve[:-1, :]
        curveLength = (((diffs ** 2).sum(1)) ** 0.5).sum()
        return curveLength

    def calculateSpaceCoordinates(self, curve):
        """
        This function returns the curvilinear abscissa parameter "s" of a curve.
        """
        s = np.zeros(len(curve))
        diffs = curve[1:] - curve[:-1]
        s[1:] = ((diffs ** 2).sum(1)) ** 0.5
        s = (np.cumsum(s))
        return s

    def leastSquareFit3D(self, curve):

        # Calculate the mean of the points, i.e. the 'center' of the cloud
        datamean = curve.mean(axis=0)

        # Do an SVD on the mean-centered data.
        uu, dd, vv = np.linalg.svd(curve - datamean, full_matrices=False)

        # Now vv[0] contains the first principal component, i.e. the direction
        # vector of the 'best fit' line in the least squares sense.

        return vv[0]

    def calculateMaxMin2D(self, curve):
        """
        This function finds the coordinates of maximum and minimum values of a 2D curve.
        """

        x = curve[:, 0]
        y = curve[:, 1]
        z = curve[:, 2]

        dx = np.gradient(x)
        dy = np.gradient(y)
        dz = np.gradient(z)

        # Since these are numerical derivatives, it is very likely that they are never exactly equal to zero
        # for this reason we look for points in which the derivative of the next point changes its sign compared to
        # the previous point (we do the same thing for each axis)
        sign = np.sign(dx)
        sign[sign == 0] = -1
        critx = np.where(np.diff(sign))[0]  # coordinates
        crit_ptsx = x[critx]  # values

        sign = np.sign(dy)
        sign[sign == 0] = -1
        crity = np.where(np.diff(sign))[0]
        crit_ptsy = y[crity]

        sign = np.sign(dz)
        sign[sign == 0] = -1
        critz = np.where(np.diff(sign))[0]
        crit_ptsz = z[critz]

        return critx, crit_ptsx, crity, crit_ptsy, critz, crit_ptsz

    def calculateInflectionPoints2D(self, curve):
        """
        This function finds the inflection points of a 2D curve (coordinates and values).
        """

        x = curve[:, 0]
        y = curve[:, 1]
        z = curve[:, 2]

        dx = np.gradient(x)
        dy = np.gradient(y)
        dz = np.gradient(z)

        dxx = np.gradient(dx)
        dyy = np.gradient(dy)
        dzz = np.gradient(dz)

        # Since these are numerical derivatives, it is very likely that they are never exactly equal to zero
        # for this reason we look for points in which the derivative of the next point changes its sign compared to
        # the previous point (we do the same thing for each axis)
        sign = np.sign(dxx)
        sign[sign == 0] = -1
        inflx = np.where(np.diff(sign))[0]
        infl_ptsx = x[inflx]

        sign = np.sign(dyy)
        sign[sign == 0] = -1
        infly = np.where(np.diff(sign))[0]
        infl_ptsy = y[infly]

        sign = np.sign(dzz)
        sign[sign == 0] = -1
        inflz = np.where(np.diff(sign))[0]
        infl_ptsz = z[inflz]

        return inflx, infl_ptsx, infly, infl_ptsy, inflz, infl_ptsz

    def calculateAngles(self, curve, s, infl, crit):
        """
        This function calculates the angle between inflection points which has as its vertex a point of
        local maximum or minimum.
        """

        fl = np.hstack((0, infl, len(curve) - 1))
        curve = (np.vstack((s, curve))).T

        angle_deg = []
        for i in range(0, len(fl) - 1):
            for j in range(0, len(crit)):
                if crit[j] > fl[i] and crit[j] < fl[i + 1]:
                    v1 = curve[fl[i]] - curve[crit[j]]
                    v2 = curve[fl[i + 1]] - curve[crit[j]]
                    dotProd = np.dot(v1, v2)
                    d1 = np.sqrt(sum((v1 ** 2)))
                    d2 = np.sqrt(sum((v2 ** 2)))
                    angle_rad = np.arccos(dotProd / (d1 * d2))
                    angle_deg.append(np.degrees(angle_rad))

        return 180 - np.asarray(angle_deg)

    def createPlot2D(self, curve, savePath):
        """
        This function create a 2D plot of the centerline?s trajectories (x vs s, y vs s, z vs s).
        """
        try:
            import matplotlib
            import matplotlib.pyplot as plt
            matplotlib.use("Agg")
        except ModuleNotFoundError:
            if slicer.util.confirmOkCancelDisplay(
                    "This module requires 'matplotlib' Python package. Click OK to install it now."):
                slicer.util.pip_install("matplotlib")
                import matplotlib
                import matplotlib.pyplot as plt
                matplotlib.use("Agg")

        # Use curvilinear abscissa
        s = self.calculateSpaceCoordinates(curve)

        # Create a subplot
        matplotlib.rcParams['toolbar'] = 'None'
        fig1, ax = plt.subplots(3, 1, constrained_layout=True)
        ax[0].plot(s, curve[:, 0])
        ax[0].set_xlabel("Length (mm)")
        ax[0].set_ylabel("X (mm)")
        ax[0].grid(True)
        ax[0].set_aspect('equal', 'box')
        ax[1].plot(s, curve[:, 1])
        ax[1].set_xlabel("Length (mm)")
        ax[1].set_ylabel("Y (mm)")
        ax[1].grid(True)
        ax[1].set_aspect('equal', 'box')
        ax[2].plot(s, curve[:, 2])
        ax[2].set_xlabel("Length (mm)")
        ax[2].set_ylabel("Z (mm)")
        ax[2].grid(True)
        ax[2].set_aspect('equal', 'box')
        fig1.suptitle("2D Trajectories", fontsize=14)

        fig1.savefig(os.path.join(savePath, "Trajectories.png"), dpi=100)
        plt.close()

    def saveData(self, propertiesArray, nSegments, savePath):

        reportPath = os.path.join(savePath, "Indices_report.csv")
        if not os.path.exists(reportPath):
            firstRow = ['ID', 'Length', 'TI', 'ICM', 'SOAM', 'avgAbsCurv', 'RMSCurv', 'SqrDerCur', 'avgAbsTorsion',
                        '2DTX', '2DTY', '2DTZ', 'CorT45', 'CorT90']
            # firstRow is repeated for each segment
            completeFirstRow = []
            for column in firstRow:
                completeFirstRow.append(column)
            for i in range(1, nSegments):  # Add a column for each segment
                newCols = [f"{column}{i}" for column in firstRow[1:]]
                completeFirstRow.extend(newCols)
            with open(reportPath, 'w', newline='') as f_object:
                writer_object = writer(f_object)
                writer_object.writerow(completeFirstRow)
                f_object.close()

        with open(reportPath, 'a', newline='') as f_object:
            writer_object = writer(f_object)
            writer_object.writerow(propertiesArray)
            f_object.close()

    def createPdf(self, savePath, examID):

        reportTablePath = os.path.join(savePath, examID + "_reportTable.csv")
        reportPath = os.path.join(savePath, examID + "_report.pdf")

        try:
            from fpdf import FPDF
        except ModuleNotFoundError:
            if slicer.util.confirmOkCancelDisplay(
                    "This module requires 'fpdf2' Python package. Click OK to install it now."):
                slicer.util.pip_install("fpdf2")
                from fpdf import FPDF

        pdf = FPDF(orientation='L')
        pdf.add_page()
        ch = 8  # cell height
        pdf.set_font('helvetica', 'B', size=18)
        pdf.cell(w=0, h=20, txt="Coronary tortuosity report", ln=1)
        pdf.ln(ch)

        pdf.set_font('helvetica', size=12)
        pdf.cell(w=0, h=20, txt="Tortuosity Indices", ln=1)
        pdf.set_font('helvetica', size=8)

        with open(reportTablePath, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                for item in row:
                    pdf.cell(20, 8, str(item), border=1)
                pdf.ln(ch)

        pdf.ln(ch)
        pdf.set_font('helvetica', size=12)
        pdf.multi_cell(w=0, h=5, txt='Classification')
        pdf.ln(2)
        pdf.set_font('helvetica', size=10)
        pdf.multi_cell(w=500, h=5, txt='This bar indicates the value of the avg curvature index (blue dot) with'
                                       ' respect to the distribution of cases')
        pdf.ln(ch)
        colormapsPath = os.path.join(savePath, "Colormaps.png")
        pdf.image(colormapsPath, x=10, y=None, w=250, h=0, type='PNG')

        pdf.add_page(orientation='L')

        pdf.ln(ch)
        trajPath = os.path.join(savePath, "Trajectories.png")
        pdf.image(trajPath, x=pdf.w / 4, y=pdf.h / 4, w=150, h=0, type='PNG')

        orientations = ["From left to right and from right to left", "From posterior to anterior and from anterior "
                        "to posterior", "From inferior to superior and from superior to inferior"]
        i = 0
        for axis in range(0, 6):
            if axis % 2 == 0:
                pdf.add_page(orientation='L')
                pdf.cell(w=500, h=5, txt=orientations[i])
                i += 1
            viewPath = os.path.join(savePath, "view3D_axis_" + str(axis) + ".png")
            pdf.image(viewPath, x=pdf.w / 6, y=None, w=150, h=0, type='PNG')

        pdf.output(reportPath)

    def createColorPlot(self, index, savePath):

        if index >= 0.25:  # avg abs curvature
            index = 0.248

        import matplotlib.pyplot as plt
        new_cmap = self.truncateColormap(plt.get_cmap('hsv'), 0.375, 0.0)

        index_range = (np.arange(0, 0.25, 0.001))
        fig, ax = plt.subplots(figsize=(25, 0.75))
        x = index_range
        extent = [x[0] - (x[1] - x[0]) / 2., x[-1] + (x[1] - x[0]) / 2., 0, 1]
        ax.imshow(x[np.newaxis, :], cmap=new_cmap, aspect="auto", extent=extent)
        ax.plot(index, 0.5, 'bo', markersize=12)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlim(extent[0], extent[1])

        fig.savefig(savePath)
        plt.close()

    def truncateColormap(self, cmap, minval=0.0, maxval=1.0, n=100):

        import matplotlib.colors as colors
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    def capture3DWiews(self, savePath, centerlineModelNode):
        """
        Capture 3D view in png file with transparent background
        axis = 0 --> left to right
        axis = 1 --> right to left
        axis = 2 --> posterior to anterior
        axis = 3 --> anterior to posterior
        axis = 4 --> inferior to superior
        axis = 5 --> superior to inferior
        """

        # Set up centerline coloring by curvature
        centerlineModelNode.GetDisplayNode().SetActiveScalar("Curvature", vtk.vtkAssignAttribute.POINT_DATA)
        # centerlineModelNode.GetDisplayNode().SetActiveScalar("Torsion", vtk.vtkAssignAttribute.POINT_DATA)
        centerlineModelNode.GetDisplayNode().SetAndObserveColorNodeID(
            "vtkMRMLColorTableNodeFileColdToHotRainbow.txt")
        centerlineModelNode.GetDisplayNode().SetScalarVisibility(True)
        centerlineModelNode.GetDisplayNode().SetScalarRangeFlag(0)
        centerlineModelNode.GetDisplayNode().SetScalarRange(0.1, 0.5)  # TODO: set optimal scalar range
        colorLegendDisplayNode = slicer.modules.colors.logic().AddDefaultColorLegendDisplayNode(centerlineModelNode)
        colorLegendDisplayNode.SetLabelFormat("%4.2f mm-1")
        colorLegendDisplayNode.SetTitleText("Curvature")

        # Set up view layout and content
        slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUp3DView)

        # Set background to black (required for transparent background)
        view = slicer.app.layoutManager().threeDWidget(0).threeDView()
        view.mrmlViewNode().SetBackgroundColor(0, 0, 0)
        view.mrmlViewNode().SetBackgroundColor2(0, 0, 0)

        # Hide box and axes labels
        view.mrmlViewNode().SetBoxVisible(0)
        view.mrmlViewNode().SetAxisLabelsVisible(0)

        # Hide other model, markups and segmentation nodes
        modelNodes = slicer.util.getNodesByClass("vtkMRMLModelNode")
        for model in modelNodes:
            model.GetDisplayNode().SetVisibility(False)
        segmentationNodes = slicer.util.getNodesByClass("vtkMRMLSegmentationNode")
        for segmentation in segmentationNodes:
            segmentation.GetDisplayNode().SetVisibility(False)
        markupsNodes = slicer.util.getNodesByClass("vtkMRMLMarkupsFiducialNode")
        for markup in markupsNodes:
            markup.GetDisplayNode().SetVisibility(False)

        # Display only the selected centerline
        centerlineModelNode.GetDisplayNode().SetVisibility(True)

        view.resetFocalPoint()
        view.resetCamera()
        for axis in range(0, 6):
            view.rotateToViewAxis(axis)
            # Update the view
            view.forceRender()
            # Capture RGBA image
            renderWindow = view.renderWindow()
            renderWindow.SetAlphaBitPlanes(1)
            wti = vtk.vtkWindowToImageFilter()
            wti.SetInputBufferTypeToRGBA()
            wti.SetInput(renderWindow)
            wti.Update()
            writer = vtk.vtkPNGWriter()
            file_name = "view3D_axis_" + str(axis) + ".png"
            writer.SetFileName(os.path.join(savePath, file_name))
            writer.SetInputConnection(wti.GetOutputPort())
            writer.SetCompressionLevel(1)
            writer.Update()
            writer.Write()

        # Restore original background color and setup
        view.mrmlViewNode().SetBackgroundColor(0.7568627450980392, 0.7647058823529411, 0.9098039215686274)
        view.mrmlViewNode().SetBackgroundColor2(0.4549019607843137, 0.47058823529411764, 0.7450980392156863)
        view.mrmlViewNode().SetBoxVisible(1)
        view.mrmlViewNode().SetAxisLabelsVisible(1)
        slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)

        view.resetFocalPoint()
        view.resetCamera()
        # Update the view
        view.forceRender()

        # Restore other model, markups and segmentation nodes
        modelNodes = slicer.util.getNodesByClass("vtkMRMLModelNode")
        for model in modelNodes:
            model.GetDisplayNode().SetVisibility(True)
        segmentationNodes = slicer.util.getNodesByClass("vtkMRMLSegmentationNode")
        for segmentation in segmentationNodes:
            segmentation.GetDisplayNode().SetVisibility(True)
        markupsNodes = slicer.util.getNodesByClass("vtkMRMLMarkupsFiducialNode")
        for markup in markupsNodes:
            markup.GetDisplayNode().SetVisibility(True)

        # Restore scalar visibility
        centerlineModelNode.GetDisplayNode().SetScalarVisibility(False)

        # Remove legend
        slicer.mrmlScene.RemoveNode(colorLegendDisplayNode)

    def switchModule(self, moduleName):
        """
        Switch to another module
        """
        pluginHandlerSingleton = slicer.qSlicerSubjectHierarchyPluginHandler.instance()
        pluginHandlerSingleton.pluginByName('Default').switchToModule(moduleName)


#
# TortuosityEvaluatorTest
#


class TortuosityEvaluatorTest(ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """Do whatever is needed to reset the state - typically a scene clear will be enough."""
        slicer.mrmlScene.Clear()

    def runTest(self):
        """Run as few or as many tests as needed here."""
        self.setUp()
        self.test_TortuosityEvaluator1()

    def test_TortuosityEvaluator1(self):
        """Ideally you should have several levels of tests.  At the lowest level
        tests should exercise the functionality of the logic with different inputs
        (both valid and invalid).  At higher levels your tests should emulate the
        way the user would interact with your code and confirm that it still works
        the way you intended.
        One of the most important features of the tests is that it should alert other
        developers when their changes will have an impact on the behavior of your
        module.  For example, if a developer removes a feature that you depend on,
        your test should break so they know that the feature is needed.
        """

        self.delayDisplay("Starting the test")

        # Get/create input data

        import SampleData

        registerSampleData()
        inputVolume = SampleData.downloadSample("TortuosityEvaluator1")
        self.delayDisplay("Loaded test data set")

        inputScalarRange = inputVolume.GetImageData().GetScalarRange()
        self.assertEqual(inputScalarRange[0], 0)
        self.assertEqual(inputScalarRange[1], 695)

        outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
        threshold = 100

        # Test the module logic

        logic = TortuosityEvaluatorLogic()

        # Test algorithm with non-inverted threshold
        logic.process(inputVolume, outputVolume, threshold, True)
        outputScalarRange = outputVolume.GetImageData().GetScalarRange()
        self.assertEqual(outputScalarRange[0], inputScalarRange[0])
        self.assertEqual(outputScalarRange[1], threshold)

        # Test algorithm with inverted threshold
        logic.process(inputVolume, outputVolume, threshold, False)
        outputScalarRange = outputVolume.GetImageData().GetScalarRange()
        self.assertEqual(outputScalarRange[0], inputScalarRange[0])
        self.assertEqual(outputScalarRange[1], inputScalarRange[1])

        self.delayDisplay("Test passed")
