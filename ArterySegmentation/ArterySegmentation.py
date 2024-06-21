# slicer imports
import traceback
import vtk, qt, ctk, slicer
from typing import Annotated, Optional
from slicer.i18n import tr as _
from slicer.i18n import translate
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
from slicer.parameterNodeWrapper import (
    parameterNodeWrapper,
    WithinRange,
)
from slicer import (
    vtkMRMLScalarVolumeNode,
    vtkMRMLMarkupsFiducialNode,
    vtkMRMLModelNode,
    vtkMRMLLabelMapVolumeNode
)


#
# ArterySegmentation
#


class ArterySegmentation(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("Artery Segmentation")
        self.parent.categories = [translate("qSlicerAbstractCoreModule", "SlicerCATE")]
        self.parent.dependencies = []  # TODO: add here list of module names that this module requires
        self.parent.contributors = ["Michela Ferrari (University of Pavia & Fondazione IRCCS Policlino San Matteo)"]
        # TODO: update with short description of the module and a link to online module documentation
        # _() function marks text as translatable to other languages
        self.parent.helpText = _(""" Documentation is available at: 
            <a href="https://github.com/organization/projectname#ArterySegmentation">module documentation</a>.
        """)
        # TODO: replace with organization, grant and thanks
        self.parent.acknowledgementText = _("""""")


#
# ArterySegmentationParameterNode
#


@parameterNodeWrapper
class ArterySegmentationParameterNode:
    """
    The parameters needed by module.

    inputVolume - Input volume for segmentation.
    vesselnessVolume - Filtered volume (from preprocessing).
    segmentationSeeds - Seed points for initialization and evolution.
    segmentationStoppers - Stoppers for initialization and evolution (only with Fast Marching method).
    segmentationLabelmap - Output binary labelmap.
    segmentationModel - Output segmentation model.
    segmentationInflation - Weight assigned to model inflation.
    segmentationCurvature - Weight assigned to surface regularization.
    segmentationAttraction - Weight that regulates the attraction of the surface of the image gradient modulus ridges.
    segmentationIterations - Number of deformation steps the model will perform.
    """

    inputVolume: vtkMRMLScalarVolumeNode
    vesselnessVolume: vtkMRMLScalarVolumeNode
    segmentationSeeds: vtkMRMLMarkupsFiducialNode
    segmentationStoppers: vtkMRMLMarkupsFiducialNode
    segmentationLabelmap: vtkMRMLLabelMapVolumeNode
    segmentationModel: vtkMRMLModelNode
    segmentationInflation: Annotated[float, WithinRange(-100.00, 100.00)] = 0.00
    segmentationCurvature: Annotated[float, WithinRange(-100.00, 100.00)] = 0.00
    segmentationAttraction: Annotated[float, WithinRange(-100.00, 100.00)] = 100.00
    segmentationIterations: Annotated[int, WithinRange(0, 100)] = 10


#
# ArterySegmentationWidget
#


class ArterySegmentationWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
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

    def setup(self) -> None:
        """Called when the user opens the module the first time and the widget is initialized."""
        ScriptedLoadableModuleWidget.setup(self)

        # Load widget from .ui file (created by Qt Designer).
        # Additional widgets can be instantiated manually and added to self.layout.
        uiWidget = slicer.util.loadUI(self.resourcePath("UI/ArterySegmentation.ui"))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)

        # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
        # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
        # "setMRMLScene(vtkMRMLScene*)" slot.
        uiWidget.setMRMLScene(slicer.mrmlScene)

        # Create logic class. Logic implements all computations that should be possible to run
        # in batch mode, without a graphical user interface.
        self.logic = ArterySegmentationLogic()

        # Connections

        # These connections ensure that we update parameter node when scene is closed
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

        # Buttons
        self.ui.segmentationPreviewButton.connect('clicked(bool)', self.onSegmentationPreviewButton)
        self.ui.segmentationApplyButton.connect('clicked(bool)', self.onSegmentationApplyButton)
        self.ui.changeModuleButton1.connect('clicked(bool)', self.onChangeModuleButton1)
        self.ui.changeModuleButton2.connect('clicked(bool)', self.onChangeModuleButton2)

        # Other connections
        self.ui.segmentationInputSelector.connect('currentNodeChanged(vtkMRMLNode*)', self.updateThresholdRange)
        self.ui.segmentationVesselnessSelector.connect('currentNodeChanged(vtkMRMLNode*)', self.updateThresholdRange)
        self.ui.segmentationThresholdSlider.connect('valuesChanged(double,double)', self.onThresholdSliderChanged)
        self.ui.fastMarchingRadio.connect('clicked(bool)', self.onInitializationRadio)
        self.ui.collidingFrontsRadio.connect('clicked(bool)', self.onInitializationRadio)

        # Collapse advanced steps
        self.ui.advancedSegmentationCollapsibleButton.collapsed = True

        # Make sure parameter node is initialized (needed for module reload)
        self.initializeParameterNode()

        # Initialize radio buttons
        self.onInitializationRadio()

        # Make sure module dependencies are installed
        self.importRequiredModules()

    def cleanup(self) -> None:
        """Called when the application closes and the module widget is destroyed."""
        self.removeObservers()

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
        if not self._parameterNode.inputVolume:
            firstVolumeNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
            if firstVolumeNode:
                self._parameterNode.inputVolume = firstVolumeNode

    def setParameterNode(self, inputParameterNode: Optional[ArterySegmentationParameterNode]) -> None:
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
            self.updateThresholdRange()
            self._checkCanApply()

    def _checkCanApply(self, caller=None, event=None) -> None:
        """
        Check if all the needed inputs are available before enabling "apply" buttons.
        """
        if self._parameterNode and self._parameterNode.inputVolume and \
                self._parameterNode.segmentationSeeds:
            self.ui.segmentationPreviewButton.toolTip = "Click to refresh the preview."
            self.ui.segmentationPreviewButton.enabled = True
        else:
            self.ui.segmentationPreviewButton.toolTip = "Select input volume and seeds."
            self.ui.segmentationPreviewButton.enabled = False
            self.ui.segmentationApplyButton.enabled = False

    def updateThresholdRange(self):
        """
        Update threshold range for segmentation.
        Code source:
        https://github.com/vmtk/SlicerExtension-VMTK/blob/master/LevelSetSegmentation/LevelSetSegmentation.py
        """
        if self._parameterNode is None:
            return

        # If we have a vesselnessNode, we will configure the threshold slider for it instead of the original image
        # if not, the currentNode is the input volume
        thresholdedNode = self.ui.segmentationVesselnessSelector.currentNode()
        if not thresholdedNode or thresholdedNode == 'None':
            thresholdedNode = self.ui.segmentationInputSelector.currentNode()

        if not thresholdedNode or not thresholdedNode.GetImageData():
            wasBlocked = self.ui.segmentationThresholdSlider.blockSignals(True)
            # Reset the thresholdSlider
            self.ui.segmentationThresholdSlider.minimum = 0
            self.ui.segmentationThresholdSlider.maximum = 100
            self.ui.segmentationThresholdSlider.minimumValue = 0
            self.ui.segmentationThresholdSlider.maximumValue = 100
            self.ui.segmentationThresholdSlider.blockSignals(wasBlocked)
            return

        imageData = thresholdedNode.GetImageData()
        scalarRange = imageData.GetScalarRange()
        minimumScalarValue = round(scalarRange[0], 0)
        maximumScalarValue = round(scalarRange[1], 0)

        wasBlocked = self.ui.segmentationThresholdSlider.blockSignals(True)

        self.ui.segmentationThresholdSlider.minimum = minimumScalarValue
        self.ui.segmentationThresholdSlider.maximum = maximumScalarValue

        # If the image has a small scalarRange, we have to adjust the singleStep
        if maximumScalarValue <= 10:
            self.ui.segmentationThresholdSlider.singleStep = 0.01

        if maximumScalarValue > 1:
            self.ui.segmentationThresholdSlider.minimumValue = minimumScalarValue
            self.ui.segmentationThresholdSlider.maximumValue = maximumScalarValue
        else:
            self.ui.segmentationThresholdSlider.minimumValue = minimumScalarValue + (
                    maximumScalarValue - minimumScalarValue) * 0.08
            self.ui.segmentationThresholdSlider.maximumValue = maximumScalarValue

        self.ui.segmentationThresholdSlider.blockSignals(wasBlocked)

    def onThresholdSliderChanged(self):
        """
        Define actions that need to be done when threshold slider changes.
        Code source:
        https://github.com/vmtk/SlicerExtension-VMTK/blob/master/LevelSetSegmentation/LevelSetSegmentation.py
        """
        # First, check if we have a vesselness node
        thresholdedNode = self._parameterNode.vesselnessVolume
        # If we don't have a vesselness node, check if we have an original input node
        if not thresholdedNode:
            thresholdedNode = self._parameterNode.inputVolume

        if thresholdedNode:
            displayNode = thresholdedNode.GetDisplayNode()
            if displayNode:
                slicer.util.setSliceViewerLayers(foregroundOpacity=0.6)
                displayNode.SetLowerThreshold(self.ui.segmentationThresholdSlider.minimumValue)
                displayNode.SetUpperThreshold(self.ui.segmentationThresholdSlider.maximumValue)
                displayNode.SetApplyThreshold(1)

    def resetThresholdOnDisplayNode(self):
        """
        Reset threshold on display node after some computation.
        Code source:
        https://github.com/vmtk/SlicerExtension-VMTK/blob/master/LevelSetSegmentation/LevelSetSegmentation.py
        """
        inputVolumeNode = self._parameterNode.inputVolume
        if inputVolumeNode:
            displayNode = inputVolumeNode.GetDisplayNode()
            if displayNode:
                displayNode.SetApplyThreshold(0)

        vesselnessVolumeNode = self._parameterNode.vesselnessVolume
        if vesselnessVolumeNode:
            displayNode = vesselnessVolumeNode.GetDisplayNode()
            if displayNode:
                displayNode.SetApplyThreshold(0)

    def onInitializationRadio(self):
        """
        Hide stoppers selector when not needed.
        """
        # The level set segmentation initialization method is colliding fronts by default, which does not require
        # stoppers, so stoppers selector is hidden
        if self.ui.collidingFrontsRadio.checked:
            self.ui.segmentationStopperSelector.setVisible(False)
            self.ui.segmentationStopperPlaceWidget.setVisible(False)
            self.ui.segmentationStopperLabel.setVisible(False)
        # When level set method is changed to fast marching, stopper selector is set to visible
        else:
            self.ui.segmentationStopperSelector.setVisible(True)
            self.ui.segmentationStopperPlaceWidget.setVisible(True)
            self.ui.segmentationStopperLabel.setVisible(True)

    def onSegmentationPreviewButton(self):
        """
        Run processing when user clicks preview button in segmentation section.
        """
        self.startSegmentation(True)
        self.ui.segmentationApplyButton.enabled = True

    def onSegmentationApplyButton(self):
        """
        Run processing when user clicks apply button in segmentation section.
        """
        # this is no preview, so we set the preview parameter to false
        self.startSegmentation(preview=False)

    def startSegmentation(self, preview=False):
        """
        Perform actual segmentation.
        Code adapted from:
        https://github.com/vmtk/SlicerExtension-VMTK/blob/master/LevelSetSegmentation/LevelSetSegmentation.py
        """
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            # Import LevelSetSegmentation logic from "SlicerVMTK"
            import LevelSetSegmentation
            lss = LevelSetSegmentation.LevelSetSegmentationLogic()

            # Get required inputs
            inputVolumeNode = self._parameterNode.inputVolume
            seedsNode = self._parameterNode.segmentationSeeds
            stoppersNode = self._parameterNode.segmentationStoppers
            vesselnessNode = self._parameterNode.vesselnessVolume
            labelmapNode = self._parameterNode.segmentationLabelmap
            modelNode = self._parameterNode.segmentationModel

            # Check dimensions
            if inputVolumeNode.GetImageData().GetDimensions() != vesselnessNode.GetImageData().GetDimensions():
                raise ValueError("input and filtered volumes must have the same dimensions")

            # Convert and sort input seeds (sorting is needed when a single branch with segments defined by seeds
            # is segmented)
            seeds = self.convertFiducialHierarchyToVtkIdList(seedsNode, inputVolumeNode)
            seeds.Sort()  # sort the ids in the list in ascending id order
            if seeds.GetNumberOfIds() < 2 and self.ui.collidingFrontsRadio.checked:
                # We need at least two seeds to perform segmentation
                raise ValueError("at least two seeds are required to perform segmentation")

            # Convert stoppers
            if stoppersNode:
                stoppers = self.convertFiducialHierarchyToVtkIdList(stoppersNode, inputVolumeNode)
            else:
                stoppers = vtk.vtkIdList()

            # Create labelmap node to store the binary labelmap if it's not given by the user
            if not labelmapNode:
                newLabelmapNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLLabelMapVolumeNode")
                newLabelmapNode.UnRegister(None)
                newLabelmapNode.CopyOrientation(inputVolumeNode)
                newLabelmapNode.SetName(slicer.mrmlScene.GetUniqueNameByString("SegmentationLabelmap"))
                labelmapNode = slicer.mrmlScene.AddNode(newLabelmapNode)
                labelmapNode.CreateDefaultDisplayNodes()
                self.ui.segmentationLabelmapSelector.setCurrentNode(labelmapNode)

            # Create model node to store the final model if it's not given by the user
            if not modelNode:
                newModelNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLModelNode")
                newModelNode.UnRegister(None)
                newModelNode.SetName(slicer.mrmlScene.GetUniqueNameByString("SegmentationModel"))
                modelNode = slicer.mrmlScene.AddNode(newModelNode)
                modelNode.CreateDefaultDisplayNodes()
                self.ui.segmentationModelSelector.setCurrentNode(modelNode)

            # Create input image for the initialization
            inputImage = vtk.vtkImageData()

            # Check if there's a vesselnessNode, it will be our input for the initialization
            if vesselnessNode:
                inputImage.DeepCopy(vesselnessNode.GetImageData())
            else:
                inputImage.DeepCopy(inputVolumeNode.GetImageData())  # use the original image

            # Create image for initialization
            initImageData = vtk.vtkImageData()
            # Create image for Evolution
            evolImageData = vtk.vtkImageData()

            # Get threshold values
            lowerThreshold = float(self.ui.segmentationThresholdSlider.minimumValue)
            upperThreshold = float(self.ui.segmentationThresholdSlider.maximumValue)

            # Get the selected method
            if self.ui.fastMarchingRadio.checked:
                method = 'fastmarching'
            elif self.ui.collidingFrontsRadio.checked:
                method = 'collidingfronts'

            # Perform the initialization
            initImageData.DeepCopy(lss.performInitialization(inputImage, lowerThreshold, upperThreshold, seeds,
                                                             stoppers, method))

            if not initImageData.GetPointData().GetScalars():
                # Something went wrong, the image is empty
                raise ValueError("segmentation failed, the output was empty.")

            # Check if it is a preview call
            if preview:
                # if this is a preview call, skip the evolution
                evolImageData.DeepCopy(initImageData)
            else:
                # get evolution parameters
                numberOfIterations = float(self._parameterNode.segmentationIterations)
                inflation = float(self._parameterNode.segmentationInflation)
                curvature = float(self._parameterNode.segmentationCurvature)
                attraction = float(self._parameterNode.segmentationAttraction)
                # get the selected evolution method
                if self.ui.geodesicRadio.checked:
                    levelSetsType = 'geodesic'
                elif self.ui.curvesRadio.checked:
                    levelSetsType = 'curves'

                # Perform evolution with input volume (not vesselness node)
                evolImageData.DeepCopy(lss.performEvolution(inputVolumeNode.GetImageData(),
                                                            initImageData,
                                                            int(numberOfIterations),
                                                            inflation,
                                                            curvature,
                                                            attraction,
                                                            levelSetsType))

            # Create segmentation labelmap
            labelMap = vtk.vtkImageData()
            labelMap.DeepCopy(lss.buildSimpleLabelMap(evolImageData, 5, 0))
            labelmapNode.CopyOrientation(inputVolumeNode)
            # Propagate the label map to the node
            labelmapNode.SetAndObserveImageData(labelMap)

            # Deactivate the threshold in the GUI
            self.resetThresholdOnDisplayNode()

            # Overlap the image and the result
            slicer.util.setSliceViewerLayers(background=inputVolumeNode, foreground=vesselnessNode,
                                             label=labelmapNode,
                                             foregroundOpacity=0.6 if preview else 0.1)

            # Generate 3D model using marching cubes
            model = vtk.vtkPolyData()
            # Get the ijkToRas transform for the marching cubes call
            ijkToRasMatrix = vtk.vtkMatrix4x4()
            labelmapNode.GetIJKToRASMatrix(ijkToRasMatrix)
            # Call marching cubes
            model.DeepCopy(lss.marchingCubes(evolImageData, ijkToRasMatrix, 0.0))

            # Propagate model to nodes
            modelNode.SetAndObservePolyData(model)
            modelNode.CreateDefaultDisplayNodes()
            modelDisplayNode = modelNode.GetDisplayNode()

            # Always configure the displayNode to show the model
            modelDisplayNode.SetColor(1.0, 0.55, 0.4)  # red
            modelDisplayNode.SetBackfaceCulling(0)
            modelDisplayNode.SetVisibility2D(0)
            modelDisplayNode.SetVisibility(1)
            modelDisplayNode.SetOpacity(1.0)

            # Fit slice to all sliceviewers
            slicer.app.applicationLogic().FitSliceToAll()

            # Jump all sliceViewers to the first fiducial point, if one was used
            if seedsNode:
                coordinatesRAS = [0, 0, 0]
                if isinstance(seedsNode, slicer.vtkMRMLMarkupsFiducialNode):
                    # Let's get the first children
                    seedsNode.GetNthControlPointPosition(0, coordinatesRAS)
                numberOfSliceNodes = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLSliceNode')
                for n in range(numberOfSliceNodes):
                    sliceNode = slicer.mrmlScene.GetNthNodeByClass(n, "vtkMRMLSliceNode")
                    if sliceNode:
                        sliceNode.JumpSliceByOffsetting(coordinatesRAS[0], coordinatesRAS[1],
                                                        coordinatesRAS[2])

            # Center 3D view(s) on the new model
            if coordinatesRAS:
                for d in range(slicer.app.layoutManager().threeDViewCount):
                    threeDView = slicer.app.layoutManager().threeDWidget(d).threeDView()
                    # Reset the focal point
                    threeDView.resetFocalPoint()
                    # And fly to our seed point
                    interactor = threeDView.interactor()
                    renderer = threeDView.renderWindow().GetRenderers().GetItemAsObject(0)
                    interactor.FlyTo(renderer, coordinatesRAS[0], coordinatesRAS[1], coordinatesRAS[2])

        except Exception as e:
            slicer.util.errorDisplay("Failed to compute results: "+str(e))
            traceback.print_exc()

        qt.QApplication.restoreOverrideCursor()
        slicer.util.showStatusMessage("Segmentation complete", 3000)

    def onChangeModuleButton1(self):
        self.logic.switchModule('Preprocessing')

    def onChangeModuleButton2(self):
        self.logic.switchModule('TortuosityEvaluator')

    def importRequiredModules(self):
        try:
            import LevelSetSegmentation
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

    @staticmethod
    def convertFiducialHierarchyToVtkIdList(hierarchyNode, volumeNode):
        """
        Code source:
        https://github.com/vmtk/SlicerExtension-VMTK/blob/master/LevelSetSegmentation/LevelSetSegmentation.py
        """
        outputIds = vtk.vtkIdList()

        if not hierarchyNode or not volumeNode:
            return outputIds

        if isinstance(hierarchyNode, slicer.vtkMRMLMarkupsFiducialNode) and isinstance(volumeNode,
                                                                                       slicer.vtkMRMLScalarVolumeNode):

            image = volumeNode.GetImageData()

            # now we have the children which are fiducialNodes - let's loop!
            for n in range(hierarchyNode.GetNumberOfControlPoints()):
                coordinatesRAS = [0, 0, 0]

                # grab the current coordinates
                hierarchyNode.GetNthControlPointPosition(n, coordinatesRAS)

                # convert the RAS to IJK
                coordinatesIJK = ArterySegmentationWidget.ConvertRAStoIJK(volumeNode, coordinatesRAS)

                # strip the last element since we need a 3based tuple
                coordinatesIJKlist = (
                    int(coordinatesIJK[0]), int(coordinatesIJK[1]), int(coordinatesIJK[2]))
                outputIds.InsertNextId(int(image.ComputePointId(coordinatesIJKlist)))

        # IdList was created, return it even if it might be empty
        return outputIds

    @staticmethod
    def ConvertRAStoIJK(volumeNode, rasCoordinates):
        """
        Code source:
        https://github.com/vmtk/SlicerExtension-VMTK/blob/master/LevelSetSegmentation/LevelSetSegmentation.py
        """
        rasToIjkMatrix = vtk.vtkMatrix4x4()
        volumeNode.GetRASToIJKMatrix(rasToIjkMatrix)

        # the RAS coordinates need to be 4
        if len(rasCoordinates) < 4:
            rasCoordinates.append(1)

        ijkCoordinates = rasToIjkMatrix.MultiplyPoint(rasCoordinates)

        return ijkCoordinates


#
# ArterySegmentationLogic
#


class ArterySegmentationLogic(ScriptedLoadableModuleLogic):
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
        return ArterySegmentationParameterNode(super().getParameterNode())

    def switchModule(self, moduleName):
        """
        Switch to another module
        """
        pluginHandlerSingleton = slicer.qSlicerSubjectHierarchyPluginHandler.instance()
        pluginHandlerSingleton.pluginByName('Default').switchToModule(moduleName)

