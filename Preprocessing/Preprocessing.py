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
    vtkMRMLMarkupsFiducialNode
)

# python includes
import SimpleITK as sitk
import sitkUtils


#
# Preprocessing
#


class Preprocessing(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("Preprocessing")
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
# PreprocessingParameterNode
#


@parameterNodeWrapper
class PreprocessingParameterNode:
    """
    The parameters needed by module.

    inputVolume - Input volume for preprocessing.
    vesselnessVolume - Preprocessing output volume.
    contrastSeed - Seed placed inside the vessel needed to compute the parameter gamma for Frangi filter.
    minVesselDiameter - Minimum diameter for the vessel (millimeters)
    maxVesselDiameter - Maximum diameter for the vessel (millimeters)
    vesselnessAlpha - Parameter alpha for Frangi filter.
    vesselnessBeta - Parameter beta for Frangi filter.
    filteringSteps - Number of divisions for the interval [minVesselDiameter, maxVesselDiameter].
    """

    inputVolume: vtkMRMLScalarVolumeNode
    vesselnessVolume: vtkMRMLScalarVolumeNode
    contrastSeed: vtkMRMLMarkupsFiducialNode
    minVesselDiameter: Annotated[float, WithinRange(0.00, 25.00)] = 0.60  # millimeters
    maxVesselDiameter: Annotated[float, WithinRange(0.00, 100.00)] = 6.00  # millimeters
    vesselnessAlpha: Annotated[float, WithinRange(0.0, 1.0)] = 0.3
    vesselnessBeta: Annotated[float, WithinRange(0.0, 500.0)] = 500.0
    filteringSteps: Annotated[int, WithinRange(1, 100)] = 10


#
# PreprocessingWidget
#


class PreprocessingWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
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
        uiWidget = slicer.util.loadUI(self.resourcePath("UI/Preprocessing.ui"))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)

        # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
        # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
        # "setMRMLScene(vtkMRMLScene*)" slot.
        uiWidget.setMRMLScene(slicer.mrmlScene)

        # Create logic class. Logic implements all computations that should be possible to run
        # in batch mode, without a graphical user interface.
        self.logic = PreprocessingLogic()

        # Collapse advanced steps
        self.ui.advancedFilteringCollapsibleButton.collapsed = True

        # Connections

        # These connections ensure that we update parameter node when scene is closed
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

        # Buttons
        self.ui.deleteCropButton.connect('clicked(bool)', self.onDeleteCropButton)
        self.ui.cropPreviewButton.connect('clicked(bool)', self.onCropPreviewButton)
        self.ui.filterApplyButton.connect('clicked(bool)', self.onFilterApplyButton)
        self.ui.changeModuleButton.connect('clicked(bool)', self.onChangeModuleButton)

        # Make sure parameter node is initialized (needed for module reload)
        self.initializeParameterNode()

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

    def setParameterNode(self, inputParameterNode: Optional[PreprocessingParameterNode]) -> None:
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
        if self._parameterNode and self._parameterNode.inputVolume and self._parameterNode.vesselnessVolume and \
                self._parameterNode.contrastSeed:
            self.ui.filterApplyButton.toolTip = "Compute output filtered volume."
            self.ui.cropPreviewButton.toolTip = "Create a preview for crop."
            self.ui.deleteCropButton.toolTip = "Delete current crop preview."
            self.ui.cropPreviewButton.enabled = True
            self.ui.deleteCropButton.enabled = True
            self.ui.filterApplyButton.enabled = True
        else:
            self.ui.cropPreviewButton.toolTip = "Select input and output volume nodes and place a contrast seed."
            self.ui.deleteCropButton.toolTip = "Select input and output volume nodes and place a contrast seed."
            self.ui.filterApplyButton.toolTip = "Select input and output volume nodes and place a contrast seed."
            self.ui.cropPreviewButton.enabled = False
            self.ui.deleteCropButton.enabled = False
            self.ui.filterApplyButton.enabled = False

    def importRequiredModules(self):
        try:
            import VesselnessFiltering
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

    def onCropPreviewButton(self):
        """
        Create a ROI for cropping.
        """
        try:
            # Create a small ROI (Region Of interest)
            inputVolume = self._parameterNode.inputVolume
            bounds = [0, 0, 0, 0, 0, 0]
            inputVolume.GetBounds(bounds)

            roiNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLMarkupsROINode")
            roiNode.UnRegister(None)
            roiNode.SetName(slicer.mrmlScene.GetUniqueNameByString("ROI"))

            # Center the ROI in the volume center
            roiCenter = [0, 0, 0]
            for i in range(0, 3):
                roiCenter[i] = (bounds[i * 2 + 1] + bounds[i * 2]) / 2
            roiNode.SetXYZ(roiCenter)
            #Set a ROI radius which includes the heart
            roiRadii = [75, 85, 75]  # radius in mm
            roiNode.SetRadiusXYZ(roiRadii)

            currentROINode = slicer.mrmlScene.AddNode(roiNode)
            currentROINode.CreateDefaultDisplayNodes()

        except Exception as e:
            slicer.util.errorDisplay("Failed to compute results: "+str(e))
            traceback.print_exc()

    def onDeleteCropButton(self):
        """
        Delete crop ROI.
        """
        roiNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLMarkupsROINode")
        if roiNode:
            slicer.mrmlScene.RemoveNode(roiNode)

    def onFilterApplyButton(self):
        """
        Run processing when user clicks apply button in filtering section.
        Code adapted from:
        https://github.com/vmtk/SlicerExtension-VMTK/blob/master/VesselnessFiltering/VesselnessFiltering.py
        """
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            # Apply Frangi filter to enhance vessels
            # Part of the filter logic is implemented in Vesselness Filtering Module from SlicerVMTK
            import VesselnessFiltering
            vfl = VesselnessFiltering.VesselnessFilteringLogic()

            # Input and output volumes
            filterInputVolumeNode = self._parameterNode.inputVolume
            filterOutputVolumeNode = self._parameterNode.vesselnessVolume

            contrastSeedNode = self._parameterNode.contrastSeed
            if contrastSeedNode.GetNumberOfControlPoints() < 1:
                # We need one seed to perform filtering
                raise ValueError("at least one seed is required to perform filtering.")

            costrastSeedRAS = [0, 0, 0]
            contrastSeedNode.GetNthControlPointPosition(0, costrastSeedRAS)

            # Crop volume based on ROI
            roiNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLMarkupsROINode")
            if roiNode:
                slicer.modules.cropvolume.logic().CropVoxelBased(roiNode, filterInputVolumeNode, filterInputVolumeNode)
                slicer.mrmlScene.RemoveNode(roiNode)

            # Get filter parameters
            minimumDiameterMm = float(self._parameterNode.minVesselDiameter)
            maximumDiameterMm = float(self._parameterNode.maxVesselDiameter)
            alpha = float(self._parameterNode.vesselnessAlpha)
            beta = float(self._parameterNode.vesselnessBeta)
            discretizationSteps = float(self._parameterNode.filteringSteps)

            # Perform top Hat Filtering
            TopHatFilteredVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "TopHatFilteredVolume")
            kernelRadiusX = maximumDiameterMm / 2
            kernelRadiusY = maximumDiameterMm / 2
            kernelRadiusZ = 0
            kernelType = 1  # ball
            self.logic.computeTopHatTransform(filterInputVolumeNode, TopHatFilteredVolume, kernelType, kernelRadiusX,
                                              kernelRadiusY, kernelRadiusZ)
            # Get step method
            if self.ui.equispacedRadio.checked:
                stepMethod = 0
            elif self.ui.logarithmicRadio.checked:
                stepMethod = 1

            # If the next comment is removed, filtering will be performed without top hat transform
            # TopHatFilteredVolume = filterInputVolumeNode
            ijk = vfl.getIJKFromRAS(TopHatFilteredVolume, costrastSeedRAS)
            image = TopHatFilteredVolume.GetImageData()
            diameter = int(maximumDiameterMm)
            seedValue = image.GetScalarComponentAsFloat(ijk[0], ijk[1], ijk[2], 0)
            outsideValues = [seedValue - image.GetScalarComponentAsFloat(ijk[0] + (2 * diameter), ijk[1], ijk[2], 0),  # right
                             seedValue - image.GetScalarComponentAsFloat(ijk[0] - (2 * diameter), ijk[1], ijk[2], 0),  # left
                             seedValue - image.GetScalarComponentAsFloat(ijk[0], ijk[1] + (2 * diameter), ijk[2], 0),  # top
                             seedValue - image.GetScalarComponentAsFloat(ijk[0], ijk[1] - (2 * diameter), ijk[2], 0),  # bottom
                             seedValue - image.GetScalarComponentAsFloat(ijk[0], ijk[1], ijk[2] + (2 * diameter), 0),  # front
                             seedValue - image.GetScalarComponentAsFloat(ijk[0], ijk[1], ijk[2] - (2 * diameter), 0)]  # back
            differenceValue = max(outsideValues)
            contrastMeasure = 0.25*differenceValue

            self.logic.computeVesselnessVolume(TopHatFilteredVolume, filterOutputVolumeNode, minimumDiameterMm,
                                               maximumDiameterMm, alpha, beta, contrastMeasure,
                                               int(discretizationSteps), stepMethod)

            slicer.util.setSliceViewerLayers(background=filterInputVolumeNode, foreground=filterOutputVolumeNode,
                                             foregroundOpacity=0.6)

        except Exception as e:
            slicer.util.errorDisplay("Failed to compute results: "+str(e))
            traceback.print_exc()

        qt.QApplication.restoreOverrideCursor()
        slicer.util.showStatusMessage("Filtering complete", 3000)

    def onChangeModuleButton(self):
        self.logic.switchModule('ArterySegmentation')


#
# PreprocessingLogic
#


class PreprocessingLogic(ScriptedLoadableModuleLogic):
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
        return PreprocessingParameterNode(super().getParameterNode())

    def computeTopHatTransform(self, currentVolumeNode, currentOutputVolumeNode, kernelType, radiusMmX,
                               radiusMmY, radiusMmZ):
        """
        Compute White Top-Hat transform.
        """
        inputImage = sitkUtils.PullVolumeFromSlicer(currentVolumeNode)
        filter = sitk.WhiteTopHatImageFilter()  # create the filter
        filter.SetKernelType(kernelType)
        filter.SetSafeBorder(True)

        # Filter kernel size must be specified in pixels, so whe have to convert from mm to pixels
        inputSpacing = inputImage.GetSpacing()
        radiusPixelX = int(round(radiusMmX / inputSpacing[0], 0))
        radiusPixelY = int(round(radiusMmY / inputSpacing[1], 0))
        radiusPixelZ = int(round(radiusMmZ / inputSpacing[2], 0))
        filter.SetKernelRadius([radiusPixelX, radiusPixelY, radiusPixelZ])
        outputImage = filter.Execute(inputImage)
        sitkUtils.PushVolumeToSlicer(outputImage, currentOutputVolumeNode)

    def computeVesselnessVolume(self, currentVolumeNode, currentOutputVolumeNode, minimumDiameterMm=0.5,
                                maximumDiameterMm=5, alpha=0.3, beta=500, contrastMeasure=100, discretizationSteps=10,
                                stepMethod=0):
        """
        Perform Frangi vesselness filtering.
        Code source:
        https://github.com/vmtk/SlicerExtension-VMTK/blob/master/VesselnessFiltering/VesselnessFiltering.py
        """
        # This image will later hold the inputImage
        inImage = vtk.vtkImageData()

        # Clone the inputImage
        inImage.DeepCopy(currentVolumeNode.GetImageData())
        currentOutputVolumeNode.CopyOrientation(currentVolumeNode)

        # Temporarily set spacing to allow vesselness computation performed in physical space
        inImage.SetSpacing(currentVolumeNode.GetSpacing())

        # We now compute the vesselness in RAS space, inImage has spacing and origin attached,
        # the diameters are converted to mm, we use RAS space to support anisotropic datasets
        cast = vtk.vtkImageCast()
        cast.SetInputData(inImage)
        cast.SetOutputScalarTypeToFloat()
        cast.Update()
        inImage = cast.GetOutput()

        import vtkvmtkSegmentationPython as vtkvmtkSegmentation

        v = vtkvmtkSegmentation.vtkvmtkVesselnessMeasureImageFilter()
        v.SetInputData(inImage)
        v.SetSigmaMin(minimumDiameterMm)
        v.SetSigmaMax(maximumDiameterMm)
        v.SetNumberOfSigmaSteps(discretizationSteps)
        v.SetAlpha(alpha)
        v.SetBeta(beta)
        v.SetGamma(contrastMeasure)
        v.SetSigmaStepMethod(stepMethod)
        v.Update()

        outImage = vtk.vtkImageData()
        outImage.DeepCopy(v.GetOutput())
        outImage.GetPointData().GetScalars().Modified()

        # Restore Slicer-compliant image spacing
        outImage.SetSpacing(1, 1, 1)

        # We set the outImage which has spacing 1,1,1. The ijkToRas matrix of the node will take care of that
        currentOutputVolumeNode.SetAndObserveImageData(outImage)

        # Save which volume node vesselness filtering result was saved to
        currentVolumeNode.SetAndObserveNodeReferenceID("VesselnessVolume", currentOutputVolumeNode.GetID())

    def switchModule(self, moduleName):
        """
        Switch to another module
        """
        pluginHandlerSingleton = slicer.qSlicerSubjectHierarchyPluginHandler.instance()
        pluginHandlerSingleton.pluginByName('Default').switchToModule(moduleName)
