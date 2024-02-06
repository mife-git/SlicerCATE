cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(SlicerCATE)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/SlicerCATE")
set(EXTENSION_CATEGORY "SlicerCATE")
set(EXTENSION_CONTRIBUTORS "Michela Ferrari (Fondazione IRCCS Policlinico San Matteo)")
set(EXTENSION_DESCRIPTION "Coronary Artery Tortuosity Evaluator Extension")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/SlicerCATE.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/SlicerCATE/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CoronaryArteryTortuosityEvaluator)
add_subdirectory(TortuosityEvaluator)
add_subdirectory(ArterySegmentation)
add_subdirectory(Preprocessing)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
