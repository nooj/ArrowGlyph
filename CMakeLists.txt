PROJECT("ArrowGlyph")

CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
IF(NOT PARAVIEW_LOADED)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
  SET(PARAVIEW_LOADED 1)
ELSE(NOT PARAVIEW_LOADED)
  IF (PARAVIEW_BUILD_QT_GUI)
    INCLUDE(${QT_USE_FILE})
  ENDIF (PARAVIEW_BUILD_QT_GUI)
  INCLUDE_DIRECTORIES(
    ${PARAVIEW_INCLUDE_DIRS}
    ${PARAVIEW_GUI_INCLUDE_DIRS}
    ${PARAVIEW_KWSYS_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
  )
ENDIF(NOT PARAVIEW_LOADED)

#--------------------------------------------------
# Set project include directories 
#--------------------------------------------------
INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

#--------------------------------------------------
# Define the source files wrapped by client server
#--------------------------------------------------
SET (ArrowGlyph_WRAPPED_SRCS 
  ${ArrowGlyph_SOURCE_DIR}/vtkArrowGlyphFilter.cxx
)

#--------------------------------------------------
# Define other source files used
#--------------------------------------------------
SET (ArrowGlyph_SRCS 
)

#--------------------------------------------------
# Define the filter list for the GUI
#--------------------------------------------------
SET (ArrowGlyph_CLIENT_XML 
  ${ArrowGlyph_SOURCE_DIR}/ArrowGlyph_Filters.xml
)

#--------------------------------------------------
# Define Plugin
#--------------------------------------------------
SET(LIBRARY_TARGET_NAME ArrowGlyph_plugin)
ADD_PARAVIEW_PLUGIN(
  ${LIBRARY_TARGET_NAME}
  "1.0" 

  SERVER_MANAGER_XML 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkArrowGlyphFilter.xml
  SERVER_MANAGER_SOURCES
    ${ArrowGlyph_WRAPPED_SRCS}
  SERVER_SOURCES
    ${ArrowGlyph_SRCS}  
  GUI_RESOURCE_FILES
    ${ArrowGlyph_CLIENT_XML}
)

TARGET_LINK_LIBRARIES(ArrowGlyph_plugin vtkGraphics)

#--------------------------------------------------------
# Create the configuration header.
#--------------------------------------------------------
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/vtkLocalConfigure.h.in
               ${PROJECT_BINARY_DIR}/vtkArrowGlyphConfigure.h
               @ONLY IMMEDIATE)


#--------------------------------------------------
# Install
#--------------------------------------------------
SET(INSTALL_PATH 
  "${CMAKE_INSTALL_PREFIX}/lib/paraview-${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}"
)

INSTALL(
  TARGETS 
    ArrowGlyph_plugin
  DESTINATION ${INSTALL_PATH}
)

