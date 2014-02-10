cmake_minimum_required(VERSION 2.8.9)

#-----------------------------------------------------------------------------
set(MODULE_NAME CatheterMatchRegistration) # Do not use 'project()'

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

##-----------------------------------------------------------------------------
#set(MODULE_INCLUDE_DIRECTORIES
#  )
#
#set(MODULE_SRCS
#  CatheterMatchRegistration.cxx
#  )
#
#set(MODULE_TARGET_LIBRARIES
#  ${ITK_LIBRARIES}
#  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  ADDITIONAL_SRCS
    itkCatheterMatchRigidTransformGradientDescentOptimizer.cxx
  TARGET_LIBRARIES ModuleDescriptionParser ${ITK_LIBRARIES}
  LINK_DIRECTORIES ${ModuleDescriptionParser_BINARY_DIR}
  INCLUDE_DIRECTORIES
    ${SlicerBaseCLI_SOURCE_DIR} ${SlicerBaseCLI_BINARY_DIR}
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
#  add_subdirectory(Testing)
endif()


