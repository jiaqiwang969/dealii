#----------------------------------------------------------------
# Generated CMake target import file for configuration "DebugRelease".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "deal_II.g" for configuration "DebugRelease"
set_property(TARGET deal_II.g APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUGRELEASE)
set_target_properties(deal_II.g PROPERTIES
  IMPORTED_LOCATION_DEBUGRELEASE "${_IMPORT_PREFIX}/lib/libdeal_II.g.so.10.0.0-pre"
  IMPORTED_SONAME_DEBUGRELEASE "libdeal_II.g.so.10.0.0-pre"
  )

list(APPEND _IMPORT_CHECK_TARGETS deal_II.g )
list(APPEND _IMPORT_CHECK_FILES_FOR_deal_II.g "${_IMPORT_PREFIX}/lib/libdeal_II.g.so.10.0.0-pre" )

# Import target "deal_II" for configuration "DebugRelease"
set_property(TARGET deal_II APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUGRELEASE)
set_target_properties(deal_II PROPERTIES
  IMPORTED_LOCATION_DEBUGRELEASE "${_IMPORT_PREFIX}/lib/libdeal_II.so.10.0.0-pre"
  IMPORTED_SONAME_DEBUGRELEASE "libdeal_II.so.10.0.0-pre"
  )

list(APPEND _IMPORT_CHECK_TARGETS deal_II )
list(APPEND _IMPORT_CHECK_FILES_FOR_deal_II "${_IMPORT_PREFIX}/lib/libdeal_II.so.10.0.0-pre" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
