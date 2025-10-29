if(NOT DEFINED LIDAR_EQUALIZER_EXE)
    message(FATAL_ERROR "LIDAR_EQUALIZER_EXE not defined")
endif()
if(NOT DEFINED INPUT_FILE)
    message(FATAL_ERROR "INPUT_FILE not defined")
endif()
if(NOT DEFINED OUTPUT_FILE)
    message(FATAL_ERROR "OUTPUT_FILE not defined")
endif()

get_filename_component(_out_dir "${OUTPUT_FILE}" DIRECTORY)
file(MAKE_DIRECTORY "${_out_dir}")

set(_cell_size "1.0")
if(DEFINED CELL_SIZE)
    set(_cell_size "${CELL_SIZE}")
endif()

set(_target_density "6.0")
if(DEFINED TARGET_DENSITY)
    set(_target_density "${TARGET_DENSITY}")
endif()

set(_seed "42")
if(DEFINED SEED_VALUE)
    set(_seed "${SEED_VALUE}")
endif()

set(_cmd
    "${LIDAR_EQUALIZER_EXE}"
    "${INPUT_FILE}"
    "${OUTPUT_FILE}"
    "${_cell_size}"
    "${_target_density}"
    "${_seed}"
)
if(DEFINED CLASS_SCOPE)
    list(APPEND _cmd "--class-scope=${CLASS_SCOPE}")
endif()
if(DEFINED FLAG_OVERLAP_ONLY)
    string(TOLOWER "${FLAG_OVERLAP_ONLY}" _flag_overlap_only_l)
    if(_flag_overlap_only_l MATCHES "^(0|off|false|no)$")
        list(APPEND _cmd "--flag-overlap-only=false")
    else()
        list(APPEND _cmd "--flag-overlap-only")
    endif()
endif()
if(DEFINED OVERWRITE_OVERLAP)
    string(TOLOWER "${OVERWRITE_OVERLAP}" _overwrite_overlap_l)
    if(_overwrite_overlap_l MATCHES "^(1|on|true|yes)$")
        list(APPEND _cmd "--overwrite-overlap")
    else()
        list(APPEND _cmd "--overwrite-overlap=false")
    endif()
endif()
if(DEFINED MIN_POINTS_PER_PSID)
    list(APPEND _cmd "--min-points-per-psid=${MIN_POINTS_PER_PSID}")
endif()
if(DEFINED OVERLAP_DILATE)
    list(APPEND _cmd "--overlap-dilate=${OVERLAP_DILATE}")
endif()
if(DEFINED SWATH_KEY)
    list(APPEND _cmd "--swath-key=${SWATH_KEY}")
endif()

execute_process(
    COMMAND ${_cmd}
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE _stderr
    RESULT_VARIABLE _result
)

if(NOT _result EQUAL 0)
    message(FATAL_ERROR "Equalizer command failed with code ${_result}\nstdout:\n${_stdout}\nstderr:\n${_stderr}")
endif()

string(REGEX MATCH "Input points[^\n]*" _input_line "${_stdout}")
if(NOT _input_line)
    file(REMOVE "${OUTPUT_FILE}")
    message(FATAL_ERROR "Did not find \"Input points\" line in output.\nstdout:\n${_stdout}")
endif()

string(TOLOWER "${FLAG_OVERLAP_ONLY}" _flag_overlap_only_l)
set(_overlap_only FALSE)
if(DEFINED FLAG_OVERLAP_ONLY AND NOT _flag_overlap_only_l MATCHES "^(0|off|false|no)$")
    set(_overlap_only TRUE)
endif()

if(_overlap_only)
    string(REGEX MATCH "Points flagged[^\n]*" _flagged_line "${_stdout}")
    if(NOT _flagged_line)
        file(REMOVE "${OUTPUT_FILE}")
        message(FATAL_ERROR "Did not find \"Points flagged\" line in output.\nstdout:\n${_stdout}")
    endif()
else()
    string(REGEX MATCH "Kept points[^\n]*" _kept_line "${_stdout}")
    if(NOT _kept_line)
        file(REMOVE "${OUTPUT_FILE}")
        message(FATAL_ERROR "Did not find \"Kept points\" line in output.\nstdout:\n${_stdout}")
    endif()

    string(REGEX REPLACE ".*Kept points[ ]*:[ ]*([0-9]+).*" "\\1" _kept_count "${_kept_line}")
    string(REGEX REPLACE ".*Input points[ ]*:[ ]*([0-9]+).*" "\\1" _input_count "${_input_line}")

    if("${_input_count}" STREQUAL "${_input_line}" OR "${_input_count}" STREQUAL "")
        file(REMOVE "${OUTPUT_FILE}")
        message(FATAL_ERROR "Failed to parse input point count from \"${_input_line}\"")
    endif()

    if("${_kept_count}" STREQUAL "${_kept_line}" OR "${_kept_count}" STREQUAL "")
        file(REMOVE "${OUTPUT_FILE}")
        message(FATAL_ERROR "Failed to parse kept point count from \"${_kept_line}\"")
    endif()

    math(EXPR _input_int "${_input_count}")
    math(EXPR _kept_int "${_kept_count}")

    if(_input_int LESS_EQUAL 0)
        file(REMOVE "${OUTPUT_FILE}")
        message(FATAL_ERROR "Input point count not positive (value=${_input_int})")
    endif()

    if(_kept_int LESS_EQUAL 0)
        file(REMOVE "${OUTPUT_FILE}")
        message(FATAL_ERROR "Kept point count not positive (value=${_kept_int})")
    endif()

    if(_kept_int GREATER_EQUAL _input_int)
        file(REMOVE "${OUTPUT_FILE}")
        message(FATAL_ERROR "Kept point count (${_kept_int}) is not lower than input point count (${_input_int})")
    endif()
endif()

if(NOT EXISTS "${OUTPUT_FILE}")
    message(FATAL_ERROR "Expected output file \"${OUTPUT_FILE}\" not created")
endif()

file(REMOVE "${OUTPUT_FILE}")
