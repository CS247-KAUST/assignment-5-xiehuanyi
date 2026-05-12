cmake_minimum_required (VERSION 3.21)
project (assignment_5B)

# --- 1. C++ STANDARD CONFIGURATION ---
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- 2. AUTOMATIC DEPENDENCY DOWNLOAD (FETCHCONTENT) ---
include(FetchContent)

# --- GLFW ---
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.4
)
FetchContent_MakeAvailable(glfw)

# --- GLM (Math Library) ---
FetchContent_Declare(
    glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(glm)

# --- OpenGL (cross-platform) ---
find_package(OpenGL REQUIRED)

# --- 3. INCLUDE DIRECTORIES ---
include_directories(common common/glad)

# Adding compile units to executable
add_executable(${PROJECT_NAME}
    src/glad.c
    src/CS247_prog.cpp
    src/glslprogram.cpp
    src/vboquad.cpp
)

# Pass the project source dir to the program so it can find shaders/data
# regardless of the build directory.
target_compile_definitions(${PROJECT_NAME} PRIVATE
    PROJECT_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
)

# C++ Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# --- LINK LIBRARIES ---
target_link_libraries(${PROJECT_NAME} PRIVATE glfw OpenGL::GL glm::glm)

# Cocoa frameworks needed by GLFW on macOS are linked automatically by
# the GLFW target. Nothing extra needed here.
