#ifndef CS247_PROG_H
#define CS247_PROG_H

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// framework includes
#include "glslprogram.h"
#include "vboquad.h"


////////////////
// Structures //
////////////////

// Initial window size; resized to match the dataset upon load.
const unsigned int gWindowWidth = 512;
const unsigned int gWindowHeight = 512;

int current_scalar_field;
int data_size;
bool en_arrow;
bool en_streamline;
bool en_pathline;

int sampling_rate;
float dt;
float path_time_scale;   // pathline time-advance per spatial step (decouples t from xy)

// Integration method: 0 = Euler, 1 = RK2, 2 = RK4
int integ_method;

// Colormap: 0 = off (grayscale), 1 = rainbow, 2 = cool-warm
int colormap_mode;
float blend_factor;

// Glyph appearance
bool length_by_speed;

// Integration parameters
float vec_norm;        // normalization factor: max |v| over all timesteps
float min_vec_mag;     // stop if |v|/vec_norm <= this
float max_acc_length;  // max accumulated path length (grid units)
int   max_steps;

// Rake (bonus)
int  rake_count;
bool rake_horizontal;

// Persistent seeded paths
struct StreamlineRec {
    int seed_x, seed_y;
    std::vector<glm::vec2> pts;   // grid coords
    glm::vec3 color;
};
struct PathlineRec {
    int seed_x, seed_y;
    int seed_t;
    std::vector<glm::vec2> pts;   // grid coords
    glm::vec3 color;
};
std::vector<StreamlineRec> streamlines;
std::vector<PathlineRec>   pathlines;


//////////////////////
//  Global defines  //
//////////////////////
#define TIMER_FREQUENCY_MILLIS  50

//////////////////////
// Global variables //
//////////////////////

static GLFWwindow* window;

char bmModifiers;

int clearColor;

// data handling
char* filenames[ 3 ];
bool grid_data_loaded;
bool scalar_data_loaded;
unsigned short vol_dim[ 3 ];
float* vector_array;
float* scalar_fields;
float* scalar_bounds;

GLuint scalar_field_texture;

int num_scalar_fields;
int num_timesteps;

int loaded_file;
int loaded_timestep;
float timestep;

int view_width, view_height;

GLuint displayList_idx;

int toggle_xy;

////////////////
// Prototypes //
////////////////

void drawGlyphs();
void drawStreamlines();
void drawPathlines();

void seedStreamline(int sx, int sy);
void seedPathline(int sx, int sy);
void recomputeAllStreamlines();

void loadNextTimestep( void );

void LoadData( char* base_filename );
void LoadVectorData( const char* filename );

void DownloadScalarFieldAsTexture( void );
void initGL( void );

void reset_rendering_props( void );

// Quad VBO for the textured background
VBOQuad quad;

// Shader programs: textured colormapped quad + simple line shader for overlays
GLSLProgram vectorProgram;
GLSLProgram lineProgram;
glm::mat4 model;

// Dynamic line VAO/VBO used for glyphs/streamlines/pathlines
GLuint lineVAO = 0;
GLuint lineVBO = 0;


#endif //CS247_PROG_H
