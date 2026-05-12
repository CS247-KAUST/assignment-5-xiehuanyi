// CS 247 - Scientific Visualization, KAUST
//
// Programming Assignment #5
#include <cstring>
#include "CS247_prog.h"

#ifndef PROJECT_SOURCE_DIR_STR
#define PROJECT_SOURCE_DIR_STR "."
#endif

static std::string projectPath(const char* rel) {
    return std::string(PROJECT_SOURCE_DIR_STR) + "/" + rel;
}

// cycle clear colors
static void nextClearColor()
{
    clearColor = (clearColor + 1) % 3;
    switch(clearColor)
    {
        case 0:  glClearColor(0.0f, 0.0f, 0.0f, 1.0f); break;
        case 1:  glClearColor(0.2f, 0.2f, 0.3f, 1.0f); break;
        default: glClearColor(0.7f, 0.7f, 0.7f, 1.0f); break;
    }
}

// ---------------------------------------------------------------------------
// Vector field sampling
// ---------------------------------------------------------------------------

// Read raw vector at integer grid (ix,iy) at time slice t.
static glm::vec2 vecAt(int ix, int iy, int t)
{
    int W = vol_dim[0];
    int k = 3 * (ix + iy * W) + 3 * t * data_size;
    return glm::vec2(vector_array[k], vector_array[k + 1]);
}

// Bilinear interpolation in space at fixed time slice t.
static glm::vec2 sampleBilinear(float fx, float fy, int t)
{
    int W = vol_dim[0], H = vol_dim[1];
    if (fx < 0 || fx > W - 1 || fy < 0 || fy > H - 1) return glm::vec2(0.0f);
    int x0 = (int)std::floor(fx), y0 = (int)std::floor(fy);
    int x1 = std::min(x0 + 1, W - 1);
    int y1 = std::min(y0 + 1, H - 1);
    float ax = fx - x0, ay = fy - y0;
    glm::vec2 v00 = vecAt(x0, y0, t);
    glm::vec2 v10 = vecAt(x1, y0, t);
    glm::vec2 v01 = vecAt(x0, y1, t);
    glm::vec2 v11 = vecAt(x1, y1, t);
    glm::vec2 v0 = glm::mix(v00, v10, ax);
    glm::vec2 v1 = glm::mix(v01, v11, ax);
    return glm::mix(v0, v1, ay);
}

// Trilinear interpolation: bilinear in space + linear in time.
static glm::vec2 sampleTrilinear(float fx, float fy, float ft)
{
    if (num_timesteps <= 1) return sampleBilinear(fx, fy, 0);
    if (ft < 0) ft = 0;
    if (ft > num_timesteps - 1) ft = (float)(num_timesteps - 1);
    int t0 = (int)std::floor(ft);
    int t1 = std::min(t0 + 1, num_timesteps - 1);
    float at = ft - t0;
    glm::vec2 v0 = sampleBilinear(fx, fy, t0);
    glm::vec2 v1 = sampleBilinear(fx, fy, t1);
    return glm::mix(v0, v1, at);
}

// Normalized vector (divided by global max magnitude across timesteps).
static glm::vec2 nVecBilinear(float fx, float fy, int t)  { return sampleBilinear(fx, fy, t) / vec_norm; }
static glm::vec2 nVecTrilinear(float fx, float fy, float ft) { return sampleTrilinear(fx, fy, ft) / vec_norm; }

// ---------------------------------------------------------------------------
// Coordinate helpers
// ---------------------------------------------------------------------------

static glm::vec2 gridToNDC(glm::vec2 g)
{
    float W = (float)(vol_dim[0] - 1);
    float H = (float)(vol_dim[1] - 1);
    return glm::vec2(2.0f * g.x / W - 1.0f, 2.0f * g.y / H - 1.0f);
}

static bool insideGrid(glm::vec2 p)
{
    return p.x >= 0 && p.x <= vol_dim[0] - 1 && p.y >= 0 && p.y <= vol_dim[1] - 1;
}

// Distinct color per seed (golden-ratio hue sequence).
static glm::vec3 nextSeedColor()
{
    static int n = 0;
    float h = std::fmod(n++ * 0.61803398875f + 0.13f, 1.0f);
    int   i = (int)(h * 6.0f);
    float f = h * 6.0f - i;
    float q = 1.0f - f;
    switch (i % 6) {
        case 0: return glm::vec3(1, f, 0);
        case 1: return glm::vec3(q, 1, 0);
        case 2: return glm::vec3(0, 1, f);
        case 3: return glm::vec3(0, q, 1);
        case 4: return glm::vec3(f, 0, 1);
        default: return glm::vec3(1, 0, q);
    }
}

// ---------------------------------------------------------------------------
// Streamline & Pathline integration (Euler / RK2 / RK4)
// ---------------------------------------------------------------------------

// One integration step for streamlines (fixed time slice t).
static glm::vec2 stepStream(glm::vec2 p, float dts, int t, int method)
{
    glm::vec2 k1 = nVecBilinear(p.x, p.y, t);
    if (method == 0) return p + dts * k1;
    if (method == 1) {
        glm::vec2 mid = p + 0.5f * dts * k1;
        glm::vec2 k2 = nVecBilinear(mid.x, mid.y, t);
        return p + dts * k2;
    }
    glm::vec2 a = p + 0.5f * dts * k1;
    glm::vec2 k2 = nVecBilinear(a.x, a.y, t);
    glm::vec2 b = p + 0.5f * dts * k2;
    glm::vec2 k3 = nVecBilinear(b.x, b.y, t);
    glm::vec2 c = p + dts * k3;
    glm::vec2 k4 = nVecBilinear(c.x, c.y, t);
    return p + (dts / 6.0f) * (k1 + 2.0f * k2 + 2.0f * k3 + k4);
}

// Pathline state and step (advances both space and time).
struct PathState { glm::vec2 p; float t; };

static PathState stepPath(PathState s, float dts, int method)
{
    glm::vec2 k1 = nVecTrilinear(s.p.x, s.p.y, s.t);
    float dtt = dts * path_time_scale;
    if (method == 0) {
        return { s.p + dts * k1, s.t + dtt };
    }
    if (method == 1) {
        glm::vec2 pm = s.p + 0.5f * dts * k1;
        float     tm = s.t + 0.5f * dtt;
        glm::vec2 k2 = nVecTrilinear(pm.x, pm.y, tm);
        return { s.p + dts * k2, s.t + dtt };
    }
    glm::vec2 pa = s.p + 0.5f * dts * k1;
    float     ta = s.t + 0.5f * dtt;
    glm::vec2 k2 = nVecTrilinear(pa.x, pa.y, ta);
    glm::vec2 pb = s.p + 0.5f * dts * k2;
    glm::vec2 k3 = nVecTrilinear(pb.x, pb.y, ta);
    glm::vec2 pc = s.p + dts * k3;
    glm::vec2 k4 = nVecTrilinear(pc.x, pc.y, s.t + dtt);
    return { s.p + (dts / 6.0f) * (k1 + 2.0f * k2 + 2.0f * k3 + k4), s.t + dtt };
}

// Integrate a streamline starting at (sx,sy) in the given direction.
static std::vector<glm::vec2> integrateStreamline(float sx, float sy, int t,
                                                  float dt_signed, int method)
{
    std::vector<glm::vec2> out;
    glm::vec2 p(sx, sy);
    out.push_back(p);
    float acc = 0.0f;
    for (int s = 0; s < max_steps; ++s) {
        glm::vec2 v = nVecBilinear(p.x, p.y, t);
        if (glm::length(v) < min_vec_mag) break;
        glm::vec2 pn = stepStream(p, dt_signed, t, method);
        if (!insideGrid(pn)) break;
        float seg = glm::length(pn - p);
        acc += seg;
        if (acc > max_acc_length) break;
        p = pn;
        out.push_back(p);
    }
    return out;
}

// Integrate a pathline starting at (sx,sy) at time st.
static std::vector<glm::vec2> integratePathline(float sx, float sy, float st,
                                                float dt_signed, int method)
{
    std::vector<glm::vec2> out;
    PathState s = { glm::vec2(sx, sy), st };
    out.push_back(s.p);
    float acc = 0.0f;
    for (int i = 0; i < max_steps; ++i) {
        glm::vec2 v = nVecTrilinear(s.p.x, s.p.y, s.t);
        if (glm::length(v) < min_vec_mag) break;
        PathState sn = stepPath(s, dt_signed, method);
        if (!insideGrid(sn.p)) break;
        if (sn.t < 0.0f || sn.t > num_timesteps - 1) break;
        float seg = glm::length(sn.p - s.p);
        acc += seg;
        if (acc > max_acc_length) break;
        s = sn;
        out.push_back(s.p);
    }
    return out;
}

// Combine forward + backward integration into one polyline.
static std::vector<glm::vec2> joinFwdBwd(const std::vector<glm::vec2>& fwd,
                                         const std::vector<glm::vec2>& bwd)
{
    std::vector<glm::vec2> joined;
    joined.reserve(fwd.size() + bwd.size());
    for (auto it = bwd.rbegin(); it != bwd.rend(); ++it) joined.push_back(*it);
    for (size_t i = 1; i < fwd.size(); ++i) joined.push_back(fwd[i]);
    return joined;
}

void seedStreamline(int sx, int sy)
{
    if (!grid_data_loaded) return;
    auto fwd = integrateStreamline((float)sx, (float)sy, loaded_timestep,  dt, integ_method);
    auto bwd = integrateStreamline((float)sx, (float)sy, loaded_timestep, -dt, integ_method);
    StreamlineRec r;
    r.seed_x = sx; r.seed_y = sy;
    r.pts = joinFwdBwd(fwd, bwd);
    r.color = nextSeedColor();
    streamlines.push_back(std::move(r));
    fprintf(stderr, "Streamline seeded at (%d,%d). %zu total.\n", sx, sy, streamlines.size());
}

void seedPathline(int sx, int sy)
{
    if (!grid_data_loaded) return;
    auto fwd = integratePathline((float)sx, (float)sy, (float)loaded_timestep,  dt, integ_method);
    auto bwd = integratePathline((float)sx, (float)sy, (float)loaded_timestep, -dt, integ_method);
    PathlineRec r;
    r.seed_x = sx; r.seed_y = sy; r.seed_t = loaded_timestep;
    r.pts = joinFwdBwd(fwd, bwd);
    r.color = nextSeedColor();
    pathlines.push_back(std::move(r));
    fprintf(stderr, "Pathline seeded at (%d,%d) t=%d. %zu total.\n", sx, sy, loaded_timestep, pathlines.size());
}

// Recompute every streamline against the currently loaded time slice.
void recomputeAllStreamlines()
{
    for (auto& s : streamlines) {
        auto fwd = integrateStreamline((float)s.seed_x, (float)s.seed_y, loaded_timestep,  dt, integ_method);
        auto bwd = integrateStreamline((float)s.seed_x, (float)s.seed_y, loaded_timestep, -dt, integ_method);
        s.pts = joinFwdBwd(fwd, bwd);
    }
}

// ---------------------------------------------------------------------------
// Drawing helpers
// ---------------------------------------------------------------------------

static void uploadAndDrawLines(const std::vector<glm::vec2>& ndcVerts, GLenum mode,
                               glm::vec4 color)
{
    if (ndcVerts.size() < 2) return;
    glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 ndcVerts.size() * sizeof(glm::vec2),
                 ndcVerts.data(), GL_DYNAMIC_DRAW);
    lineProgram.setUniform("lineColor", color);
    glBindVertexArray(lineVAO);
    glDrawArrays(mode, 0, (GLsizei)ndcVerts.size());
    glBindVertexArray(0);
}

void drawGlyphs()
{
    if (!grid_data_loaded || !en_arrow) return;
    int W = vol_dim[0], H = vol_dim[1];
    int stride = std::max(sampling_rate, 2);
    float baseLen = stride * 0.5f;       // grid units
    float headFrac = 0.30f;
    float headAngle = 0.45f;
    float c = std::cos(headAngle);
    float s = std::sin(headAngle);

    std::vector<glm::vec2> verts;
    verts.reserve((size_t)((W / stride) * (H / stride)) * 6);

    for (int iy = stride / 2; iy < H; iy += stride) {
        for (int ix = stride / 2; ix < W; ix += stride) {
            glm::vec2 v = vecAt(ix, iy, loaded_timestep) / vec_norm;
            float mag = glm::length(v);
            if (mag < 1e-5f) continue;
            glm::vec2 dir = v / mag;
            float L = length_by_speed ? baseLen * std::min(mag, 1.0f) : baseLen;
            if (L < 1e-3f) continue;
            glm::vec2 p0((float)ix, (float)iy);
            glm::vec2 p1 = p0 + dir * L;

            // Arrowhead segments: rotate -dir by +/- headAngle.
            glm::vec2 nd = -dir;
            glm::vec2 left  (c * nd.x - s * nd.y,  s * nd.x + c * nd.y);
            glm::vec2 right (c * nd.x + s * nd.y, -s * nd.x + c * nd.y);
            float headLen = L * headFrac;

            verts.push_back(gridToNDC(p0));
            verts.push_back(gridToNDC(p1));
            verts.push_back(gridToNDC(p1));
            verts.push_back(gridToNDC(p1 + left  * headLen));
            verts.push_back(gridToNDC(p1));
            verts.push_back(gridToNDC(p1 + right * headLen));
        }
    }
    if (verts.empty()) return;
    lineProgram.use();
    uploadAndDrawLines(verts, GL_LINES, glm::vec4(1.0f, 1.0f, 1.0f, 1.0f));
}

void drawStreamlines()
{
    if (!en_streamline || streamlines.empty()) return;
    lineProgram.use();
    for (auto& s : streamlines) {
        if (s.pts.size() < 2) continue;
        std::vector<glm::vec2> ndc;
        ndc.reserve(s.pts.size());
        for (auto& p : s.pts) ndc.push_back(gridToNDC(p));
        uploadAndDrawLines(ndc, GL_LINE_STRIP, glm::vec4(s.color, 1.0f));
    }
}

void drawPathlines()
{
    if (!en_pathline || pathlines.empty()) return;
    lineProgram.use();
    for (auto& p : pathlines) {
        if (p.pts.size() < 2) continue;
        std::vector<glm::vec2> ndc;
        ndc.reserve(p.pts.size());
        for (auto& q : p.pts) ndc.push_back(gridToNDC(q));
        uploadAndDrawLines(ndc, GL_LINE_STRIP, glm::vec4(p.color, 1.0f));
    }
}

// ---------------------------------------------------------------------------
// Callbacks
// ---------------------------------------------------------------------------

void frameBufferCallback(GLFWwindow* window, int width, int height)
{
    view_width = width;
    view_height = height;
    glViewport(0, 0, width, height);
}

static const char* INTEG_NAMES[] = { "Euler", "RK2", "RK4" };
static const char* CMAP_NAMES[]  = { "off", "rainbow", "cool-warm" };

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_RELEASE) return;

    switch (key) {
        case '1':
            toggle_xy = 0;
            LoadData(filenames[0]); loaded_file = 0;
            fprintf(stderr, "Loading %s dataset.\n", filenames[0]);
            streamlines.clear(); pathlines.clear();
            break;
        case '2':
            toggle_xy = 0;
            LoadData(filenames[1]); loaded_file = 1;
            fprintf(stderr, "Loading %s dataset.\n", filenames[1]);
            streamlines.clear(); pathlines.clear();
            break;
        case '3':
            toggle_xy = 1;
            LoadData(filenames[2]); loaded_file = 2;
            fprintf(stderr, "Loading %s dataset.\n", filenames[2]);
            streamlines.clear(); pathlines.clear();
            break;
        case '0':
            if (num_timesteps > 1) {
                loadNextTimestep();
                recomputeAllStreamlines();
                fprintf(stderr, "Timestep %d / %d.\n", loaded_timestep, num_timesteps);
            }
            break;
        case GLFW_KEY_A:
            en_arrow = !en_arrow;
            fprintf(stderr, "%s drawing arrows.\n", en_arrow ? "enabling" : "disabling");
            break;
        case GLFW_KEY_S:
            current_scalar_field = (current_scalar_field + 1) % std::max(num_scalar_fields, 1);
            DownloadScalarFieldAsTexture();
            fprintf(stderr, "Scalar field %d.\n", current_scalar_field);
            break;
        case GLFW_KEY_B:
            nextClearColor();
            break;
        case GLFW_KEY_EQUAL:
            sampling_rate = std::min(sampling_rate + 5, 100);
            fprintf(stderr, "Sampling rate -> %d.\n", sampling_rate);
            break;
        case GLFW_KEY_MINUS:
            sampling_rate = std::max(sampling_rate - 5, 5);
            fprintf(stderr, "Sampling rate -> %d.\n", sampling_rate);
            break;
        case GLFW_KEY_I:
            dt = std::min(dt + 0.05f, 5.0f);
            fprintf(stderr, "dt = %.3f\n", dt);
            recomputeAllStreamlines();
            break;
        case GLFW_KEY_K:
            dt = std::max(dt - 0.05f, 0.005f);
            fprintf(stderr, "dt = %.3f\n", dt);
            recomputeAllStreamlines();
            break;
        case GLFW_KEY_T:
            en_streamline = !en_streamline;
            fprintf(stderr, "%s drawing streamlines.\n", en_streamline ? "enabling" : "disabling");
            break;
        case GLFW_KEY_P:
            en_pathline = !en_pathline;
            fprintf(stderr, "%s drawing pathlines.\n", en_pathline ? "enabling" : "disabling");
            break;
        case GLFW_KEY_M:
            integ_method = (integ_method + 1) % 3;
            fprintf(stderr, "Integration: %s.\n", INTEG_NAMES[integ_method]);
            recomputeAllStreamlines();
            break;
        case GLFW_KEY_C:
            colormap_mode = (colormap_mode + 1) % 3;
            fprintf(stderr, "Colormap: %s.\n", CMAP_NAMES[colormap_mode]);
            break;
        case GLFW_KEY_LEFT_BRACKET:
            blend_factor = std::max(blend_factor - 0.1f, 0.0f);
            fprintf(stderr, "Blend factor: %.2f\n", blend_factor);
            break;
        case GLFW_KEY_RIGHT_BRACKET:
            blend_factor = std::min(blend_factor + 0.1f, 1.0f);
            fprintf(stderr, "Blend factor: %.2f\n", blend_factor);
            break;
        case GLFW_KEY_L:
            length_by_speed = !length_by_speed;
            fprintf(stderr, "Arrow length: %s.\n", length_by_speed ? "magnitude-based" : "constant");
            break;
        case GLFW_KEY_X:
            streamlines.clear();
            pathlines.clear();
            fprintf(stderr, "Cleared all seeds.\n");
            break;
        case GLFW_KEY_R:
            rake_count = (rake_count == 1) ? 5 : (rake_count == 5 ? 9 : 1);
            fprintf(stderr, "Rake size = %d.\n", rake_count);
            break;
        case GLFW_KEY_O:
            rake_horizontal = !rake_horizontal;
            fprintf(stderr, "Rake orientation: %s.\n", rake_horizontal ? "horizontal" : "vertical");
            break;
        case GLFW_KEY_COMMA:
            path_time_scale = std::max(path_time_scale * 0.5f, 0.001f);
            fprintf(stderr, "Pathline time scale = %.4f\n", path_time_scale);
            break;
        case GLFW_KEY_PERIOD:
            path_time_scale = std::min(path_time_scale * 2.0f, 1.0f);
            fprintf(stderr, "Pathline time scale = %.4f\n", path_time_scale);
            break;
        case GLFW_KEY_Q:
        case GLFW_KEY_ESCAPE:
            exit(0);
            break;
        default:
            fprintf(stderr,
                "\nKeyboard:\n"
                "1/2/3  load %s / %s / %s\n"
                "0      cycle timestep (recomputes streamlines)\n"
                "a      toggle arrows\n"
                "t      toggle streamlines\n"
                "p      toggle pathlines\n"
                "+/-    sampling rate\n"
                "i/k    dt\n"
                "m      integration: Euler/RK2/RK4 (current: %s)\n"
                "c      colormap: off/rainbow/cool-warm (current: %s)\n"
                "[ / ]  blend factor (current: %.2f)\n"
                "l      arrow length: const / by speed (current: %s)\n"
                "r      rake size (1/5/9, current: %d)\n"
                "o      rake orientation H/V (current: %s)\n"
                ", / .  pathline time scale (current: %.4f)\n"
                "s      cycle scalar field\n"
                "b      cycle clear color\n"
                "x      clear all seeds\n"
                "left mouse  : single seed\n"
                "right mouse : rake seed\n"
                "q/esc  quit\n",
                filenames[0], filenames[1], filenames[2],
                INTEG_NAMES[integ_method], CMAP_NAMES[colormap_mode],
                blend_factor, length_by_speed ? "by speed" : "constant",
                rake_count, rake_horizontal ? "H" : "V",
                path_time_scale);
            break;
    }
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    if (action != GLFW_PRESS) return;
    if (!grid_data_loaded) return;
    if (!en_streamline && !en_pathline) {
        fprintf(stderr, "Enable streamlines (t) or pathlines (p) to seed.\n");
        return;
    }

    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    int wx, wy;
    glfwGetWindowSize(window, &wx, &wy);
    if (wx <= 0 || wy <= 0) return;

    float gx = (float)(xpos / wx) * (vol_dim[0] - 1);
    float gy = (float)(1.0 - ypos / wy) * (vol_dim[1] - 1);
    if (gx < 0 || gx > vol_dim[0] - 1 || gy < 0 || gy > vol_dim[1] - 1) return;

    auto seedAt = [](int sx, int sy) {
        if (en_streamline) seedStreamline(sx, sy);
        if (en_pathline)   seedPathline(sx, sy);
    };

    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        seedAt((int)gx, (int)gy);
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        int N = std::max(rake_count, 1);
        if (N == 1) { seedAt((int)gx, (int)gy); return; }
        // Spread rake across 60% of the relevant dimension.
        float span = (rake_horizontal ? (vol_dim[0] - 1) : (vol_dim[1] - 1)) * 0.6f;
        float spacing = span / (N - 1);
        for (int i = 0; i < N; ++i) {
            float off = -span * 0.5f + i * spacing;
            float x = rake_horizontal ? (gx + off) : gx;
            float y = rake_horizontal ? gy        : (gy + off);
            if (x < 0 || x > vol_dim[0] - 1) continue;
            if (y < 0 || y > vol_dim[1] - 1) continue;
            seedAt((int)x, (int)y);
        }
    }
}

static void errorCallback(int error, const char* description)
{
    fprintf(stderr, "GLFW Error: %s\n", description);
}

// ---------------------------------------------------------------------------
// Data loading
// ---------------------------------------------------------------------------

void loadNextTimestep( void )
{
    loaded_timestep = (loaded_timestep + 1) % num_timesteps;
    DownloadScalarFieldAsTexture();
}

void LoadData( char* base_filename )
{
    reset_rendering_props();

    char filename[ 160 ];
    strcpy(filename, base_filename);
    strcat(filename, ".gri");

    fprintf(stderr, "loading grid file %s\n", filename);

    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Cannot open file %s for reading.\n", filename);
        return;
    }

    char header[40];
    fread(header, sizeof(char), 40, fp);
    sscanf(header, "SN4DB %hu %hu %hu %d %d %f",
           &vol_dim[0], &vol_dim[1], &vol_dim[2],
           &num_scalar_fields, &num_timesteps, &timestep);
    fclose(fp);

    fprintf(stderr, "dimensions: %d x %d x %d.\n", vol_dim[0], vol_dim[1], vol_dim[2]);
    fprintf(stderr, "info: scalar fields=%d, timesteps=%d, dt=%f.\n",
            num_scalar_fields, num_timesteps, timestep);

    loaded_timestep = 0;
    LoadVectorData(base_filename);

    glfwSetWindowSize(window, vol_dim[0], vol_dim[1]);
    grid_data_loaded = true;
}

void LoadVectorData( const char* filename )
{
    fprintf(stderr, "loading data file(s) for %s\n", filename);

    char ts_name[160];
    if (num_timesteps > 1) std::snprintf(ts_name, sizeof(ts_name), "%s.%.5d.dat", filename, 0);
    else                   std::snprintf(ts_name, sizeof(ts_name), "%s.dat", filename);

    FILE* fp = fopen(ts_name, "rb");
    if (!fp) {
        fprintf(stderr, "Cannot open file %s for reading.\n", ts_name);
        return;
    }
    fclose(fp);

    data_size = vol_dim[0] * vol_dim[1] * vol_dim[2];

    if (vector_array)  { delete[] vector_array;  vector_array  = nullptr; }
    if (scalar_fields) { delete[] scalar_fields; scalar_fields = nullptr; }
    if (scalar_bounds) { delete[] scalar_bounds; scalar_bounds = nullptr; }

    vector_array  = new float[data_size * 3 * num_timesteps];
    scalar_fields = new float[data_size * num_scalar_fields * num_timesteps];
    scalar_bounds = new float[2 * num_scalar_fields * num_timesteps];

    int num_total_fields = num_scalar_fields + 3;
    float* tmp = new float[(size_t)data_size * num_total_fields * num_timesteps];

    for (int k = 0; k < num_timesteps; ++k) {
        char times_name[160];
        if (num_timesteps > 1) std::snprintf(times_name, sizeof(times_name), "%s.%.5d.dat", filename, k);
        else                    std::snprintf(times_name, sizeof(times_name), "%s.dat", filename);
        fp = fopen(times_name, "rb");
        if (!fp) { fprintf(stderr, "Cannot open %s\n", times_name); continue; }
        fread(&tmp[(size_t)k * data_size * num_total_fields], sizeof(float),
              (size_t)data_size * num_total_fields, fp);
        fclose(fp);

        for (int i = 0; i < num_scalar_fields; ++i) {
            float min_val = std::numeric_limits<float>::infinity();
            float max_val = -std::numeric_limits<float>::infinity();
            int offset = i * data_size * num_timesteps;
            for (int j = 0; j < data_size; ++j) {
                float val = tmp[j * num_total_fields + 3 + i + (size_t)k * data_size * num_total_fields];
                scalar_fields[j + (size_t)k * data_size + offset] = val;
                if (i == 0) {
                    if (toggle_xy) {
                        vector_array[3 * j + 0 + (size_t)3 * k * data_size] = tmp[j * num_total_fields + 1 + (size_t)k * data_size * num_total_fields];
                        vector_array[3 * j + 1 + (size_t)3 * k * data_size] = tmp[j * num_total_fields + 0 + (size_t)k * data_size * num_total_fields];
                        vector_array[3 * j + 2 + (size_t)3 * k * data_size] = tmp[j * num_total_fields + 2 + (size_t)k * data_size * num_total_fields];
                    } else {
                        vector_array[3 * j + 0 + (size_t)3 * k * data_size] = tmp[j * num_total_fields + 0 + (size_t)k * data_size * num_total_fields];
                        vector_array[3 * j + 1 + (size_t)3 * k * data_size] = tmp[j * num_total_fields + 1 + (size_t)k * data_size * num_total_fields];
                        vector_array[3 * j + 2 + (size_t)3 * k * data_size] = tmp[j * num_total_fields + 2 + (size_t)k * data_size * num_total_fields];
                    }
                }
                if (val < min_val) min_val = val;
                if (val > max_val) max_val = val;
            }
            scalar_bounds[2 * i     + k * num_scalar_fields * 2] = min_val;
            scalar_bounds[2 * i + 1 + k * num_scalar_fields * 2] = max_val;
        }

        // Normalize scalar fields to [0,1] for the texture.
        for (int i = 0; i < num_scalar_fields; ++i) {
            int offset = i * data_size * num_timesteps;
            float lo = scalar_bounds[2 * i     + k * num_scalar_fields * 2];
            float hi = scalar_bounds[2 * i + 1 + k * num_scalar_fields * 2];
            if (lo < 0.0f && hi > 0.0f) {
                float scale = 0.5f / std::max(-lo, hi);
                for (int j = 0; j < data_size; ++j)
                    scalar_fields[offset + j + (size_t)k * data_size] = 0.5f + scalar_fields[offset + j + (size_t)k * data_size] * scale;
                scalar_bounds[2 * i     + k * num_scalar_fields * 2] = 0.5f + lo * scale;
                scalar_bounds[2 * i + 1 + k * num_scalar_fields * 2] = 0.5f + hi * scale;
            } else {
                float sign = hi <= 0.0f ? -1.0f : 1.0f;
                float denom = (hi - lo);
                if (std::abs(denom) < 1e-12f) denom = 1.0f;
                float scale = sign / denom;
                for (int j = 0; j < data_size; ++j)
                    scalar_fields[offset + j + (size_t)k * data_size] = (scalar_fields[offset + j + (size_t)k * data_size] - lo) * scale;
                scalar_bounds[2 * i     + k * num_scalar_fields * 2] = 0.0f;
                scalar_bounds[2 * i + 1 + k * num_scalar_fields * 2] = 1.0f;
            }
        }
    }
    delete[] tmp;

    // Compute global max vector magnitude across all timesteps for integration normalization.
    vec_norm = 1e-6f;
    for (int t = 0; t < num_timesteps; ++t) {
        for (int j = 0; j < data_size; ++j) {
            int kbase = 3 * j + (size_t)3 * t * data_size;
            float vx = vector_array[kbase + 0];
            float vy = vector_array[kbase + 1];
            float m = std::sqrt(vx * vx + vy * vy);
            if (m > vec_norm) vec_norm = m;
        }
    }
    fprintf(stderr, "max |v| = %f\n", vec_norm);

    DownloadScalarFieldAsTexture();
    scalar_data_loaded = true;
}

void DownloadScalarFieldAsTexture( void )
{
    fprintf(stderr, "downloading scalar field to 2D texture\n");

    if (scalar_field_texture == 0) {
        glGenTextures(1, &scalar_field_texture);
    }
    glBindTexture(GL_TEXTURE_2D, scalar_field_texture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    int datasize = vol_dim[0] * vol_dim[1];
    int sf_offset = (loaded_timestep + current_scalar_field * num_timesteps) * datasize;

    // Avoid 4-byte row alignment surprises with non-multiple-of-4 widths.
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_R16, vol_dim[0], vol_dim[1], 0,
                 GL_RED, GL_FLOAT, &scalar_fields[sf_offset]);
}

// ---------------------------------------------------------------------------
// GL initialization
// ---------------------------------------------------------------------------

bool initApplication(int argc, char **argv)
{
    std::string version((const char *)glGetString(GL_VERSION));
    std::stringstream stream(version);
    unsigned major = 0, minor = 0;
    char dot;
    stream >> major >> dot >> minor;
    std::cout << "OpenGL Version " << major << "." << minor << std::endl;
    if (major < 3) {
        std::cout << "Need at least OpenGL 3.x." << std::endl;
        return false;
    }
    return true;
}

void reset_rendering_props( void )
{
    num_scalar_fields = 0;
}

static void initLineGL()
{
    glGenVertexArrays(1, &lineVAO);
    glBindVertexArray(lineVAO);
    glGenBuffers(1, &lineVBO);
    glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void setup() {
    LoadData(filenames[0]);
    loaded_file = 0;

    DownloadScalarFieldAsTexture();

    // Compile & link shaders (paths derived from CMake-provided source dir).
    vectorProgram.compileShader(projectPath("shaders/vertex.vs").c_str());
    vectorProgram.compileShader(projectPath("shaders/fragment.fs").c_str());
    vectorProgram.link();

    lineProgram.compileShader(projectPath("shaders/line.vs").c_str());
    lineProgram.compileShader(projectPath("shaders/line.fs").c_str());
    lineProgram.link();

    quad.init();
    initLineGL();
}

// ---------------------------------------------------------------------------
// Render
// ---------------------------------------------------------------------------

void render() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Background scalar field with colormap.
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, scalar_field_texture);
    vectorProgram.use();
    model = glm::mat4(1.0f);
    vectorProgram.setUniform("vertexColor", glm::vec4(0.0f));
    vectorProgram.setUniform("model", model);
    vectorProgram.setUniform("txtr", 0);
    vectorProgram.setUniform("colormapMode", colormap_mode);
    vectorProgram.setUniform("blendFactor", blend_factor);
    quad.render();

    // Overlays: glyphs / streamlines / pathlines.
    drawGlyphs();
    drawStreamlines();
    drawPathlines();
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

int main(int argc, char** argv)
{
    view_width = 0;
    view_height = 0;
    toggle_xy = 0;

    en_arrow = false;
    en_streamline = false;
    en_pathline = false;
    sampling_rate = 15;
    dt = 0.5f;
    path_time_scale = 0.05f;

    integ_method = 1;          // RK2 by default
    colormap_mode = 0;
    blend_factor = 0.7f;
    length_by_speed = true;

    min_vec_mag = 1e-4f;
    max_acc_length = 4096.0f;
    max_steps = 5000;

    rake_count = 1;
    rake_horizontal = false;

    reset_rendering_props();

    vector_array = nullptr;
    scalar_fields = nullptr;
    scalar_bounds = nullptr;
    scalar_field_texture = 0;
    grid_data_loaded = false;
    scalar_data_loaded = false;
    current_scalar_field = 0;
    clearColor = 0;

    static char path0[256], path1[256], path2[256];
    std::snprintf(path0, sizeof(path0), "%s/data/block/c_block",          PROJECT_SOURCE_DIR_STR);
    std::snprintf(path1, sizeof(path1), "%s/data/tube/tube",              PROJECT_SOURCE_DIR_STR);
    std::snprintf(path2, sizeof(path2), "%s/data/hurricane/hurricane_p_tc", PROJECT_SOURCE_DIR_STR);
    filenames[0] = path0;
    filenames[1] = path1;
    filenames[2] = path2;

    glfwSetErrorCallback(errorCallback);

    if (!glfwInit()) exit(EXIT_FAILURE);

    // Request GL 4.1 core (max supported on macOS).
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(gWindowWidth, gWindowHeight,
                              "AMCS/CS247 Scientific Visualization", nullptr, nullptr);
    if (!window) { glfwTerminate(); exit(EXIT_FAILURE); }

    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetFramebufferSizeCallback(window, frameBufferCallback);

    glfwMakeContextCurrent(window);

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    view_width = width;
    view_height = height;

    gladLoadGL();

    if (!initApplication(argc, argv)) { glfwTerminate(); exit(EXIT_FAILURE); }

    setup();

    // Print menu once on launch.
    keyCallback(window, GLFW_KEY_BACKSLASH, 0, GLFW_PRESS, 0);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    while (!glfwWindowShouldClose(window))
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        render();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return EXIT_SUCCESS;
}
