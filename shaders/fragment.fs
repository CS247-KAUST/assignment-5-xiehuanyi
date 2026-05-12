#version 410 core

in vec2 texCoord;

uniform sampler2D txtr;
uniform vec4 vertexColor;
uniform int colormapMode;     // 0 = off (grayscale), 1 = rainbow, 2 = cool-warm
uniform float blendFactor;    // 0 = pure grayscale, 1 = full color

out vec4 fragColor;

vec3 rainbow(float t) {
    t = clamp(t, 0.0, 1.0);
    if (t < 0.25)      return mix(vec3(0.0, 0.0, 1.0), vec3(0.0, 1.0, 1.0), t / 0.25);
    else if (t < 0.50) return mix(vec3(0.0, 1.0, 1.0), vec3(0.0, 1.0, 0.0), (t - 0.25) / 0.25);
    else if (t < 0.75) return mix(vec3(0.0, 1.0, 0.0), vec3(1.0, 1.0, 0.0), (t - 0.50) / 0.25);
    else               return mix(vec3(1.0, 1.0, 0.0), vec3(1.0, 0.0, 0.0), (t - 0.75) / 0.25);
}

vec3 coolWarm(float t) {
    t = clamp(t, 0.0, 1.0);
    vec3 cool   = vec3(0.230, 0.299, 0.754);
    vec3 middle = vec3(0.865, 0.865, 0.865);
    vec3 warm   = vec3(0.706, 0.016, 0.150);
    if (t < 0.5) return mix(cool, middle, t / 0.5);
    else         return mix(middle, warm, (t - 0.5) / 0.5);
}

void main() {
    float scalar = texture(txtr, texCoord).r;
    vec3 gray = vec3(scalar);
    vec3 mapped = gray;
    if (colormapMode == 1) mapped = rainbow(scalar);
    else if (colormapMode == 2) mapped = coolWarm(scalar);
    vec3 color = mix(gray, mapped, clamp(blendFactor, 0.0, 1.0));
    fragColor = vec4(color, 1.0) + vertexColor;
}
