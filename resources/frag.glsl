#version 120

uniform vec3 lightPos;
uniform vec3 lightPos1;
uniform float intensity_1;
uniform vec3 lightPos2;
uniform float intensity_2;
uniform vec3 ka;
uniform vec3 kd;
uniform vec3 ks;
uniform float s;

//varying vec3 color; // passed from the vertex shader
varying vec3 normal;
varying vec3 v_position;

void main()
{
vec3 n = normalize(normal);
vec3 e = vec3(0-v_position.x,0-v_position.y,0-v_position.z);//eye vector
e = normalize(e);

vec3 lv =vec3(lightPos1.x-v_position.x,lightPos1.y-v_position.y,lightPos1.z-v_position.z);//light vector
lv = normalize(lv);
vec3 cd = kd * max(0,dot(lv,n));
vec3 h =vec3(e.x+lv.x,e.y+lv.y,e.z+lv.z);
h = normalize(h);
vec3 cs = ks*pow(max(0,dot(h,n)),s);
vec3 c1 = ka + cd + cs;

vec3 lv2 =vec3(lightPos2.x-v_position.x,lightPos2.y-v_position.y,lightPos2.z-v_position.z);//light vector
lv2 = normalize(lv2);
vec3 cd2 = kd * max(0,dot(lv2,n));
vec3 h2 =vec3(e.x+lv2.x,e.y+lv2.y,e.z+lv2.z);
h2 = normalize(h2);
vec3 cs2 = ks*pow(max(0,dot(h2,n)),s);
vec3 c2 = ka + cd2 + cs2;

vec3 c= c1*intensity_1+c2*intensity_2;

gl_FragColor = vec4(c.r, c.g, c.b, 1.0);

}
