#version 120

uniform mat4 P;
uniform mat4 MV;

attribute vec4 aPos; // in object space
attribute vec3 aNor; // in object space

//varying vec3 color; // Pass to fragment shader
varying vec3 normal;
varying vec3 v_position;

void main()
{
	gl_Position = P * MV * aPos;
	vec4 n = vec4(aNor,0.0);
	n = MV * n;
	normal = vec3(n.x,n.y,n.z);

	vec4 p = aPos;
	p = MV * p;
	v_position = vec3(p.x,p.y,p.z);


}
