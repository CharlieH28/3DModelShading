#version 330 core
#extension GL_ARB_explicit_attrib_location : require
#extension GL_ARB_explicit_uniform_location : require

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;

uniform float As;
uniform float Ds;
uniform float Ss;
uniform float Ts;

//layout (location = 1) vec3 uColor;

layout(location = 1) uniform sampler2D u_texture;


// viewspace data (this must match the output of the fragment shader)
in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} f_in;

// framebuffer output
out vec4 fb_color;


vec3 lightColor = vec3(1, 1, 1);
//vec3 lightDir = lightDir(0.25,0.25,-1);

void main() {
	// calculate lighting (hack)
	vec3 eye = normalize(-f_in.position);
	float light = abs(dot(normalize(f_in.normal), eye));

	vec3 color = mix(uColor / 4, uColor,light);
	vec3 text_color = vec3(texture(u_texture, f_in.textureCoord));

	// load shader 
	float ambientStrength = 0.2;
	vec3 ambient = As * lightColor;

	vec3 norm = normalize(f_in.normal);
	//lightDire = normalize(-lightDirection);

	float diff = max(dot(norm, eye), 0.0);
	vec3 diffuse = diff * lightColor*Ds;

	float specularStrength = 0.8;
	vec3 reflectDir = reflect(-eye, norm);
	vec3 viewDir = normalize(-f_in.position);

	float spec = pow(max(dot(eye, reflectDir),0.0), 32);
	vec3 specular = Ss * spec * lightColor;


	vec3 tex = (text_color*Ts);
	vec3 result = (ambient+diffuse+specular) * color+tex;

	// output to the frambuffer
	fb_color = vec4(result, 1);
}