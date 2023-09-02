#version 330 core
#extension GL_ARB_explicit_attrib_location : require
#extension GL_ARB_explicit_uniform_location : require

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;
uniform vec3 lightPos;
uniform vec3 viewPos;


uniform float As;
uniform float Ds;
uniform float Ss;
uniform float Ts;
uniform float Cs;
uniform float ONs;

uniform float roughness;
uniform float material;


//layout (location = 1) vec3 uColor;

uniform sampler2D shadowMap;

layout (location = 1) uniform sampler2D u_texture;
layout (location = 5) uniform sampler2D normal_texture;

// viewspace data (this must match the output of the fragment shader)

in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
	mat3 TBN;
} f_in;

// framebuffer output
out vec4 fb_color;
/*
      float mirror = 0.0;
      float air = 1.0;
      float bubble = 1.1;
      float ice = 1.31;
      float water = 1.33;
      float glass = 1.5;

      float standard = 2.0;
      float steel = 2.5;

*/
//float material = 1.5;

vec3 lightColor = vec3(1, 1, 1);
//vec3 lightPos = vec3(-0.8,-0.8,0.8);

void main() {

	//NORMAL
	vec3 normal = vec3(texture(normal_texture, f_in.textureCoord)).xyz;
	normal = normal * 2.0 - 1.0;   
	normal = normalize(normal);

	// calculate lighting direction and the view direction
	//vec3 viewDir = normalize(-f_in.position);
	//vec3 lightDir = normalize((lightPos));

	vec3 viewDir = f_in.TBN * normalize(-f_in.position);
	vec3 lightDir = f_in.TBN * normalize((lightPos));

	float light = abs(dot(normalize(f_in.normal), lightDir));

	vec3 color = mix(uColor / 4, uColor, light);


	vec3 text_color = vec3(texture(u_texture, f_in.textureCoord));

	// PHONG
	float ambientStrength = 0.2;
	vec3 ambient = As * lightColor;

	vec3 norm = normalize(f_in.normal);


	float diff = max(dot(normal, lightDir), 0.0);
	vec3 diffuse = diff * lightColor * Ds;

	// SPEC
	float specularStrength = 0.8;
	vec3 reflectDir = reflect(-lightDir, normal);

	float shininess = 50;

	float spec = pow(max(dot(lightDir, reflectDir),0.0), shininess);
	vec3 specular = Ss * spec * lightColor;

	// COOK

	float PI = 3.1415926535897932384626433832795;
	float NdotV = clamp(dot(norm, viewDir), 0.0, 1.0);
	float NdotL = clamp(dot(norm, lightDir), 0.0, 1.0);

	vec3 H = normalize(lightDir + viewDir);
	float NdotH = clamp(dot(norm, H), 0.0, 1.0);
	float VdotH = max(0, dot(lightDir,H));

	//Fresnel
	float Kr = pow((1.0 - material) / (1.0 + material), 2.0);
	float F = Kr + (1.0 - Kr) * pow((1.0 - VdotH), 5.0);

	// Geometric Attenuation Factor
	float roughness2 = pow(clamp(roughness, 0.01, 0.99), 2.0);
    float g1 = (NdotL * material) / (NdotL + sqrt(roughness2 + (1.0 - roughness2) * pow(NdotL, 2.0)));
    float g2 = (NdotV * material) / (NdotV + sqrt(roughness2 + (1.0 - roughness2) * pow(NdotV, 2.0)));
    float G = g1 * g2;

	// Normal Distribution Function
	float NH2 = pow(NdotH, 2.0);
    roughness2 = pow(clamp(roughness, 0.01, 0.99), 2.0);

    float denom = NH2 * roughness2 + (1.0 - NH2);
    float D = roughness2 / (PI * pow(denom, 2.0));


	vec3 cookTorrance = Cs * (lightColor * (F*G*D) / ((PI * NdotV)));

	// OREN NAYAR
	float angleVN = acos(NdotV);
	float angleLN = acos(diff);

    float alpha = max(angleVN, angleLN);
    float beta = min(angleVN, angleLN);
    float gamma = cos(angleVN - angleLN);

    float roughness_sq = roughness * roughness;

    float A = 1.0 - 0.5 * (roughness_sq / (roughness_sq + 0.57));
    float B = 0.45 * (roughness_sq / (roughness_sq + 0.09));
    float C = sin(alpha) * tan(beta);

    float OrenNayar = ONs*diff * (A + (B * max(0.0, gamma) * C));
	
	//RESULT
	vec3 tex = (text_color);
	
	//if (Ts > 0)
	//	color = tex
	   

	vec3 result = (ambient+diffuse+cookTorrance+specular+OrenNayar) *  text_color;

	// output to the framebuffer
	fb_color = vec4(result, 1);
}