#version 330 core
out vec4 FragColor;

in VS_OUT {
    vec3 FragPos;
    vec3 Normal;
    vec2 TexCoords;
    vec4 FragPosLightSpace;
} fs_in;

uniform sampler2D shadowMap;
uniform vec3 viewPos;

uniform float As;
uniform float Ds;
uniform float Ss;
uniform float Cs;
uniform float Ts;
uniform float roughness;

uniform vec3 fixedLightPos;
uniform vec3 uColor;

float ShadowCalculation(vec4 fragPosLightSpace)
{
    vec3 projCoords = fragPosLightSpace.xyz;
    projCoords = projCoords * 0.5 + 0.5;
    float closestDepth = texture(shadowMap, projCoords.xy).r; 
    float currentDepth = projCoords.z;
    vec3 lightDir = normalize(fixedLightPos - fs_in.FragPos);
    float bias = max(0.05 * (1.0 - dot(fs_in.Normal, lightDir)), 0.005); 
    float shadow = currentDepth - bias > closestDepth ? 0.7 : 0.0;  
    return shadow;
}

void main()
{
    //vec3 color = vec3(0.3);

    vec3 lightDir = normalize(fixedLightPos - fs_in.FragPos);
    float light = abs(dot(normalize(fs_in.Normal), lightDir));

    vec3 color = mix(uColor / 4, uColor, light);

    vec3 viewDir = normalize(viewPos - fs_in.FragPos);
    vec3 halfwayDir = normalize(lightDir + viewDir);

    float diff = max(dot(fs_in.Normal, lightDir), 0.0);
    float diffuse = Ds * diff;

    float spec = pow(max(dot(fs_in.Normal, halfwayDir), 0.0), roughness);
    float specular = Ss * spec;

    float shadow = ShadowCalculation(fs_in.FragPosLightSpace);

    color = color * As + (1.0 - shadow) * (diffuse + specular);

    FragColor = vec4(color, 1.0);
}
