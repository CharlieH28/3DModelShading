
#version 330 core


/// uniform data
uniform mat4 uModelViewMatrix;
uniform mat4 lightSpaceMatrix;


layout (location = 0) in vec3 aPos;

void main() {

	//gl_position = lightSpaceMatrix * uModelViewMatrix * vec4(aPos,1.0);
	gl_Position = lightSpaceMatrix * uModelViewMatrix  * vec4(aPos, 1.0);
}