
#pragma once

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"
#include "cgra/cgra_image.hpp"
#include "skeleton_model.hpp"


// Basic model that holds the shader, mesh and transform for drawing.
// Can be copied and modified for adding in extra information for drawing
// including textures for texture mapping etc.
struct basic_model {
	GLuint shader = 0;


	cgra::gl_mesh mesh;
	glm::vec3 color{1,1,0};

	// Phong shading --> TO BE CHANGED
	float ambient = 0.2;
	float diffuse = 0.5;
	float specular = 0.8;
	float cook = 0.1;
	float oren_nayar = 0.1;
	float roughness = 0.1;
	float material = 1.5;
	float texture_strength = 0;
	glm::vec3 lightdirection = glm::vec3(-0.8, -0.8, 0.8);
	glm::vec3 viewdirection = glm::vec3(-0.8, -0.8, 0.8);

	glm::mat4 modelTransform{1.0};
	GLuint m_texture;
	cgra::rgba_image m_texture_data;
	cgra::rgba_image m_normal_texture_data;

	GLuint m_texture_normal;

	void draw(const glm::mat4 &view, const glm::mat4 proj);
	void plane_draw(const glm::mat4& view, const glm::mat4 proj);
	void depthMapping(const glm::mat4& view, const glm::mat4 proj, const glm::mat4 lightSpaceMatrix,GLuint shadowMapTexture, unsigned int depthProg);
	void drawWithShadow(const glm::mat4& view, const glm::mat4 proj,  const glm::mat4 lightSpaceMatrix, unsigned int shadowMapTexture, unsigned int secondPassProg);
};

struct basic_sphere {
	int radius = 5;
	int slices = 20;
	int stacks = 20;

	void draw();
};

// Main application class
//
class Application {
private:
	// window
	glm::vec2 m_windowsize;
	GLFWwindow *m_window;


	glm::vec3 lightdirection = glm::vec3(-0.8, -0.8, 0.8);


	// oribital camera
	int width, height;
	glm::mat4 view, proj;
	float m_pitch = .86;
	float m_yaw = -.86;
	float m_distance = 20;
	bool isShadow = false;

	// last input
	bool m_leftMouseDown = false;
	glm::vec2 m_mousePosition;

	// drawing flags
	bool m_show_axis = false;
	bool m_show_grid = false;
	bool m_showWireframe = false;
	bool m_attachToLight = false;
	bool rebuild_models = false;

	// geometry
	basic_model m_model;
	basic_model m_plane;

	//Asssignment
	int numLatitudeLines = 30;
	int numLongitudeLines = 30;
	float radius = 1.0f;

	// CubeSphere
	int subdivisions = 0;

	//Torus
	float outerRadius = 5.0; // outer radius
	float innerRadius = 1.0; // inner radius

	int outerSubdivisions = 15; // outer subdivisions
	int innerSubdivisions = 15; // inner subdivisions


	//SHADOW
	const unsigned int SHADOW_WIDTH = 5000, SHADOW_HEIGHT = 5000;

	unsigned int depthMapFBO;
	unsigned int depthMap;
	GLuint shadowMapProg;
	GLuint secondPassProg;
public:
	// setup
	Application(GLFWwindow *);

	// disable copy constructors (for safety)
	Application(const Application&) = delete;
	Application& operator=(const Application&) = delete;

	//CORE
	cgra::mesh_builder drawUVSphere();

	//COMPLETION
	cgra::mesh_builder drawSquareSphere();

	//CHALLENGE
	cgra::mesh_builder drawTorus();

	cgra::mesh_builder drawPlane();
	// SHADERS
	GLuint shader();

	// rendering callbacks (every frame)
	void render();
	void renderGUI();

	// input callbacks
	void cursorPosCallback(double xpos, double ypos);
	void mouseButtonCallback(int button, int action, int mods);
	void scrollCallback(double xoffset, double yoffset);
	void keyCallback(int key, int scancode, int action, int mods);
	void charCallback(unsigned int c);



	//SHADOW
	glm::mat4 lightSpaceMat();
	void shadowMapPass();
	void shadowRender();
	void newRender();
	void newRenderSecond();
	void lightingPass();
};

