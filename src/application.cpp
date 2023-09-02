
// std
#include <iostream>
#include <string>
#include <chrono>

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

// project
#include "application.hpp"
#include "cgra/cgra_geometry.hpp"
#include "cgra/cgra_gui.hpp"
#include "cgra/cgra_image.hpp"
#include "cgra/cgra_shader.hpp"
#include "cgra/cgra_wavefront.hpp"


using namespace std;
using namespace cgra;
using namespace glm;

#define M_PI 3.1415926535897932384626433832795

void basic_model::draw(const glm::mat4 &view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;

	mat4 lightView = lookAt(lightdirection, -lightdirection, vec3(0, 1, 0));


	mat4 lightview = lightView * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));

	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));
	glUniform3fv(glGetUniformLocation(shader, "lightPos"), 1, value_ptr(lightdirection));
	glUniform3fv(glGetUniformLocation(shader, "viewPos"), 1, value_ptr(viewdirection));

	glUniform1f(glGetUniformLocation(shader, "As"), ambient);
	glUniform1f(glGetUniformLocation(shader, "Ds"), diffuse);
	glUniform1f(glGetUniformLocation(shader, "Ss"), specular);
	glUniform1f(glGetUniformLocation(shader, "Cs"), cook);
	glUniform1f(glGetUniformLocation(shader, "ONs"), oren_nayar);

	glUniform1f(glGetUniformLocation(shader, "Ts"), texture_strength);
	glUniform1f(glGetUniformLocation(shader, "roughness"), roughness);
	glUniform1f(glGetUniformLocation(shader, "material"), material);

	glUniform1i(1, 4);
	glUniform1i(5, 3);

	mesh.draw(); // draw
}
void basic_model::plane_draw(const glm::mat4& view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));
	glUniform3fv(glGetUniformLocation(shader, "lightPos"), 1, value_ptr(lightdirection));

	glUniform1f(glGetUniformLocation(shader, "As"), ambient);
	glUniform1f(glGetUniformLocation(shader, "Ds"), diffuse);
	glUniform1f(glGetUniformLocation(shader, "Ss"), specular);
	glUniform1f(glGetUniformLocation(shader, "Cs"), cook);
	glUniform1f(glGetUniformLocation(shader, "ONs"), oren_nayar);

	glUniform1f(glGetUniformLocation(shader, "Ts"), texture_strength);
	glUniform1f(glGetUniformLocation(shader, "roughness"), roughness);
	glUniform1f(glGetUniformLocation(shader, "material"), material);


	glUniform1i(3, 4);
	//glUniform1i(2, shader);
	mesh.draw(); // draw
}

void basic_model::depthMapping(const glm::mat4& view, const glm::mat4 proj, const glm::mat4 lightSpaceMatrix, GLuint shadowMapTexture, unsigned int depthProg) {
	glUseProgram(depthProg);
	glm::mat4 modelview = view * modelTransform;
	glUniformMatrix4fv(glGetUniformLocation(depthProg, "lightSpaceMatrix"), 1, GL_FALSE, glm::value_ptr(lightSpaceMatrix));
	glUniformMatrix4fv(glGetUniformLocation(depthProg, "uModelViewMatrix"), 1, GL_FALSE, glm::value_ptr(modelview));
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, shadowMapTexture);
	mesh.draw();
}


void  basic_model::drawWithShadow(const glm::mat4& view, const glm::mat4 proj,const glm::mat4 lightSpaceMatrix, unsigned int shadowMapTexture, unsigned int secondPassProg)
{
	glUseProgram(secondPassProg);
	glm::mat4 modelview = view * modelTransform;

	//vertex data
	glUniformMatrix4fv(glGetUniformLocation(secondPassProg, "uProjectionMatrix"), 1, GL_FALSE, glm::value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(secondPassProg, "uModelViewMatrix"), 1, GL_FALSE, glm::value_ptr(modelview));
	glUniformMatrix4fv(glGetUniformLocation(secondPassProg, "uViewMatrix"), 1, GL_FALSE, glm::value_ptr(view));
	glUniformMatrix4fv(glGetUniformLocation(secondPassProg, "lightSpaceMatrix"), 1, GL_FALSE, glm::value_ptr(lightSpaceMatrix));

	//fragment data
	glUniform1i(glGetUniformLocation(secondPassProg, "shadowMap"), 0);
	glUniform3fv(glGetUniformLocation(secondPassProg, "lightPos"), 1, glm::value_ptr(lightdirection));
	glUniform3fv(glGetUniformLocation(secondPassProg, "viewPos"), 1, glm::value_ptr(view));
	glUniform3fv(glGetUniformLocation(secondPassProg, "uColor"), 1, value_ptr(color));

	GLint ambientUniformLocation = glGetUniformLocation(secondPassProg, "As");
	glUniform1f(ambientUniformLocation, ambient);
	GLint diffuseUniformLocation = glGetUniformLocation(secondPassProg, "Ds");
	glUniform1f(diffuseUniformLocation, diffuse);
	GLint specularUniformLocation = glGetUniformLocation(secondPassProg, "Ss");
	glUniform1f(specularUniformLocation, specular);



	glUniformMatrix4fv(glGetUniformLocation(secondPassProg, "lightSpaceMatrix"), 1, GL_FALSE, glm::value_ptr(lightSpaceMatrix));
	glBindTexture(GL_TEXTURE_2D, shadowMapTexture);

	mesh.draw();
}

mesh_builder Application::drawUVSphere() {
	// One vertex at every latitude-longitude intersection,
	mesh_builder mb;
	mesh_vertex mv;


	float latStep =	M_PI / numLatitudeLines;
	float longStep = 2 * M_PI / numLongitudeLines;

	float x, y, z;
	float nx, ny, nz;
	float u, v;

	float uStep = 1.0f / numLongitudeLines;
	float vStep = 1.0f / numLatitudeLines;
	float uCoord = 0;
	float vCoord = 0;

	for (int i = 0; i < numLatitudeLines + 1; i++) {
		for (int j = 0; j < numLongitudeLines + 1; j++) {
			
			// Positions
			x = radius * sin(latStep * i) * cos(longStep * j);
			y = radius * sin(latStep * i) * sin(longStep * j);
			z = radius * cos(latStep * i);
			vec3 positions = vec3(x, y, z);

			//	UVs
			u = uCoord;
			v = vCoord;
			vec2 uv = vec2(u, v);
			
			
			// Add the vertex to the sphere mesh
			mb.push_vertex({ positions, positions, uv}); 

			uCoord += uStep;
		}

		uCoord = 0.0f;
		vCoord += vStep;
	}

	// Add indices for the sphere faces
	for (int i = 0; i < numLatitudeLines + 1; i++) {
		for (int j = 0; j < numLongitudeLines + 1; j++) {
			
			int p1 = i * numLongitudeLines + j;
			int p2 = i * numLongitudeLines + j + 1;

			int p3 = (i + 1) * numLongitudeLines + j + 1;
			int p4 = (i + 1) * numLongitudeLines + j;

			// Triangle 1
			mb.push_index(p1);
			mb.push_index(p2);
			mb.push_index(p3);

			// Triangle 2
			mb.push_index(p1);
			mb.push_index(p3);
			mb.push_index(p4);
		}
	}
	return mb;
}
mesh_builder Application::drawSquareSphere() {

	mesh_builder mb;
	mesh_vertex mv;
	std::vector<vec3> positions;
	std::vector<vec3> normals;
	std::vector<vec2> uvs;

	std::vector<vec3> positiveXvec;

	std::vector<unsigned int> indices;
	const float DEG2RAD = acos(-1) / 180.0f;

	vec3 normLong; // normal of longitudinal plane rotating along Y-axis
	vec3 normLat; // normal of latitudinal plane rotating along Z-axis
	vec3 normIntersect; // direction vector intersecting 2 planes, n1 x n2

	float a1;           // longitudinal angle along Y-axis
	float a2;           // latitudinal angle along Z-axis


	// compute the number of vertices per row, 2^n + 1
	int pointsPerRow = (int)pow(2, subdivisions) + 1;

	// rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
	for (unsigned int i = 0; i < pointsPerRow; ++i)
	{
			// normal for latitudinal plane
			// if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
			// therefore, it is rotating (0,1,0) vector by latitude angle a2
			a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
			normLat = vec3(-sin(a2), cos(a2), 0);

			// rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
			for (unsigned int j = 0; j < pointsPerRow; ++j)
			{
				// normal for longitudinal plane
				// if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
				// therefore, it is rotating (0,0,-1) vector by longitude angle a1
				a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
				normLong = vec3(-sin(a1), 0, -cos(a1));
				
				// find direction vector of intersected line, n1 x n2
				normIntersect.x = normLong[1] * normLat[2] - normLong[2] * normLat[1];
				normIntersect.y = normLong[2] * normLat[0] - normLong[0] * normLat[2];
				normIntersect.z = normLong[0] * normLat[1] - normLong[1] * normLat[0];

				// normalize direction vector
				float scale = 1 / sqrt(normIntersect.x * normIntersect.x + normIntersect.y * normIntersect.y + normIntersect.z * normIntersect.z);
				normIntersect *= scale;

				// add a vertex into array
				positiveXvec.push_back(normIntersect);
			}
		}
		float s, t;
		int k = 0, k1, k2;

		// build +X face
		for (unsigned int i = 0; i < pointsPerRow; ++i)
		{
			k1 = i * pointsPerRow;     // index for curr row
			k2 = k1 + pointsPerRow;    // index for next row
			t = (float)i / (pointsPerRow - 1);

			for (unsigned int j = 0; j < pointsPerRow; ++j, k += 1, ++k1, ++k2)
			{

				s = (float)j / (pointsPerRow - 1);

				// vec version
				positions.push_back(positiveXvec[k]*radius);
				normals.push_back(vec3(positiveXvec[k]));
				uvs.push_back(vec2(s, t));

				// add indices
				if (i < (pointsPerRow - 1) && j < (pointsPerRow - 1))
				{
					// first triangle
					indices.push_back(k1);
					indices.push_back(k2);
					indices.push_back(k1 + 1);
					
					// sec triangle
					indices.push_back(k1+1);
					indices.push_back(k2);
					indices.push_back(k2 + 1);
				}
			}
		}
		
		//array size and index for building next face
		unsigned int startIndex;                    // starting index for next face
		int vertexSize = (int)positions.size();      // vertex array size of +X face
		int indexSize = (int)indices.size();        // index array size of +X face

		// build -X face by negating x and z
		startIndex = positions.size();
		for (int i = 0, j = 0; i < vertexSize; i += 1, j += 1)
		{
			positions.push_back(vec3(-positions[i].x, positions[i].y, -positions[i].z));
			uvs.push_back(vec2(uvs[j].s, uvs[j].t));
			normals.push_back(vec3(-normals[i].x, normals[i].y, -normals[i].z));
		}
		for (int i = 0; i < indexSize; ++i)
		{
			indices.push_back(startIndex + indices[i]);
		}

		//build +Y face by swapping x=>y, y=>-z, z=>-x
		startIndex = positions.size();
		for (int i = 0, j = 0; i < vertexSize; i += 1, j += 1)
		{
			positions.push_back(vec3(-positions[i].z, positions[i].x, -positions[i].y));
			uvs.push_back(vec2(uvs[j].s, uvs[j].t));
			normals.push_back(vec3(-normals[i].z, normals[i].x, -normals[i].y));
		}
		for (int i = 0; i < indexSize; ++i)
		{
			indices.push_back(startIndex + indices[i]);
		}

		// build -Y face by swapping x=>-y, y=>z, z=>-x
		startIndex = positions.size();
		for (int i = 0, j = 0; i < vertexSize; i += 1, j += 1)
		{
			positions.push_back(vec3(-positions[i].z, -positions[i].x, positions[i].y));
			uvs.push_back(vec2(uvs[j].s, uvs[j].t));
			normals.push_back(vec3(-normals[i].z, -normals[i].x, normals[i].y));
		}
		for (int i = 0; i < indexSize; ++i)
		{
			indices.push_back(startIndex + indices[i]);
		}

		// build +Z face by swapping x=>z, z=>-x
		startIndex = positions.size();
		for (int i = 0, j = 0; i < vertexSize; i += 1, j += 1)
		{
			positions.push_back(vec3(-positions[i].z, positions[i].y, positions[i].x));
			uvs.push_back(vec2(uvs[j].s, uvs[j].t));
			normals.push_back(vec3(-normals[i].z, normals[i].y, normals[i].x));
		}
		for (int i = 0; i < indexSize; ++i)
		{
			indices.push_back(startIndex + indices[i]);
		}

		// build -Z face by swapping x=>-z, z=>x
		startIndex = positions.size();
		for (int i = 0, j = 0; i < vertexSize; i += 1, j += 1)
		{
			positions.push_back(vec3(positions[i].z, positions[i].y, -positions[i].x));
			uvs.push_back(vec2(uvs[j].s, uvs[j].t));
			normals.push_back(vec3(normals[i].z, normals[i].y, -normals[i].x));
		}
		for (int i = 0; i < indexSize; ++i)
		{
			indices.push_back(startIndex + indices[i]);
		}

		
	//// create mesh 
	for (unsigned int i = 0; i < indices.size(); ++i) {
		
			mb.push_index(i);
			mb.push_vertex(mesh_vertex{
				positions[indices[i]],
				normals[indices[i]],
				uvs[indices[i]]
				});
	}

	return mb;
}

mesh_builder Application::drawTorus() {
	// intialise geometry
	mesh_builder mb;
	mesh_vertex mv;
	std::vector<vec3> positions;
	std::vector<vec3> normals;
	std::vector<vec2> uvs;
	std::vector<unsigned int> indices;
	
	// create torus geometry
	float du = 2 * M_PI / outerSubdivisions;
	float dv = 2 * M_PI / innerSubdivisions;

	// Torus Geometry Calculations
	for (int i = 0; i < outerSubdivisions; i++) {

		float u = i * du;
		for (int j = 0; j <= innerSubdivisions; j++) {

			float v = (j % innerSubdivisions) * dv;
			for (int k = 0; k < 2; k++){
			
				float uu = u + k*du;
				// compute vertex
				float x = (outerRadius + innerRadius * cos(v)) * cos(uu);
				float y = (outerRadius + innerRadius * cos(v)) * sin(uu);
				float z = innerRadius * sin(v);

				// add vertex
				positions.push_back(vec3(x, y, z));

				// compute normal 
				float nx = cos(v) * cos(uu);
				float ny = cos(v) * sin(uu);
				float nz = sin(v);

				// add normal 
				normals.push_back(vec3(nx, ny, nz));

				// compute texture coords
				float tx = uu / (2 * M_PI);
				float ty = v / (2 * M_PI);
				
				// add tex coords
				uvs.push_back(vec2(tx, ty));

			}

			v += dv;
		}
	}

	// Torus Indices assigning
	for (unsigned int i = 2; i < positions.size(); ++i) {
		indices.push_back(i-2);
		indices.push_back(i-1);
		indices.push_back(i);
	}
	

	// create mesh from indices
	for (unsigned int i = 0; i < indices.size(); ++i) {


		mb.push_index(i);
		mb.push_vertex(mesh_vertex{
			positions[indices[i]],
			normals[indices[i]],
			uvs[indices[i]]
			});
	}

	return mb;
}
mesh_builder Application::drawPlane() {
	mesh_builder mb;
	mesh_vertex mv;
	std::vector<vec3> positions;
	std::vector<vec3> normals;
	std::vector<vec2> uvs;

	std::vector<unsigned int> indices;
	int subdivisions = 1; // You can adjust this value to control the plane's subdivisions

	int size = 20;
	int pointsPerRow = (int)pow(2, subdivisions) + 1;

	float u, v;

	// Generate vertices and indices for the plane
	for (int i = 0; i <= subdivisions; ++i) {
		float y = static_cast<float>(i) / static_cast<float>(subdivisions) - 0.5f ;
		//v = (float)i / (pointsPerRow - 1);
		for (int j = 0; j <= subdivisions; ++j) {
			float x = static_cast<float>(j) / static_cast<float>(subdivisions) - 0.5f;

			// Calculate vertex position
			positions.push_back(vec3(x, 0.0f, y) * vec3(size) +vec3(0,-2,0));

			// Calculate vertex normal (for a simple plane, this will be the same for all vertices)
			normals.push_back(vec3(0, 1, 0.0f));

			// Calculate vertex texture coordinates
			float u = (float)i / subdivisions ;
			float v = (float)j / subdivisions;
			//u = (float)j / (pointsPerRow - 1);
			std::cout << u << "\n";
			std::cout << v << "\n";
			uvs.push_back(vec2(u, v));
		}
	}

	// Generate indices for the plane
	for (int i = 0; i < subdivisions; ++i) {

		for (int j = 0; j < subdivisions; ++j) {
			int v0 = i * (subdivisions + 1) + j;
			int v1 = v0 + 1;
			int v2 = (i + 1) * (subdivisions + 1) + j;
			int v3 = v2 + 1;

			// First triangle
			indices.push_back(v0);
			indices.push_back(v2);
			indices.push_back(v1);

			// Second triangle
			indices.push_back(v1);
			indices.push_back(v2);
			indices.push_back(v3);
		}
	}

	// Create mesh from indices
	for (unsigned int i = 0; i < indices.size(); ++i) {
		mb.push_index(i);
		mb.push_vertex(mesh_vertex{
			positions[indices[i]],
			normals[indices[i]],
			uvs[indices[i]]
			});
	}

	//lighting 
	m_plane.ambient = 0.7;
	m_plane.color = { 1,1,1 };

	return mb;
}


Application::Application(GLFWwindow *window) : m_window(window) {
	

	m_model.m_texture_data = rgba_image((CGRA_SRCDIR + std::string("//res//textures//Texture.png")));
	m_model.m_normal_texture_data = rgba_image((CGRA_SRCDIR + std::string("//res//textures//NormalMap.png")));

	m_model.m_texture = m_model.m_texture_data.uploadTexture();
	m_model.m_texture_normal = m_model.m_normal_texture_data.uploadNormalTexture();

	std::cout <<"gluint: " << m_model.m_texture << "\n";

	m_plane.m_texture_data = rgba_image((CGRA_SRCDIR + std::string("//res//textures//checkerboard.jpg")));
	
	m_plane.m_texture = m_plane.m_texture_data.uploadCheckersTexture();
	//m_plane.shader = 1;

	std::cout << "pgluint: " << m_plane.m_texture << "\n";

	shader_builder sb;
    sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//cook_torrance_vert.glsl"));
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//cook_torrance_frag.glsl"));
	GLuint shader = sb.build();

	shader_builder plane_sb;
	plane_sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//oren_nayar_vert.glsl"));
	plane_sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//oren_nayar_frag.glsl"));
	GLuint plane_shader = plane_sb.build();



	m_model.shader = shader;
	m_model.mesh = drawUVSphere().build();


	m_plane.shader = plane_shader;
	m_plane.mesh = drawPlane().build();




	// Shadow
	shader_builder DepthMap;
	DepthMap.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//depth_map_vert.glsl"));
	DepthMap.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//depth_map_frag.glsl"));

	shader_builder secondPass;
	secondPass.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//shadow_vert.glsl"));
	secondPass.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//shadow_frag.glsl"));

	shadowMapProg = DepthMap.build();
	secondPassProg = secondPass.build();

	glGenFramebuffers(1, &depthMapFBO);
	glGenTextures(1, &depthMap);
	glBindTexture(GL_TEXTURE_2D, depthMap);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthMap, 0);

	// Set the color buffer of the framebuffer to none
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);

	// Check if the framebuffer is complete and report an error if it's not
	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
		std::cout << "Error: Framebuffer is not complete!" << std::endl;
	}

	// Unbind the framebuffer to avoid accidental changes
	glBindFramebuffer(GL_FRAMEBUFFER, 0);


}

GLuint Application::shader() {
	shader_builder sb;
	sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//cook_torrance_vert.glsl"));
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//cook_torrance_frag.glsl"));
	GLuint shader = sb.build();

	return shader;
}


void Application::render() {

	// retrieve the window hieght
	glfwGetFramebufferSize(m_window, &width, &height);

	m_windowsize = vec2(width, height); // update window size
	glViewport(0, 0, width, height); // set the viewport to draw to the entire window

	// clear the back-buffer
	glClearColor(0.3f, 0.3f, 0.4f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// enable flags for normal/forward rendering
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	// projection matrix
	proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);

	// view matrix
	view = translate(mat4(1), vec3(0, 0, -m_distance))
		* rotate(mat4(1), m_pitch, vec3(1, 0, 0))
		* rotate(mat4(1), m_yaw, vec3(0, 1, 0));
	
	
	vec3 viewPos;
	

	m_model.lightdirection = lightdirection;
	m_plane.lightdirection = lightdirection;

	m_model.viewdirection = vec3(m_pitch, m_yaw, m_distance);
	m_plane.viewdirection = vec3(m_pitch, m_yaw, m_distance);

	if (m_attachToLight) {
		view = lookAt(lightdirection, -lightdirection, vec3(0, 1, 0));
	}

	

	// helpful draw options
	if (m_show_grid) drawGrid(view, proj);
	if (m_show_axis) drawAxis(view, proj);
	glPolygonMode(GL_FRONT_AND_BACK, (m_showWireframe) ? GL_LINE : GL_FILL);

	if (isShadow) {
		shadowRender();
	}
	else { 
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		m_model.draw(view, proj);
	}

}


void Application::renderGUI() {

	// setup window
	ImGui::SetNextWindowPos(ImVec2(5, 5), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(300, 200), ImGuiSetCond_Once);
	ImGui::Begin("Options", 0);

	// display current camera parameters
	ImGui::Text("Application %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::SliderFloat("Pitch", &m_pitch, -pi<float>() / 2, pi<float>() / 2, "%.2f");
	ImGui::SliderFloat("Yaw", &m_yaw, -pi<float>(), pi<float>(), "%.2f");
	ImGui::SliderFloat("Distance", &m_distance, 0, 100, "%.2f", 2.0f);

	// helpful drawing options
	ImGui::Checkbox("Show axis", &m_show_axis);
	ImGui::SameLine();
	ImGui::Checkbox("Show grid", &m_show_grid);
	ImGui::Checkbox("Wireframe", &m_showWireframe);
	ImGui::Checkbox("LightView", &m_attachToLight);
	ImGui::Checkbox("Show shadow scene", &isShadow);

	ImGui::SameLine();
	if (ImGui::Button("Screenshot")) rgba_image::screenshot(true);

	
	ImGui::Separator();

	// example of how to use input boxes
	static float exampleInput;
	if (ImGui::InputFloat("example input", &exampleInput)) {
		cout << "example input changed to " << exampleInput << endl;
	}

	ImGui::SliderInt("lattitude", &numLatitudeLines, 5, 50);
	ImGui::SliderInt("longitude", &numLongitudeLines, 5, 50);
	ImGui::SliderFloat("Radius", &radius, 1, 10, "%.2f", 2.0f);
	ImGui::Separator();
	ImGui::SliderInt("Subdivisions", &subdivisions, 0, 5);
	ImGui::Separator();
	ImGui::SliderFloat("Outer Radius", &outerRadius, 1, 10, "%.2f", 2.0f);
	ImGui::SliderFloat("Inner Radius", &innerRadius, 1, 10, "%.2f", 2.0f);
	ImGui::SliderInt("Outer Subdivisions", &outerSubdivisions, 5, 50);
	ImGui::SliderInt("Inner Subdivisions", &innerSubdivisions, 5, 50);
	
	// Geometry
	static int geometry_index = 0;
	if (ImGui::Combo("Geometry", &geometry_index, "Sphere\0Cube\0Torus\0", 3)) {
		switch (geometry_index) {		
		case 0: 
			
			m_model.shader = shader();
			m_model.mesh = drawUVSphere().build();
			break;
		case 1:
			m_model.shader = shader();
			m_model.mesh = drawSquareSphere().build();
			break;
		case 2: 
			m_model.shader = shader();
			m_model.mesh = drawTorus().build();
			break;
		}
	}

	ImGui::Separator();
	ImGui::SliderFloat3("Light Position", value_ptr(lightdirection), -10, 10, "%.1f");
	ImGui::SliderFloat3("Model Color", value_ptr(m_model.color), 0, 1, "%.2f");
	ImGui::Separator();
	ImGui::SliderFloat("Ambient", &m_model.ambient, 0, 1, "%.1f");
	ImGui::SliderFloat("Diffuse", &m_model.diffuse, 0, 1, "%.1f");
	ImGui::SliderFloat("Specular", &m_model.specular, 0, 1, "%.1f");
	ImGui::Separator();
	ImGui::SliderFloat("Cook", &m_model.cook, 0, 1, "%.1f");
	ImGui::SliderFloat("Oren Nayar", &m_model.oren_nayar, 0, 1, "%.1f");
	ImGui::SliderFloat("Roughness", &m_model.roughness, 0, 1, "%.1f");
	ImGui::SliderFloat("Material", &m_model.material, 0, 3, "%.1f");
	ImGui::Separator();
	ImGui::SliderFloat("Texture", &m_model.texture_strength, 0, 1, "%.1f");
	// finish creating window
	ImGui::End();
}


void Application::cursorPosCallback(double xpos, double ypos) {
	if (m_leftMouseDown) {
		vec2 whsize = m_windowsize / 2.0f;

		// clamp the pitch to [-pi/2, pi/2]
		m_pitch += float(acos(glm::clamp((m_mousePosition.y - whsize.y) / whsize.y, -1.0f, 1.0f))
			- acos(glm::clamp((float(ypos) - whsize.y) / whsize.y, -1.0f, 1.0f)));
		m_pitch = float(glm::clamp(m_pitch, -pi<float>() / 2, pi<float>() / 2));

		// wrap the yaw to [-pi, pi]
		m_yaw += float(acos(glm::clamp((m_mousePosition.x - whsize.x) / whsize.x, -1.0f, 1.0f))
			- acos(glm::clamp((float(xpos) - whsize.x) / whsize.x, -1.0f, 1.0f)));
		if (m_yaw > pi<float>()) m_yaw -= float(2 * pi<float>());
		else if (m_yaw < -pi<float>()) m_yaw += float(2 * pi<float>());
	}

	// updated mouse position
	m_mousePosition = vec2(xpos, ypos);
}


void Application::mouseButtonCallback(int button, int action, int mods) {
	(void)mods; // currently un-used

	// capture is left-mouse down
	if (button == GLFW_MOUSE_BUTTON_LEFT)
		m_leftMouseDown = (action == GLFW_PRESS); // only other option is GLFW_RELEASE
}


void Application::scrollCallback(double xoffset, double yoffset) {
	(void)xoffset; // currently un-used
	m_distance *= pow(1.1f, -yoffset);
}


void Application::keyCallback(int key, int scancode, int action, int mods) {
	(void)key, (void)scancode, (void)action, (void)mods; // currently un-used
}


void Application::charCallback(unsigned int c) {
	(void)c; // currently un-used
}



// SHADOW STUFF

mat4 Application::lightSpaceMat() {

	float near_plane = -50, far_plane = 30;

	glm::vec3 lightPosition = -lightdirection;
	glm::vec3 lightInvDir = -glm::normalize(lightPosition);
	mat4 depthProjectionMatrix = ortho<float>(-50, 50, -50, 50, near_plane, far_plane);
	mat4 depthViewMatrix = lookAt(lightPosition, vec3(0, 0, 0), vec3(0, 1, 0));
	mat4 depthModelMatrix = mat4(1.0);
	glm::mat4 lightSpaceMatrix = depthProjectionMatrix * depthViewMatrix * depthModelMatrix;


	return lightSpaceMatrix;
}

void Application::shadowMapPass() {



	// Check what shadowProg
	//glUseProgram(shadowMapProg);


	glm::mat4 lightSpaceMatrix = lightSpaceMat();
	glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
	glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
	glClear(GL_DEPTH_BUFFER_BIT);


	// Possibly models??
	m_model.depthMapping(view, proj, lightSpaceMatrix, depthMap, shadowMapProg);
	m_plane.depthMapping(view, proj, lightSpaceMatrix, depthMap, shadowMapProg);

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void Application::shadowRender() {

	// clear bufferd
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glUseProgram(secondPassProg);

	shadowMapPass();
	lightingPass();
}

// LIGHTING PASS
void Application::lightingPass() {

	glViewport(0, 0, width, height);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBindTexture(GL_TEXTURE_2D, depthMap);
	newRenderSecond();
}

void Application::newRender() {
	glm::mat4 lightSpaceMatrix = lightSpaceMat();
	glUseProgram(shadowMapProg);
	GLint lightSpaceMatrixLocation = glGetUniformLocation(shadowMapProg, "lightSpaceMatrix");
	glUniformMatrix4fv(lightSpaceMatrixLocation, 1, GL_FALSE, glm::value_ptr(lightSpaceMatrix));
	glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
	glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
}

void Application::newRenderSecond() {


	glUseProgram(secondPassProg);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, depthMap);

	GLint shadowMapLoc = glGetUniformLocation(secondPassProg, "shadowMap");
	glUniform1i(shadowMapLoc, 0);
	glUniform3fv(glGetUniformLocation(secondPassProg, "viewPos"), 1, value_ptr(view));
	mat4 lightSpaceMatrix = lightSpaceMat();

	m_plane.drawWithShadow(view, proj, lightSpaceMatrix, depthMap, secondPassProg);
	m_model.drawWithShadow(view, proj, lightSpaceMatrix, depthMap, secondPassProg);
}








