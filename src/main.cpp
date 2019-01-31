#include "rmpch.h"

#ifndef _GLIBCXX_USE_NANOSLEEP
#define _GLIBCXX_USE_NANOSLEEP
#endif
#include <thread>

#include "Scene.h"
#include <Eigen/Core>

#define RECORD_VIDEO 0

using namespace std;
using namespace Eigen;

bool keyToggles[256] = {false}; // only for English keyboards!

GLFWwindow *window; // Main application window
string RESOURCE_DIR = ""; // Where the resources are loaded from

shared_ptr<Camera> camera;
shared_ptr<Program> prog;
shared_ptr<Program> progSimple;
shared_ptr<Program> progSoft;
shared_ptr<Scene> scene;

char pixels[4 * 1920 * 1080];
int steps = 0;

static void error_callback(int error, const char *description)
{
	cerr << description << endl;
}

static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
}

static void char_callback(GLFWwindow *window, unsigned int key)
{
	keyToggles[key] = !keyToggles[key];
	switch(key) {
		case 'h':
			scene->step();
			break;
		case 'r':
			scene->reset();
			break;
	}
}

static void cursor_position_callback(GLFWwindow* window, double xmouse, double ymouse)
{
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if(state == GLFW_PRESS) {
		camera->mouseMoved(xmouse, ymouse);
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// Get the current mouse position.
	double xmouse, ymouse;
	glfwGetCursorPos(window, &xmouse, &ymouse);
	// Get current window size.
	int width, height;
	glfwGetWindowSize(window, &width, &height);
	if(action == GLFW_PRESS) {
		bool shift = mods & GLFW_MOD_SHIFT;
		bool ctrl  = mods & GLFW_MOD_CONTROL;
		bool alt   = mods & GLFW_MOD_ALT;
		camera->mouseClicked(xmouse, ymouse, shift, ctrl, alt);
	}
}

static void init()
{
	GLSL::checkVersion();
	
	// Set background color
	glClearColor(0.6f, 0.6f, 0.6f, 1.0f);
	// Enable z-buffer test
	glEnable(GL_DEPTH_TEST);
	// Enable alpha blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	progSimple = make_shared<Program>();
	progSimple->setShaderNames(RESOURCE_DIR + "simple_vert.glsl", RESOURCE_DIR + "simple_frag.glsl");
	progSimple->setVerbose(true); // Set this to true when debugging.
	progSimple->init();
	progSimple->addUniform("P");
	progSimple->addUniform("MV");
	//progSimple->setVerbose(false);

	prog = make_shared<Program>();
	prog->setVerbose(true); // Set this to true when debugging.
	prog->setShaderNames(RESOURCE_DIR + "vert.glsl", RESOURCE_DIR + "frag.glsl");
	prog->init();
	prog->addUniform("P");
	prog->addUniform("MV");
	prog->addAttribute("aPos");
	prog->addAttribute("aNor");
	prog->addUniform("intensity_1");
	prog->addUniform("lightPos1");
	prog->addUniform("lightPos2");
	prog->addUniform("intensity_2");
	prog->addUniform("ka");
	prog->addUniform("kd");
	prog->addUniform("ks");
	prog->addUniform("s");
	//prog->setVerbose(false);

	progSoft = make_shared<Program>();
	progSoft->setVerbose(true); // Set this to true when debugging.
	progSoft->setShaderNames(RESOURCE_DIR + "phong_vert.glsl", RESOURCE_DIR + "phong_frag.glsl");
	progSoft->init();
	progSoft->addUniform("P");
	progSoft->addUniform("MV");
	progSoft->addUniform("kdFront");
	progSoft->addUniform("kdBack");
	progSoft->addAttribute("aPos");
	progSoft->addAttribute("aNor");
	//progSoft->setVerbose(false);


	camera = make_shared<Camera>();
	camera->setInitDistance(140.5f);

	scene = make_shared<Scene>();
	scene->load(RESOURCE_DIR);
	//scene->tare();
	scene->init();
	
	// If there were any OpenGL errors, this will print something.
	// You can intersperse this line in your code to find the exact location
	// of your OpenGL error.
	GLSL::checkError(GET_FILE_LINE);
}

void render()
{
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height);
	
	// Use the window size for camera.
	glfwGetWindowSize(window, &width, &height);
	camera->setAspect((float)width/(float)height);
	
	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(keyToggles[(unsigned)'c']) {
		glEnable(GL_CULL_FACE);
	} else {
		glDisable(GL_CULL_FACE);
	}
	if(keyToggles[(unsigned)'l']) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	
	auto P = make_shared<MatrixStack>();
	auto MV = make_shared<MatrixStack>();
	
	// Apply camera transforms
	P->pushMatrix();
	camera->applyProjectionMatrix(P);
	MV->pushMatrix();
	camera->applyViewMatrix(MV);

	// Draw grid
	progSimple->bind();
	glUniformMatrix4fv(progSimple->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(progSimple->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glLineWidth(2.0f);
	float x0 = -10.0f;
	float x1 = 10.0f;
	float z0 = -10.0f;
	float z1 = 10.0f;
	int gridSize = 10;
	glLineWidth(1.0f);
	glBegin(GL_LINES);
	for(int i = 1; i < gridSize; ++i) {
		if(i == gridSize/2) {
			glColor3f(0.1f, 0.1f, 0.1f);
		} else {
			glColor3f(0.8f, 0.8f, 0.8f);
		}
		float x = x0 + i / (float)gridSize * (x1 - x0);
		glVertex3f(x, 0.0f, z0);
		glVertex3f(x, 0.0f, z1);
	}
	for(int i = 1; i < gridSize; ++i) {
		if(i == gridSize/2) {
			glColor3f(0.1f, 0.1f, 0.1f);
		} else {
			glColor3f(0.8f, 0.8f, 0.8f);
		}
		float z = z0 + i / (float)gridSize * (z1 - z0);
		glVertex3f(x0, 0.0f, z);
		glVertex3f(x1, 0.0f, z);
	}
	glEnd();
	glColor3f(0.4f, 0.4f, 0.4f);
	glBegin(GL_LINE_LOOP);
	glVertex3f(x0, 0.0f, z0);
	glVertex3f(x1, 0.0f, z0);
	glVertex3f(x1, 0.0f, z1);
	glVertex3f(x0, 0.0f, z1);
	glEnd();
	progSimple->unbind();

	// Draw scene
	//prog->bind();
	scene->draw(MV, prog, progSimple,progSoft, P);
	
	//////////////////////////////////////////////////////
	// Cleanup
	//////////////////////////////////////////////////////
	
	// Pop stacks
	MV->popMatrix();
	P->popMatrix();
	
	GLSL::checkError(GET_FILE_LINE);
}

void stepperFunc()
{
	
	while(true) {
		
		if(keyToggles[(unsigned)' ']) {
			
			scene->step();
			steps += 1;

		}
		this_thread::sleep_for(chrono::microseconds(1));
	}
}

int main(int argc, char **argv)
{
	Eigen::initParallel();
	//omp_set_num_threads(1);
	Eigen::setNbThreads(1);


	if(argc < 2) {
		cout << "Please specify the resource directory." << endl;
		return 0;
	}
	RESOURCE_DIR = argv[1] + string("/");
	
	// Set error callback.
	glfwSetErrorCallback(error_callback);
	// Initialize the library.
	if(!glfwInit()) {
		return -1;
	}
	// Create a windowed mode window and its OpenGL context.
	window = glfwCreateWindow(1920, 1080, "Reduced Coordinate Rigid Body FEM Simulation", NULL, NULL);
	if(!window) {
		glfwTerminate();
		return -1;
	}
	// Make the window's context current.
	glfwMakeContextCurrent(window);
	// Initialize GLEW.
	glewExperimental = true;
	if(glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return -1;
	}
	glGetError(); // A bug in glewInit() causes an error that we can safely ignore.
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	// Set vsync.
	glfwSwapInterval(1);
	// Set keyboard callback.
	glfwSetKeyCallback(window, key_callback);
	// Set char callback.
	glfwSetCharCallback(window, char_callback);
	// Set cursor position callback.
	glfwSetCursorPosCallback(window, cursor_position_callback);
	// Set mouse button callback.
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	// Initialize scene.
	init();
	// Start simulation thread.
	
	thread stepperThread(stepperFunc);

	//// Loop until the user closes the window.
	while(!glfwWindowShouldClose(window)) {
			
		// Render scene.
		render();
		// Swap front and back buffers.
		glfwSwapBuffers(window);
		// Poll for and process events.
		glfwPollEvents();
	}
	// Quit program.
	stepperThread.detach();
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
