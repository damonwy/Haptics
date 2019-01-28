#pragma once

#include <iostream>
#include <fstream>
#include <functional>

#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <omp.h>

#ifdef _MEX_
#include "mex.h"
#endif

#define _USE_MATH_DEFINES
#include <cmath> 

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>

#include <string>
#include <cstring>
#include <cassert>
#include <json.hpp>

#include "SE3.h"
#include "MatlabDebug.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"
#include "Camera.h"
#include "GLSL.h"
