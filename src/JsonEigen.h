#pragma once
#ifndef MUSCLEMASS_SRC_JSONEIGEN_H_
#define MUSCLEMASS_SRC_JSONEIGEN_H_

#include <Eigen/Core>
#include <json.hpp>


using json = nlohmann::json;

namespace Eigen {
	void to_json(json &j, const Vector3d &v);

	void to_json(json &j, const Matrix3d &m);

	void to_json(json &j, const Matrix4d &m4);

	void from_json(const json &j, Vector3d &v);

	void from_json(const json &j, Matrix3d &m);

	void from_json(const json &j, Matrix4d &m4);
}

#endif