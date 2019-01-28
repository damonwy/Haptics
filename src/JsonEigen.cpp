#include "rmpch.h"
#include "JsonEigen.h"

namespace Eigen {
	void to_json(json &j, const Vector3d &v) {
		std::vector<double> data{ v(0), v(1), v(2) };
		j = data;
	}

	void to_json(json &j, const Matrix3d &m) {
		std::vector<std::string> rows(3);

		for (int i = 0; i < 3; ++i) {
			auto &row = rows[i];
			for (int j = 0; j < 3; ++j) row += std::to_string(m(i, j)) + " ";
			row += std::to_string(m(i, 3));
		}

		j.push_back(rows[0]);
		j.push_back(rows[1]);
		j.push_back(rows[2]);
	}

	void to_json(json &j, const Matrix4d &m4) {
		std::vector<std::string> rows(4);

		for (int i = 0; i < 4; ++i) {
			auto &row = rows[i];
			for (int j = 0; j < 4; ++j) row += std::to_string(m4(i, j)) + " ";
			row += std::to_string(m4(i, 4));
		}

		j.push_back(rows[0]);
		j.push_back(rows[1]);
		j.push_back(rows[2]);
		j.push_back(rows[3]);
	}

	void from_json(const json &j, Vector3d &v) {
		std::vector<double> data = j.get<std::vector<double>>();
		v(0) = data[0];
		v(1) = data[1];
		v(2) = data[2];
	}

	void from_json(const json &j, Matrix3d &m) {
		const std::vector<std::string> rows = j.get<std::vector<std::string>>();
		assert(rows.size() == 3);

		uint32_t count = 0;
		for (const auto &r : rows) {
			std::stringstream ss(r);
			ss >> m(count, 0) >> m(count, 1) >> m(count, 2);
			count += 1;
		}
	}


	void from_json(const json &j, Matrix4d &m4) {
		const std::vector<std::string> rows = j.get<std::vector<std::string>>();
		assert(rows.size() == 4);

		uint32_t count = 0;
		for (const auto &r : rows) {
			std::stringstream ss(r);
			ss >> m4(count, 0) >> m4(count, 1) >> m4(count, 2) >> m4(count, 3);
			count += 1;
		}
	}
}