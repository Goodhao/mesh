#pragma once
#include "Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>


class MeshParameterization {
public:
	MeshParameterization(MeshKernel::SurfaceMesh& _mesh) :
		mesh(_mesh) {
	};
	void GetBoundaryPoints();
	void GetPlane();
	void GetProjetion();
	void UpdateBoundary();
	void Parameterization();
	void Execute();
private:
	MeshKernel::SurfaceMesh& mesh;
	std::vector<MeshKernel::Vertex> boundary_points, boundary_projection;
	std::map<int, int> resorted_id, reverse_id;
	double _a, _b, _c, _d; // 平面 ax+by+cz+d=0
	int n; // 内点个数
};