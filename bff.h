#pragma once
#include "Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>


class bff {
public:
	bff(MeshKernel::SurfaceMesh& _mesh) :
		mesh(_mesh) {
	};
	void Init();
	double CotanAngles(MeshKernel::EdgeHandle eh);
	void Parameterization();
	void InitCurvature();
	std::vector<MeshKernel::EdgeHandle> AroundEdge(MeshKernel::VertexHandle vh);
	MeshKernel::EdgeHandle NextEdge(MeshKernel::VertexHandle vh);
	void Laplace();
	void Execute();
private:
	MeshKernel::SurfaceMesh& mesh;
	int n, bn, in;
	std::vector<MeshKernel::VertexHandle> points;
	std::map<int, int> to_new, to_old;
	std::map<int, bool> vis;
	Eigen::MatrixXd A, Aii, Aib, Abi, Abb;
	Eigen::MatrixXd K, k;
};