#pragma once
#include "Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <set>
#include <algorithm>


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
	void fill_holes();
	void Laplace();
	void Execute();
	double Angle(MeshKernel::VertexHandle vh, MeshKernel::VertexHandle vh1, MeshKernel::VertexHandle vh2);
private:
	MeshKernel::SurfaceMesh& mesh;
	int n, bn, in;
	std::vector<MeshKernel::VertexHandle> points;
	std::unordered_map<int, int> to_new, to_old;
	std::unordered_map<int, bool> vis;
	std::vector<std::set<MeshKernel::EdgeHandle> > loops;
	std::set<MeshKernel::VertexHandle> fake_points;
	std::map<std::tuple<int, int, int>, int> points_to_faces;
	int longest_loop_idx;
	Eigen::SparseMatrix<double> A, Aii, Aib, Abi, Abb;
	Eigen::MatrixXd K, k;
};