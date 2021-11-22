#include "meshparameterization.h"

void MeshParameterization::Execute() {
	printf("Mesh Parameterization is starting\n");
	GetBoundaryPoints();
	GetPlane();
	GetProjetion();
	UpdateBoundary();
	Parameterization();
	printf("Mesh Parameterization success\n");
}

void MeshParameterization::GetBoundaryPoints() {
	int idx = 0;
	for (auto vp : mesh.allvertices()) {
		if (!mesh.isOnBoundary(vp.first)) {
			reverse_id[vp.first] = idx;
			resorted_id[idx++] = vp.first;
		}
	}
	n = idx;
	for (auto vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			boundary_points.push_back(vp.second);
			reverse_id[vp.first] = idx;
			resorted_id[idx++] = vp.first;
		}
	}
}

void MeshParameterization::GetPlane() {
	Eigen::MatrixXf A(3,3);
	Eigen::VectorXf b(3);
	for (auto p : boundary_points) {
		Eigen::MatrixXf tmp(3, 3);
		for (int i = 0; i < 3; i++) {
			tmp(i, 0) = p.x();
			tmp(i, 1) = p.y();
			tmp(i, 2) = p.z();
		}
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (i == 0) tmp(i, j) *= p.x();
				else if (i == 1) tmp(i, j) *= p.y();
				else if (i == 2) tmp(i, j) *= p.z();
			}
		}
		A += tmp;
		b(0) += -p.x(), b(1) += -p.y(), b(2) += -p.z();
	}
	Eigen::MatrixXf res = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	_a = res(0, 0);
	_b = res(1, 0);
	_c = res(2, 0);
	_d = 1;
}

void MeshParameterization::GetProjetion() {
	std::vector<MeshKernel::Vertex> v = boundary_points;
	assert(_d == 1);
	assert(_a * _a + _b * _b + _c * _c > 0);
	for (auto &p : v) {
		double t = _a * p.x() + _b * p.y() + _c * p.z() + _d;
		t /= _a * _a + _b * _b + _c * _c;
		//std::cout << t << std::endl;
		p.x() -= _a * t;
		p.y() -= _b * t;
		p.z() -= _c * t;
	}
	boundary_projection = v;
}

void MeshParameterization::UpdateBoundary() {
	int idx = 0;
	for (auto vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			mesh.vertices(vp.first) = boundary_projection[idx++];
		}
	}
	/*
	for (auto vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			std::cout << mesh.vertices(vp.first).x() << " " << mesh.vertices(vp.first).y()
				<< " " << mesh.vertices(vp.first).z() << std::endl;
		}
	}*/
}

void MeshParameterization::Parameterization() {
	int nn = mesh.VertexSize();
	Eigen::MatrixXf lambda = Eigen::MatrixXf::Zero(nn, nn);
	for (int i = 0; i < nn; i++) {
		for (int j = 0; j < nn; j++) {
			if (i == j) lambda(i, j) = 0.0;
			else lambda(i, j) = 1.0 / (int)mesh.NeighborVh(MeshKernel::VertexHandle(resorted_id[i])).size();
		}
	}
	Eigen::MatrixXf A = Eigen::MatrixXf::Zero(n, n);
	for (int i = 0; i < n; i++) {
		A(i, i) = 1.0;
	}
	for (int i = 0; i < n; i++) {
		for (auto idx : mesh.NeighborVh(MeshKernel::VertexHandle(resorted_id[i]))) {
			if (!mesh.isOnBoundary(idx)) {
				A(i, reverse_id[idx]) = -lambda(i, reverse_id[idx]);
			}
		}
	}
	Eigen::VectorXf bx = Eigen::VectorXf::Zero(n);
	Eigen::VectorXf by = Eigen::VectorXf::Zero(n);
	Eigen::VectorXf bz = Eigen::VectorXf::Zero(n);
	Eigen::VectorXf solx = Eigen::VectorXf::Zero(n);
	Eigen::VectorXf soly = Eigen::VectorXf::Zero(n);
	Eigen::VectorXf solz = Eigen::VectorXf::Zero(n);
	for (int i = 0; i < n; i++) {
		for (auto idx : mesh.NeighborVh(MeshKernel::VertexHandle(resorted_id[i]))) {
			if (mesh.isOnBoundary(idx)) {
				int j = reverse_id[idx];
				bx(i) += lambda(i, j) * mesh.vertices(idx).x();
				by(i) += lambda(i, j) * mesh.vertices(idx).y();
				bz(i) += lambda(i, j) * mesh.vertices(idx).z();
			}
		}
	}
	solx = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bx);
	soly = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(by);
	solz = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bz);
	for (int i = 0; i < n; i++) {
		mesh.vertices(MeshKernel::VertexHandle(resorted_id[i])).setPosition(solx[i], soly[i], solz[i]);
	}
}
