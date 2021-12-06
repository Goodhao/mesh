#include "bff.h"
//#include <fstream>

#define PI acos(-1)

void bff::Execute() {
	printf("boundary first flattening is starting\n");
	//printf("detecting and filling holes\n");
	fill_holes();
	//printf("holes filled\n");
	Init();
	Laplace();
	InitCurvature();
	Parameterization();
	printf("boundary first flattening success\n");
}

void bff::fill_holes() {
	MeshKernel::VertexHandle s = MeshKernel::VertexHandle(0);
	MeshKernel::VertexHandle t = MeshKernel::VertexHandle(0);
	vis.clear();
	int idx = 0;
	double max_len = -1e6;
	longest_loop_idx = 0;
	for (auto& ep : mesh.alledges()) {
		auto eh = ep.first;
		if (mesh.isOnBoundary(eh) && !vis[eh]) {
			s = ep.second.vh1();
			double loop_len = 0;
			std::set<MeshKernel::EdgeHandle> emp;
			loops.push_back(emp);
			while (1) {
				auto eh = NextEdge(s);
				if ((int)eh > mesh.EdgeSize()) break;
				t = mesh.NeighborVhFromEdge(s, eh);
				loops[idx].insert(eh);
				vis[eh] = 1;
				loop_len += (mesh.vertices(s) - mesh.vertices(t)).norm();
				s = t;
			}
			if (loop_len > max_len) {
				longest_loop_idx = idx;
				max_len = loop_len;
			}
			idx++;
		}
	}
	/*
	std::cout << idx << std::endl;
	for (int i = 0; i < idx; i++) {
		std::cout << i << ": ";
		for (auto x : loops[i]) std::cout << mesh.edges(x).vh1() << " ";
		std::cout << std::endl;
	}
	std::cout << longest_loop_idx << std::endl;*/
	MeshKernel::Vertex cc;
	for (int i = 0; i < idx; i++) if (i != longest_loop_idx) {
		cc = MeshKernel::Vertex(0, 0, 0);
		for (auto x : loops[i]) {
			cc = cc + mesh.vertices(mesh.edges(x).vh1());
		}
		cc = cc / (int)loops[i].size();
		MeshKernel::VertexHandle cch = mesh.AddVertex(cc);
		fake_points.insert(cch);
		for (auto x : loops[i]) {
			mesh.AddFace(std::vector<MeshKernel::VertexHandle>{cch, mesh.edges(x).vh1(), mesh.edges(x).vh2()});
		}
	}
	vis.clear();
}

MeshKernel::EdgeHandle bff::NextEdge(MeshKernel::VertexHandle vh) {
	for (auto& eh : mesh.NeighborEh(vh)) {
		if (mesh.isOnBoundary(eh)) {
			if (!vis[eh]) {
				return eh;
			}
		}
	}
	return MeshKernel::EdgeHandle(mesh.EdgeSize() + 10);
}

double bff::Angle(MeshKernel::VertexHandle vh, MeshKernel::VertexHandle vh1, MeshKernel::VertexHandle vh2) {
	auto d1 = mesh.vertices(vh1) - mesh.vertices(vh);
	auto d2 = mesh.vertices(vh2) - mesh.vertices(vh);
	return acos(d1 * d2 / (d1.norm() * d2.norm()));
}

std::vector<MeshKernel::EdgeHandle> bff::AroundEdge(MeshKernel::VertexHandle vh) {
	std::set<int> s;
	std::map<int, MeshKernel::EdgeHandle> vtoe;
	std::vector<MeshKernel::EdgeHandle> res;
	int start = -1;
	MeshKernel::VertexHandle cur, nxt;
	for (auto& eh : mesh.NeighborEh(vh)) {
		auto vh2 = mesh.NeighborVhFromEdge(vh, eh);
		if (start == -1 || mesh.isOnBoundary(eh)) start = vh2;
		s.insert(vh2);
		vtoe[vh2] = eh;
	}
	cur = MeshKernel::VertexHandle(start);
	while (1) {
		res.push_back(vtoe[cur]);
		s.erase(cur);
		bool ok = 0;
		for (auto& eh : mesh.NeighborEh(cur)) {
			auto p = mesh.NeighborVhFromEdge(cur, eh);
			if (s.count(p)) {
				if (points_to_faces.count(std::make_tuple(int(vh), int(cur), int(p)))) {
					nxt = p;
					ok = 1;
				}
			}
		}
		if (!ok) break;
		cur = nxt;
		if (cur == MeshKernel::VertexHandle(start)) break;
	}
	return res;
}

void bff::Parameterization() {
	// ��һ�� ���㹲������u
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
	solver.compute(Aii);
	// Eigen::MatrixXd ui = Aii.inverse() * -K;
	Eigen::MatrixXd ui = solver.solve(-K);
	Eigen::MatrixXd ub = Eigen::MatrixXd::Zero(bn, 1);
	Eigen::MatrixXd u(n, 1);
	u << ui,
		ub;
	// for (int i = 0; i < in; i++) std::cout << ui(i, 0) << std::endl;
	// �ڶ��� �����±߽��ĽǶ�k_tilde
	Eigen::MatrixXd h = -Abi * ui;
	Eigen::MatrixXd k_tilde = k - h;
	// ������������߽�gamma
	MeshKernel::VertexHandle s = MeshKernel::VertexHandle(0), t = MeshKernel::VertexHandle(0);
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			s = vp.first;
			break;
		}
	}
	MeshKernel::VertexHandle s2 = s; // record

	// ΢������ʹ�ñ߽�պ�
	Eigen::SparseMatrix<double> N(bn, bn), N_inverse(bn, bn);
	Eigen::SparseMatrix<double> T(2, bn);
	Eigen::MatrixXd len = Eigen::MatrixXd::Zero(bn, 1), len_tilde = Eigen::MatrixXd::Zero(bn, 1);

	double angle = 0;
	vis.clear();
	while (1) {
		auto eh = NextEdge(s);
		if ((int)eh > mesh.EdgeSize()) break;
		t = mesh.NeighborVhFromEdge(s, eh);
		vis[eh] = 1;
		angle += k_tilde(to_new[s] - in, 0);
		T.coeffRef(0, to_new[s] - in) = cos(angle), T.coeffRef(1, to_new[s] - in) = sin(angle);
		len(to_new[s] - in, 0) = exp(0.5 * (u(to_new[s], 0) + u(to_new[t], 0))) *
			(mesh.vertices(mesh.edges(eh).vh1()) - mesh.vertices(mesh.edges(eh).vh2())).norm();
		N.coeffRef(to_new[s] - in, to_new[s] - in) = 1 / len(to_new[s] - in, 0);
		N_inverse.coeffRef(to_new[s] - in, to_new[s] - in) = len(to_new[s] - in, 0);
		s = t;
	}
	len_tilde = len - (N_inverse * T.transpose() * Eigen::MatrixXd(T * N_inverse * T.transpose()).inverse() * T) * len;
	vis.clear();
	angle = 0;
	s = s2;
	std::map<int, MeshKernel::Vertex> gamma;
	gamma[s] = MeshKernel::Vertex(0, 0, 0);
	while (1) {
		auto eh = NextEdge(s);
		if ((int)eh > mesh.EdgeSize()) break;
		t = mesh.NeighborVhFromEdge(s, eh);
		vis[eh] = 1;
		angle += k_tilde(to_new[s] - in, 0);
		gamma[t] = gamma[s] + MeshKernel::Vertex(T.coeffRef(0, to_new[s] - in), T.coeffRef(1, to_new[s] - in), 0) * len_tilde(to_new[s] - in, 0);
		s = t;
	}

	//std::cout << "bn:" << bn << " " << gamma.size() << std::endl;
	//std::cout << " total angle:" << angle * 180 / PI << std::endl;
	Eigen::MatrixXd gamma_re = Eigen::MatrixXd::Zero(bn, 1);
	Eigen::MatrixXd gamma_im = Eigen::MatrixXd::Zero(bn, 1);
	for (auto x : gamma) {
		gamma_re(to_new[x.first] - in, 0) = x.second.x();
		gamma_im(to_new[x.first] - in, 0) = x.second.y();
	}
	/*
	std::ofstream fout("data.txt");
	fout << " total angle:" << angle * 180 / PI << std::endl;
	for (auto x : gamma) {
		fout << x.second.x() << " " << x.second.y() << " " << x.second.z() << std::endl;
	}*/

	// ���Ĳ�������gammaΪf
	//Eigen::MatrixXd ai = -Aii.inverse() * Aib * gamma_re;
	//Eigen::MatrixXd bi = -Aii.inverse() * Aib * gamma_im;
	Eigen::MatrixXd ai = solver.solve(-Aib * gamma_re);
	Eigen::MatrixXd bi = solver.solve(-Aib * gamma_im);
	Eigen::MatrixXd a(n, 1);
	a << ai,
		gamma_re;
	Eigen::MatrixXd b(n, 1);
	b << bi,
		gamma_im;
	for (auto x : fake_points) mesh.DeleteVertex(x);
	for (int i = 0; i < n; i++) {
		auto vh = MeshKernel::VertexHandle(to_old[i]);
		mesh.vertices(vh).setPosition(a(i, 0), b(i, 0), 0.0);
	}
}

void bff::InitCurvature() {
	K = Eigen::MatrixXd::Zero(in, 1);
	k = Eigen::MatrixXd::Zero(bn, 1);
	for (auto& vp : mesh.allvertices()) {
		double defect = 0;
		bool first = 1;
		MeshKernel::Vertex last;
		if (!mesh.isOnBoundary(vp.first)) {
			for (auto& eh : AroundEdge(vp.first)) {
				MeshKernel::Vertex cur;
				auto vh2 = mesh.NeighborVhFromEdge(vp.first, eh);
				cur = mesh.vertices(vh2) - vp.second;
				if (!first) {
					defect += acos(cur * last / (cur.norm() * last.norm()));
				}
				last = cur;
				first = 0;
			}
			auto vh2 = mesh.NeighborVhFromEdge(vp.first, *AroundEdge(vp.first).begin());
			auto cur = mesh.vertices(vh2) - vp.second;
			defect += acos(cur * last / (cur.norm() * last.norm()));
			defect = 2 * PI - defect;
			K(to_new[vp.first], 0) = defect;
		}
		else {
			for (auto& eh : AroundEdge(vp.first)) {
				MeshKernel::Vertex cur;
				auto vh2 = mesh.NeighborVhFromEdge(vp.first, eh);
				cur = mesh.vertices(vh2) - vp.second;
				if (!first) {
					defect += acos(cur * last / (cur.norm() * last.norm()));
				}
				last = cur;
				first = 0;
			}
			defect = PI - defect;
			k(to_new[vp.first] - in, 0) = defect;
		}
	}
}

void bff::Init() {
	n = mesh.VertexSize();
	for (auto& vp : mesh.allvertices()) {
		if (!mesh.isOnBoundary(vp.first)) {
			points.push_back(vp.first);
			to_new[vp.first] = points.size() - 1;
			to_old[points.size() - 1] = vp.first;
		}
	}
	in = points.size();
	bn = n - in;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			points.push_back(vp.first);
			to_new[vp.first] = points.size() - 1;
			to_old[points.size() - 1] = vp.first;
		}
	}
	for (auto& fp : mesh.allfaces()) {
		auto &t = fp.second.getVertexHandle();
		std::vector<int> tmp;
		for (auto x : t) tmp.push_back(x);
		assert((int)tmp.size() == 3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) if (i != j) {
				for (int k = 0; k < 3; k++) if (k != i && k != j) {
					points_to_faces[std::make_tuple(tmp[i], tmp[j], tmp[k])] = fp.first;
				}
			}
		}
	}
}

double bff::CotanAngles(MeshKernel::EdgeHandle eh) {
	double cot = 0;
	for (auto& fh : mesh.NeighborFh(eh)) {
		for (auto& vh : mesh.faces(fh).getVertexHandle()) {
			auto vh1 = mesh.edges(eh).vh1(), vh2 = mesh.edges(eh).vh2();
			if (vh != vh1 && vh != vh2) {
				auto v = mesh.vertices(vh), v1 = mesh.vertices(vh1), v2 = mesh.vertices(vh2);
				cot += 1.0 / tan(acos((v - v1) * (v - v2) / (v - v1).norm() / (v - v2).norm()));
			}
		}
	}
	return cot;
}

void bff::Laplace() {
	A.resize(n, n);
	Aii.resize(in, in);
	Aib.resize(in, bn);
	Abi.resize(bn, in);
	Abb.resize(bn, bn);
	for (auto& ep : mesh.alledges()) {
		auto& e = ep.second;
		auto eh = ep.first;
		auto vh1 = to_new[e.vh1()], vh2 = to_new[e.vh2()];
		A.coeffRef(vh1, vh2) = A.coeffRef(vh2, vh1) = -0.5 * CotanAngles(eh);
		A.coeffRef(vh1, vh1) += -A.coeffRef(vh1, vh2);
		A.coeffRef(vh2, vh2) += -A.coeffRef(vh1, vh2);
	}
	Aii = A.block(0, 0, in, in);
	Aib = A.block(0, in, in, bn);
	Abi = A.block(in, 0, bn, in);
	Abb = A.block(in, in, bn, bn);
}