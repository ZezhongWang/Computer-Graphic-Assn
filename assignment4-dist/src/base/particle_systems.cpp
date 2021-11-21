#include "particle_systems.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>

using namespace std;
using namespace FW;

namespace {

	inline Vec3f fGravity(float mass) {
		return Vec3f(0, -9.8f * mass, 0);
	}

	// force acting on particle at pos1 due to spring attached to pos2 at the other end
	inline Vec3f fSpring(const Vec3f& pos1, const Vec3f& pos2, float k, float rest_length) {
		// YOUR CODE HERE (R2)
		auto d = pos1 - pos2;
		return -k * (d.length() - rest_length) * d / d.length();
	}

	inline Vec3f fDrag(const Vec3f& v, float k) {
		return -k * v;
	}

} // namespace

void SimpleSystem::reset() {
	current_state_ = State(1, Vec3f(0, radius_, 0));
}

State SimpleSystem::evalF(const State& state) const {
	State f(1, Vec3f(-state[0].y, state[0].x, 0));
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H
// using the implicit Euler method, the simple system should converge towards origin -- as opposed to the explicit Euler, which diverges outwards from the origin.
void SimpleSystem::evalJ(const State&, SparseMatrix& result, bool initial) const {
	if (initial) {
		result.coeffRef(1, 0) = 1.0f;
		result.coeffRef(0, 1) = -1.0f;
	}
}
#endif

Points SimpleSystem::getPoints() {
	return Points(1, current_state_[0]);
}

Lines SimpleSystem::getLines() {
	static const auto n_lines = 50u;
	auto l = Lines(n_lines * 2);
	const auto angle_incr = 2 * FW_PI / n_lines;
	for (auto i = 0u; i < n_lines; ++i) {
		l[2 * i] = l[2 * i + 1] =
			Vec3f(radius_ * FW::sin(angle_incr * i), radius_ * FW::cos(angle_incr * i), 0);
	}
	rotate(l.begin(), l.begin() + 1, l.end());
	return l;
}

void SpringSystem::reset() {
	const auto start_pos = Vec3f(0.1f, -0.5f, 0.0f);
	const auto spring_k = 30.0f;
	const auto rest_length = 0.5f;
	current_state_ = State(4);
	// YOUR CODE HERE (R2)
	// Set the initial state for a particle system with one particle fixed
	// at origin and another particle hanging off the first one with a spring.
	// Place the second particle initially at start_pos.
	spring_ = Spring(0, 1, spring_k, rest_length);
	current_state_[0] = Vec3f(0, 0, 0);
	current_state_[1] = Vec3f(0, 0, 0);		// start velocity
	current_state_[2] = start_pos;
	current_state_[3] = Vec3f(0, 0, 0);		// start velocity
}

State SpringSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	State f(4);
	// YOUR CODE HERE (R2)
	// Return a derivative for the system as if it was in state "state".
	// You can use the fGravity, fDrag and fSpring helper functions for the forces.
	f[0] = current_state_[0];
	f[1] = 0;		// fixed
	f[2] = current_state_[3];
	auto fG = fGravity(mass);
	auto fD = fDrag(f[2], drag_k);
	auto fS = fSpring(current_state_[spring_.i2 * 2], current_state_[spring_.i1 * 2], spring_.k, spring_.rlen);
	f[3] = (fG + fD + fS) / mass;
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

// This is a very useful read for the Jacobians of the spring forces. It deals with spring damping as well, we don't do that -- our drag is simply a linear damping of velocity (that results in some constants in the Jacobian).
// http://blog.mmacklin.com/2012/05/04/implicitsprings/

void SpringSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	// EXTRA: Evaluate the Jacobian into the 'result' matrix here. Only the free end of the spring should have any nonzero values related to it.
}
#endif

Points SpringSystem::getPoints() {
	auto p = Points(2);
	p[0] = current_state_[0]; p[1] = current_state_[2];
	return p;
}

Lines SpringSystem::getLines() {
	auto l = Lines(2);
	l[0] = current_state_[0]; l[1] = current_state_[2];
	return l;
}

void PendulumSystem::reset() {
	const auto spring_k = 1000.0f;
	const auto start_point = Vec3f(0);
	const auto end_point = Vec3f(0.05, -1.5, 0);
	current_state_ = State(2 * n_);
	// YOUR CODE HERE (R4)
	// Set the initial state for a pendulum system with n_ particles
	// connected with springs into a chain from start_point to end_point with uniform intervals.
	// The rest length of each spring is its length in this initial configuration.
	for (int i = 0; i < n_ - 1; ++i) {
		Spring spring = Spring(i, i+1, spring_k, 1.5f / (n_ - 1));
		springs_.push_back(spring);
		current_state_[pos(i)] = Vec3f(end_point.x / (n_ - 1) * i, end_point.y / (n_ - 1) * i, 0);
		current_state_[vel(i)] = Vec3f(0, 0, 0);
	}
	current_state_[pos(n_ - 1)] = Vec3f(end_point.x, end_point.y, 0);
	current_state_[vel(n_ - 1)] = Vec3f(0, 0, 0);
}

State PendulumSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 0.5f;
	auto f = State(2 * n_);
	// YOUR CODE HERE (R4)
	// As in R2, return a derivative of the system state "state".
	f[0] = f[1] = 0;
	for (int i = 1; i < n_; ++i) {
		f[2 * i] = current_state_[vel(i)];
		// add gravity & drag
		f[2 * i + 1] = (fGravity(mass) + fDrag(current_state_[vel(i)], drag_k)) / mass;
	}

	for (auto& spring : springs_) {
		// add resultant force
		int i1 = spring.i1;
		int i2 = spring.i2;
		auto fS = fSpring(current_state_[pos(i1)], current_state_[pos(i2)], spring.k, spring.rlen) / mass;
		f[2 * i1 + 1] += fS;
		f[2 * i2 + 1] += -fS;
	}
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void PendulumSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {

	const auto drag_k = 0.5f;
	const auto mass = 0.5f;

	// EXTRA: Evaluate the Jacobian here. Each spring has an effect on four blocks of the matrix -- both of the positions of the endpoints will have an effect on both of the velocities of the endpoints.
}
#endif


Points PendulumSystem::getPoints() {
	auto p = Points(n_);
	for (auto i = 0u; i < n_; ++i) {
		p[i] = current_state_[i * 2];
	}
	return p;
}

Lines PendulumSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(current_state_[2 * s.i1]);
		l.push_back(current_state_[2 * s.i2]);
	}
	return l;
}

void ClothSystem::reset() {
	const auto spring_k = 300.0f;
	const auto width = 1.5f, height = 1.5f; // width and height of the whole grid
	current_state_ = State(2 * x_ * y_);
	// YOUR CODE HERE (R5)
	// Construct a particle system with a x_ * y_ grid of particles,
	// connected with a variety of springs as described in the handout:
	// structural springs, shear springs and flex springs.
	// structural springs
	for (int i = 0; i < y_; ++i) {
		for (int j = 0; j < x_; ++j) {
			// horizontal
			int idx = i * x_ + j;
			if (j < x_ - 1) {
				Spring spring = Spring(idx, idx + 1, spring_k, width / (x_ - 1));
				springs_.push_back(spring);
			}
			// vertical
			if (i < y_ - 1) {
				Spring spring = Spring(idx, idx + x_, spring_k, height / (y_ - 1));
				springs_.push_back(spring);
			}
			current_state_[pos(i, j)] = Vec3f(width / (x_ - 1) * j, height / (y_ - 1) * i, 0);
			current_state_[vel(i, j)] = Vec3f(0, 0, 0);
		}
	}
}

State ClothSystem::evalF(const State& state) const {
	const auto drag_k = 0.08f;
	const auto n = x_ * y_;
	static const auto mass = 0.025f;
	auto f = State(2 * n);
	// YOUR CODE HERE (R5)
	// This will be much like in R2 and R4.
	for (int i = 0; i < y_; ++i) {
		for (int j = 0; j < x_; ++j) {
			// add gravity & drag
			int idx = i * x_ + j;
			f[2 * idx] = current_state_[vel(i, j)];
			f[2 * idx + 1] = (fGravity(mass) + fDrag(current_state_[vel(i, j)], drag_k)) / mass;
		}
	}

	for (auto& spring : springs_) {
		// add resultant force
		int i1 = spring.i1;
		int i2 = spring.i2;
		auto fS = fSpring(current_state_[2 * i1], current_state_[2 * i2], spring.k, spring.rlen) / mass;
		f[2 * i1 + 1] += fS;
		f[2 * i2 + 1] += -fS;
	}

	// holding points
	f[0] = f[1] = 0;
	f[2 * (x_ - 1)] = 0;// f[2 * (x_ - 1) + 1] = 0;

	for (int x = 0; x < x_; ++x) {
		for (int y = 0; y < y_; ++y) {
			int idx = y * x_ + x;
			if (f[2 * idx + 1].isZero()) {
				std::cout << "is Zero: " << idx <<"\t";
			}
		}
	}
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void ClothSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.08f;
	static const auto mass = 0.025f;

	// EXTRA: Evaluate the Jacobian here. The code is more or less the same as for the pendulum.
}

#endif

Points ClothSystem::getPoints() {
	auto n = x_ * y_;
	auto p = Points(n);
	for (auto i = 0u; i < n; ++i) {
		p[i] = current_state_[2 * i];
	}
	return p;
}

Lines ClothSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(current_state_[2 * s.i1]);
		l.push_back(current_state_[2 * s.i2]);
	}
	return l;
}
State FluidSystem::evalF(const State&) const {
	return State();
}

