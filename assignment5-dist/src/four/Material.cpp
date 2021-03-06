#include "material.hpp"

#include "lights.hpp"
#include "VecUtils.h"

using namespace FW;

Vec3f PhongMaterial::shade(const Ray &ray, const Hit &hit, 
		const Vec3f &dir_to_light, 
		const Vec3f &incident_intensity,
		bool shade_back) const
{
	// YOUR CODE HERE (R4)
	// Ambient light was already dealt with in R1; implement the rest
	// of the Phong reflectance model (diffuse and specular) here.
	// Start with diffuse and add specular when diffuse is working.
	// NOTE: if shade_back flag is set,
	// you should treat the material as two-sided, i.e., flip the
	// normal if the ray hits the surface from behind. Otherwise
	// you should return zero for hits coming from behind the surface.
	// Remember, when computing the specular lobe, you shouldn't add
	// anything if the light is below the local horizon!
	Vec3f answer = Vec3f(0.0f);
	Vec3f hit_point = ray.pointAtParameter(hit.t);
	Vec3f norm = hit.normal;
	Vec3f pointToEye = (ray.origin - hit_point).normalized();

	if (shade_back) {
		if (dot(pointToEye, norm) < 0.0f) {
			norm = -1.0f * norm;
		}
	}
	if (dot(dir_to_light, hit.normal) < 0.0f) {
		return answer;
	}
	// diffuse
	answer += incident_intensity * FW::clamp(dot(dir_to_light, norm), 0.0f, 1.0f) * diffuse_color(hit_point);
	// specular
	Vec3f r_vec = - pointToEye - 2 * dot(-pointToEye, norm) * norm;
	r_vec.normalize();
	answer += incident_intensity * specular_color_ * FW::pow(FW::clamp(dot(r_vec, dir_to_light), 0.0f, 1.0f), exponent_);
	return answer;
}

Vec3f ProceduralMaterial::diffuse_color(const Vec3f& point) const {
	Vec3f a1 = m1_->diffuse_color(point);
	Vec3f a2 = m2_->diffuse_color(point);
	Vec3f pt = VecUtils::transformPoint(matrix_, point);
	float v = interpolation(pt);
	return a1 * v + a2 * (1 - v);
}

Vec3f ProceduralMaterial::reflective_color(const Vec3f& point) const {
	Vec3f a1 = m1_->reflective_color(point);
	Vec3f a2 = m2_->reflective_color(point);
	Vec3f pt = VecUtils::transformPoint(matrix_, point);
	float v = interpolation(pt);
	return a1 * v + a2 * (1 - v);
}

Vec3f ProceduralMaterial::transparent_color( const Vec3f& point) const {
	Vec3f a1 = m1_->transparent_color(point);
	Vec3f a2 = m2_->transparent_color(point);
	Vec3f pt = VecUtils::transformPoint(matrix_, point);
	float v = interpolation(pt);
	return a1 * v + a2 * (1 - v);
}

float ProceduralMaterial::refraction_index(const Vec3f& point) const {
	float a1 = m1_->refraction_index(point);
	float a2 = m2_->refraction_index(point);
	Vec3f pt = VecUtils::transformPoint(matrix_, point);
	float v = interpolation(pt);
	return a1 * v + a2 * (1 - v);
}

Vec3f ProceduralMaterial::shade(const Ray &ray, const Hit &hit, 
		const Vec3f &dirToLight, 
		const Vec3f &lightColor,
		bool shade_back) const {
	Vec3f pt = ray.pointAtParameter(hit.t);
	Vec3f a1 = m1_->shade(ray,hit,dirToLight,lightColor, shade_back);
	Vec3f a2 = m2_->shade(ray,hit,dirToLight,lightColor, shade_back);
	Vec3f pt2 = VecUtils::transformPoint(matrix_, pt);
	float v = interpolation(pt2);
	Vec3f answer = a1 * v + a2 * (1 - v);
	return answer;
}

float Checkerboard::interpolation(const Vec3f& point) const {
	int count = 1;
	count += int(point[0]) % 2;
	count += int(point[1]) % 2;
	count += int(point[2]) % 2;
	if(point[0] < 0) count ++;
	if(point[1] < 0) count ++;
	if(point[2] < 0) count ++;
	return float(count & 0x01);
}

