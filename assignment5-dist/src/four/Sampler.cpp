#include "Sampler.h"

#include <cassert>

Sampler::Sampler(int nSamples) :
	m_nSamples(nSamples)
{}

Sampler* Sampler::constructSampler(Args::SamplePatternType t, int numsamples)
{
	if( t == Args::Pattern_Uniform ) {
		return new UniformSampler(numsamples);
	} else if ( t == Args::Pattern_Regular ) {
		return new RegularSampler(numsamples);
	} else if ( t == Args::Pattern_Jittered ) {
		return new JitteredSampler(numsamples);
	} else {
		assert(false && "Bad sampler type");
		return nullptr;
	}
}

UniformSampler::UniformSampler(int nSamples) :    
	Sampler(nSamples)
{}

Vec2f UniformSampler::getSamplePosition(int i) {
	// make the compiler shut up by referencing the variable
	(void)i;	
	// YOUR CODE HERE (R9)
	// Return a uniformly distributed random 2-vector within the unit square [0,1]^2

	return Vec2f((rand() * 1.0) / (RAND_MAX * 1.0), (rand() * 1.0) / (RAND_MAX * 1.0));
}

RegularSampler::RegularSampler(int nSamples) :
	Sampler(nSamples)
{
	// test that it's a perfect square
	int dim = (int)sqrtf(float(nSamples));
	assert(nSamples == dim * dim);
	m_dim = dim;
}

Vec2f RegularSampler::getSamplePosition(int n) {
	// YOUR CODE HERE (R9)
	// Return a sample through the center of the Nth subpixel.
	// The starter code only supports one sample per pixel.
	int r = n / m_dim;
	int c = n % m_dim;
	float step = 1.0f / m_dim;
	float offset = step / 2.0f;
	return Vec2f(offset + (c * step), offset + (r * step));
}

JitteredSampler::JitteredSampler(int nSamples) :
	Sampler(nSamples)
{
	// test that it's a perfect square
	int dim = (int)sqrtf(float(nSamples));
	assert(nSamples == dim * dim);
	m_dim = dim;
}

Vec2f JitteredSampler::getSamplePosition(int n) {
	// YOUR CODE HERE (R9)
	// Return a randomly generated sample through Nth subpixel.
	int row = n / m_dim;
	int col = n % m_dim;
	float step = 1.0f / m_dim;
	float offset = step / 2.0f;
	float x = (rand() * 1.0) / (RAND_MAX * 1.0);
	return Vec2f(x * step + (col * step), x * step + (row * step));
}

