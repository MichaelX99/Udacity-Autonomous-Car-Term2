/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 1000;
	default_random_engine gen;

	for (int i = 0; i < num_particles; i++)
	{
		normal_distribution<double> dist_x(x, std[0]);
		normal_distribution<double> dist_y(y, std[1]);
		normal_distribution<double> dist_theta(theta, std[2]);

		Particle temp;
		temp.id = i;
		temp.x = dist_x(gen);
		temp.y = dist_y(gen);
		temp.theta = dist_theta(gen);
		temp.weight = 1;

		particles.push_back(temp);
		weights.push_back(1);

		is_initialized = true;
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	if (fabs(yaw_rate) < .0001)
	{
		yaw_rate = 1;
	}

	double coef = velocity / yaw_rate;
	double move = yaw_rate * delta_t;

	double x, y, theta, new_x, new_y, new_theta;
	for (int i = 0; i < num_particles; i++)
	{
		x = particles[i].x;
		y = particles[i].y;
		theta = particles[i].theta;

		normal_distribution<double> dist_x(0, std_pos[0]);
		normal_distribution<double> dist_y(0, std_pos[1]);
		normal_distribution<double> dist_theta(0, std_pos[2]);

		new_x = x + coef * (sin(theta + move) - sin(theta));
		new_y = y + coef * (cos(theta) - cos(theta + move));
		new_theta = theta + move;

		particles[i].x = new_x + dist_x(gen);
		particles[i].y = new_y + dist_y(gen);
		particles[i].theta = new_theta + dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	double best_distance, check;
	int best_id;
	for (int i = 0; i < observations.size(); i++)
	{
		best_distance = 100000000000;
		best_id = 0;
		for (int j = 0; j < predicted.size(); j++)
		{
			check = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
			if (check < best_distance)
			{
				best_distance = check;
				best_id = predicted[j].id;
			}
		}
		observations[i].id = best_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	double sigma = std_landmark[0] * std_landmark[1];
	double coef = 1.0 / (6.28 * sigma);
	double xsqr = pow(std_landmark[0], 2);
	double ysqr = pow(std_landmark[1], 2);


	// Preallocation
	vector<LandmarkObs> transformed_observations, nearest_landmark;
	double x, y, theta, curr_dist, weight_init, pred_x, pred_y, first, second;
	Map::single_landmark_s current_landmark;
	LandmarkObs temp_observation, temp_landmark;

	for (int i = 0; i < particles.size(); i++)
	{

			for (int j = 0; j < observations.size(); j++)
			{
				temp_observation.x = particles[i].x + observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta);
				temp_observation.y = particles[i].y + observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta);
				temp_observation.id = observations[j].id;

				transformed_observations.push_back(temp_observation);
			}

			// Find nearest map landmark
			for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
			{
					current_landmark = map_landmarks.landmark_list[j];

					curr_dist = dist(particles[i].x, particles[i].y, current_landmark.x_f, current_landmark.y_f);

					if (curr_dist < sensor_range)
					{
							temp_landmark = {current_landmark.id_i, current_landmark.x_f, current_landmark.y_f};
							nearest_landmark.push_back(temp_landmark);

					}
			}

			// Call data association function
			dataAssociation(nearest_landmark, transformed_observations);

			weight_init = 1.0;
			
			// Update weights
			for (int k = 0; k < transformed_observations.size(); k++)
			{
					for (int j = 0; j < nearest_landmark.size(); j++)
					{
							if (nearest_landmark[j].id == transformed_observations[k].id)
							{
									pred_x = nearest_landmark[j].x;
									pred_y = nearest_landmark[j].y;
									break;
							}
					}

					first = -pow((transformed_observations[k].x - pred_x), 2)/(2*xsqr);
					second = -pow(transformed_observations[k].y - pred_y, 2)/(2*ysqr);

					weight_init = weight_init * coef * exp(first + second);

			}

			particles[i].weight = weight_init;
			weights[i] = weight_init;
			transformed_observations.clear();
			nearest_landmark.clear();
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// Found here:
	// https://discussions.udacity.com/t/resampling-algorithm-using-resampling-wheel/241313/12
	random_device rd;
  mt19937 gen(rd());
  discrete_distribution<> d(weights.begin(), weights.end());
	vector<Particle> new_particles;
	int i;
  for (i=0; i < num_particles; i++)
	{
  	new_particles.push_back(particles[d(gen)]);
  }
	particles.clear();
	for (i = 0; i < num_particles; i++)
	{
		particles.push_back(new_particles[i]);
	}
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
