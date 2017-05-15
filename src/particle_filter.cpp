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

#include "particle_filter.h"

//Output log files for each step
//#define ENABLE_LOG

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  num_particles = 10;

	// Add random Gaussian noise to each particle.
  //
  default_random_engine gen;
  double std_x, std_y, std_theta;
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];
  
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  
  weights = vector<double>();
  
  particles = vector<Particle>();
  
  for(int i=0;i<num_particles;i++)
  {
    Particle p;
    p.id = i;
    p.theta = dist_theta(gen);
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.weight = 1.0f;
    particles.push_back(p);
    weights.push_back(p.weight);
  }

  write("log_init.txt");
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  default_random_engine gen;
  double std_x, std_y, std_theta;
  std_x = std_pos[0];
  std_y = std_pos[1];
  std_theta = std_pos[2];
  
  //Generate 0 mean distributions to add to each particle dimension
  normal_distribution<double> dist_x(0, std_x);
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> dist_theta(0, std_theta);

  for(auto&& p : particles)
  {
    if(fabs(yaw_rate)>.0001)
    {
      p.x = p.x + (velocity/yaw_rate)*(sin(p.theta+yaw_rate*delta_t)-sin(p.theta));
      p.y = p.y + (velocity/yaw_rate)*(cos(p.theta)-cos(p.theta+yaw_rate*delta_t));
      p.theta = p.theta + yaw_rate*delta_t;
      
      p.x = p.x + dist_x(gen);
      p.y = p.y + dist_y(gen);
      p.theta = p.theta + dist_theta(gen);
    }
    else
    {
      p.x += velocity*delta_t*cos(p.theta)+dist_x(gen);
      p.y += + velocity*delta_t*sin(p.theta)+dist_y(gen);
      p.theta += dist_theta(gen);
    }
  }

  write("log_predict.txt");
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  
  

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

  //pi/2 1.5708 offset because 0 angle on coordinate system with 0 angle vehicle results in pi/2 measurements
  //for targets in front of the vehicle. 
  //ignore the above comment, it works as x = x*cos - y*sin + xp, y = x*sin+y*cos + yp
  double offset = 0;
  double weight_sum = 0;
  //for each particle, translate observations into map locations
  for(auto&& p : particles)
  {
    std::vector<LandmarkObs> map_observations = std::vector<LandmarkObs>();
    for(auto const obs : observations)
    {
      LandmarkObs l;
      l.x = p.x + obs.x*cos(p.theta-offset)-obs.y*sin(p.theta-offset);
      l.y = p.y + obs.x*sin(p.theta-offset)+obs.y*cos(p.theta-offset);
    //correlate the measurement to a landmark on the map.  The Map contains a list of all the landmarks
      l.id = -1;
      double lastDist = sensor_range + 1;
      for(auto const lm : map_landmarks.landmark_list)
      {
        double d = dist(lm.x_f,lm.y_f,l.x,l.y);
        if(d<lastDist)
        {
          lastDist = d;
          l.id=lm.id_i;
        }
      }
      map_observations.push_back(l);
    }

    //use the multi-variate gaussian dist to determine the liklihood of the match
    //update the weight with that
    double scaler = 1/(2*3.1416*std_landmark[0]*std_landmark[1]);
    double product = 1;
    for(auto&& obs : map_observations)
    {
      double xval = pow((obs.x-map_landmarks.landmark_list[obs.id-1].x_f),2)/(2*std_landmark[0]*std_landmark[0]);
      double yval = pow((obs.y-map_landmarks.landmark_list[obs.id-1].y_f),2)/(2*std_landmark[1]*std_landmark[1]);

      double w = scaler * exp(-1*(xval+yval));
      product = product * w;
    }
    //set particle weight to product
    p.weight = product;
    //synch weights array with particle weights
    weights[p.id] = product;
    
    weight_sum += product;
  }
  
  if(weight_sum>0)
  {
    for(auto&& w : weights)
    {
      w = w/weight_sum;
    }
  }
  
  write("log_updateWeights.txt");
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  default_random_engine gen;
  discrete_distribution<> dd_weights(weights.begin(),weights.end());
  
  vector<Particle> particles2(num_particles);
  
  for(int i =0; i<num_particles;i++)
  {
    int indx = dd_weights(gen);
    particles2[i] = particles[indx];
    particles2[i].id = i;
  }
  
  particles = particles2;

  write("log_resample.txt");
}

void ParticleFilter::write(std::string filename) {
  
#ifdef ENABLE_LOG
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
  dataFile<<"Particles\n";
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
  dataFile<<"Weights\n";
  for (int i = 0; i <num_particles; i++)
  {
    dataFile<<weights[i]<<"\n";
  }
	dataFile.close();
#endif
}
