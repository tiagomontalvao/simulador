#include <cstdio>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
using namespace std;

// global coisos
random_device rd;
mt19937 mt;

// global constants
const double p1 = 0.3;
const double p2 = 0.1;
const double p3 = 0.3;
const double p = 0.3;

// distributions
normal_distribution<double> normal;
exponential_distribution<double> expo;
poisson_distribution<int> poisson;
uniform_real_distribution<double> unif(0, 1);

int getPacketSize() {
  double sample = unif(mt);
  if (sample < p1) return 64;
  if (sample < p1 + 448*p/1436) return 64 + int(1436*(sample-p1)/p);
  if (sample < p1 + p2 + 448*p/1436) return 512;
  if (sample < 1 - p3) return 512 + int(1436*(sample-(p1+p2+448*p/1436))/p);
  return 1500;
}

int main() {
  mt = mt19937(rd());
  normal_distribution<double> normal;
  map<int,int> hist;
  for (int i = 0; i < 50000; ++i) hist[getPacketSize()]++;
  for (auto x: hist) printf("%d: %d\n", x.first, x.second);
  return 0;
}
