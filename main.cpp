#include <cstdio>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <string>
#include <vector>

using namespace std;

// global coisos
random_device rd;
mt19937 mt(rd());

// global constants
constexpr double p1 = 0.3;
constexpr double p2 = 0.1;
constexpr double p3 = 0.3;
constexpr double p = 0.3;
constexpr int NO_CHANNELS = 30;
constexpr int VOICE_PACKET_SIZE = 64;

// distributions
uniform_real_distribution<double> unif(0, 1);
exponential_distribution<double> time_between_data_packets;
exponential_distribution<double> time_between_voice_packet_groups(1.0/650);
geometric_distribution<int> no_voice_packets(1.0/22);

int getPacketSize() {
  double sample = unif(mt);
  if (sample < p1) return 64;
  if (sample < p1 + 448*p/1436) return 64 + int(1436*(sample-p1)/p);
  if (sample < p1 + p2 + 448*p/1436) return 512;
  if (sample < 1 - p3) return 512 + int(1436*(sample-(p1+p2+448*p/1436))/p);
  return 1500;
}

enum Packet_type {DATA, VOICE};

struct Event {
  uint64_t packet_size;
  uint64_t arrival;
  uint64_t channel;
  uint64_t vpackets;
  uint64_t idx;
  Packet_type p_type;

  // TODO: fix later
  bool operator<(const Event &rhs) const {
    if (p_type != rhs.p_type)
      return p_type == VOICE;
    else
      return arrival < rhs.arrival;
  }

  // Constructor for data events
  Event(double sim_time): p_type(DATA) {
    packet_size = getPacketSize();
    arrival = sim_time + time_between_data_packets(mt);
  }

  Event(uint64_t channel, double sim_time): channel(channel), p_type(VOICE) {
    packet_size = VOICE_PACKET_SIZE;
    arrival = sim_time + time_between_voice_packet_groups(mt);
    vpackets = no_voice_packets(mt);
    idx = 0;
  }
};

void run_simulation(int round_size, double rho) {
  double sim_time = 0;
  double mean = 755/2e6;
  priority_queue<Event> events;

  time_between_data_packets = exponential_distribution<double>(rho/mean);

  events.push(Event(sim_time));
  for (int i = 0; i < NO_CHANNELS; ++i) {
    events.push(Event(i, sim_time));
  }

  int round_idx = 0;

  while (round_idx++ < round_size) {
    Event e = events.top();
    events.pop();

    sim_time = e.arrival + e.size / ;

    // process event

    // create new event of the same type
  }
}

int main(int argc, char **argv) {
  int rounds, round_size;

  for (int i = 0; i < rounds; ++i) {
    run_simulation(round_size, 0.1);
  }

  return 0;
}
