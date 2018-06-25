#include <cstdio>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <string>
#include <vector>
#include <list>

using namespace std;

// global coisos
random_device rd;
mt19937 mt(rd());

// global constants
constexpr double p1 = 0.3;
constexpr double p2 = 0.1;
constexpr double p3 = 0.3;
constexpr double p = 0.3;
constexpr double SPEED = 2e3; // bytes per ms
constexpr double MEAN_DATA_PACKET_SIZE = 755;
constexpr double MEAN_DATA_SERVICE_TIME = MEAN_DATA_PACKET_SIZE / SPEED;
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
  int packet_size;
  int arrival;
  int channel;
  int vpackets;
  int idx;
  int round_idx;
  Packet_type p_type;

  // TODO: fix later (uh... fixed?)
  bool operator<(const Event &rhs) const {
      return arrival < rhs.arrival;
  }

  // Constructor for data events
  Event(double sim_time): p_type(DATA) {
    packet_size = getPacketSize();
    arrival = sim_time + time_between_data_packets(mt);
  }

  // Constructor for voice events
  Event(int channel, double sim_time): channel(channel), p_type(VOICE) {
    packet_size = VOICE_PACKET_SIZE;
    arrival = sim_time + time_between_voice_packet_groups(mt);
    vpackets = no_voice_packets(mt);
    idx = 0;
  }
};

void run_simulation(int transient_phase_size, int rounds, int round_size, double rho, bool interrupt) {
  double sim_time = 0;
  list<Event> events;
  time_between_data_packets = exponential_distribution<double>(rho/MEAN_DATA_SERVICE_TIME);


  // Initialize events queue
  events.emplace_back(sim_time);
  for (int i = 0; i < NO_CHANNELS; ++i) {
    events.emplace_back(i, sim_time);
  }

  int events_processed = 0;

  while (events_processed < round_size * rounds) {
    Event &e = events.front();

    // without interrupt
    double next_sim_time = sim_time + e.packet_size / SPEED;

    // process event
    if (e.p_type == DATA) {
      auto new_event = events.front();
      for (auto &event: events) {
        if (event.p_type == VOICE && event.arrival < next_sim_time) {
          swap(e, event);
          // new_event = event;
          break;
        }
      }


		  events.emplace_back(e.arrival);
  	} else {

  	}

    // create new event of the same type
    //~ events.pop();
    events.sort();
  }
}

int main(int argc, char **argv) {
  int transient_phase_size, rounds, round_size;

  for (double rho = 0.1; rho <= 0.71; rho += 0.1) {
    run_simulation(transient_phase_size, rounds, round_size, rho, false);
    run_simulation(transient_phase_size, rounds, round_size, rho, true);
  }

  return 0;
}
