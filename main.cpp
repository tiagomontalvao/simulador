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

// Simulation state variables
double sim_time;
int Nq1;
int Nq2;

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
enum Event_type {ARRIVAL, DISPATCH, INTERRUPTION};

struct Event {
	Event_type type;
	Packet_type packet_type;
	double arrival;
	int packet_size;
	int channel;
	int vpackets;
	int vgroup_idx;
	int round_idx;

	// TODO: fix later (uh... fixed?)
	bool operator<(const Event &rhs) const {
			return arrival < rhs.arrival;
	}

	// Constructor for data arrival events
	Event(double sim_time): packet_type(DATA) {
		packet_size = getPacketSize();
		arrival = sim_time + time_between_data_packets(mt);
	}

	// Constructor for voice arrival events
	Event(int channel, double sim_time): channel(channel), packet_type(VOICE) {
		packet_size = VOICE_PACKET_SIZE;
		arrival = sim_time + time_between_voice_packet_groups(mt);
		vpackets = no_voice_packets(mt);
		vgroup_idx = 0;
	}
};

void events_queue_init(list<Event> &events) {
	events.emplace_back(sim_time);
	
	for (int i = 0; i < NO_CHANNELS; ++i) {
		events.emplace_back(i, sim_time);
	}

	events.sort();
}

void dbg_show_queue(list<Event> &events) {
	cout << "{";
	for (auto &e: events) {
		cout << "  (" << (e.packet_type == DATA ? "Data" : "Voice") << ", " << e.packet_size << " B";
		if (e.packet_type == VOICE)
			cout << " on #" << e.channel;

		cout << " at " << e.arrival << " ms)\n";
	}
	cout << "}" << endl;
}

void run_simulation(int transient_phase_size, int rounds, int round_size, double rho, bool interrupt) {
	sim_time = 0;
	list<Event> events;
	time_between_data_packets = exponential_distribution<double>(rho/MEAN_DATA_SERVICE_TIME);

	// Initialize events queue
	events_queue_init(events);

	int events_processed = 0;
	int n = 10000;
	// dbg_show_queue(events);
	
	while (events_processed < n) {
		Event &e = events.front();

		// without interrupt
		double next_sim_time = sim_time + e.packet_size / SPEED;

		// process event
		auto next_voice = events.end();
		if (e.packet_type == DATA) {
			for (auto event = events.begin(); event != events.end(); ++event) {
				if (event->packet_type == VOICE && event->arrival < next_sim_time) {
					next_voice = event;
					// new_event = event;
					break;
				}
			}

			if (next_voice != events.end()) {
				dbg_show_queue(events);
				events.push_front(*next_voice);
				events.erase(next_voice);
				dbg_show_queue(events);
			}
		// events.emplace_back(e.arrival);
		} else {

		}

		// create new event of the same type
		//~ events.pop();
		// events.sort();
		// create_next_arrival();


		events_processed++;
	}
}

int main(int argc, char **argv) {
	int transient_phase_size, rounds, round_size;

	transient_phase_size = rounds = round_size = 1;
	for (double rho = 0.1; rho <= 0.71; rho += 0.1) {
		run_simulation(transient_phase_size, rounds, round_size, rho, false);
		run_simulation(transient_phase_size, rounds, round_size, rho, true);
	}

	return 0;
}
