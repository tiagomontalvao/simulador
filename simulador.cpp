#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <string>

using namespace std;

/*========================================
=            TYPE DEFINITIONS            =
========================================*/

enum Packet_type {DATA, VOICE};

enum Event_type {ARRIVAL, DISPATCH};

struct Event {
    // Every event should have a <type> and the point in time <t> at which it happens
    Event_type type;
    double t;
    int id;

    // DISPATCH events will keep the time of arrival
    double arrival_t;

    // When interrupted, data packets will record the time at which they returned to the queue
    double queue_t;

    // Variables used for both data and voice packets
    Packet_type packet_type;
    int packet_size;

    // Variables specific to voice packets
    int channel;
    int vgroup_size;
    int vgroup_idx;

    // Member functions prototypes
    bool operator<(const Event &rhs) const;
    Event();
    Event(int channel);
};

bool Event::operator<(const Event &rhs) const {
    if (t != rhs.t) return t < rhs.t;
    if (type != rhs.type) return type == DISPATCH;
    if (packet_type != rhs.packet_type) return packet_type == VOICE;
    return id < rhs.id;
}

/*=========================================
=            PROBLEM CONSTANTS            =
=========================================*/

constexpr double p1 = 0.3;
constexpr double p2 = 0.1;
constexpr double p3 = 0.3;
constexpr double p = 0.3;
constexpr double SPEED = 2e3/8; // bytes per ms
constexpr double MEAN_DATA_PACKET_SIZE = 755; // bytes
constexpr int NO_CHANNELS = 30;
constexpr int VOICE_PACKET_SIZE = 64; // bytes
constexpr double MEAN_DATA_SERVICE_TIME = MEAN_DATA_PACKET_SIZE / SPEED; // ms
constexpr double VOICE_SERVICE_TIME = VOICE_PACKET_SIZE / SPEED; // ms
constexpr double TIME_BETWEEN_VOICE_PACKETS = 16; // ms

/*==================================================
=            SIMULATION STATE VARIABLES            =
==================================================*/

double sim_t;
int packets_processed;
set<Event> event_queue;
deque<Event> data_queue;
deque<Event> voice_queue;
Event event_on_server;
bool idle;
int id_counter;

// stats

double T1, W1, X1, Nq1;
double T2, W2, X2, Nq2;

double tmpX1;
double tmpW1, tmpW2;

double total_time;
int data_packets_processed;
int voice_packets_processed;
int last_t1, last_t2;

/*===========================================
=            FUNCTION PROTOTYPES            =
===========================================*/

/*----------  RANDOM NUMBER GENERATION SETUP  ----------*/
int get_packet_size();
/*----------  INITIALIZATION  ----------*/
void simulation_state_init(double rho);
void event_queue_init();
/*----------  MAIN  ----------*/
/*----------  SIMULATION LOGIC  ----------*/
void run_simulation(int transient_phase_size, int rounds, int round_size, double rho, bool interrupt);
void handle_arrival(Event &cur_event, bool interrupt);
void handle_data_arrival(Event &cur_event);
void handle_voice_arrival(Event &cur_event, bool interrupt);
void handle_dispatch(Event &cur_event);
void handle_data_dispatch(Event &cur_event);
void handle_voice_dispatch(Event &cur_event);
void serve_next_packet();
void enter_the_server(const Event &event, bool force=false);
void unschedule_data_dispatch();
/*----------  EVENT CREATION  ----------*/
Event make_arrival_from(const Event &event);
Event make_dispatch_from(const Event &event);
/*----------  DEBUG & LOG  ----------*/
void log (const Event &event, bool verbose=true, const string &prefix = "");
void dbg_show_queue(bool verbose=true);

/*======================================================
=            RANDOM NUMBER GENERATION SETUP            =
======================================================*/

random_device rd;
mt19937 mt(rd());

// distributions
uniform_real_distribution<double> unif(0, 1);
exponential_distribution<double> time_between_data_packets;
exponential_distribution<double> time_between_voice_groups(1.0/650);
geometric_distribution<int> voice_group_size(1.0/22);

int get_packet_size() {
    double sample = unif(mt);
    if (sample < p1) return 64;
    if (sample < p1 + 448*p/1436) return 64 + round(1436*(sample-p1)/p);
    if (sample < p1 + p2 + 448*p/1436) return 512;
    if (sample < 1 - p3) return 512 + round(1436*(sample-(p1+p2+448*p/1436))/p);
    return 1500;
}

/*============================
=            MAIN            =
============================*/

int main(int argc, char **argv) {
    int transient_phase_size, rounds, round_size;
    ios_base::sync_with_stdio(false);

    transient_phase_size = 500;
    rounds = 1000;
    round_size = 500;

    for (double rho = 0.1; rho <= 0.71; rho += 0.1) {
        run_simulation(transient_phase_size, rounds, round_size, rho, false);
        run_simulation(transient_phase_size, rounds, round_size, rho, true);
    }

    return 0;
}

/*======================================
=            INITIALIZATION            =
======================================*/

void simulation_state_init(double rho) {
    sim_t = 0;
    event_queue.clear();
    data_queue.clear();
    voice_queue.clear();
    idle = true;
    time_between_data_packets = exponential_distribution<double>(rho/MEAN_DATA_SERVICE_TIME);
    id_counter = 0;
}

void statistics_init() {
    T1 = W1 = X1 = Nq1 = 0;
    T2 = W2 = X2 = Nq2 = 0;

    tmpX1 = 0;
    tmpW1 = tmpW2 = 0;

    packets_processed = 0;
    total_time = 0;
    data_packets_processed = voice_packets_processed = 0;

    // Last time the data/voice queue changed size
    last_t1 = last_t2 = 0;
}

void event_queue_init() {
    assert(event_queue.empty());
    event_queue.emplace();

    for (int i = 0; i < NO_CHANNELS; ++i) {
        event_queue.emplace(i);
    }
}

/*========================================
=            SIMULATION LOGIC            =
========================================*/

void run_simulation(int transient_phase_size, int rounds, int round_size, double rho, bool interrupt) {
    // Simulation state initialization
    simulation_state_init(rho);

    // Initialize incremental sums and counters
    statistics_init();

    // Initialize event queue
    event_queue_init();

    int n = rounds * round_size;

    while (packets_processed < n) {
        // Get next event and then remove it from event queue
        Event cur_event = *event_queue.begin();
        event_queue.erase(event_queue.begin());

        // Update simulation time
        sim_t = cur_event.t;

        // Deal with the event
        if (cur_event.type == ARRIVAL)
            handle_arrival(cur_event, interrupt);
        else
            handle_dispatch(cur_event);
    }

    cout << defaultfloat;
    cout << (interrupt ? "Com" : "Sem") << " interrupção. ρ₁ = " << rho << endl;
    cout << fixed << setprecision(6);
    cout << "E[X₁]:  " << X1/data_packets_processed << endl;
    cout << "E[W₁]:  " << W1/data_packets_processed << endl;
    cout << "E[T₁]:  " << T1/data_packets_processed << endl;
    cout << "E[Nq₁]: " << Nq1/total_time << endl;
    // cout << "λ₁:     " << (Nq1/total_time) / (W1/data_packets_processed) << endl;
    // cout << "ρ₁:     " << (Nq1/total_time) / (W1/data_packets_processed) * (X1/data_packets_processed) << endl;
    // cout << "E[Nq1] = lambda1 * E[W1] = " << (rho/MEAN_DATA_SERVICE_TIME) * (W1/data_packets_processed) << endl;
    cout << endl;
    cout << "E[X₂]:  " << X2/voice_packets_processed << endl;
    cout << "E[W₂]:  " << W2/voice_packets_processed << endl;
    cout << "E[T₂]:  " << T2/voice_packets_processed << endl;
    cout << "E[Nq₂]: " << Nq2/total_time << endl;
    // cout << "λ₂:     " << (Nq2/total_time) / (W2/voice_packets_processed) << endl;
    // cout << "ρ₂:     " << (Nq2/total_time) / (W2/voice_packets_processed) * (X2/voice_packets_processed) << endl;
    cout << endl << defaultfloat;

}

void handle_arrival(Event &cur_event, bool interrupt) {
    assert(cur_event.type == ARRIVAL);

    if (cur_event.packet_type == DATA) {
        handle_data_arrival(cur_event);
    } else {
        handle_voice_arrival(cur_event, interrupt);
    }

    // Before the current arrival event vanishes,
    // a new arrival event is created from it
    event_queue.insert(make_arrival_from(cur_event));
}

void handle_data_arrival(Event &cur_event) {
    assert(cur_event.type == ARRIVAL && cur_event.packet_type == DATA);

    Nq1 += (sim_t - last_t1) * data_queue.size();
    last_t1 = sim_t;

    if (!idle) {
        cur_event.queue_t = sim_t;
        data_queue.push_back(cur_event);
    } else {
        enter_the_server(cur_event);
    }
}

void handle_voice_arrival(Event &cur_event, bool interrupt) {
    assert(cur_event.type == ARRIVAL && cur_event.packet_type == VOICE);

    Nq2 += (sim_t - last_t2) * voice_queue.size();
    last_t2 = sim_t;

    if (!idle && interrupt && event_on_server.packet_type == DATA) {
        enter_the_server(cur_event, true);
    } else if (idle){
        enter_the_server(cur_event);
    } else {
        voice_queue.push_back(cur_event);
    }
}

void handle_dispatch(Event &cur_event) {
    assert(cur_event.type == DISPATCH);

    packets_processed++;
    total_time = sim_t;

    if (cur_event.packet_type == DATA)
        handle_data_dispatch(cur_event);
    else
        handle_voice_dispatch(cur_event);

    serve_next_packet();
}

void handle_data_dispatch(Event &cur_event) {
    assert(cur_event.type == DISPATCH && cur_event.packet_type == DATA);

    data_packets_processed++;

    X1 += cur_event.packet_size / SPEED + tmpX1;
    tmpX1 = 0;

    W1 += tmpW1;
    tmpW1 = 0;

    T1 += sim_t - cur_event.arrival_t;

    Nq1 += (sim_t - last_t1) * data_queue.size();
    last_t1 = sim_t;
}

void handle_voice_dispatch(Event &cur_event) {
    assert(cur_event.type == DISPATCH && cur_event.packet_type == VOICE);
    voice_packets_processed++;

    X2 += cur_event.packet_size / SPEED;

    W2 += tmpW2;
    tmpW2 = 0;

    T2 += sim_t - cur_event.arrival_t;

    Nq2 += (sim_t - last_t2) * voice_queue.size();
    last_t2 = sim_t;
}

inline void serve_next_packet() {
    // Server is set to idle until a packet manages to get in
    idle = true;

    // Serve next packet waiting in line (voice packets first)
    if (!voice_queue.empty()) {
        enter_the_server(voice_queue.front());
        tmpW2 += sim_t - voice_queue.front().arrival_t;
        voice_queue.pop_front();
    } else if (!data_queue.empty()) {
        enter_the_server(data_queue.front());
        tmpW1 += sim_t - data_queue.front().queue_t;
        data_queue.pop_front();
    }
}

inline void enter_the_server(const Event &event, bool force) {
    assert(event.type == ARRIVAL);

    // Voice packets may enter the server by force when interruption is on
    if (force) {
        assert(!idle && event.packet_type == VOICE && event.type == ARRIVAL);
        unschedule_data_dispatch();
    }

    // Schedule the dispatch of the entering event
    event_queue.insert(make_dispatch_from(event));

    // Update server state after a packet enters
    event_on_server = event;
    idle = false;
}

void unschedule_data_dispatch() {
    // Find the DISPATCH event in the event queue
    auto it = find_if(event_queue.begin(), event_queue.end(),
        [] (const Event &e) {
            return e.type == DISPATCH;
        });

    assert(it != event_queue.end() && it->type == DISPATCH && it->packet_type == DATA);

    // Update the partial service time for the interrupted data packet
    tmpX1 += it->packet_size / SPEED - (it->t - sim_t);

    // Also update the product Nq1 * time for the interrupted data packet,
    // since its return to the queue will increase the queue size
    Nq1 += (sim_t - last_t1) * data_queue.size();
    last_t1 = sim_t;

    // Return interrupted data packet to the front of the queue.
    event_on_server.queue_t = sim_t;
    data_queue.push_front(event_on_server);

    // Remove it from the event queue
    event_queue.erase(it);
}

/*======================================
=            EVENT CREATION            =
======================================*/

// Constructs first data arrival event (assumes simulation clock is at 0)
Event::Event(): type(ARRIVAL) {
    arrival_t = t = time_between_data_packets(mt);
    packet_type = DATA;
    packet_size = get_packet_size();
    id = id_counter++;
}

// Constructs first voice arrival event (assumes simulation clock is at 0)
Event::Event(int channel): type(ARRIVAL), channel(channel) {
    arrival_t = t = time_between_voice_groups(mt);
    packet_type = VOICE;
    packet_size = VOICE_PACKET_SIZE;
    vgroup_size = voice_group_size(mt) + 1;
    vgroup_idx = 0;
    id = id_counter++;
}

Event make_arrival_from(const Event &event) {
    assert(event.type == ARRIVAL);
    Event new_event(event);

    if (event.packet_type == VOICE) {
        // At least 16ms must pass, whatever the case is
        new_event.t += TIME_BETWEEN_VOICE_PACKETS;

        // If this isn't the last event in the voice group,
        // the next arrival happens in the same group...
        if (event.vgroup_idx+1 < event.vgroup_size) {
            new_event.vgroup_idx++;
        } else {
            // ...If not, the first arrival in a voice group starts
            // after an exponentially distributed silence period
            new_event.t += time_between_voice_groups(mt);
            new_event.vgroup_size = voice_group_size(mt) + 1;
            new_event.vgroup_idx = 0;
        }

    } else {
        // A new data arrival happens at an exponential
        // time after the last arrival
        new_event.t += time_between_data_packets(mt);
        new_event.packet_size = get_packet_size();
    }

    // Backup arrival time
    new_event.arrival_t = new_event.t;

    // Every event is assigned a unique id
    new_event.id = id_counter++;

    return new_event;
}

Event make_dispatch_from(const Event &event) {
    assert(event.type == ARRIVAL);
    Event new_event(event);

    new_event.type = DISPATCH;

    if (event.packet_type == VOICE) {
        new_event.t = sim_t + VOICE_SERVICE_TIME;
    } else {
        new_event.t = sim_t + event.packet_size / SPEED;
    }

    return new_event;
}

/*===================================
=            DEBUG & LOG            =
===================================*/

void log (const Event &e, bool verbose, const string &prefix) {
    std::cout << fixed;
    std::cout << setprecision(2);

    if (!verbose) {
        cout << (e.type == ARRIVAL ? "+" : "-");
        cout << (e.packet_type == DATA ? "D    " : "V");
        if (e.packet_type == VOICE) {
            cout << "#" << left << setw(2) << e.channel << " ";
        }

        return;
    }

    cout << prefix;

    cout << setw(10) << e.t << ": ";
    cout << (e.type == ARRIVAL ? "+" : "-");
    cout << (e.packet_type == DATA ? "D " : "V ");
    cout << left << setw(7) << e.id;
    cout << right << setw(5) << e.packet_size << "B";

    if (e.packet_type == VOICE) {
        cout << setw(3) << " #" << setw(2) << e.channel << " ";
        cout << e.vgroup_idx+1 << "/" << e.vgroup_size;
    }

    cout << endl;
}

void dbg_show_queue(bool verbose) {
    cout << endl;
    if (!verbose) cout << "{ ";
    int k = 0;
    for (const auto &e: event_queue) {
        if (!verbose && k == 16) cout << endl << "  ";
        log(e, verbose);
        k++;
    }
    if (!verbose) cout << "}";
    cout << "\ndata queue size : " << data_queue.size() << endl;
    cout << "voice queue size: " << voice_queue.size() << endl;

    cout << endl;
}
