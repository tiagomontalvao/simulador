#include <algorithm>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <string>

// #define NDEBUG
#include <cassert>

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
    // And yes, this is useless for VOICE events
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
    bool operator==(const Event &rhs) const;
    Event();
    Event(int channel);

};

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
/*----------  SIMULATION LOGIC  ----------*/
void run_simulation(int transient_phase_size, int rounds, int round_size, double rho, bool interrupt);
void handle_arrival(Event &cur_event, bool interrupt);
void handle_data_arrival(Event &cur_event);
void handle_voice_arrival(Event &cur_event, bool interrupt);
void handle_dispatch(Event &cur_event);
void handle_data_dispatch(Event &cur_event);
void handle_voice_dispatch(Event &cur_event);
void serve_next_packet();
void enter_the_server(const Event &cur_event, bool force=false);
void replace_data_dispatch(const Event &cur_event);
/*----------  EVENT CREATION  ----------*/
Event make_next_data_arrival(const Event &event);
Event make_data_dispatch(const Event &event);
Event make_next_voice_arrival(const Event &event);
Event make_voice_dispatch(const Event &event);
/*----------  MISC  ----------*/
bool is_dispatch(const Event &event);
/*----------  DEBUG & LOG  ----------*/
void log (const Event &event, bool verbose=true, const string &prefix = "");
void dbg_show_queue(bool verbose=true);

/*============================
=            MAIN            =
============================*/

int main(int argc, char **argv) {
    int transient_phase_size, rounds, round_size;
    ios_base::sync_with_stdio(false);

    transient_phase_size = 500;
    rounds = 100;
    round_size = 50000;

    for (double rho = 0.4; rho <= 0.71; rho += 0.1) {
        run_simulation(transient_phase_size, rounds, round_size, rho, true);
        break;
        run_simulation(transient_phase_size, rounds, round_size, rho, true);
    }

    return 0;
}

/*======================================================
=            RANDOM NUMBER GENERATION SETUP            =
======================================================*/

random_device rd;
mt19937 mt(rd());

// distributions
uniform_real_distribution<double> unif(0, 1);
exponential_distribution<double> time_between_data_packets;
exponential_distribution<double> time_between_voice_packet_groups(1.0/650);
geometric_distribution<int> voice_group_size(1.0/22);

int get_packet_size() {
    double sample = unif(mt);
    if (sample < p1) return 64;
    if (sample < p1 + 448*p/1436) return 64 + round(1436*(sample-p1)/p);
    if (sample < p1 + p2 + 448*p/1436) return 512;
    if (sample < 1 - p3) return 512 + round(1436*(sample-(p1+p2+448*p/1436))/p);
    return 1500;
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
    last_t1 = last_t2 = 0;
}

void event_queue_init() {
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
        Event cur_event = *event_queue.begin();
        event_queue.erase(event_queue.begin());
        
        sim_t = cur_event.t;

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
        event_queue.insert(make_next_data_arrival(cur_event));
    } else {
        handle_voice_arrival(cur_event, interrupt);
        event_queue.insert(make_next_voice_arrival(cur_event));
    }
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
    if (!voice_queue.empty()) {
        enter_the_server(voice_queue.front());

        tmpW2 += sim_t - voice_queue.front().arrival_t;

        voice_queue.pop_front();
    } else if (!data_queue.empty()) {
        enter_the_server(data_queue.front());

        tmpW1 += sim_t - data_queue.front().queue_t;

        data_queue.pop_front();
    } else {
        idle = true;
    }
}

inline void enter_the_server(const Event &cur_event, bool force) {
    assert(cur_event.type == ARRIVAL);

    if (cur_event.packet_type == DATA) {
        event_queue.insert(make_data_dispatch(cur_event));
    } else if (!force) {
        event_queue.insert(make_voice_dispatch(cur_event));
    } else {
        replace_data_dispatch(cur_event);
    }

    event_on_server = cur_event;
    idle = false;
}

void replace_data_dispatch(const Event &cur_event) {
    assert(cur_event.type == ARRIVAL && cur_event.packet_type == VOICE);

    auto it = find_if(event_queue.begin(), event_queue.end(), is_dispatch);
    assert(it != event_queue.end() && it->type == DISPATCH && it->packet_type == DATA);

    tmpX1 += it->packet_size / SPEED - (it->t - sim_t);

    event_queue.erase(it);
    event_queue.insert(make_voice_dispatch(cur_event));

    Nq1 += (sim_t - last_t1) * data_queue.size();
    last_t1 = sim_t;

    event_on_server.queue_t = sim_t;
    data_queue.push_front(event_on_server);
}

/*======================================
=            EVENT CREATION            =
======================================*/

// Constructs first data arrival event (assumes simulation clock is at 0)
Event::Event(): type(ARRIVAL) {
    t = time_between_data_packets(mt);
    queue_t = arrival_t = t;
    packet_type = DATA;
    packet_size = get_packet_size();
    id = id_counter++;
}

// Constructs first voice arrival event (assumes simulation clock is at 0)
Event::Event(int channel): type(ARRIVAL), channel(channel) {
    t = time_between_voice_packet_groups(mt);
    packet_type = VOICE;
    packet_size = VOICE_PACKET_SIZE;
    vgroup_size = voice_group_size(mt) + 1;

    arrival_t = t;
    vgroup_idx = 0;
    id = id_counter++;
}

Event make_next_data_arrival(const Event &event) {
    Event new_event(event);
    new_event.t += time_between_data_packets(mt);
    new_event.queue_t = new_event.arrival_t = new_event.t;
    new_event.packet_size = get_packet_size();
    new_event.id = id_counter++;

    return new_event;
}

Event make_data_dispatch(const Event &event) {
    Event new_event(event);
    new_event.type = DISPATCH;
    new_event.t = sim_t + event.packet_size / SPEED;

    return new_event;
}

Event make_next_voice_arrival(const Event &event) {
    Event new_event(event);
    new_event.t += TIME_BETWEEN_VOICE_PACKETS;
    new_event.arrival_t = new_event.t;
    new_event.id = id_counter++;

    if (event.vgroup_idx+1 < event.vgroup_size) {
        new_event.vgroup_idx++;
        return new_event;
    }

    new_event.t += time_between_voice_packet_groups(mt);
    new_event.vgroup_size = voice_group_size(mt) + 1;

    new_event.arrival_t = new_event.t;
    new_event.vgroup_idx = 0;

    return new_event;
}

Event make_voice_dispatch(const Event &event) {
    Event new_event(event);
    new_event.type = DISPATCH;
    new_event.t = sim_t + VOICE_SERVICE_TIME;

    return new_event;
}

/*============================
=            MISC            =
============================*/

bool Event::operator<(const Event &rhs) const {
    return t == rhs.t ? (packet_type == VOICE && rhs.packet_type == DATA) : t < rhs.t;
}

bool Event::operator==(const Event &rhs) const {
    return id == rhs.id;
}

bool is_dispatch(const Event &event) {
    return event.type == DISPATCH;
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

    cout << (e.type == ARRIVAL ? "+ " : "- ");
    cout << (e.packet_type == DATA ? "D " : "V ");
    cout << setw(10) << e.t << "ms ";
    cout << setw(5) << e.packet_size << "B";

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
