#include <algorithm>
#include <iomanip>
#include <iostream>
#include <list>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace std;

/*========================================
=            TYPE DEFINITIONS            =
========================================*/

enum Packet_type {DATA, VOICE};

enum Event_type {ARRIVAL, DISPATCH};

struct Event {
    // Every event should have a <type> and the point in time <t> at which it happens.
    Event_type type;
    double t;

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
list<Event> event_queue;
list<Event> data_queue;
list<Event> voice_queue;
Event *event_on_server;
double last_voice_event_t;
double Nq2;
double L;
long long data_packets;
double total_time;
int packets_processed;
bool idle;
double X2;

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
/*----------  EVENT CREATION  ----------*/
Event make_next_data_arrival(const Event &event);
Event make_data_dispatch(const Event &event);
Event make_next_voice_arrival(const Event &event);
Event make_voice_dispatch(const Event &event);
/*----------  MISC  ----------*/
bool is_dispatch(const Event &event);
/*----------  DEBUG & LOG  ----------*/
void log (Event &event, const string &prefix = "");
void dbg_show_queue();

/*============================
=            MAIN            =
============================*/

int main(int argc, char **argv) {
    int transient_phase_size, rounds, round_size;

    transient_phase_size = rounds = round_size = 1;

    run_simulation(transient_phase_size, rounds, round_size, 0.7, false);
    return 0;

    for (double rho = 0.7; rho <= 0.71; rho += 0.1) {
        run_simulation(transient_phase_size, rounds, round_size, rho, false);
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
geometric_distribution<int> no_voice_packets(1.0/22);

int get_packet_size() {
    double sample = unif(mt);
    if (sample < p1) return 64;
    if (sample < p1 + 448*p/1436) return 64 + int(1436*(sample-p1)/p);
    if (sample < p1 + p2 + 448*p/1436) return 512;
    if (sample < 1 - p3) return 512 + int(1436*(sample-(p1+p2+448*p/1436))/p);
    return 1500;
}

/*======================================
=            INITIALIZATION            =
======================================*/

void simulation_state_init(double rho) {
    sim_t = 0;
    Nq2 = 0;
    L = 0;
    X2 = 0;
    data_packets = 0;
    total_time = 0;
    last_voice_event_t = 0;
    event_queue.clear();
    data_queue.clear();
    voice_queue.clear();
    event_on_server = nullptr;
    idle = true;
    time_between_data_packets = exponential_distribution<double>(rho/MEAN_DATA_SERVICE_TIME);
}

void event_queue_init() {
    event_queue.emplace_back();
    
    for (int i = 0; i < NO_CHANNELS; ++i) {
        event_queue.emplace_back(i);
    }

    event_queue.sort();
}

/*========================================
=            SIMULATION LOGIC            =
========================================*/

void run_simulation(int transient_phase_size, int rounds, int round_size, double rho, bool interrupt) {
    // Simulation state initialization
    simulation_state_init(rho);

    // Initialize event queue
    event_queue_init();
    // dbg_show_queue();

    packets_processed = 0;
    int n = 500000;

    while (packets_processed < n) {
        Event cur_event = event_queue.front();
        sim_t = cur_event.t;

        if (cur_event.type == ARRIVAL)
            handle_arrival(cur_event, interrupt);
        else 
            handle_dispatch(cur_event);

        event_queue.pop_front();
        event_queue.sort();
    }

    std::cout << setprecision(12);
    cout << L / data_packets << endl;
}

void handle_arrival(Event &cur_event, bool interrupt) {
    if (cur_event.packet_type == DATA)
        handle_data_arrival(cur_event);
    else
        handle_voice_arrival(cur_event, interrupt);
}

void handle_data_arrival(Event &cur_event) {
    // cout << "data" << endl;

    L += cur_event.packet_size;
    data_packets++;

    if (event_on_server) {
        // cout << "busy" << endl;
        data_queue.push_back(cur_event);
    } else {
        // cout << "idle" << endl;
        // event_queue.emplace_back(DISPATCH, DATA, sim_t + cur_event.packet_size / SPEED);
        event_queue.push_back(make_data_dispatch(cur_event));
        event_on_server = &cur_event;
    }

    // event_queue.emplace_back(sim_t);
    event_queue.push_back(make_next_data_arrival(cur_event));

    // cout << data_queue.size() << endl;
}

void handle_voice_arrival(Event &cur_event, bool interrupt) {
    if (!voice_queue.empty()) {
        Nq2 += (cur_event.t - last_voice_event_t) * (double)voice_queue.size();
        // cout << cur_event.t << " - " << last_voice_event_t << " = " << cur_event.t - last_voice_event_t << endl;
    }
    last_voice_event_t = cur_event.t;
    total_time = cur_event.t;

    if (event_on_server) {
        // cout << "busy" << endl;
        if (event_on_server->packet_type == DATA) {
            // cout << " with data" << endl;
            if (interrupt){
                // cout << "interruption" << endl;
                data_queue.push_front(*event_on_server);
                // event_queue.remove_if(is_dispatch);
                auto it = find_if(event_queue.begin(), event_queue.end(), is_dispatch);
                if (it != event_queue.end())
                    event_queue.erase(it);
                // event_queue.emplace_back(DISPATCH, VOICE, sim_t + VOICE_SERVICE_TIME);
                event_queue.push_back(make_voice_dispatch(cur_event));
                event_on_server = &cur_event;
            } else {
                // cout << "no interruption" << endl;
                voice_queue.push_back(cur_event);
            }

        } else {
            // cout << "with data" << endl;
            voice_queue.push_back(cur_event);
        }
    } else {
        // cout << "idle" << endl;
        // event_queue.emplace_back(DISPATCH, VOICE, sim_t + VOICE_SERVICE_TIME);
        event_queue.push_back(make_voice_dispatch(cur_event));
        event_on_server = &cur_event;
    }

    event_queue.push_back(make_next_voice_arrival(cur_event));

}

void handle_dispatch(Event &cur_event) {
    packets_processed++;

    if (cur_event.packet_type == DATA)
        handle_data_dispatch(cur_event);
    else
        handle_voice_dispatch(cur_event);
}

void handle_data_dispatch(Event &cur_event) {
    if (!voice_queue.empty()) {
        // cout << "voice waiting" << endl;
        // event_queue.emplace_back(DISPATCH, VOICE, sim_t + VOICE_SERVICE_TIME);
        event_queue.push_back(make_voice_dispatch(voice_queue.front()));
        event_on_server = &voice_queue.front();
        voice_queue.pop_front();
    } else if (!data_queue.empty()) {
        // cout << "data waiting" << endl;
        // event_queue.emplace_back(DISPATCH, DATA, sim_t + data_queue.front().packet_size / SPEED);
        event_queue.push_back(make_data_dispatch(data_queue.front()));
        event_on_server = &data_queue.front();
        data_queue.pop_front();
    } else {
        // cout << "idle" << endl;
        event_on_server = nullptr;
    }

    // cout << data_queue.size() << endl;
}

void handle_voice_dispatch(Event &cur_event) {
    // cout << "voice" << endl;
    if (!voice_queue.empty()) {
        Nq2 += (cur_event.t - last_voice_event_t) * voice_queue.size();
        // cout << cur_event.t << " - " << last_voice_event_t << " = " << cur_event.t - last_voice_event_t << endl;
    }
    last_voice_event_t = cur_event.t;
    total_time = cur_event.t;

    if (!voice_queue.empty()) {
        // cout << "voice waiting" << endl;
        // event_queue.emplace_back(DISPATCH, VOICE, sim_t + VOICE_SERVICE_TIME);
        event_queue.push_back(make_voice_dispatch(voice_queue.front()));
        event_on_server = &voice_queue.front();
        voice_queue.pop_front();
    } else if (!data_queue.empty()) {
        // cout << "data waiting" << endl;
        // event_queue.emplace_back(DISPATCH, DATA, sim_t + data_queue.front().packet_size / SPEED);
        event_queue.push_back(make_data_dispatch(data_queue.front()));
        event_on_server = &data_queue.front();
        data_queue.pop_front();
    } else {
        // cout << "idle" << endl;
        event_on_server = nullptr;
    }
}

/*======================================
=            EVENT CREATION            =
======================================*/

// Constructs first data arrival event (assumes simulation clock is at 0)
Event::Event(): type(ARRIVAL) {
    t = time_between_data_packets(mt);
    packet_type = DATA;
    packet_size = get_packet_size();
}

// Constructs first voice arrival event (assumes simulation clock is at 0)
Event::Event(int channel): type(ARRIVAL), channel(channel) {
    t = time_between_voice_packet_groups(mt);
    packet_type = VOICE;
    packet_size = VOICE_PACKET_SIZE;
    vgroup_size = no_voice_packets(mt);
    while (!vgroup_size) {
        t += time_between_voice_packet_groups(mt);
        vgroup_size = no_voice_packets(mt);
    }

    vgroup_idx = 0;
}

Event make_next_data_arrival(const Event &event) {
    Event new_event(event);
    new_event.t = event.t + time_between_data_packets(mt);
    new_event.packet_size = get_packet_size();

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

    if (event.vgroup_idx+1 < event.vgroup_size) {
        new_event.vgroup_idx++;
        return new_event;
    }

    new_event.t += time_between_voice_packet_groups(mt);
    new_event.vgroup_size = no_voice_packets(mt);
    while (!new_event.vgroup_size) {
        new_event.t += time_between_voice_packet_groups(mt);
        new_event.vgroup_size = no_voice_packets(mt);
    }

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
    return t < rhs.t;
}

bool is_dispatch(const Event &event) {
    return event.type == DISPATCH;
}

/*===================================
=            DEBUG & LOG            =
===================================*/

void log (Event &e, const string &prefix) {
    std::cout << fixed;
    std::cout << setprecision(2);

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

void dbg_show_queue() {
    cout << endl;
    for (auto &e: event_queue) {
        log(e);
    }
    cout << endl;
}
