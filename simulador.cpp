#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <deque>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

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

    // A unique id for each arrival. Used to guarantee std::set won't ignore events
    int id;

    // Identifies the round the event belongs to
    int round_id;

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

struct Metrics {
    // The round this set of metrics belongs to
    int round_id;

    // Number of packets to process
    int target;

    // The instant data collection started
    // < 0 if it hasn't started
    double start_t;

    // The instant data collection stopped
    // < 0 if it hasn't stopped
    double end_t;

    double T1, W1, X1, Nq1;
    double T2, W2, X2, Nq2, Delta, Delta_sq, Vdelta;

    double tmpX1;
    double tmpW1;

    int packets_processed;
    int data_packets_processed;
    int voice_packets_processed;
    int delta_intervals;
    int last_t1, last_t2;

    Metrics() {}

    Metrics(int rnd, int tgt);
    void init();
    void update_Nq_sum(const Event &event);
    void update_W_sum(const Event &event);
    void update_on_dispatch(const Event &cur_event);
    void update_on_interruption(const Event &event);
    Metrics operator+(const Metrics &rhs) const;
    Metrics operator/(int rhs) const;

private:
    bool should_collect(const Event &event);
    void compute_estimators();
};

Metrics::Metrics(int rnd, int tgt):
    round_id(rnd), target(tgt), start_t(-1), end_t(-1),
    T1(0), W1(0), X1(0), Nq1(0), T2(0), W2(0), X2(0), Nq2(0),
    Delta(0), Delta_sq(0), Vdelta(0), tmpX1(0), tmpW1(0),
    packets_processed(0), data_packets_processed(0),
    voice_packets_processed(0), delta_intervals(0), last_t1(0), last_t2(0) {}

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
set<Event> event_queue;
deque<Event> data_queue;
deque<Event> voice_queue;
Event event_on_server;
bool idle;
int id_counter;
int cur_round;
double last_voice_dispatch_t[NO_CHANNELS];
vector<Metrics> round_metrics;

/*===========================================
=            FUNCTION PROTOTYPES            =
===========================================*/

/*----------  RANDOM NUMBER GENERATION SETUP  ----------*/
int get_packet_size();
/*----------  INITIALIZATION  ----------*/
void simulation_state_init(double rho);
void statistics_init();
void event_queue_init();
/*----------  MAIN  ----------*/
/*----------  SIMULATION LOGIC  ----------*/
void run_simulation(int warmup_period, int rounds, int round_size, double rho, bool interrupt);
void handle_arrival(Event &cur_event, bool interrupt);
void handle_data_arrival(Event &cur_event);
void handle_voice_arrival(Event &cur_event, bool interrupt);
void serve_next_packet();
void enter_the_server(const Event &event, bool force=false);
void unschedule_data_dispatch();
/*----------  EVENT CREATION  ----------*/
Event make_arrival_from(const Event &event);
Event make_dispatch_from(const Event &event);
/*----------  STATISTICS  ----------*/
// See class Metrics definition
/*----------  DEBUG & LOG  ----------*/
void log (const Event &event, bool verbose=true, const string &prefix = "");
void dbg_show_queue(bool verbose=true);

/*======================================================
=            RANDOM NUMBER GENERATION SETUP            =
======================================================*/

random_device rd;
mt19937 mt(rd());

// distributions
uniform_real_distribution<double> unif(0.0, 1.0);
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

int main(int argc, char *argv[]) {
    int warmup_period, rounds, round_size, interrupt;
    double rho;
    ios_base::sync_with_stdio(false);

    if (argc == 6) {
        sscanf(argv[1], "%d", &warmup_period);
        sscanf(argv[2], "%d", &rounds);
        sscanf(argv[3], "%d", &round_size);
        sscanf(argv[4], "%lf", &rho);
        sscanf(argv[5], "%d", &interrupt);

        run_simulation(warmup_period, rounds, round_size, rho, interrupt);
    } else {
        warmup_period = 5e2;
        rounds = 50;
        round_size = 1e5;
        rho = 0.1;

        run_simulation(warmup_period, rounds, round_size, rho, false);
    }

    
    // for (double rho = 0.1; rho <= 0.71; rho += 0.1) {
        // run_simulation(warmup_period, rounds, round_size, rho, true);
    // }

    return 0;
}

/*======================================
=            INITIALIZATION            =
======================================*/

void simulation_state_init(int warmup_period, int rounds, int round_size, double rho) {
    sim_t = 0;
    event_queue.clear();
    data_queue.clear();
    voice_queue.clear();
    idle = true;
    time_between_data_packets = exponential_distribution<double>(rho/MEAN_DATA_SERVICE_TIME);
    id_counter = 0;

    // position 'rounds' is for the warmup period
    cur_round = warmup_period ? rounds : 0;

    round_metrics.resize(rounds+1); // rounds + warmup
    for (int i = 0; i < rounds; ++i) {
        round_metrics[i] = Metrics(i, round_size);
    }
    round_metrics[rounds] = Metrics(rounds, warmup_period);

    if (warmup_period)
        round_metrics[rounds].start_t = 0;
    else
        round_metrics[0].start_t = 0;
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

#define PRINT(X) cout << #X << " " << X << " ";

void run_simulation(int warmup_period, int rounds, int round_size, double rho, bool interrupt) {
    // Simulation state initialization
    simulation_state_init(warmup_period, rounds, round_size, rho);

    // Initialize event queue
    event_queue_init();

    // Probably useless
    int n = 0;

    while (round_metrics[rounds-1].end_t < 0) {
        // Get next event and then remove it from event queue
        Event cur_event = *event_queue.begin();
        event_queue.erase(event_queue.begin());

        // Update simulation time
        sim_t = cur_event.t;

        // Deal with the event
        if (cur_event.type == ARRIVAL) {
            // New arrivals belong to the current round
            cur_event.round_id = cur_round;

            round_metrics[cur_round].update_Nq_sum(cur_event);
            handle_arrival(cur_event, interrupt);
        } else {
            round_metrics[cur_round].update_on_dispatch(cur_event);
            
            if (cur_event.packet_type == VOICE)
                last_voice_dispatch_t[cur_event.channel] = sim_t;
            
            serve_next_packet();
            n++;
        }

        // Go to the next round, if the current one is finished
        if (round_metrics[cur_round].end_t >= 0) {
            cur_round = (cur_round == rounds) ? 0 : cur_round+1;
            round_metrics[cur_round].init();
        }
    }

    // #ifdef PYTHON_OUT
    for (auto rm: round_metrics) {
        if (rm.round_id == rounds) continue;
        PRINT(rm.W1);
        PRINT(rm.X1);
        PRINT(rm.T1);
        PRINT(rm.Nq1);
        PRINT(rm.W2);
        PRINT(rm.X2);
        PRINT(rm.T2);
        PRINT(rm.Nq2);
        PRINT(rm.Delta);
        PRINT(rm.Vdelta);
        cout << endl;
    }
    // #endif

    // auto print_metrics = [&rounds] (const Metrics &rm) {
    //     cout << fixed << setprecision(6);
    //     cout << "E[W1]:    " << setw(15) << left << (rm.W1 ? rm.W1 : 0);
    //     cout << "E[W2]:    " << (rm.W2 ? rm.W2 : 0) << endl;
    //     cout << "E[X1]:    " << setw(15) << left << (rm.X1 ? rm.X1 : 0);
    //     cout << "E[X2]:    " << (rm.X2 ? rm.X2 : 0) << endl;
    //     cout << "E[T1]:    " << setw(15) << left << (rm.T1 ? rm.T1 : 0);
    //     cout << "E[T2]:    " << (rm.T2 ? rm.T2 : 0) << endl;
    //     cout << "E[Nq1]:   " << setw(15) << left << (rm.Nq1 ? rm.Nq1 : 0);
    //     cout << "E[Nq2]:   " << (rm.Nq2 ? rm.Nq2 : 0) << endl;
    //     cout << defaultfloat;
    // };

    // auto print_round = [&rounds, &print_metrics] (const Metrics &rm, bool warmup=false) {
    //     cout << defaultfloat;
    //     if (warmup)
    //         cout << "Warm-up round";
    //     else
    //         cout << "Round " << rm.round_id+1 << "/" << rounds;
    //     cout << fixed << setprecision(2);
    //     cout << " (" << rm.start_t << " --> " << rm.end_t << "):\n";
    //     print_metrics(rm);
    //     // cout << "#:        " << rm.packets_processed << endl;
    // };

    // print_round(round_metrics[rounds], true);
    // for (auto rm: round_metrics) {
    //     if (rm.round_id == rounds) continue;
    //     print_round(rm);
    // }
    // cout << endl;

    // cout << fixed << setprecision(2);
    // cout << (interrupt ? "Com" : "Sem") << " interrupção. Rho1 = " << rho << ":" << endl;
    // print_metrics(accumulate(round_metrics.begin(), round_metrics.end()-1, Metrics(0, 0)) / rounds);
    // cout << n << " pacotes processados" << endl << endl;

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

    if (!idle) {
        data_queue.push_back(cur_event);
    } else {
        enter_the_server(cur_event);
    }
}

void handle_voice_arrival(Event &cur_event, bool interrupt) {
    assert(cur_event.type == ARRIVAL && cur_event.packet_type == VOICE);

    if (!idle && interrupt && event_on_server.packet_type == DATA) {
        enter_the_server(cur_event, true);
    } else if (idle){
        enter_the_server(cur_event);
    } else {
        voice_queue.push_back(cur_event);
    }
}

void serve_next_packet() {
    // Server is set to idle until a packet manages to get in
    idle = true;

    // Serve next packet waiting in line (voice packets first)
    if (!voice_queue.empty()) {
        enter_the_server(voice_queue.front());
        voice_queue.pop_front();
    } else if (!data_queue.empty()) {
        enter_the_server(data_queue.front());
        data_queue.pop_front();
    }

    // If someone left the queue and entered the
    // server, update its waiting time
    if (!idle)
        round_metrics[cur_round].update_W_sum(event_on_server);
}

void enter_the_server(const Event &event, bool force) {
    assert(event.type == ARRIVAL);

    // Voice packets may enter the server by force when interruption is on
    // but the data packet on the server must be returned to the queue
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

    round_metrics[cur_round].update_on_interruption(*it);

    // Update time at which the interrupted packet entered the queue
    event_on_server.queue_t = sim_t;
    // Return interrupted data packet to the front of the queue
    data_queue.push_front(event_on_server);

    // Remove it from the event queue
    event_queue.erase(it);
}

/*======================================
=            EVENT CREATION            =
======================================*/

// Constructs first data arrival event (assumes simulation clock is at 0)
Event::Event(): type(ARRIVAL) {
    queue_t = t = time_between_data_packets(mt);
    packet_type = DATA;
    packet_size = get_packet_size();
    id = id_counter++;
}

// Constructs first voice arrival event (assumes simulation clock is at 0)
Event::Event(int channel): type(ARRIVAL), channel(channel) {
    t = time_between_voice_groups(mt);
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
        new_event.queue_t = new_event.t;
    }

    // Every arrival event is assigned a unique id
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

/*==================================
=            STATISTICS            =
==================================*/

Metrics Metrics::operator+(const Metrics &b) const {
    auto a = Metrics(*this);
    a.T1 += b.T1;
    a.W1 += b.W1;
    a.X1 += b.X1;
    a.Nq1 += b.Nq1;
    a.T2 += b.T2;
    a.W2 += b.W2;
    a.X2 += b.X2;
    a.Nq2 += b.Nq2;
    a.Delta += b.Delta;
    a.Vdelta += b.Vdelta;

    return a;
}

Metrics Metrics::operator/(int b) const {
    auto a = Metrics(*this);
    a.T1 /= b;
    a.W1 /= b;
    a.X1 /= b;
    a.Nq1 /= b;
    a.T2 /= b;
    a.W2 /= b;
    a.X2 /= b;
    a.Nq2 /= b;
    a.Delta /= b;
    a.Vdelta /= b;

    return a;
}

void Metrics::init() {
    if (start_t >= 0) return;
    start_t = sim_t;
    last_t1 = last_t2 = sim_t;
}

bool Metrics::should_collect(const Event &event) {
    // Returns true if:
    // 1) Metrics for this round are not finished
    // 2) These are not metrics for the warm-up period
    // 3) The event belongs to this round
    return end_t < 0 && round_id != (int)round_metrics.size()-1 && round_id == event.round_id;
}

void Metrics::compute_estimators() {
    double total_time = end_t - start_t;

    if (data_packets_processed) {
        T1 /= data_packets_processed;
        W1 /= data_packets_processed;
        X1 /= data_packets_processed;
    } else {
        T1 = W1 = X1 = 0;
    }
    Nq1 /= total_time;

    if (delta_intervals > 1) {
        Vdelta = Delta_sq/(delta_intervals-1) - (Delta*Delta) / ((double)delta_intervals * (delta_intervals-1));
        Delta /= delta_intervals;
    } else {
        Delta = Vdelta = 0;
    }

    if (voice_packets_processed) {
        T2 /= voice_packets_processed;
        W2 /= voice_packets_processed;
        X2 /= voice_packets_processed;
    } else {
        T2 = W2 = X2 = 0;
    }
    Nq2 /= total_time;
}

// Updates W cumulative sums
// Should be called only when a packet leaves the queue
// and enters the server.
void Metrics::update_W_sum(const Event& event) {
    if (!should_collect(event)) return;

    if (event.packet_type == DATA) {
        // Because data packets may never leave the system
        // during the round, a temporary variable is used
        tmpW1 += sim_t - event.queue_t;
    } else {
        // If a voice packet enters the server,
        // it cannot be interrupted, thus W2 can
        // be directly written to
        W2 += sim_t - event.t;
    }
}

// Updates the cumulative product Nq x t
// Should be called on every arrival, dispatch and interruption.
void Metrics::update_Nq_sum(const Event& event) {
    if (!should_collect(event)) return;

    if (event.packet_type == DATA) {
        Nq1 += (sim_t - last_t1) * data_queue.size();
        last_t1 = sim_t;
    } else {
        Nq2 += (sim_t - last_t2) * voice_queue.size();
        last_t2 = sim_t;
    }
}

// Updates sums on dispatches.
void Metrics::update_on_dispatch(const Event& cur_event) {
    assert(cur_event.type == DISPATCH && cur_event.t == sim_t);

    if (!should_collect(cur_event)) {
        // In case this is the warm-up period,
        // count every packet so the round can finish
        if (round_id == (int)round_metrics.size()-1) {
            packets_processed++;
            if (packets_processed == target) {
                end_t = sim_t;
            }
        }
        return;
    }

    if (cur_event.packet_type == DATA) {
        data_packets_processed++;

        X1 += cur_event.packet_size / SPEED + tmpX1;
        tmpX1 = 0;

        W1 += tmpW1;
        tmpW1 = 0;

        T1 += sim_t - event_on_server.t;

    } else {
        voice_packets_processed++;

        X2 += VOICE_SERVICE_TIME;

        T2 += sim_t - event_on_server.t;

        if (cur_event.vgroup_idx > 0) {
            double Dj = sim_t - last_voice_dispatch_t[cur_event.channel];
            Delta += Dj;
            Delta_sq += Dj * Dj;
            delta_intervals++;
        }
    }

    update_Nq_sum(cur_event);

    packets_processed++;

    if (packets_processed == target) {
        end_t = sim_t;
        compute_estimators();
    }

}

// Interrupted data packets get 
void Metrics::update_on_interruption(const Event& event) {
    assert(event.type == DISPATCH && event.packet_type == DATA);

    if (!should_collect(event)) return;

    // Update the partial service time for the interrupted data packet
    tmpX1 += event.packet_size / SPEED - (event.t - sim_t);

    // Also update the product Nq1 * time for the interrupted data packet,
    // since its return to the queue will increase the queue size
    update_Nq_sum(event);
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
