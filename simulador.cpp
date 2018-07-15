#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "boost/math/distributions/students_t.hpp"

using namespace std;
using namespace boost;

/*========================================
=            TYPE DEFINITIONS            =
========================================*/

enum Packet_type {DATA, VOICE};

enum Event_type {ARRIVAL, DEPARTURE};

typedef unordered_map<string, double> Metrics;

struct Event {
    // Every event should have a <type> and the point in time <t> at which it happens
    Event_type type;
    double t;

    // A unique id for each arrival. Used to guarantee std::set won't ignore events
    int id;

    // Identifies the round the event belongs to
    int round_id;

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
    if (type != rhs.type) return type == DEPARTURE;
    if (packet_type != rhs.packet_type) return packet_type == VOICE;
    return id < rhs.id;
}

struct RoundMetrics {
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

    Metrics sum;
    Metrics sum_of_squares;
    Metrics param;

    double tmpX1;

    int packets_processed;
    int data_packets_processed;
    int voice_packets_processed;
    int delta_intervals;
    double last_t1, last_t2;

    RoundMetrics() {}

    RoundMetrics(int rnd, int tgt);
    void init();
    void update_Nq_sum(const Event &event);
    void update_on_departure(const Event &cur_event);
    void update_on_interruption(const Event &event);

private:
    bool should_collect(const Event &event);
    void compute_estimators();
};

RoundMetrics::RoundMetrics(int rnd, int tgt):
    round_id(rnd), target(tgt), start_t(-1), end_t(-1),
    tmpX1(0), packets_processed(0), data_packets_processed(0),
    voice_packets_processed(0), delta_intervals(0), last_t1(0), last_t2(0) {}

struct FinalMetrics {
    int num_samples;

    Metrics sum;
    Metrics sum_of_squares;
    Metrics mean;
    Metrics ci_halfwidth;
    Metrics precision;

    FinalMetrics();
    void add_metrics(const RoundMetrics &round_metrics);
    bool has_enough_precision(bool interrupt);
    void show();
};

FinalMetrics::FinalMetrics(): num_samples(0) {}

/*=========================================
=            PROBLEM CONSTANTS            =
=========================================*/

constexpr double p1 = 0.3;
constexpr double p2 = 0.1;
constexpr double p3 = 0.3;
constexpr double p = 0.3;
constexpr double SPEED = 2e3/8; // bytes per ms
constexpr double MEAN_DATA_PACKET_SIZE = 755; // bytes
constexpr int NUM_CHANNELS = 30;
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
double last_voice_departure_t[NUM_CHANNELS];
vector<RoundMetrics> round_metrics;
FinalMetrics results;

/*===========================================
=            FUNCTION PROTOTYPES            =
===========================================*/

/*----------  RANDOM NUMBER GENERATION SETUP  ----------*/
int get_packet_size();
/*----------  INITIALIZATION  ----------*/
void simulation_state_init(int warmup_period, int round_size, double rho);
void event_queue_init();
/*----------  MAIN  ----------*/
/*----------  SIMULATION LOGIC  ----------*/
void run_simulation(int warmup_period, int round_size, double rho, bool interrupt);
bool simulation_should_stop(bool interrupt);
void handle_arrival(Event &cur_event, bool interrupt);
void handle_data_arrival(Event &cur_event);
void handle_voice_arrival(Event &cur_event, bool interrupt);
void serve_next_packet();
void enter_the_server(const Event &event, bool force=false);
void unschedule_data_departure();
/*----------  EVENT CREATION  ----------*/
Event make_arrival_from(const Event &event);
Event make_departure_from(const Event &event);
/*----------  STATISTICS  ----------*/
pair<double, double> get_ci(double mu, double var, int n);
double variance(double sumX, double sumXsq, int n);
// See class RoundMetrics and FinalMetrics declarations

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
    ios_base::sync_with_stdio(false);

    const string mt_state_filename = "mt_state.sav";
    ifstream mt_state_r(mt_state_filename, ios::in | ios::binary);

    if (mt_state_r.good()) {
        mt_state_r >> mt;
        mt_state_r.close();
    }

    int warmup_period, round_size, interrupt;
    double rho;

    if (argc == 5) {
        sscanf(argv[1], "%d", &warmup_period);
        sscanf(argv[2], "%d", &round_size);
        sscanf(argv[3], "%lf", &rho);
        sscanf(argv[4], "%d", &interrupt);

        run_simulation(warmup_period, round_size, rho, interrupt);
    } else {
        warmup_period = 10000;
        round_size = 15000;
        
        for (int i = 1; i < 8; ++i) {
            run_simulation(warmup_period, round_size, i/10.0, false);
        }
        cout << endl;

        for (int i = 1; i < 8; ++i) {
            run_simulation(warmup_period, round_size, i/10.0, true);
        }
    }

    ofstream mt_state_w(mt_state_filename, ios::out | ios::binary);
    mt_state_w << mt;
    mt_state_w.close();

    return 0;
}

/*======================================
=            INITIALIZATION            =
======================================*/

void simulation_state_init(int warmup_period, int round_size, double rho) {
    sim_t = 0;
    event_queue.clear();
    data_queue.clear();
    voice_queue.clear();
    idle = true;
    time_between_data_packets = exponential_distribution<double>(rho/MEAN_DATA_SERVICE_TIME);
    id_counter = 0;
    cur_round = 0;

    results = FinalMetrics();

    round_metrics.clear();
    if (warmup_period)
        round_metrics.emplace_back(-1, warmup_period);
    else
        round_metrics.emplace_back(0, round_size);

    round_metrics[0].init();
}

void event_queue_init() {
    assert(event_queue.empty());
    event_queue.emplace();

    for (int i = 0; i < NUM_CHANNELS; ++i) {
        event_queue.emplace(i);
    }
}

/*========================================
=            SIMULATION LOGIC            =
========================================*/

void run_simulation(int warmup_period, int round_size, double rho, bool interrupt) {
    // Simulation state initialization
    simulation_state_init(warmup_period, round_size, rho);

    // Initialize event queue
    event_queue_init();

    while (true) {
        // Get next event and then remove it from event queue
        Event cur_event = *event_queue.begin();
        event_queue.erase(event_queue.begin());

        // Update simulation time
        sim_t = cur_event.t;

        // Deal with the event
        if (cur_event.type == ARRIVAL) {
            round_metrics[cur_round].update_Nq_sum(cur_event);
            handle_arrival(cur_event, interrupt);
        } else {
            round_metrics[cur_round].update_on_departure(cur_event);
            
            if (cur_event.packet_type == VOICE) {
                last_voice_departure_t[cur_event.channel] = sim_t;
            }
            
            serve_next_packet();
        }

        // Go to the next round, if the current one is finished
        if (round_metrics[cur_round].end_t >= 0) {
            // If this is the warm-up period ending,
            // start the actual first run
            if (round_metrics[cur_round].round_id == -1) {
                round_metrics[0] = RoundMetrics(0, round_size);
            } else {
                // Add this run's measurements to final results
                results.add_metrics(round_metrics[cur_round]);

                // Stop simulation if the required precision has been achieved
                if (results.has_enough_precision(interrupt))
                    break;

                // Create a new round if simulation didn't stop
                round_metrics.emplace_back(cur_round++, round_size);
            }

            // Initialize next round
            round_metrics[cur_round].init();
        }
    }
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
}

void enter_the_server(const Event &event, bool force) {
    assert(event.type == ARRIVAL);

    // Voice packets may enter the server by force when interruption is on
    // but the data packet on the server must be returned to the queue
    if (force) {
        assert(!idle && event.packet_type == VOICE && event.type == ARRIVAL);
        unschedule_data_departure();
    }

    // Schedule the departure of the entering event
    event_queue.insert(make_departure_from(event));

    // Update server state after a packet enters
    event_on_server = event;
    idle = false;
}

void unschedule_data_departure() {
    // Find the DEPARTURE event in the event queue
    auto it = find_if(event_queue.begin(), event_queue.end(),
        [] (const Event &e) {
            return e.type == DEPARTURE;
        });

    assert(it != event_queue.end() && it->type == DEPARTURE && it->packet_type == DATA);

    round_metrics[cur_round].update_on_interruption(*it);

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
    t = time_between_data_packets(mt);
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
    }

    // Every arrival event is assigned a unique id
    new_event.id = id_counter++;

    return new_event;
}

Event make_departure_from(const Event &event) {
    assert(event.type == ARRIVAL);
    Event new_event(event);

    new_event.type = DEPARTURE;

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

double variance(double sum, double sum_of_squares, int num_samples) {
    return sum_of_squares/(num_samples-1) - (sum*sum)/((double) num_samples * (num_samples-1));
}

void RoundMetrics::init() {
    if (start_t >= 0) return;
    start_t = sim_t;
    last_t1 = last_t2 = sim_t;
}

bool RoundMetrics::should_collect(const Event &event) {
    // Returns true if:
    // 1) Metrics for this round are not finished
    // 2) These are not metrics for the warm-up period
    // 3) The event belongs to this round
    return end_t < 0 && round_id != -1; // && round_id == event.round_id;
}

// Once the round is finished, the parameters estimators
// must be evaluated
void RoundMetrics::compute_estimators() {
    double total_time = end_t - start_t;

    if (data_packets_processed) {
        param["ET1"] = sum["T1"] / data_packets_processed;
        param["EW1"] = sum["W1"] / data_packets_processed;
        param["EX1"] = sum["X1"] / data_packets_processed;
    } else {
        param["ET1"] = param["EW1"] = param["EX1"] = 0;
    }
    param["ENq1"] = sum["Nq1"] / total_time;

    if (voice_packets_processed) {
        param["ET2"] = sum["T2"] / voice_packets_processed;
        param["EW2"] = sum["W2"] / voice_packets_processed;
    } else {
        param["ET2"] = param["EW2"] = 0;
    }
    param["ENq2"] = sum["Nq2"] / total_time;

    if (delta_intervals > 0) {
        param["EDelta"] = sum["Delta"] / delta_intervals;
    } else {
        param["EDelta"] = 0;
    }

    if (delta_intervals > 1) {
        param["VDelta"] = variance(sum["Delta"], sum_of_squares["Delta"], delta_intervals);
    } else {
        param["VDelta"] = 0;
    }

}

// Updates the cumulative product Nq x Δt
// Should be called on every arrival, departure and interruption
void RoundMetrics::update_Nq_sum(const Event& event) {
    if (event.packet_type == DATA) {
        sum["Nq1"] += (sim_t - last_t1) * data_queue.size();
        last_t1 = sim_t;
    } else {
        sum["Nq2"] += (sim_t - last_t2) * voice_queue.size();
        last_t2 = sim_t;
    }
}

// Updates sums on departures.
void RoundMetrics::update_on_departure(const Event& cur_event) {
    assert(cur_event.type == DEPARTURE && cur_event.t == sim_t);
    
    update_Nq_sum(cur_event);

    if (!should_collect(cur_event)) {
        // In case this is the warm-up period, count every 
        // packet so that the warm-up period can finish
        if (round_id < 0) {
            packets_processed++;
            if (packets_processed == target)
                end_t = sim_t;
        }
        return;
    }

    if (cur_event.packet_type == DATA) {
        data_packets_processed++;

        sum["X1"] += cur_event.packet_size / SPEED + tmpX1;
        tmpX1 = 0;

        sum["T1"] += sim_t - event_on_server.t;

        sum["W1"] += sim_t - event_on_server.t - (cur_event.packet_size / SPEED + tmpX1);

    } else {
        voice_packets_processed++;

        sum["T2"] += sim_t - event_on_server.t;

        sum["W2"] += sim_t - event_on_server.t - VOICE_SERVICE_TIME;

        // If the voice packet is not the first in its group,
        // compute the interval Δj between this and the last departure.
        if (cur_event.vgroup_idx > 0) {
            double Dj = sim_t - last_voice_departure_t[cur_event.channel];
            sum["Delta"] += Dj;
            sum_of_squares["Delta"] += Dj * Dj;
            delta_intervals++;
        }
    }

    packets_processed++;

    if (packets_processed == target) {
        end_t = sim_t;
        compute_estimators();
    }

}

// Interrupted data packets get partial service times and change
// data queue size
void RoundMetrics::update_on_interruption(const Event& event) {
    assert(event.type == DEPARTURE && event.packet_type == DATA);

    // Update Nq1 * time for the interrupted data packet,
    // since its return to the queue will increase the queue size
    update_Nq_sum(event);

    if (!should_collect(event)) return;

    // Update the partial service time for the interrupted data packet
    tmpX1 += event.packet_size / SPEED - (event.t - sim_t);

}

void FinalMetrics::add_metrics(const RoundMetrics &round_metrics) {
    assert(round_metrics.end_t >= 0);
    num_samples++;

    for (auto const &it: round_metrics.param) {
        sum[it.first] += it.second;
        sum_of_squares[it.first] += it.second * it.second;

        #ifdef PYTHON_SCRIPT
        cerr << it.first << " " << it.second << " ";
        #endif

        // if (it.first[it.first.size()-1] != '1') {
        //     cout << fixed << setprecision(5);
        //     cout << it.first << " " << sum[it.first]/num_samples << " ";
        // }
    }

    // cout << endl;

    #ifdef PYTHON_SCRIPT
    cerr << endl;
    #endif

}

bool FinalMetrics::has_enough_precision(bool interrupt) {
    if (num_samples < 8) return false;
    
    // Compute confidence intervals for each parameter
    math::students_t t_dist(num_samples-1);
    double t_ppf = quantile(complement(t_dist, 0.1/2));
    for (auto const &it: sum) {
        auto const &param = it.first;
        double var = variance(sum[param], sum_of_squares[param], num_samples);
        
        mean[param] = sum[param] / num_samples;
        ci_halfwidth[param] = t_ppf * sqrt(var/num_samples);
        precision[param] = ci_halfwidth[param] / mean[param];
    }

    // If there are no interruptions, check if precisions
    // for the data channel parameters are okay (except E[X1])
    if (!interrupt) {
        if (precision["EW1"] > 0.05 || precision["ET1"] > 0.05 || precision["ENq1"] > 0.05) {
            return false;
        }
    }

    // Check if precisions for the voice channels (and E[X1]) parameters are okay
    if (precision["EW2"] <= 0.05 && precision["ET2"] <= 0.05 && precision["ENq2"] <= 0.05 && precision["EX1"] <= 0.05) {
        show();
        return true;
    }

    return false;
}

void FinalMetrics::show() {
    #define show_param(x) cout << mean[x] - ci_halfwidth[x] << "," << mean[x] << "," << mean[x] + ci_halfwidth[x] << ",";

    cout << fixed << setprecision(4);
    show_param("ET1");
    show_param("EW1");
    show_param("EX1");
    show_param("ENq1");
    show_param("ET2");
    show_param("EW2");
    show_param("ENq2");
    show_param("EDelta");
    show_param("VDelta");
    cout << defaultfloat << endl;
}
