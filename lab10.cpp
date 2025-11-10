#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include <cmath>

enum ParticleType
{
    EMPTY = 0,
    CO_AD = 1,
    O_AD = 2,
    O_V = 3,
};

const int L = 20;
const int N = L * L;
const double TOTAL_TIME = 300.0;

const double k01 = 1e4;
const double P_CO = 0.3;
const double k1 = P_CO * k01;
const double k_min1 = 3e2;
const double k2 = 2.5e3;
const double k3 = 2.5e4;
const double k4 = 0.11;
const double k5 = 6.5e-3;
const double k6_7_8 = 1e5;

struct SystemState
{
    std::vector<std::vector<int>> lattice;
    double time;
    int step;
    int empty_count;
    int co_ad_count;
    int o_ad_count;
    int o_v_count;
};

enum EventType
{
    CREATE_CO_AD = 1,
    CREATE_CO = -1,
    CREATE_O_AD = 2,
    CREATE_CO_2 = 3,
    CREATE_O_V = 4,
    RECOVERY_CO_2 = 5,
    MOVE_CO_AD = 6,
    MOVE_O_AD = 7,
    MOVE_O_V = 8
};

struct Event
{
    int eventType;
    double rate;
    std::vector<std::pair<int, int>> pos;
};

std::mt19937 rng;
std::uniform_real_distribution<double> uniform(0.0, 1.0);

void initializeLattice(SystemState& system, bool empty_lattice)
{
    system.lattice.resize(L, std::vector<int>(L, EMPTY));
    system.time = 0.0;
    system.step = 0;
    system.empty_count = 0;
    system.co_ad_count = 0;
    system.o_ad_count = 0;
    system.o_v_count = 0;
    
    if (!empty_lattice)
    {
        int target_o = static_cast<int>(0.4 * N);
        int target_ov = static_cast<int>(0.4 * N);
        
        for (int i = 0; i < target_o; i++)
        {
            int pos_x;
            int pos_y;
            do
            {
                pos_x = static_cast<int>(uniform(rng) * L);
                pos_y = static_cast<int>(uniform(rng) * L);
            }
            while (system.lattice[pos_x][pos_y] != EMPTY);
            system.lattice[pos_x][pos_y] = O_AD;
            system.o_ad_count++;
        }
        
        for (int i = 0; i < target_ov; i++) {
            int pos_x;
            int pos_y;
            do
            {
                pos_x = static_cast<int>(uniform(rng) * L);
                pos_y = static_cast<int>(uniform(rng) * L);
            } 
            while (system.lattice[pos_x][pos_y] != EMPTY);
            system.lattice[pos_x][pos_y] = O_V;
            system.o_v_count++;
        }
    }
}

std::vector<std::pair<int, int>> getNeighbors(int x, int y)
{
    std::vector<std::pair<int, int>> neighbors;
    neighbors.push_back(std::pair((x + 1) % L, y));
    neighbors.push_back(std::pair(x - 1 >= 0 ? x - 1 : L - 1, y));
    neighbors.push_back(std::pair(x, (y + 1) % L));
    neighbors.push_back(std::pair(x, y - 1 >= 0 ? y - 1 : L - 1));
    return neighbors;
}

std::vector<Event> calculateEvents(const SystemState& system)
{
    std::vector<Event> events;
    
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            int curr_cell = system.lattice[x][y];
            // 4. [O] -> [O]v
            if (curr_cell == O_AD)
            {
                events.push_back({CREATE_O_V, k4, {{x, y}}});
            }
            // 1+. CO + * -> [CO]
            else if (curr_cell == EMPTY) events.push_back({CREATE_CO_AD, k1, {{x, y}}});
            // 1-. [CO] -> CO + *
            else if (curr_cell == CO_AD) events.push_back({CREATE_CO, k_min1, {{x, y}}});

            auto neighbours = getNeighbors(x, y);
            for (int i = 0; i < neighbours.size(); i++)
            {
                int n_x = neighbours[i].first;
                int n_y = neighbours[i].second;
                int n_cell = system.lattice[n_x][n_y];

                if (n_cell == EMPTY)
                {
                    // 2. O2 + 2* -> 2[O]
                    if (curr_cell == EMPTY) events.push_back({CREATE_O_AD, k2, {{x, y}, {n_x, n_y}}});
                    // 6. Миграция [CO]
                    else if (curr_cell == CO_AD) events.push_back({MOVE_CO_AD, k6_7_8, {{x, y}, {n_x, n_y}}});
                    // 7. Миграция [O]
                    else if (curr_cell == O_AD) events.push_back({MOVE_O_AD, k6_7_8, {{x, y}, {n_x, n_y}}});
                    // 8. Миграция [O]v
                    else if (curr_cell == O_V) events.push_back({MOVE_O_V, k6_7_8, {{x, y}, {n_x, n_y}}});
                }

                if (curr_cell == CO_AD)
                {
                    // 3. [CO] + [O] -> 2* + CO2
                    if (n_cell == O_AD) events.push_back({CREATE_CO_2, k3, {{x, y}, {n_x, n_y}}});
                    // 5. [CO] + [O]v -> 2* + CO2
                    else if (n_cell == O_V) events.push_back({RECOVERY_CO_2, k5, {{x, y}, {n_x, n_y}}});
                }
            }
        }
    }
    
    return events;
}

void executeEvent(SystemState& system, const Event& event) {
    auto eventType = event.eventType;
    double rate = event.rate;
    auto pos = event.pos;
    int pos_x = pos[0].first;
    int pos_y = pos[0].second;
    int pos2_x, pos2_y;
    if (pos.size() == 2)
    {
        pos2_x = pos[1].first;
        pos2_y = pos[1].second;
    }
    
    switch (eventType) {
        case CREATE_CO_AD:
            system.lattice[pos_x][pos_y] = CO_AD;
            system.co_ad_count++;
            break;
        case CREATE_CO:
            system.lattice[pos_x][pos_y] = EMPTY;
            system.co_ad_count--;
            break;
        case CREATE_O_AD:
            system.lattice[pos_x][pos_y] = O_AD;
            system.lattice[pos2_x][pos2_y] = O_AD;
            system.o_ad_count += 2;
            break;
        case CREATE_CO_2:
            system.lattice[pos_x][pos_y] = EMPTY;
            system.lattice[pos2_x][pos2_y] = EMPTY;
            system.o_ad_count--;
            system.co_ad_count--;
            break;
        case CREATE_O_V:
            system.lattice[pos_x][pos_y] = O_V;
            system.o_ad_count--;
            system.o_v_count++;
            break;
        case RECOVERY_CO_2:
            system.lattice[pos_x][pos_y] = EMPTY;
            system.lattice[pos2_x][pos2_y] = EMPTY;
            system.o_v_count--;
            system.co_ad_count--;
            break;
        case MOVE_CO_AD:
            system.lattice[pos_x][pos_y] = EMPTY;
            system.lattice[pos2_x][pos2_y] = CO_AD;
            break;
        case MOVE_O_AD:
            system.lattice[pos_x][pos_y] = EMPTY;
            system.lattice[pos2_x][pos2_y] = O_AD;
            break;
        case MOVE_O_V:
            system.lattice[pos_x][pos_y] = EMPTY;
            system.lattice[pos2_x][pos2_y] = O_V;
            break;
    }
}

std::string formatDuration(std::chrono::steady_clock::duration duration) {
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
    
    std::ostringstream oss;
    oss << std::setw(2) << std::setfill('0') << hours.count() << ":"
        << std::setw(2) << std::setfill('0') << minutes.count() << ":"
        << std::setw(2) << std::setfill('0') << seconds.count();
    return oss.str();
}

int main() {
    rng.seed(std::chrono::steady_clock::now().time_since_epoch().count());
    SystemState system;
    bool empty_lattice = false;
    initializeLattice(system, empty_lattice);
    std::ofstream datafile("simulation_data_" + std::to_string(L) + "_" + std::to_string(TOTAL_TIME) + "_" + std::to_string(empty_lattice) + "_" + std::to_string(P_CO) + ".csv");
    datafile << "step,time,empty_count,co_ad_count,o_ad_count,o_v_count\n";  
    auto start_time = std::chrono::steady_clock::now();
    std::cout << "Start" << std::endl;
    std::cout << "Size: " << L << "x" << L << std::endl;
    std::cout << "Total time: " << TOTAL_TIME << " s" << std::endl;
    
    while (system.time < TOTAL_TIME) {
        auto events = calculateEvents(system);

        if (events.size() == 0)
        {
            std::cout << "No events?\n";
            break;
        }
        
        double total_rate = 0.0;
        for (const auto& event : events) {
            total_rate += event.rate;
        }
        
        double random_value = uniform(rng) * total_rate;
        double sum_rate = 0.0;
        for (const auto& event : events)
        {
            sum_rate += event.rate;
            if (random_value <= sum_rate)
            {
                int empty_count = N - system.co_ad_count - system.o_ad_count - system.o_v_count;
                datafile << system.step << "," << system.time << ","
                        << empty_count << "," << system.co_ad_count << "," 
                        << system.o_ad_count << "," << system.o_v_count << "\n";
                datafile.flush();

                if (system.step % 10000 == 0)
                {
                    auto current_time = std::chrono::steady_clock::now();
                    auto elapsed = current_time - start_time;
                    
                    std::cout << "[" << formatDuration(elapsed) << "] "
                        << "Step: " << std::setw(10) << system.step 
                        << ", Time: " << std::setw(10) << system.time << " s"
                        << ", Empty: " << std::setw(10) << empty_count
                        << ", CO: " << std::setw(10) << system.co_ad_count 
                        << ", O: " << std::setw(10) << system.o_ad_count
                        << ", Ov: " << std::setw(10) << system.o_v_count
                        << std::endl;
                }

                executeEvent(system, event);
                system.time += log(event.rate) / total_rate;
                system.step++;
                break;
            }
        }
    }
    
    int empty_count = N - system.co_ad_count - system.o_ad_count - system.o_v_count;
    datafile << system.step << "," << system.time << ","
                << empty_count << "," << system.co_ad_count << "," 
                << system.o_ad_count << "," << system.o_v_count << "\n";
    
    datafile.close();
    
    auto end_time = std::chrono::steady_clock::now();
    auto total_duration = end_time - start_time;
    
    std::cout << "[" << formatDuration(total_duration) << "] "
                     << "Step: " << system.step 
                     << ", Time: " << system.time << " s"
                     << ", Empty: " << empty_count
                     << ", CO: " << system.co_ad_count 
                     << ", O: " << system.o_ad_count
                     << ", Ov: " << system.o_v_count << std::endl;
    
    return 0;
}