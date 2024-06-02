#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <chrono>

#include "callback.h"
#include "col_gen.h"
#include "comb.h"
#include "cvrpsep/strcomb.h"
#include "glm.h"
#include "gurobi_c++.h"
#include "load.h"
#include "mstar.h"
#include "rcc_sep.h"
#include "utils.h"
#include <unordered_set>

struct DistanceLocationPair {
    int distance;
    int location;

    DistanceLocationPair(int d, int i) : distance(d), location(i) {}
};

class CVRP {
private:
    

public:
    // Instance Properties
    std::string name;
    std::vector<std::vector<int>> cost;
    std::vector<int> demand;
    int capacity;
    std::vector<std::pair<double, double>> coords;
    int dimension;
    int upper_bound;
    std::vector<std::vector<int>> ub_routes;
    std::set<int> vertices;
    std::set<int> customers;
    std::set<std::tuple<int, int>> edges;
    std::set<std::tuple<int, int>> asym_edges;
    std::set<std::tuple<int, int>> edges_subset;
    std::set<std::tuple<int, int>> edges_not_subset;
    std::set<std::tuple<int, int>> starting_edges;
    std::set<std::tuple<int, int>> ub_edges;
    
    // Model Properties
    GRBModel* model;
    std::vector<std::vector<GRBVar>> x;
    std::vector<std::vector<std::vector<GRBVar>>> f;

    // Algorithmic Properties
    int inner_iterations = 0;
    int col_gen_iters = 0;
    int variables_eliminated = 0;

    int total_rcc = 0;
    int current_rcc = 0;

    int total_combs = 0;
    int current_combs = 0;

    int total_mstars = 0;
    int current_mstars = 0;

    int total_glm = 0;
    int current_glm = 0;

    std::vector<std::vector<GRBColumn>> var_columns;

    bool column_generation = false;
    bool flow_variables = false;

    // Time Properties
    double total_time = 0.0;
    double setup_time = 0.0;
    double starting_subset_time = 0.0;
    double LP_solve_time = 0.0;
    double separation_time = 0.0;
    double add_vars_time = 0.0;
    double remove_constrs_time = 0.0;
    double eliminate_vars_time = 0.0;

    double RCC_separation_time = 0.0;
    double RCC_heuristic_time = 0.0;
    double RCC_exact_time = 0.0;
    double RCC_shrink_time = 0.0;
    double RCC_exact_solve_time = 0.0;
    double RCC_convert_shrunk_cuts_time = 0.0;
    double add_RCC_time = 0.0;
    double determine_RCC_time = 0.0;

    double SC_separation_time = 0.0;
    double SC_heuristic_time = 0.0;
    double SC_exact_time = 0.0;
    double SC_exact_solve_time = 0.0;
    double add_SC_time = 0.0;
    double determine_SC_time = 0.0;

    double MSTAR_separation_time = 0.0;
    double MSTAR_heuristic_time = 0.0;
    double MSTAR_exact_time = 0.0;
    double MSTAR_exact_solve_time = 0.0;
    double add_MSTAR_time = 0.0;
    double determine_MSTAR_time = 0.0;

    double GLM_total_separation_time = 0.0;
    double GLM_separation_time = 0.0;
    double add_GLM_time = 0.0;
    double determine_GLM_time = 0.0;

    // Callbacks
    my_rc_callback rc_callback;
    mycallback callback;


    // Initializes instance
    CVRP(const std::string& filename) {
        // Load Data
        InstanceData inst_data = load_instance(filename);

        // Instance Properties
        name = filename;
        cost = inst_data.distances;
        demand = inst_data.demand;
        capacity = inst_data.capacity;
        coords = inst_data.coords;
        dimension = inst_data.dimension;
        upper_bound = inst_data.upper_bound;
        ub_routes = inst_data.routes;

        // Populate vertices, customers, and edges
        for (int i = 0; i < dimension; ++i) {
            vertices.insert(i);
            if (i > 0) customers.insert(i);
        }
        for (int i : vertices) {
            for (int j : customers) {
                if (i < j) edges.insert({ i, j });
            }
        }

        edges_subset = edges;

        // Model properties
        model = nullptr;
        x.resize(dimension);
        var_columns.resize(dimension);
        for (int i = 0; i < dimension; ++i) {
            x[i].resize(dimension); // Reserve space for the inner vector
            var_columns[i].resize(dimension);
            for (int j = i + 1; j < dimension; j++) {
                var_columns[i][j] = GRBColumn();
            }
        }
        //f.clear();
    }

    
    // Computes the intial sets of edges
    void set_starting_edges(bool include_bk_sol = true, int no_of_shortest_edges = 15) {
        auto start = std::chrono::high_resolution_clock::now();

        std::vector<std::set<int>> route_locations(ub_routes.size() + 1);
        column_generation = true;
        edges_subset.clear();
        for (int i = 0; i < ub_routes.size(); i++) {
            int start_edge = 0;
            for (int j = 0; j < ub_routes[i].size(); j++) {
                route_locations[i].insert(ub_routes[i][j]);
                int end_edge = ub_routes[i][j];
                if (start_edge < end_edge) edges_subset.insert({ start_edge, end_edge });
                else edges_subset.insert({ end_edge, start_edge });
                start_edge = end_edge;
            }
            edges_subset.insert({ 0, start_edge });
        }

        ub_edges = edges_subset;

        for (int i = 0; i < dimension; i++) {
            std::vector<DistanceLocationPair> distances;
            for (int j = 0; j < dimension; j++) {
                if (j < i) distances.push_back(DistanceLocationPair(cost[j][i], j));
                else if (i < j) distances.push_back(DistanceLocationPair(cost[i][j], j));
            }

            // Sort distances in ascending order
            std::sort(distances.begin(), distances.end(), [](const DistanceLocationPair& a, const DistanceLocationPair& b) {
                return a.distance < b.distance;
                });

            for (int k = 0; k < no_of_shortest_edges ; k++) {
                int j = distances[k].location;
                if (j < i) edges_subset.insert({ j, i });
                else if (i < j) edges_subset.insert({ i, j });
            }
        }

        std::set_difference(edges.begin(), edges.end(), edges_subset.begin(), edges_subset.end(),
            std::inserter(edges_not_subset, edges_not_subset.begin()));

        starting_edges = edges_subset;
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
        starting_subset_time += duration.count();
    }


    // Sets up the base CVRP model (TI)
    void setup_base_model(bool relax = true, bool silent = true, GRBEnv env = GRBEnv()) {
        auto start = std::chrono::high_resolution_clock::now();
        // Create Model
        model = new GRBModel(env);

        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);
        // Add Variables
        double ub;
        double obj;
        char vtype;
        std::string name;
        for (auto& edge : edges_subset) {
            int i = std::get<0>(edge);
            int j = std::get<1>(edge);
            name = "x[" + itos(i) + "][" + itos(j) + "]";
            ub = (i == 0) ? 2.0 : 1.0; // UB = 2 if i = 0 else 1
            obj = (i < j) ? cost[i][j] : cost[j][i];
            if (relax) {
                vtype = GRB_CONTINUOUS;
            }
            else {
                vtype = (i == 0) ? GRB_INTEGER : GRB_BINARY;
            }
            
            x[i][j] = model->addVar(0.0, ub, obj, vtype, name);
        }

        // Add Constraints
        for (int k : customers) {
            GRBLinExpr expr;
            std::set<int> k_set;
            k_set.insert(k);
            std::set<std::tuple<int, int, double>> edges_not_constraint;
            for (const auto& edge : delta_edges(edges, vertices, k_set)) {
                int i = std::get<0>(edge);
                int j = std::get<1>(edge);
                if (edges_subset.count({ i, j }) == 1) expr += x[i][j];
                else edges_not_constraint.insert({ i, j, 1.0 });
            }
            GRBConstr constr = model->addConstr(expr, GRB_EQUAL, 2, "y_" + itos(k));
            for (const auto& var : edges_not_constraint) {
                int i = std::get<0>(var);
                int j = std::get<1>(var);
                double coeff = std::get<2>(var);
                var_columns[i][j].addTerm(coeff, constr);
            }
            expr.clear(); // Clear expr for next iteration
            k_set.clear();
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
        setup_time += duration.count();
    }


    // Adds all x varaibles to the model
    void add_all_x_variables() {
        // Adds all x variables to the model
        add_all_variables(edges_not_subset, var_columns, model, x, edges_subset, cost, add_vars_time);
    }


    // Adds flow variables and constraints to the model
    void add_flow_variables() {
        // Adds all x and flow variables to the model - Upper bounds are set to for 
        for (auto& edge : edges) {
            asym_edges.insert(edge);
            asym_edges.insert({ std::get<1>(edge), std::get<0>(edge) });
        }
        f.resize(dimension);
        for (int i = 0; i < dimension; ++i) {
            f[i].resize(dimension); // Reserve space for the inner vector
            for (int j = 0; j < dimension; j++) {
                f[i][j].resize(dimension);
            }
        }


        flow_variables = true;
        add_all_variables(edges_not_subset, var_columns, model, x, edges_subset, cost, add_vars_time);

        for (const auto& edge : asym_edges) {
            int i = std::get<0>(edge);
            int j = std::get<1>(edge);
            double ub = (i == 0) ? 2.0 : 1.0;
            ub = (column_generation) ? 0.0 : ub;
            
            for (int k : customers) {
                std::string name = "f[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
                f[i][j][k] = model->addVar(0.0, ub, NULL, GRB_CONTINUOUS, name);
            }
        }

        if (column_generation) {
            for (auto& edge : edges_not_subset) {
                int i = std::get<0>(edge);
                int j = std::get<1>(edge);
                x[i][j].set(GRB_DoubleAttr_UB, 0.0);
            }

            for (auto& edge : edges_subset) {
                int i = std::get<0>(edge);
                int j = std::get<1>(edge);
                for (int k : customers) {
                    double ub = (i == 0) ? 2.0 : 1.0;
                    f[i][j][k].set(GRB_DoubleAttr_UB, ub);
                    f[j][i][k].set(GRB_DoubleAttr_UB, ub);
                }
            }
        }

        model->update();

        std::vector<std::vector<GRBLinExpr>> flow_minus(dimension);
        std::vector<std::vector<GRBLinExpr>> flow_plus(dimension);
        
        for (int l : vertices) {
            flow_minus[l].resize(dimension);
            flow_plus[l].resize(dimension);
            for (int k : customers) {
                GRBLinExpr in;
                GRBLinExpr out;
                for (const auto& edge : asym_edges) {
                    int i = std::get<0>(edge);
                    int j = std::get<1>(edge);
                    if (i == l) out += f[i][j][k];
                    if (j == l) in += f[i][j][k];
                }
                flow_minus[l][k] = in;
                flow_plus[l][k] = out;
            }
        }

        for (int k : customers) {
            model->addConstr(flow_plus[0][k] == 2);
            model->addConstr(flow_minus[k][k] == 2);
            model->addConstr(flow_minus[0][k] == 0);
            model->addConstr(flow_plus[k][k] == 0);
            for (int l : customers) {
                if (l != k) {
                    model->addConstr(flow_plus[l][k] == flow_minus[l][k]);
                    model->addConstr(flow_plus[l][k] == flow_plus[k][l]);
                }
            }

            for (auto& edge : edges) {
                int i = std::get<0>(edge);
                int j = std::get<1>(edge);
                model->addConstr(f[i][j][k] + f[j][i][k] <= x[i][j]);
            }
        }

        for (auto& edge : edges) {
            GRBLinExpr Lexpr;
            int i = std::get<0>(edge);
            int j = std::get<1>(edge);
            for (int k : customers) {
                if (k != i && k != j) Lexpr += demand[k] * (f[i][j][k] + f[j][i][k]);
            }
            model->addConstr(Lexpr <= (capacity - demand[i] - demand[j]) * x[i][j]);
        }
    }


    // Makes the model integral
    void make_integral() {
        // Makes the model integral
        solve_model();

        for (const auto& edge : edges) {
            int i = std::get<0>(edge);
            int j = std::get<1>(edge);
            if (i > 0) x[i][j].set(GRB_CharAttr_VType, GRB_BINARY);
        }
        model->getEnv().set(GRB_IntParam_OutputFlag, 1);
    }


    // Adds callback only separating RC inequalities
    void add_rc_callback(float eps_for_int, float eps_for_rcc_viol, CnstrMgrPointer MyOldCutsCMP, CnstrMgrPointer MyCutsCMP, int rcc_implementation, bool callback_in_nodes) {
        model->set(GRB_IntParam_LazyConstraints, 1);
        model->update();
        rc_callback = my_rc_callback(model, vertices, customers, edges, edges_subset, x, eps_for_int, eps_for_rcc_viol, dimension, demand, capacity,
             MyOldCutsCMP, MyCutsCMP, total_rcc, current_rcc, rcc_implementation, var_columns, determine_RCC_time, callback_in_nodes);
        model->setCallback(&rc_callback);
    }


    // Adds callback separating all inequalities
    void add_separation_callback(float eps_for_int, float eps_for_rcc_viol, CnstrMgrPointer MyOldCutsCMP, CnstrMgrPointer MyCutsCMP, int rcc_implementation) {
        model->set(GRB_IntParam_LazyConstraints, 1);
        model->update();
        callback = mycallback(model, vertices, customers, edges, edges_subset, x, eps_for_int, eps_for_rcc_viol, dimension, demand, capacity,
            MyOldCutsCMP, MyCutsCMP, total_rcc, current_rcc, rcc_implementation, var_columns, determine_RCC_time);
        model->setCallback(&callback);
    }


    // Optimizes the current model
    void solve_model() {
        auto start = std::chrono::high_resolution_clock::now();
        model->optimize();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
        LP_solve_time += duration.count();
    }


    // Retrieves all variables with positive solution value
    std::vector<std::tuple<int, int, double>> get_solution(double eps = 0.001) {
        std::vector<std::tuple<int, int, double>> solution;
        for (auto& edge : edges_subset) {
            int i = std::get<0>(edge);
            int j = std::get<1>(edge);
            double sol_val = x[i][j].get(GRB_DoubleAttr_X);
            if (sol_val >= eps) {
                solution.push_back(std::make_tuple(i, j, sol_val));
            }
        }
        return solution;
    }


    // Retrieves the objective value
    double get_objective_value() {
        return model->get(GRB_DoubleAttr_ObjVal);
    }

    
    // Separates RC inequalities
    std::pair<bool, double> separate_rcc(
        float eps_for_int,
        float eps_for_viol,
        float eps_for_early_exact_termination,
        float type2_max,
        float type3_min,
        float type3_max,
        int gcd,
        bool heuristic,
        bool exact,
        int rcc_implementation,
        bool shrink_rcc,
        CnstrMgrPointer MyOldCutsCMP,
        CnstrMgrPointer MyCutsCMP,
        GRBEnv env
    ){
        auto start_sep = std::chrono::high_resolution_clock::now();

        double max_violation = 0.0;
        std::vector<RC_Inequality> new_cuts;
        std::set<int> shrunk_customers;
        std::vector<int> shrunk_demand;
        ReachPtr shrunk_mapping = nullptr;
        std::vector<std::tuple<int, int, double>> shrunk_solution;
        std::set<int> shrunk_vertices;        

        std::vector<std::tuple<int, int, double>> solution = get_solution(eps_for_int);

        if (heuristic) {
            auto start_heur = std::chrono::high_resolution_clock::now();

            separate_rcc_heuristically(
                new_cuts, max_violation, solution, dimension, demand, capacity,
                eps_for_int, eps_for_viol, MyOldCutsCMP, MyCutsCMP);

            auto stop_heur = std::chrono::high_resolution_clock::now();
            auto duration_heur = std::chrono::duration_cast<std::chrono::duration<double>>(stop_heur - start_heur);
            RCC_heuristic_time += duration_heur.count();
        }
        if (exact && max_violation <= 0.0f) {
            auto start_exact = std::chrono::high_resolution_clock::now();

            if (shrink_rcc) {
                ReachInitMem(&shrunk_mapping, dimension + 1);
                shrink_rcc_graph(
                    shrunk_mapping, shrunk_demand, shrunk_vertices, shrunk_customers, shrunk_solution,
                    solution, dimension, demand, capacity, RCC_shrink_time, MyOldCutsCMP);

                int shrunk_gcd = capacity;
                for (int i = 1; i < demand.size(); ++i) {
                    shrunk_gcd = std::gcd(gcd, demand[i]); // Update gcd with the gcd of current element and previous gcd
                }

                separate_rcc_exactly(
                    new_cuts, max_violation,
                    shrunk_solution, shrunk_customers, capacity, shrunk_demand,
                    eps_for_early_exact_termination, shrunk_gcd, eps_for_int, eps_for_viol, RCC_exact_solve_time, env, shrunk_mapping);
                convert_cuts_rcc(new_cuts, shrunk_mapping, RCC_convert_shrunk_cuts_time);
            }
            else {
                separate_rcc_exactly(new_cuts, max_violation,
                    solution, customers, capacity, demand,
                    eps_for_early_exact_termination, gcd, eps_for_int, eps_for_viol, RCC_exact_solve_time, env, shrunk_mapping);
            }

            auto stop_exact = std::chrono::high_resolution_clock::now();
            auto duration_exact = std::chrono::duration_cast<std::chrono::duration<double>>(stop_exact - start_exact);
            RCC_exact_time += duration_exact.count();
        }
        auto start_add = std::chrono::high_resolution_clock::now();

        add_rc_cuts(new_cuts, model, x, dimension, total_rcc, current_rcc, vertices, customers,edges,
            edges_subset, var_columns, type2_max, type3_min, type3_max, rcc_implementation, eps_for_viol, determine_RCC_time, MyCutsCMP);

        auto stop_add = std::chrono::high_resolution_clock::now();
        auto duration_add = std::chrono::duration_cast<std::chrono::duration<double>>(stop_add - start_add);
        add_RCC_time += duration_add.count();

        auto stop_sep = std::chrono::high_resolution_clock::now();
        auto duration_sep = std::chrono::duration_cast<std::chrono::duration<double>>(stop_sep - start_sep);
        RCC_separation_time += duration_sep.count();
        return { (MyCutsCMP->Size > 0), max_violation};
    }



    // Separates SC inequalities
    std::pair<bool, double> separate_strengthened_comb_cuts(
        float eps_for_int,
        float eps_for_viol,
        float eps_for_early_termination,
        int gcd,
        int Qmin,
        bool heuristic,
        int max_exact_teeth,
        bool shrink,
        CnstrMgrPointer MyOldCutsCMP,
        CnstrMgrPointer MyCutsCMP,
        GRBEnv env
    ) {
        auto start_sep = std::chrono::high_resolution_clock::now();

        double max_violation = 0.0;
        std::vector<SC_Inequality> new_cuts;
        int no_new_cuts;
        std::vector<std::tuple<int, int, double>> shrunk_solution;
        std::set<int> shrunk_vertices;
        std::set<int> shrunk_customers;
        std::vector<int> shrunk_demand;
        ReachPtr shrunk_mapping = nullptr;
        int shrunk_graph_size;
        
        std::vector<std::tuple<int, int, double>> solution = get_solution(eps_for_int);
        if (heuristic) {
            auto start_heur = std::chrono::high_resolution_clock::now();

            separate_comb_heuristically(max_violation, new_cuts, solution, dimension, demand, capacity, Qmin, MyOldCutsCMP, MyCutsCMP);

            auto stop_heur = std::chrono::high_resolution_clock::now();
            auto duration_heur = std::chrono::duration_cast<std::chrono::duration<double>>(stop_heur - start_heur);
            SC_heuristic_time += duration_heur.count();
        }
        if (max_exact_teeth >= 2 && max_violation <= 0.0f) {
            auto start_exact = std::chrono::high_resolution_clock::now();

            if (shrink) {
                ReachInitMem(&shrunk_mapping, dimension + 1);

                shrink_comb_graph(
                    shrunk_mapping, shrunk_graph_size, shrunk_demand, shrunk_vertices, shrunk_customers, shrunk_solution,
                    solution, dimension, demand, capacity, Qmin
                );

                int shrunk_gcd = capacity;
                for (int i = 1; i < demand.size(); ++i) {
                    shrunk_gcd = std::gcd(gcd, demand[i]); // Update gcd with the gcd of current element and previous gcd
                }

                separate_sc_exactly(max_violation, new_cuts, shrunk_mapping,
                    shrunk_solution, edges_subset, shrunk_vertices, shrunk_customers, capacity, shrunk_demand, dimension, eps_for_early_termination, shrunk_gcd, max_exact_teeth, eps_for_int, eps_for_viol, SC_exact_solve_time, env);
                convert_cuts_comb(new_cuts, shrunk_mapping, shrunk_graph_size, dimension);
            }
            else {
                separate_sc_exactly(max_violation, new_cuts, shrunk_mapping, solution, edges_subset, vertices, customers, capacity,
                    demand, dimension, eps_for_early_termination, gcd, max_exact_teeth, eps_for_int, eps_for_viol, SC_exact_solve_time, env);
            }

            auto stop_exact = std::chrono::high_resolution_clock::now();
            auto duration_exact = std::chrono::duration_cast<std::chrono::duration<double>>(stop_exact - start_exact);
            SC_exact_time += duration_exact.count();
        }
        auto start_add = std::chrono::high_resolution_clock::now();

        add_strenghtened_comb_cut(no_new_cuts, new_cuts, model, x, dimension, total_combs, current_combs,
            vertices, edges, edges_subset, var_columns, eps_for_viol, determine_SC_time, MyCutsCMP);

        auto stop_add = std::chrono::high_resolution_clock::now();
        auto duration_add = std::chrono::duration_cast<std::chrono::duration<double>>(stop_add - start_add);
        add_SC_time += duration_add.count();

        auto stop_sep = std::chrono::high_resolution_clock::now();
        auto duration_sep = std::chrono::duration_cast<std::chrono::duration<double>>(stop_sep - start_sep);
        SC_separation_time += duration_sep.count();
        return { (no_new_cuts > 0), max_violation };
    }


    // Separates HPM inequalities
    std::pair<bool, double> separate_mstars(
        float eps_for_int,
        float eps_for_viol,
        float eps_for_early_exact_termination,
        float type2_max,
        float type3_min,
        float type3_max,
        int gcd,
        bool heuristic,
        bool exact,
        bool shrink_graph,
        CnstrMgrPointer MyOldCutsCMP,
        CnstrMgrPointer MyCutsCMP,
        GRBEnv env
    ) {
        auto start_sep = std::chrono::high_resolution_clock::now();

        double max_violation = 0.0;
        std::vector<MSTAR_Inequality> new_cuts;
        int no_new_cuts = 0;
        std::set<int> shrunk_customers;
        std::vector<int> shrunk_demand;
        ReachPtr shrunk_mapping = nullptr;
        std::vector<std::tuple<int, int, double>> shrunk_solution;
        std::set<int> shrunk_vertices;

        std::vector<std::tuple<int, int, double>> solution = get_solution(eps_for_int);

        if (heuristic && true) {
            auto start_heur = std::chrono::high_resolution_clock::now();

            separate_mstar_heuristically(max_violation, new_cuts, solution, dimension, demand, capacity, MyOldCutsCMP, MyCutsCMP);

            auto stop_heur = std::chrono::high_resolution_clock::now();
            auto duration_heur = std::chrono::duration_cast<std::chrono::duration<double>>(stop_heur - start_heur);
            MSTAR_heuristic_time += duration_heur.count();
        }
        if (exact && max_violation <= 0.0f) {
            auto start_exact = std::chrono::high_resolution_clock::now();

            if (shrink_graph) {
                ReachInitMem(&shrunk_mapping, dimension + 1);
                shrink_rcc_graph(
                    shrunk_mapping, shrunk_demand, shrunk_vertices, shrunk_customers, shrunk_solution,
                    solution, dimension, demand, capacity, RCC_shrink_time, MyOldCutsCMP);

                int shrunk_gcd = capacity;
                for (int i = 1; i < demand.size(); ++i) {
                    shrunk_gcd = std::gcd(gcd, demand[i]); // Update gcd with the gcd of current element and previous gcd
                }

                separate_mstar_exactly(
                    max_violation, new_cuts, shrunk_mapping,
                    shrunk_solution, shrunk_vertices, shrunk_customers, capacity, shrunk_demand, dimension,
                    eps_for_early_exact_termination, shrunk_gcd, eps_for_int, eps_for_viol, MSTAR_exact_solve_time, env);
                convert_cuts_mstar(new_cuts, shrunk_mapping);
            }
            else {
                separate_mstar_exactly(
                    max_violation, new_cuts, shrunk_mapping,
                    solution, vertices, customers, capacity, demand, dimension,
                    eps_for_early_exact_termination, gcd, eps_for_int, eps_for_viol, MSTAR_exact_solve_time, env);
            }

            auto stop_exact = std::chrono::high_resolution_clock::now();
            auto duration_exact = std::chrono::duration_cast<std::chrono::duration<double>>(stop_exact - start_exact);
            MSTAR_exact_time += duration_exact.count();
        }
        auto start_add = std::chrono::high_resolution_clock::now();

        add_mstars(no_new_cuts, new_cuts, model, x, dimension, total_mstars, current_mstars, vertices, customers, edges,
            edges_subset, var_columns, type2_max, type3_min, type3_max, eps_for_viol, determine_MSTAR_time, MyCutsCMP);

        auto stop_add = std::chrono::high_resolution_clock::now();
        auto duration_add = std::chrono::duration_cast<std::chrono::duration<double>>(stop_add - start_add);
        add_MSTAR_time += duration_add.count();

        auto stop_sep = std::chrono::high_resolution_clock::now();
        auto duration_sep = std::chrono::duration_cast<std::chrono::duration<double>>(stop_sep - start_sep);
        MSTAR_separation_time += duration_sep.count();
        return { (no_new_cuts > 0), max_violation };
    }


    // Separates GLM inequalities
    std::pair<bool, double> separate_glm(
        float eps_for_int,
        float eps_for_viol,
        float type1_max
    ) {
        auto start_sep = std::chrono::high_resolution_clock::now();

        double max_violation = 0.0;
        std::vector<GLM_Inequality> new_cuts;
        int no_new_cuts = 0;

        std::vector<std::tuple<int, int, double>> solution = get_solution(eps_for_int);

        auto start_heur = std::chrono::high_resolution_clock::now();
        separate_glm_CVRPSEP(max_violation, new_cuts, solution, dimension, demand, capacity);
        auto stop_heur = std::chrono::high_resolution_clock::now();
        auto duration_heur = std::chrono::duration_cast<std::chrono::duration<double>>(stop_heur - start_heur);
        GLM_separation_time += duration_heur.count();

        auto start_add = std::chrono::high_resolution_clock::now();
        add_glm(no_new_cuts, new_cuts, model, x, dimension, total_glm, current_glm, vertices, customers, demand, capacity, edges,
            edges_subset, var_columns, type1_max, eps_for_viol, determine_GLM_time);
        auto stop_add = std::chrono::high_resolution_clock::now();
        auto duration_add = std::chrono::duration_cast<std::chrono::duration<double>>(stop_add - start_add);
        add_GLM_time += duration_add.count();

        auto stop_sep = std::chrono::high_resolution_clock::now();
        auto duration_sep = std::chrono::duration_cast<std::chrono::duration<double>>(stop_sep - start_sep);
        GLM_total_separation_time += duration_sep.count();
        return { (no_new_cuts > 0), max_violation };
    }


    // Main separation routine that separates all types if inequalities and adds them to the model
    std::tuple<bool, bool> separate_cuts(
        GRBEnv env, CnstrMgrPointer MyCutsCMP, CnstrMgrPointer MyOldCutsCMP,
        auto& start_time,
        auto& last_print_time,
        double Time_Limit,
        int current_sep_iter,
        float eps_for_int,
        float slack_eps,
        int remove_during_sep,
        float eps_for_rcc_viol,
        float eps_for_sc_violation,
        float eps_for_early_sc_termination,
        bool rcc_heu,
        bool rcc_exact,
        float eps_for_early_exact_termination,
        int rcc_implementation,
        bool iterate,
        bool shrink_rcc,
        float& rcc2_max,
        float& rcc3_min,
        float& rcc3_max,
        int& gcd,
        bool comb_heur,
        int max_exact_sc_teeth,
        bool shrink_comb,
        int Qmin,
        bool mstar_heur,
        bool mstar_exact,
        bool& shrink_mstar,
        float& eps_for_mstar_violation,
        float& eps_for_early_mstar_termination,
        bool glm_sep,
        float& eps_for_glm_violation, bool stop_when_tailoff
    ) {
        bool continue_separation;
        int beginning_cuts = total_rcc + total_combs + total_mstars + total_glm;
        bool time_limit_reached = false;
        std::pair<bool, double> new_cuts;

        int no_of_extra_cut_types = 3;
        int current_extra_cut_type = 0;
        
        while (true) {
            double max_rcc_violation = 0.0;
            double max_sc_violation = 0.0;
            double max_mstar_violation = 0.0;
            double max_glm_violation = 0.0;
            continue_separation = false;
            inner_iterations++;
            solve_model();

            auto current_time = std::chrono::steady_clock::now();
            double total_time = std::chrono::duration_cast<std::chrono::duration<double>>(current_time - start_time).count();
            if (Time_Limit > 0 && total_time >= Time_Limit) {
                time_limit_reached = true;
                break;
            }

            //double total_time = std::chrono::duration_cast<std::chrono::duration<double>>(current_time - start_time).count();
            //std::cout << std::fixed << std::setprecision(2) << get_objective_value() / upper_bound * 100 << ";";
            //std::cout << std::fixed << std::setprecision(2) << total_time << "\n";
            
            // Check if 5 seconds have passed since the last print
            if (current_sep_iter == inner_iterations - 1) {
                std::cout << "CG Iteration " << col_gen_iters << ": ";
                std::cout << std::fixed << std::setprecision(2) << "Obj - " << get_objective_value() << " ";
                std::cout << std::fixed << std::setprecision(0) << "Time - " << total_time << "s\n";
                //std::cout << total_rcc << " " << total_combs << " " << total_mstars << std::endl;
            }
            else if (std::chrono::duration_cast<std::chrono::seconds>(current_time - last_print_time).count() >= 5) {
                std::cout << std::fixed << std::setprecision(2) << "Obj - " << get_objective_value() << " ";
                std::cout << std::fixed << std::setprecision(2) << get_objective_value() / upper_bound * 100 << ";";
                std::cout << std::fixed << std::setprecision(0) << "Time - " << total_time << "s\n";
                last_print_time = current_time; // Update the last print time
            }

            if (remove_during_sep > 0 && inner_iterations % remove_during_sep == 0) {
                remove_constraints_by_slack(slack_eps, dimension, model, edges_not_subset, var_columns, current_rcc, current_combs, current_mstars, current_glm, remove_constrs_time);
            }

            if (rcc_heu || rcc_exact) {
                new_cuts = separate_rcc(eps_for_int, eps_for_rcc_viol, eps_for_early_exact_termination,
                    rcc2_max, rcc3_min, rcc3_max, gcd,
                    rcc_heu, rcc_exact, rcc_implementation, shrink_rcc,
                    MyOldCutsCMP, MyCutsCMP, env);
                continue_separation = new_cuts.first;
                max_rcc_violation = new_cuts.second;
                if (stop_when_tailoff) continue_separation = max_rcc_violation >= 0.2f;
                //std::cout << std::setprecision(2) << "RCC " << max_rcc_violation << std::endl;
            }

            if (max_rcc_violation < 0.2f && (glm_sep || mstar_heur || mstar_exact || comb_heur || max_exact_sc_teeth > 0)) {
                int cut_types_separated = 0;
                while (cut_types_separated < no_of_extra_cut_types) {
                    if (current_extra_cut_type == 0) {
                        cut_types_separated++;
                        current_extra_cut_type++;
                        current_extra_cut_type %= no_of_extra_cut_types;
                        if (glm_sep) {
                            new_cuts = separate_glm(eps_for_int, eps_for_glm_violation, rcc2_max);
                            max_glm_violation = new_cuts.second;
                            continue_separation = (continue_separation || max_glm_violation >= 0.05f);
                        }
                        //std::cout << std::setprecision(2) << max_glm_violation << std::endl;
                    }
                    if (max_glm_violation >= 0.05f) break;

                    if (current_extra_cut_type == 1) {
                        cut_types_separated++;
                        current_extra_cut_type++;
                        current_extra_cut_type %= no_of_extra_cut_types;
                        if (mstar_heur || mstar_exact) {
                            new_cuts = separate_mstars(
                                eps_for_int, eps_for_mstar_violation, eps_for_early_mstar_termination,
                                rcc2_max, rcc3_min, rcc3_max, gcd,
                                mstar_heur, mstar_exact, shrink_mstar,
                                MyOldCutsCMP, MyCutsCMP, env);
                            max_mstar_violation = new_cuts.second;
                            continue_separation = (continue_separation || max_mstar_violation >= 0.05f);
                        }
                        //std::cout << std::setprecision(2) << "MSTAR " << max_mstar_violation << std::endl;
                    }
                    if (max_mstar_violation >= 0.05f) break;

                    if (current_extra_cut_type == 2) {
                        cut_types_separated++;
                        current_extra_cut_type++;
                        current_extra_cut_type %= no_of_extra_cut_types;
                        if (comb_heur || max_exact_sc_teeth > 0) {
                            new_cuts = separate_strengthened_comb_cuts(
                                eps_for_int, eps_for_sc_violation, eps_for_early_sc_termination,
                                gcd, Qmin, comb_heur, max_exact_sc_teeth, shrink_comb,
                                MyOldCutsCMP, MyCutsCMP, env);
                            max_sc_violation = new_cuts.second;
                            continue_separation = (continue_separation || max_sc_violation >= 0.1);
                        }
                        //std::cout << std::setprecision(2) << "SC " << max_mstar_violation << std::endl;
                    }                    
                    if (max_sc_violation >= 0.1f) break;
                }
            }
            if (!continue_separation) break;
            //std::cout << total_rcc << " " << total_combs << " " << total_mstars << " " << total_glm << std::endl;

            for (int i = 0; i < MyCutsCMP->Size; i++)
            {
                CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
            }
            MyCutsCMP->Size = 0;

            if (!iterate) break;
        }
        return { (total_rcc + total_combs + total_mstars + total_glm> beginning_cuts), time_limit_reached };
    }


    // Main routine for solving the (LP relaxation) of an instance
    bool solve_instance(
        GRBEnv env, CnstrMgrPointer MyCutsCMP, CnstrMgrPointer MyOldCutsCMP,
        auto& start_time, double time_limit,
        bool exact_if_terminate,
        bool iterate_sep,
        float eps_for_int,
        int col_gen_type, int no_col_gen_edges, float eps_for_col_gen,
        float slack_eps,
        int remove_during_sep,
        float eps_for_elim,
        int rcc_implementation, bool rcc_heu, bool rcc_exact, bool shrink_rcc, float eps_for_rcc_viol, float eps_for_early_rcc_termination,
        bool comb_heur, int max_exact_sc_teeth, bool shrink_comb, float eps_for_sc_violation, float eps_for_early_sc_termination, bool sc_only_when_terminate,
        bool mstar_heur, bool mstar_exact, bool shrink_mstar, float eps_for_mstar_violation, float eps_for_early_mstar_termination,
        bool mstar_only_when_terminate, bool glm_sep, float eps_for_glm_violation, bool glm_only_when_terminate
    ) {
        auto start = std::chrono::high_resolution_clock::now();
        auto last_print_time = start;

        // Set standard values for non-necesaary parameters
        if (!column_generation || !iterate_sep) exact_if_terminate = false;
        if (eps_for_elim > 0) col_gen_type = 0;
        if (col_gen_type < 2) no_col_gen_edges = 0;
        if (no_col_gen_edges == -1) no_col_gen_edges = dimension;
        if (slack_eps <= 0.0f) remove_during_sep = 0;

        bool new_cuts;
        bool new_vars = false;
        bool time_limit_reached;
        std::tuple<bool, bool> cuts;

        int gcd = capacity;
        for (int i = 1; i < demand.size(); ++i) {
            gcd = std::gcd(gcd, demand[i]); // Update gcd with the gcd of current element and previous gcd
        }

        float rcc2_max = int((dimension + 1) / 2.0f);
        float sqrt_val = float(std::sqrt(dimension * dimension - 10 * dimension + 9));
        float rcc3_min = (-sqrt_val + 3 * dimension - 3) / 4;
        float rcc3_max = (sqrt_val + 3 * dimension - 3) / 4;
              
        int Qmin = *std::min_element(demand.begin() + 1, demand.end());


        int times_exact = 0;
        while (true) {
            col_gen_iters++;
            int current_sep_iter = inner_iterations;
            if (!exact_if_terminate) {
                cuts = separate_cuts(env, MyCutsCMP, MyOldCutsCMP, start_time, last_print_time,
                    time_limit, current_sep_iter,
                    eps_for_int, slack_eps, remove_during_sep, eps_for_rcc_viol,
                    eps_for_sc_violation, eps_for_early_sc_termination,
                    rcc_heu, rcc_exact, eps_for_early_rcc_termination,
                    rcc_implementation, iterate_sep, shrink_rcc, rcc2_max,
                    rcc3_min, rcc3_max, gcd, comb_heur,
                    max_exact_sc_teeth, shrink_comb, Qmin,
                    mstar_heur, mstar_exact, shrink_mstar, eps_for_mstar_violation, eps_for_early_mstar_termination,
                    glm_sep, eps_for_glm_violation, false);
            }
            else {
                cuts = separate_cuts(env, MyCutsCMP, MyOldCutsCMP, start_time, last_print_time,
                    time_limit, current_sep_iter,
                    eps_for_int, slack_eps, remove_during_sep, eps_for_rcc_viol,
                    eps_for_sc_violation, eps_for_early_sc_termination, 
                    rcc_heu, false, eps_for_early_rcc_termination,
                    rcc_implementation, iterate_sep, shrink_rcc, rcc2_max,
                    rcc3_min, rcc3_max, gcd, comb_heur && !sc_only_when_terminate,
                    0, shrink_comb, Qmin,
                    mstar_heur && !mstar_only_when_terminate, false, shrink_mstar, eps_for_mstar_violation, eps_for_early_mstar_termination,
                    glm_sep && !glm_only_when_terminate, eps_for_glm_violation, false);
            }

            new_cuts = std::get<0>(cuts);
            time_limit_reached = std::get<1>(cuts);
            if (time_limit_reached) break;

            if (column_generation) {
                if (!iterate_sep || slack_eps > 0) solve_model();
                if (slack_eps > 0) {
                    remove_constraints_by_slack(slack_eps, dimension, model, edges_not_subset, var_columns, current_rcc, current_combs, current_mstars, current_glm, remove_constrs_time);
                }
                if (flow_variables) {
                    new_vars = add_variables_with_flow(eps_for_col_gen, edges_subset, edges_not_subset, asym_edges, x, f, customers, add_vars_time);
                }
                else {
                    new_vars = add_variables(
                        col_gen_type, no_col_gen_edges, eps_for_col_gen,
                        edges_not_subset, var_columns, model, x, edges_subset, cost, add_vars_time);
                }
            }

            if (!new_vars) {
                if (iterate_sep || !new_cuts) {
                    if (((mstar_heur || mstar_exact) && mstar_only_when_terminate)
                            || ((comb_heur || max_exact_sc_teeth > 0) && sc_only_when_terminate)
                            || (glm_sep && glm_only_when_terminate)
                            ) {
                        std::cout << "OTHER sep" << std::endl;
                        cuts = separate_cuts(env, MyCutsCMP, MyOldCutsCMP, start_time, last_print_time,
                            time_limit, current_sep_iter,
                            eps_for_int, slack_eps, remove_during_sep, eps_for_rcc_viol,
                            eps_for_sc_violation, eps_for_early_sc_termination,
                            rcc_heu, false, eps_for_early_rcc_termination,
                            rcc_implementation, iterate_sep, shrink_rcc, rcc2_max,
                            rcc3_min, rcc3_max, gcd, comb_heur,
                            0, shrink_comb, Qmin,
                            mstar_heur, false, shrink_mstar, eps_for_mstar_violation, eps_for_early_mstar_termination,
                            glm_sep, eps_for_glm_violation, true);
                        new_cuts = std::get<0>(cuts);
                        time_limit_reached = std::get<1>(cuts);
                    }
                    if (exact_if_terminate && !new_cuts) {
                        times_exact++;
                        std::cout << "Exact" << std::endl;
                        cuts = separate_cuts(env, MyCutsCMP, MyOldCutsCMP, start_time, last_print_time,
                            time_limit, current_sep_iter,
                            eps_for_int, slack_eps, remove_during_sep, eps_for_rcc_viol,
                            eps_for_sc_violation, eps_for_early_sc_termination, 
                            rcc_heu, rcc_exact, eps_for_early_rcc_termination,
                            rcc_implementation, iterate_sep, shrink_rcc, rcc2_max,
                            rcc3_min, rcc3_max, gcd, comb_heur && !sc_only_when_terminate,
                            max_exact_sc_teeth * (1 - sc_only_when_terminate), shrink_comb, Qmin,
                            mstar_heur && !mstar_only_when_terminate, mstar_exact && !mstar_only_when_terminate,
                            shrink_mstar, eps_for_mstar_violation, eps_for_early_mstar_termination,
                            glm_sep && !glm_only_when_terminate, eps_for_glm_violation, false);
                        new_cuts = std::get<0>(cuts);
                        time_limit_reached = std::get<1>(cuts);
                    }
                    if (!new_cuts || !exact_if_terminate || time_limit_reached) break;
                }
            }
            if (eps_for_elim > 0) {
                eliminate_variables(model, x, edges, edges_subset, edges_not_subset, cost, var_columns, upper_bound,
                    eps_for_elim, eps_for_int, variables_eliminated, eliminate_vars_time);
            }

            auto current_time = std::chrono::steady_clock::now();
            double total_time = std::chrono::duration_cast<std::chrono::duration<double>>(current_time - start_time).count();
            if (time_limit > 0 && total_time >= time_limit) {
                time_limit_reached = true;
                break;
            }
        }

        if (time_limit_reached) {
            std::cout << "Time Limit has been reached\n";
            add_variables(
                0, 0, eps_for_col_gen,
                edges_not_subset, var_columns, model, x, edges_subset, cost, add_vars_time);
        }


        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
        separation_time += duration.count();

        return time_limit_reached;
    }


    // Prints the times of several stages of the solution process
    void print_times() {
        std::cout << std::fixed << std::setprecision(2) << "Total time: " << total_time << " seconds\n";

        std::cout << std::fixed << std::setprecision(2) << "  Setup time: " << setup_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "  Starting Subset time: " << starting_subset_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "  Separation time: " << separation_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "    LP solve time: " << LP_solve_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "    Remove Constrs time: " << remove_constrs_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "    Add vars time: " << add_vars_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "    Elim vars time: " << eliminate_vars_time << " seconds\n";

        std::cout << std::fixed << std::setprecision(2) << "    RCC separation time: " << RCC_separation_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      RCC heuristic time: " << RCC_heuristic_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      RCC exact time: " << RCC_exact_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "        Shrink time: " << RCC_shrink_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "        Solve time: " << RCC_exact_solve_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "        Convert Cuts time: " << RCC_convert_shrunk_cuts_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      Add RCC time: " << add_RCC_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "        Determine RCC time: " << determine_RCC_time << " seconds\n";

        std::cout << std::fixed << std::setprecision(2) << "    SC separation time: " << SC_separation_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      SC heuristic time: " << SC_heuristic_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      SC exact time: " << SC_exact_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "        Solve time: " << SC_exact_solve_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      Add SC time: " << add_SC_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "        Determine SC time: " << determine_SC_time << " seconds\n";

        std::cout << std::fixed << std::setprecision(2) << "    MSTAR separation time: " << MSTAR_separation_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      MSTAR heuristic time: " << MSTAR_heuristic_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      MSTAR exact time: " << MSTAR_exact_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "        Solve time: " << MSTAR_exact_solve_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      Add MSTAR time: " << add_MSTAR_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "        Determine MSTAR time: " << determine_MSTAR_time << " seconds\n";

        std::cout << std::fixed << std::setprecision(2) << "    GLM separation time: " << GLM_separation_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      GLM exact time: " << GLM_separation_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "      Add GLM time: " << add_GLM_time << " seconds\n";
        std::cout << std::fixed << std::setprecision(2) << "        Determine GLM time: " << determine_GLM_time << " seconds\n";
    }


    // Writes the optimal solution to a file
    void write_opt_solution() {
        std::ofstream myfile;
        myfile.open("Results/Optimal Solutions/" + name + ".sol");
        std::vector< std::tuple<int, int, double>> sol = get_solution();

        std::set<std::tuple<int, int, double>> variables_left(sol.begin(), sol.end());
        int route = 1;
        while (variables_left.size() > 0) {
            myfile << "Route #" << route << ": ";
            int current_city = 0;
            int tour_length = 0;
            while (current_city != 0 || tour_length == 0) {
                for (auto edge : variables_left) {
                    int i = std::get<0>(edge);
                    int j = std::get<1>(edge);
                    double sol_val = std::get<2>(edge);
                    if (i == current_city || j == current_city) {
                        if (sol_val == 2.0) {
                            current_city == 0;
                            myfile << j << " ";
                        }
                        else current_city = (i == current_city) ? j : i;

                        if (current_city > 0) myfile << current_city << " ";
                        tour_length++;
                        variables_left.erase(edge);
                        break;
                    }
                }
            }
            myfile << "\n";
            route++;
        }
        myfile << "Cost " << int(model->get(GRB_DoubleAttr_ObjVal));
        myfile.close();
    }
};

