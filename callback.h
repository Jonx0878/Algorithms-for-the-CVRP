#pragma once

#include "gurobi_c++.h"
#include "rcc_sep.h"
#include "glm.h"
#include "comb.h"
#include "mstar.h"


class mycallback : public GRBCallback {
public:
    GRBModel* model;
    std::set<int> vertices;
    std::set<int> customers;
    std::set<std::tuple<int, int>> edges;
    std::set<std::tuple<int, int>> edges_subset;
    std::vector<std::vector<GRBVar>> x;
    float eps_for_integrality;
    float eps_for_rcc_violation;
    int dimension;
    std::vector<int> demand;
    int capacity;
    CnstrMgrPointer MyOldCutsCMP;
    CnstrMgrPointer MyCutsCMP;

    int total_rcc;
    int current_rcc;

    int rcc_implementation;
    float type2_max;
    float type3_min;
    float type3_max;

    std::vector<std::vector<GRBColumn>> var_columns;
    double determine_RCC_time;

    int Qmin;

    mycallback(
        GRBModel* m,
        std::set<int>& vert,
        std::set<int>& cust,
        std::set<std::tuple<int, int>>& e,
        std::set<std::tuple<int, int>>& e_subset,
        std::vector<std::vector<GRBVar>>& vars,
        float& eps_for_int,
        float& eps_for_rcc_viol,
        int& dim,
        std::vector<int>& dem,
        int& cap,
        CnstrMgrPointer OldCMP,
        CnstrMgrPointer NewCMP,
        int& total_rc,
        int& current_rc,
        int& rc_implementation,
        std::vector<std::vector<GRBColumn>>& columns,
        double& determine_RC_time
    ) {
        model = m;
        vertices = vert;
        customers = cust;
        edges = e;
        edges_subset = e_subset;
        x = vars;
        eps_for_integrality = eps_for_int;
        eps_for_rcc_violation = eps_for_rcc_viol;
        dimension = dim;
        demand = dem;
        capacity = cap;
        MyOldCutsCMP = OldCMP;
        MyCutsCMP = NewCMP;

        total_rcc = total_rc;
        current_rcc = current_rc;

        rcc_implementation = rc_implementation;

        var_columns = columns;
        determine_RCC_time = determine_RC_time;

        type2_max = int((dimension + 1) / 2.0f);
        float sqrt_val = float(std::sqrt(dimension * dimension - 10 * dimension + 9));
        type3_min = (-sqrt_val + 3 * dimension - 3) / 4;
        type3_max = (sqrt_val + 3 * dimension - 3) / 4;

        Qmin = *std::min_element(demand.begin() + 1, demand.end());
    }

    mycallback() {};

protected:
    // Callback separating only RC inequalities at incumbents whilst separating all inequalities once at each node
    void callback() {
        if (where == GRB_CB_MIPNODE) {
            if (getIntInfo(GRB_CB_MIPNODE_STATUS) == 2) {
                std::vector<std::tuple<int, int, double>> solution;
                for (auto& edge : edges_subset) {
                    int i = std::get<0>(edge);
                    int j = std::get<1>(edge);
                    double sol_val = getNodeRel(x[i][j]);
                    if (sol_val >= eps_for_integrality) {
                        solution.push_back(std::make_tuple(i, j, sol_val));
                    }
                }
                double max_rcc_viol = add_RC(solution);
                add_GLM(solution);
                add_HPM(solution);
                add_SC(solution);
            }
        } 
        if (where == GRB_CB_MIPSOL) {
            std::vector<std::tuple<int, int, double>> solution;
            for (auto& edge : edges_subset) {
                int i = std::get<0>(edge);
                int j = std::get<1>(edge);
                double sol_val = getSolution(x[i][j]);
                if (sol_val >= eps_for_integrality) {
                    solution.push_back(std::make_tuple(i, j, sol_val));
                }
            }
            double max_rcc_viol = add_RC(solution);

        }

        for (int i = 0; i < MyCutsCMP->Size; i++)
        {
            CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
        }
        MyCutsCMP->Size = 0;
    }

    double add_RC(std::vector<std::tuple<int, int, double>> solution) {
        double max_violation = 0.0;
        std::vector<RC_Inequality> new_cuts;
        separate_rcc_heuristically(new_cuts, max_violation, solution, dimension, demand, capacity,
            eps_for_integrality, eps_for_rcc_violation, MyOldCutsCMP, MyCutsCMP);

        for (const auto& cut : new_cuts) {
            GRBLinExpr lhs;
            char sense;
            int rhs = cut.rhs;
            std::set<int> set = cut.set;
            bool check_violation = cut.check_violation;
            bool add_to_cnstrmgr = cut.add_to_cnstrmgr;
            std::set<std::tuple<int, int, double>> edges_not_subset;

            auto start = std::chrono::high_resolution_clock::now();

            determine_rcc(lhs, sense, rhs, set, x, vertices, customers,
                edges, edges_subset, edges_not_subset,
                type2_max, type3_min, type3_max, rcc_implementation);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
            determine_RCC_time += duration.count();

            addLazy(lhs, sense, rhs);
        }
        return max_violation;
    }


    void add_GLM(std::vector<std::tuple<int, int, double>> solution) {
        double max_violation = 0.0;
        std::vector<GLM_Inequality> new_cuts;
        separate_glm_CVRPSEP(max_violation, new_cuts, solution, dimension, demand, capacity);

        for (const auto& glm : new_cuts) {
            GRBLinExpr lhs;
            char sense = GRB_LESS_EQUAL;
            int rhs = 0;
            bool check_violation = glm.check_violation;
            std::set<std::tuple<int, int, double>> edges_not_subset;

            auto start = std::chrono::high_resolution_clock::now();

            determine_glm(lhs, sense, rhs, glm, x, vertices, customers,
                edges, edges_subset, edges_not_subset, demand, capacity, type2_max);
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
            //time += duration.count();

            addLazy(lhs, sense, rhs);
        }
    }


    void add_HPM(std::vector<std::tuple<int, int, double>> solution) {
        double max_violation = 0.0;
        std::vector<MSTAR_Inequality> new_cuts;
        separate_mstar_heuristically(max_violation, new_cuts, solution, dimension, demand, capacity, MyOldCutsCMP, MyCutsCMP);

        for (const auto& mstar : new_cuts) {
            GRBLinExpr lhs;
            char sense = GRB_LESS_EQUAL;
            int rhs = 0;
            bool check_violation = mstar.check_violation;
            bool add_to_cnstrmgr = mstar.from_exact;
            std::set<std::tuple<int, int, double>> edges_not_subset;

            auto start = std::chrono::high_resolution_clock::now();

            determine_mstar(lhs, sense, rhs, mstar, x, vertices, customers,
                edges, edges_subset, edges_not_subset,
                type2_max, type3_min, type3_max);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
            //time += duration.count();

            addLazy(lhs, sense, rhs);
        }
    }


    void add_SC(std::vector<std::tuple<int, int, double>> solution) {
        double max_violation = 0.0;
        std::vector<SC_Inequality> new_cuts;
        separate_comb_heuristically(max_violation, new_cuts, solution, dimension, demand, capacity, Qmin, MyOldCutsCMP, MyCutsCMP);

        for (const auto& cut : new_cuts) {
            GRBLinExpr lhs;
            char sense = GRB_GREATER_EQUAL;
            int rhs = cut.rhs;
            std::vector<std::set<int>> teeth = cut.teeth;
            std::set<std::tuple<int, int, double>> edges_not_subset;

            auto start = std::chrono::high_resolution_clock::now();

            determine_strengthened_comb_cut(lhs, teeth, model, x, vertices, edges, edges_subset, edges_not_subset);
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
            //time += duration.count();

            addLazy(lhs, sense, rhs);
        }
    }
};



class my_rc_callback : public GRBCallback {
public:
    GRBModel* model;
    std::set<int> vertices;
    std::set<int> customers;
    std::set<std::tuple<int, int>> edges;
    std::set<std::tuple<int, int>> edges_subset;
    std::vector<std::vector<GRBVar>> x;
    float eps_for_integrality;
    float eps_for_rcc_violation;
    int dimension;
    std::vector<int> demand;
    int capacity;
    CnstrMgrPointer MyOldCutsCMP;
    CnstrMgrPointer MyCutsCMP;

    int total_rcc;
    int current_rcc;

    int rcc_implementation;
    float type2_max;
    float type3_min;
    float type3_max;

    std::vector<std::vector<GRBColumn>> var_columns;
    double determine_RCC_time;

    bool mipnode_callback;

    my_rc_callback(
        GRBModel* m,
        std::set<int>& vert,
        std::set<int>& cust,
        std::set<std::tuple<int, int>>& e,
        std::set<std::tuple<int, int>>& e_subset,
        std::vector<std::vector<GRBVar>>& vars,
        float& eps_for_int,
        float& eps_for_rcc_viol,
        int& dim,
        std::vector<int>& dem,
        int& cap,
        CnstrMgrPointer OldCMP,
        CnstrMgrPointer NewCMP,
        int& total_rc,
        int& current_rc,
        int& rc_implementation,
        std::vector<std::vector<GRBColumn>>& columns,
        double& determine_RC_time,
        bool callback_in_nodes
    ) {
        model = m;
        vertices = vert;
        customers = cust;
        edges = e;
        edges_subset = e_subset;
        x = vars;
        eps_for_integrality = eps_for_int;
        eps_for_rcc_violation = eps_for_rcc_viol;
        dimension = dim;
        demand = dem;
        capacity = cap;
        MyOldCutsCMP = OldCMP;
        MyCutsCMP = NewCMP;

        total_rcc = total_rc;
        current_rcc = current_rc;

        rcc_implementation = rc_implementation;

        var_columns = columns;
        determine_RCC_time = determine_RC_time;

        type2_max = int((dimension + 1) / 2.0f);
        float sqrt_val = float(std::sqrt(dimension * dimension - 10 * dimension + 9));
        type3_min = (-sqrt_val + 3 * dimension - 3) / 4;
        type3_max = (sqrt_val + 3 * dimension - 3) / 4;

        mipnode_callback = callback_in_nodes;
    }

    my_rc_callback() {};

protected:
    // Callback separating only RC inequalities at incumbents and possibly at nodes
    void callback() {
        if ((where == GRB_CB_MIPNODE) && mipnode_callback) {
            if (getIntInfo(GRB_CB_MIPNODE_STATUS) == 2) {
                std::vector<std::tuple<int, int, double>> solution;
                for (auto& edge : edges_subset) {
                    int i = std::get<0>(edge);
                    int j = std::get<1>(edge);
                    double sol_val = getNodeRel(x[i][j]);
                    if (sol_val >= eps_for_integrality) {
                        solution.push_back(std::make_tuple(i, j, sol_val));
                    }
                }
                double max_rcc_viol = add_RC(solution);
            }
        }
        if (where == GRB_CB_MIPSOL) {
            std::vector<std::tuple<int, int, double>> solution;
            for (auto& edge : edges_subset) {
                int i = std::get<0>(edge);
                int j = std::get<1>(edge);
                double sol_val = getSolution(x[i][j]);
                if (sol_val >= eps_for_integrality) {
                    solution.push_back(std::make_tuple(i, j, sol_val));
                }
            }
            double max_rcc_viol = add_RC(solution);
        }

        for (int i = 0; i < MyCutsCMP->Size; i++)
        {
            CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
        }
        MyCutsCMP->Size = 0;
    }

    double add_RC(std::vector<std::tuple<int, int, double>> solution) {
        double max_violation = 0.0;
        std::vector<RC_Inequality> new_cuts;
        separate_rcc_heuristically(new_cuts, max_violation, solution, dimension, demand, capacity,
            eps_for_integrality, eps_for_rcc_violation, MyOldCutsCMP, MyCutsCMP);

        for (const auto& cut : new_cuts) {
            GRBLinExpr lhs;
            char sense;
            int rhs = cut.rhs;
            std::set<int> set = cut.set;
            bool check_violation = cut.check_violation;
            bool add_to_cnstrmgr = cut.add_to_cnstrmgr;
            std::set<std::tuple<int, int, double>> edges_not_subset;

            auto start = std::chrono::high_resolution_clock::now();

            determine_rcc(lhs, sense, rhs, set, x, vertices, customers,
                edges, edges_subset, edges_not_subset,
                type2_max, type3_min, type3_max, rcc_implementation);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
            determine_RCC_time += duration.count();

            //total_rcc++;
            //current_rcc++;
            addLazy(lhs, sense, rhs);
        }
        return max_violation;
    }
};
