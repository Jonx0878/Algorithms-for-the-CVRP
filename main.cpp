#include "comb.h"
#include "cvrp.h"
#include "glm.h"
#include "mstar.h"
#include "pre_tests.h"
#include "col_gen.h"

//General Settings
const double EPS_FOR_INTEGRALITY = 0.001f;
const bool EXACT_IF_TERMINATE = true; // To not apply exact separation during inner loops unless it terminates otherwise 
//- True only makes sense if both exact and heur separation - Has no effect if no column generation
const bool ITERATE_SEP = true; //To iterate separation algorithm - Otherwise it separates 1 round of cuts
const double TIME_LIMIT = 72000; // Time Limit in seconds - Set to 0 for no time limit

// Column Generation settings
const int COL_GEN_TYPE = 2; // 0 - Add all edges and resolve and repeat; 1 - Add all edges; 2 - Add m first edges; 3 - Add m best edges; - Defaults to 0 if eliminating variables
const int NO_COL_GEN_EDGES = -1; // Only has effect if COL_GEN_TYPE is 2 or 3 - Set to -1 one to be |V|
const float EPS_FOR_COL_GEN = 0.0001f;
const float SLACK_EPSILON = 0.01f; // Set to 0 to not remove constraints
const int REMOVE_DURING_SEP = 10; // How many separation iterations to run before removing constraints - 0 to not remove
const float EPS_FOR_ELIM = 0.00f; // To eliminate variables during column generation - Set to 0 to not eliminate variables

//RCC settings
const int RCC_IMPLEMENTATION = 0; //0 - By set_size; 1 - type 2 or 3 by # Vars; 2 - all three types by # vars
const bool HEURISTIC_RCC_SEPARATION = true;
const bool EXACT_RCC_SEPARATION = true;
const float EPS_FOR_EARLY_RCC_TERMINATION = 0.0f; //Set to 0 for always finding most violated cut
const double EPS_FOR_RCC_VIOLATION = 0.001f;
const bool SHRINK_RCC_GRAPH = true; // Only for exact separation

//SC settings
const bool HEURISTIC_SC_SEPARATION = true;
const int MAX_EXACT_SC_TEETH = 0; // Set to 0 for no exact separation
const float EPS_FOR_EARLY_SC_TERMINATION = 0.0f; //Set to 0 for always finding most violated cut
const double EPS_FOR_SC_VIOLATION = 0.001f;
const bool SHRINK_SC_GRAPH = true; // Only for exact separation
const bool SC_ONLY_WHEN_TERMINATE = false;

//MSTAR settings
const bool HEURISTIC_MSTAR_SEPARATION = true;
const bool EXACT_MSTAR_SEPARATION = false;
const float EPS_FOR_EARLY_MSTAR_TERMINATION = 0.0f; //Set to 0 for always finding most violated cut
const double EPS_FOR_MSTAR_VIOLATION = 0.01f;
const bool SHRINK_MSTAR_GRAPH = false; // Only for exact separation
const bool MSTAR_ONLY_WHEN_TERMINATE = true;

//GLM settings
const bool GLM_SEPARATION = true;
const double EPS_FOR_GLM_VIOLATION = 0.01f;
const bool GLM_ONLY_WHEN_TERMINATE = false;


void solve_cvrp(std::string filename, bool stop_after_lp, bool save_solutions) {
    // Setup
    CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;
    CMGR_CreateCMgr(&MyCutsCMP, 100);
    CMGR_CreateCMgr(&MyOldCutsCMP, 100); /* Contains no cuts initially */
    GRBEnv ENV = GRBEnv();

    CVRP instance = CVRP(filename);
    instance.set_starting_edges();

    auto start = std::chrono::high_resolution_clock::now();
    // Create Base Model
    instance.setup_base_model(true, true, ENV);


    // Solve LP Relaxation
    instance.model->optimize();

    instance.solve_instance(
        ENV, MyCutsCMP,
        MyOldCutsCMP,
        start,
        TIME_LIMIT,
        EXACT_IF_TERMINATE,
        ITERATE_SEP,
        EPS_FOR_INTEGRALITY,
        COL_GEN_TYPE,
        NO_COL_GEN_EDGES,
        EPS_FOR_COL_GEN,
        SLACK_EPSILON,
        REMOVE_DURING_SEP,
        EPS_FOR_ELIM,
        RCC_IMPLEMENTATION,
        HEURISTIC_RCC_SEPARATION,
        EXACT_RCC_SEPARATION,
        SHRINK_RCC_GRAPH,
        EPS_FOR_RCC_VIOLATION,
        EPS_FOR_EARLY_RCC_TERMINATION,
        HEURISTIC_SC_SEPARATION,
        MAX_EXACT_SC_TEETH,
        SHRINK_SC_GRAPH,
        EPS_FOR_SC_VIOLATION,
        EPS_FOR_EARLY_SC_TERMINATION,
        SC_ONLY_WHEN_TERMINATE,
        HEURISTIC_MSTAR_SEPARATION,
        EXACT_MSTAR_SEPARATION,
        SHRINK_MSTAR_GRAPH,
        EPS_FOR_MSTAR_VIOLATION,
        EPS_FOR_EARLY_MSTAR_TERMINATION,
        MSTAR_ONLY_WHEN_TERMINATE,
        GLM_SEPARATION,
        EPS_FOR_GLM_VIOLATION,
        GLM_ONLY_WHEN_TERMINATE);

    // Write Lp Model
    if (save_solutions) {
        instance.model->write("Results/LP_models/" + filename + ".mps");
    }

    // Solve instance to optimality
    if (!stop_after_lp) {
        instance.add_all_x_variables();
        instance.make_integral();
        instance.add_rc_callback(EPS_FOR_INTEGRALITY, EPS_FOR_RCC_VIOLATION, MyOldCutsCMP, MyCutsCMP, RCC_IMPLEMENTATION, true);


        auto current_time = std::chrono::high_resolution_clock::now();
        int time_left = TIME_LIMIT - int(std::chrono::duration_cast<std::chrono::duration<double>>(current_time - start).count());

        instance.model->set(GRB_DoubleParam_TimeLimit, time_left);
        instance.model->set(GRB_IntParam_OutputFlag, 1);

        instance.model->optimize();
        if (save_solutions) instance.write_opt_solution();
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
    instance.total_time = duration.count();

    // Prints
    std::cout << std::fixed << std::setprecision(2) << "Solve time: " << duration.count() << " seconds\n";
    std::cout << std::fixed << std::setprecision(2) << "Objective: " << instance.get_objective_value() << std::endl;
    std::cout << std::fixed << std::setprecision(2) << "Final Bound: " << instance.get_objective_value() / instance.upper_bound * 100 << "%\n";
    std::cout << "# RCC: " << instance.total_rcc << std::endl;
    std::cout << "# SC: " << instance.total_combs << std::endl;
    std::cout << "# MSTAR: " << instance.total_mstars << std::endl;
    std::cout << "# GLM: " << instance.total_glm << std::endl;
    std::cout << "# Sep Iter: " << instance.inner_iterations << std::endl;
    std::cout << "# CG Iter: " << instance.col_gen_iters << std::endl;
    std::cout << "# Vars elim: " << instance.variables_eliminated << std::endl;

    instance.print_times();
}


int main() {
    bool stop_after_lp = true;
    bool save_solutions = false;

    solve_cvrp("ORTEC-n242-k12", stop_after_lp, save_solutions);
    
    return 0;
}
