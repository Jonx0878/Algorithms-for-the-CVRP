#pragma once

#include "cvrpsep/cnstrmgr.h"
#include "cvrpsep/combsep.h"
#include "cvrpsep/strcomb.h"
#include "gurobi_c++.h"
#include "utils.h"
#include <map>
#include <set>
#include <vector>
#include <chrono>

struct SC_Inequality {
	std::vector<std::set<int>> teeth;
	int rhs;
	bool check_violation;
	bool add_to_cnstrmgr;
};


// Shrinks the grpah before separation of SC inequalities
void shrink_comb_graph(
	ReachPtr CompsRPtr,
	int& NoOfComponents,
	std::vector<int>& shrunk_demand,
	std::set<int>& shrunk_vertices,
	std::set<int>& shrunk_customers,
	std::vector<std::tuple<int, int, double>>& shrunk_solution,
	const std::vector<std::tuple<int, int, double>>& solution,
	const int& dimension,
	std::vector<int>& demand,
	const int& capacity,
	const int& QMin
) {
	char UseDepotMatch;
	int i, j, DemandSum, MinVehicles, NoOfEdges, MaxNoOfCuts, NoOfCustomers;
	int* CompNr, *Demand, *EdgeTail, *EdgeHead, *SuperDemand;
	double CAP, XVal, XSum;
	double* CompBoundary, *EdgeX;
	double** XMatrix, **SMatrix;
	ReachPtr SupportPtr;

	// Initialize Values
	NoOfCustomers = dimension - 1;
	CAP = static_cast<double>(capacity);
	Demand = &demand[0];
	MaxNoOfCuts = int(std::fmin(NoOfCustomers, 100));

	/* Allocate memory for the three vectors EdgeTail, EdgeHead, and EdgeX */
	EdgeTail = new int[solution.size() + 2];
	EdgeHead = new int[solution.size() + 2];
	EdgeX = new double[solution.size() + 2];
	EdgeHead[0] = -1;
	EdgeTail[0] = -1;
	EdgeX[0] = -1.0;
	int k = 1; // Initialize k to 1
	for (const auto& sol_tuple : solution) {
		int i = std::get<0>(sol_tuple);
		if (i == 0) {
			i = dimension;
		}
		int j = std::get<1>(sol_tuple);
		double value = std::get<2>(sol_tuple);
		EdgeHead[k] = i;
		EdgeTail[k] = j;
		EdgeX[k] = value;
		k++;
	}
	NoOfEdges = k - 1; // Adjust NoOfEdges calculation

	ReachInitMem(&SupportPtr, NoOfCustomers + 1);

	XMatrix = MemGetDM(NoOfCustomers + 2, NoOfCustomers + 2);
	for (i = 1; i <= NoOfCustomers + 1; i++) {
		for (j = 1; j <= NoOfCustomers + 1; j++) {
			XMatrix[i][j] = 0.0;
		}
	}

	for (i = 1; i <= NoOfEdges; i++)
	{
		ReachAddForwArc(SupportPtr, EdgeTail[i], EdgeHead[i]);
		ReachAddForwArc(SupportPtr, EdgeHead[i], EdgeTail[i]);

		XMatrix[EdgeTail[i]][EdgeHead[i]] = EdgeX[i];
		XMatrix[EdgeHead[i]][EdgeTail[i]] = EdgeX[i];
	}

	CompNr = MemGetIV(NoOfCustomers + 2);
	SuperDemand = MemGetIV(NoOfCustomers + 2);
	CompBoundary = MemGetDV(NoOfCustomers + 2);

	SMatrix = MemGetDM(NoOfCustomers + 2, NoOfCustomers + 2);

	DemandSum = 0;
	for (i = 1; i <= NoOfCustomers; i++) DemandSum += Demand[i];

	MinVehicles = 0;
	while ((MinVehicles * CAP) < DemandSum) MinVehicles++;

	XSum = 0.0;
	i = NoOfCustomers + 1;
	for (k = 1; k <= SupportPtr->LP[i].CFN; k++)
	{
		j = SupportPtr->LP[i].FAL[k];
		XVal = XMatrix[i][j];
		XSum += XVal;
	}

	if (XSum <= (2.0 * MinVehicles) + 0.01)
		UseDepotMatch = 1;
	else
		UseDepotMatch = 0;

	STRCOMB_Shrink(SupportPtr,
		NoOfCustomers,
		Demand,
		QMin,
		UseDepotMatch,
		XMatrix,
		CompNr,
		SuperDemand,
		CompBoundary,
		SMatrix,
		CompsRPtr,
		&NoOfComponents);


	//Retrieve Shrunk Graph
	shrunk_demand.resize(NoOfComponents);

	shrunk_vertices.insert(0);
	for (i = 1; i < NoOfComponents; i++) {
		shrunk_vertices.insert(i);
		shrunk_customers.insert(i);
		shrunk_demand[i] = SuperDemand[i];

		for (j = i+1; j <= NoOfComponents; j++) {
			if (SMatrix[i][j] >= 0.001f) {
				if (j == NoOfComponents) {
					shrunk_solution.push_back({ 0, i, SMatrix[i][j] });
				}
				else {
					shrunk_solution.push_back({ i, j, SMatrix[i][j] });
				}
			}
		}
	}
}


// Converts sets obtained from shrunk graph to base graph cuts for SC inequalities
void convert_cuts_comb(
	std::vector<SC_Inequality>& cuts_found,
	ReachPtr graph_mapping,
	const int& shrunk_size,
	const int& instance_size
) {
	std::vector<std::set<int>> teeth;
	int rhs;
	bool check_viol;
	bool from_exact;
	std::set<int> new_tooth;
	std::vector<std::set<int>> new_teeth;

	for (int i = 0; i < cuts_found.size(); i++) {
		teeth = cuts_found[i].teeth;
		rhs = cuts_found[i].rhs;
		check_viol = cuts_found[i].check_violation;
		from_exact = cuts_found[i].add_to_cnstrmgr;

		if (from_exact) {
			new_teeth.resize(teeth.size());
			for (int k = 0; k < teeth.size(); k++) {
				for (int j : teeth[k]) {
					if (j == 0) j = shrunk_size;
					for (int l = 1; l <= graph_mapping->LP[j].CFN; l++)
					{
						int v = graph_mapping->LP[j].FAL[l];
						if (v == instance_size) v = 0;
						new_tooth.insert(v);
					}
				}
				new_teeth[k] = new_tooth;
				new_tooth.clear();
			}
			cuts_found[i] = SC_Inequality( new_teeth, rhs, check_viol, from_exact );
			new_tooth.clear();
		}
	}
}


// Determines the lhs, sense, and rhs of an SC inequality
void determine_strengthened_comb_cut(
	GRBLinExpr& lhs,
	const std::vector<std::set<int>>& teeth,
	GRBModel* model,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<int>& vertices,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset
) {
	std::map<std::tuple<int, int>, double> edge_count_map;

	for (int k = 0; k < teeth.size(); k++) {
		for (const auto& edge : delta_edges(edges, vertices, teeth[k])) {
			edge_count_map[edge]++;
		}
	}
	for (const auto& [edge, coeff] : edge_count_map) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += coeff * x[i][j];
		else edges_not_subset.insert({ i, j, coeff });
	}
}


// Adds an SC inequality to the model
void add_strenghtened_comb_cut(
	int& no_new_cuts,
	const std::vector<SC_Inequality>& cuts_found,
	GRBModel* model,
	const std::vector<std::vector<GRBVar>>& x,
	const int& dimension,
	int& total_combs,
	int& current_combs,
	const std::set<int>& vertices,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::vector<std::vector<GRBColumn>>& var_columns,
	const float& eps_for_viol,
	double& time, CnstrMgrPointer MyCutsCMP
) {
	for (const auto& cut : cuts_found) {
		GRBLinExpr lhs;
		char sense = GRB_GREATER_EQUAL;
		int rhs = cut.rhs;
		std::vector<std::set<int>> teeth = cut.teeth;
		bool check_violation = cut.check_violation;
		bool add_to_cnstrmgr = cut.add_to_cnstrmgr;
		std::set<std::tuple<int, int, double>> edges_not_subset;

		auto start = std::chrono::high_resolution_clock::now();

		determine_strengthened_comb_cut(lhs, teeth, model, x, vertices, edges, edges_subset, edges_not_subset);

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
		time += duration.count();

		if (check_violation) {
			//std::cout << lhs.getValue() << sense << rhs - eps_for_viol << std::endl;
			if (lhs.getValue() <= rhs - eps_for_viol) goto addcut;
		}
		else {
		addcut:
			total_combs++;
			current_combs++;
			std::string name = "SC" + itos(total_combs);
			GRBConstr constr = model->addConstr(lhs, sense, rhs, name);
			for (const auto& var : edges_not_subset) {
				int i = std::get<0>(var);
				int j = std::get<1>(var);
				double coeff = std::get<2>(var);
				var_columns[i][j].addTerm(coeff, constr);
			}
			no_new_cuts++;

			if (add_to_cnstrmgr) {
				int no_of_teeth = int(teeth.size()) - 1;
				// Create Handle
				std::vector<int> handle_list(teeth[0].size() + 1);
				int index = 1;
				for (int i : teeth[0]) {
					handle_list[index] = i;
					index++;
				}

				// Create Teeth
				int teeth_size = no_of_teeth;
				std::vector<int> teeth_list(no_of_teeth*(dimension));
				for (int t = 1; t <= no_of_teeth; t++) {
					teeth_list[t] = teeth_size + 1; // Teeth t starts at index teeth_list[t]
					for (int i : teeth[t]) {
						teeth_size++;
						teeth_list[teeth_size] = i;
					}
				}

				CMGR_AddExtCnstr(
					MyCutsCMP,
					CMGR_CT_STR_COMB,
					no_of_teeth, //No of teeth
					teeth[0].size(), // Handle Size
					&handle_list[0], // Handle
					teeth_size, // Teeth size
					&teeth_list[0], // Teeth
					rhs
				);
			}
		}
	}
}


// Separate SC inequalities using CVRPSEP
void separate_comb_heuristically(
	double& max_violation,
	std::vector<SC_Inequality>& new_cuts,
	const std::vector<std::tuple<int, int, double>>& solution,
	const int& dimension,
	std::vector<int>& demand,
	const int& capacity,
	const int& QMin,
	CnstrMgrPointer MyOldCutsCMP,
	CnstrMgrPointer MyCutsCMP
) 
{
	int NoOfCustomers = dimension - 1;
	int* Demand = &demand[0];
	int NoOfEdges, MaxNoOfCuts;
	int* EdgeTail, * EdgeHead;
	double* EdgeX;

	/* Allocate memory for the three vectors EdgeTail, EdgeHead, and EdgeX */
	EdgeTail = new int[solution.size() + 2];
	EdgeHead = new int[solution.size() + 2];
	EdgeX = new double[solution.size() + 2];
	EdgeHead[0] = -1;
	EdgeTail[0] = -1;
	EdgeX[0] = -1.0;
	int k = 1; // Initialize k to 1
	for (const auto& sol_tuple : solution) {
		int i = std::get<0>(sol_tuple);
		if (i == 0) {
			i = dimension;
		}
		int j = std::get<1>(sol_tuple);
		double value = std::get<2>(sol_tuple);
		EdgeHead[k] = i;
		EdgeTail[k] = j;
		EdgeX[k] = value;
		k++;
	}
	NoOfEdges = k - 1; // Adjust NoOfEdges calculation

	MaxNoOfCuts = int(std::fmin(NoOfCustomers, 100));

	COMBSEP_SeparateCombs(
		NoOfCustomers,
		Demand,
		capacity,
		QMin,
		NoOfEdges,
		EdgeTail,
		EdgeHead,
		EdgeX,
		MaxNoOfCuts,
		&max_violation,
		MyCutsCMP);

	//Retrieve Cuts

	/* Allocate memory for n+2 rows and TMAX+1 columns in InTooth, */
	/* where TMAX is an upper bound on the number of teeth in any */
	/* of the separated strengthened comb inequalities. */
	std::vector<std::set<int>> teeth;
	int NoOfTeeth, MinIdx, MaxIdx, j, RHS;
	int TMAX = NoOfCustomers;
	char** InTooth = new char* [NoOfCustomers + 2];
	for (int i = 0; i < NoOfCustomers + 2; ++i) {
		InTooth[i] = new char[TMAX + 1];
	}

	for (int i = 0; i < MyCutsCMP->Size; i++) {
		if (MyCutsCMP->CPL[i]->CType == CMGR_CT_STR_COMB)
		{
			NoOfTeeth = MyCutsCMP->CPL[i]->Key;
			teeth.resize(NoOfTeeth + 1);

			for (int Node = 1; Node <= NoOfCustomers + 1; Node++)
				for (int Tooth = 0; Tooth <= NoOfTeeth; Tooth++)
					InTooth[Node][Tooth] = 0;
			for (int k = 1; k <= MyCutsCMP->CPL[i]->IntListSize; k++)
			{
				j = MyCutsCMP->CPL[i]->IntList[k];
				InTooth[j][0] = 1; /* Node j is in the handle */
				teeth[0].insert(j);
			}

			for (int t = 1; t <= NoOfTeeth; t++)
			{
				MinIdx = MyCutsCMP->CPL[i]->ExtList[t];
				if (t == NoOfTeeth)
					MaxIdx = MyCutsCMP->CPL[i]->ExtListSize;
				else
					MaxIdx = MyCutsCMP->CPL[i]->ExtList[t + 1] - 1;
				for (k = MinIdx; k <= MaxIdx; k++)
				{
					j = MyCutsCMP->CPL[i]->ExtList[k];
					InTooth[j][t] = 1; /* Node j is in tooth t */
					j = (j > NoOfCustomers) ? 0 : j;
					teeth[t].insert(j);
				}
			}

			/* Now InTooth[j][t] == 1 if and only if */
			/* node j is in tooth number t. The depot */
			/* is node number NoOfCustomers+1, and the */
			/* handle is represented as tooth number 0. */
			RHS = int(MyCutsCMP->CPL[i]->RHS);

			new_cuts.emplace_back(SC_Inequality(teeth, RHS, (max_violation < 0.1f), false));
			teeth.clear();
		}
	}
}


// Sets up the model used for exact separation of SC inequalities
GRBModel* setup_sc_model(
	std::vector<std::vector<GRBVar>>& tau,
	std::vector<std::vector<GRBVar>>& mu,
	const std::vector<std::tuple<int, int, double>>& solution,
	ReachPtr shrunk_mapping,
	const std::set<std::tuple<int, int>>& edges_subset,
	const std::set<int>& customers,
	const std::set<int>& vertices,
	const int& capacity,
	const std::vector<int>& demand, const int& gcd, const int& teeth,
	int& standard_dim,
	GRBEnv env
) {
	// Create Model
	GRBModel* model = new GRBModel(env);
	std::string name;

	// Setup variables
	int no_of_locations = int(vertices.size());

	tau.resize(no_of_locations); // For i in T_k
	std::vector<std::vector<std::vector<GRBVar>>> eta(no_of_locations); // For i in T_k AND T_l
	std::vector<std::vector<GRBVar>> iota(no_of_locations); // For i in T_k \ H
	mu.resize(3); // No of vehicles for teeth k
	GRBVar alpha; // For odd RHS
	std::vector<std::vector<GRBVar>> beta(teeth + 1); // 1 if T_i AND T_j subset H must hold
	std::vector<std::vector<std::vector<GRBVar>>> rho(no_of_locations); // Linearisation of eta_i^0k * eta_i^0l


	// Resize each inner vector
	for (int j = 0; j <= 2; j++) {
		mu[j].resize(teeth + 1);
	}

	for (int i : vertices) {
		tau[i].resize(teeth + 1);
		eta[i].resize(teeth + 1);
		iota[i].resize(teeth + 1);
		rho[i].resize(teeth + 1);

		for (int k = 0; k < teeth; k++) {
			eta[i][k].resize(teeth + 1);
			rho[i][k].resize(teeth + 1);
		}
	}

	for (int k = 1; k < teeth; k++) {
		beta[k].resize(teeth + 1);
	}

	// Change parameters
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);
	model->set(GRB_DoubleParam_Cutoff, 0.0);
	//model->set(GRB_IntParam_NonConvex, 1);

	// Add Variables

	double max_vehicles = 0;
	for (int i : customers) {
		max_vehicles += demand[i];
	}
	alpha = model->addVar(ceil((3 * teeth - 1) / 2.0f), floor((3*teeth*max_vehicles - 1)/2.0f), NULL, GRB_INTEGER, "alpha");

	for (int i : customers) {
		tau[i][0] = model->addVar(0, 1, NULL, GRB_BINARY, "tau[" + std::to_string(i) + ", " + std::to_string(0) + "]");

		for (int l = 1; l <= teeth; l++) {
			eta[i][0][l] = model->addVar(0, 1, NULL, GRB_CONTINUOUS, "eta[" + std::to_string(i) + ", " + std::to_string(0) + ", " + std::to_string(l) + "]");
		}

		iota[i][0] = model->addVar(0, 1, NULL, GRB_CONTINUOUS, "iota[" + std::to_string(i) + ", " + std::to_string(0) + "]");
	}

	for (int k = 1; k <= teeth; k++) {
		for (int i : vertices) {
			tau[i][k] = model->addVar(0, 1, NULL, GRB_BINARY, "tau[" + std::to_string(i) + ", " + std::to_string(k) + "]");

			for (int l = k + 1; l <= teeth; l++) {
				eta[i][k][l] = model->addVar(0, 1, NULL, GRB_CONTINUOUS, "eta[" + std::to_string(i) + ", " + std::to_string(k) + ", " + std::to_string(l) + "]");
				rho[i][k][l] = model->addVar(0, 1, NULL, GRB_CONTINUOUS, "rho[" + std::to_string(i) + ", " + std::to_string(k) + ", " + std::to_string(l) + "]");
			}

			iota[i][k] = model->addVar(0, 1, NULL, GRB_CONTINUOUS, "iota[" + std::to_string(i) + ", " + std::to_string(k) + "]");
		}

		for (int j = 0; j <= 2; j++) {
			mu[j][k] = model->addVar(1, max_vehicles, NULL, GRB_INTEGER, "mu[" + std::to_string(j) + ", " + std::to_string(k) + "]");
		}

		for (int l = k + 1; l <= teeth; l++) {
			beta[k][l] = model->addVar(0, 1, NULL, GRB_CONTINUOUS, "beta[" + std::to_string(k) + ", " + std::to_string(l) + "]");
		}
	}

	mu[0][1].set(GRB_CharAttr_VType, GRB_CONTINUOUS);

	model->update();

	// Add Constraints
	GRBLinExpr Lexpr;
	GRBQuadExpr Qexpr;

	for (int k = 1; k <= teeth; k++) {
		for (int j = 0; j <= 2; j++) {
			Lexpr += mu[j][k];
		}
	}
	model->addConstr(Lexpr, GRB_EQUAL, 2 * alpha + 1); // Odd RHS
	Lexpr.clear();

	// Vehicle Constraints
	for (int k = 1; k <= teeth; k++) {
		for (int i : customers) {
			Lexpr += demand[i] * eta[i][0][k];
		}
		model->addConstr(capacity * (mu[0][k] - 1) + gcd, GRB_LESS_EQUAL, Lexpr);
		model->addConstr(capacity * (mu[0][k]), GRB_GREATER_EQUAL, Lexpr);
		Lexpr.clear();


		for (int i : customers) {
			Qexpr += demand[i] * (iota[i][k] + tau[0][k] - 2 * iota[i][k] * tau[0][k]);
		}
		model->addQConstr(capacity * (mu[1][k] - 1) + gcd, GRB_LESS_EQUAL, Qexpr);
		model->addQConstr(capacity * (mu[1][k]), GRB_GREATER_EQUAL, Qexpr);
		Qexpr.clear();

		for (int i : customers) {
			Qexpr += demand[i] * (tau[i][k] + tau[0][k] - 2 * tau[i][k] * tau[0][k]);
		}
		model->addQConstr(capacity * (mu[2][k] - 1) + gcd, GRB_LESS_EQUAL, Qexpr);
		model->addQConstr(capacity * (mu[2][k]), GRB_GREATER_EQUAL, Qexpr);
		Qexpr.clear();
	}

	// Handle Bounds
	for (int i : customers) {
		Lexpr += shrunk_mapping->LP[i].CFN * tau[i][0];
	}
	model->addConstr(Lexpr, GRB_GREATER_EQUAL, 1);
	model->addConstr(Lexpr, GRB_LESS_EQUAL, standard_dim - 2);
	Lexpr.clear();

	for (int k = 1; k <= teeth; k++) {
		// Teeth Bounds
		shrunk_mapping->LP[vertices.size()].CFN * tau[0][k];
		for (int i : customers) {
			Lexpr += shrunk_mapping->LP[i].CFN * tau[i][k];
		}
		model->addConstr(Lexpr, GRB_GREATER_EQUAL, 2);
		model->addConstr(Lexpr, GRB_LESS_EQUAL, standard_dim - 1);
		Lexpr.clear();

		//Handle-Tooth intersect
		for (int i : customers) {
			Lexpr += eta[i][0][k];
		}
		model->addConstr(Lexpr, GRB_GREATER_EQUAL, 1);
		Lexpr.clear();

		//Tooth \ handle
		Lexpr += tau[0][k];
		for (int i : customers) {
			Lexpr += iota[i][k];
		}
		model->addConstr(Lexpr, GRB_GREATER_EQUAL, 1);
		Lexpr.clear();

		// Internal Tooth OR
		for (int l = k + 1; l <= teeth; l++) {
			for (int i : customers) {
				model->addConstr(beta[k][l] >= rho[i][k][l]);
				model->addConstr(tau[i][0] - eta[i][k][l] >= beta[k][l] - 1);
			}
			model->addConstr(eta[0][k][l] <= 1 - beta[k][l]);
		}
	}

	// Intersect variables - Eta
	for (int i : customers) {
		for (int l = 1; l <= teeth; l++) {
			model->addConstr(eta[i][0][l] >= tau[i][0] + tau[i][l] - 1);
			model->addConstr(eta[i][0][l] <= tau[i][0]);
			model->addConstr(eta[i][0][l] <= tau[i][l]);
		}
	}
	for (int i : vertices) {
		for (int k = 1; k < teeth; k++) {
			for (int l = k + 1; l <= teeth; l++) {
				model->addConstr(eta[i][k][l] >= tau[i][k] + tau[i][l] - 1);
				model->addConstr(eta[i][k][l] <= tau[i][k]);
				model->addConstr(eta[i][k][l] <= tau[i][l]);
			}
		}
	}

	// Rho
	for (int i : customers) {
		for (int k = 1; k < teeth; k++) {
			for (int l = k + 1; l <= teeth; l++) {
				model->addConstr(rho[i][k][l] >= eta[i][0][k] + eta[i][0][l] - 1);
				model->addConstr(rho[i][k][l] <= eta[i][0][k]);
				model->addConstr(rho[i][k][l] <= eta[i][0][l]);
			}
		}
	}

	// Subtract variables - Iota
	for (int i : customers) {
		for (int k = 1; k <= teeth; k++) {
			model->addConstr(iota[i][k] >= tau[i][k] - tau[i][0]);
			model->addConstr(iota[i][k] <= tau[i][k]);
			model->addConstr(iota[i][k] <= 1 - tau[i][0]);
		}
	}


	// Set Objective
	GRBQuadExpr obj = 2*alpha + 2;

	for (int k = 0; k <= teeth; k++) {
		for (const auto& edge : solution) {
			int i = std::get<0>(edge);
			int j = std::get<1>(edge);
			double value = std::get<2>(edge);
			if (k == 0 && i == 0) {
				obj -= value * tau[j][k];
			}
			else {
				obj -= value * (tau[i][k] + tau[j][k] - 2 * tau[i][k] * tau[j][k]);
			}
		}
	}

	model->setObjective(obj, GRB_MAXIMIZE);

	return model;
}


// Separates SC inequalities exactly
void separate_sc_exactly(
	double& max_violation,
	std::vector<SC_Inequality>& new_cuts,
	ReachPtr shrunk_mapping,
	const std::vector<std::tuple<int, int, double>>& solution,
	const std::set<std::tuple<int, int>>& edges_subset,
	const std::set<int>& vertices,
	const std::set<int>& customers,
	const int& capacity,
	const std::vector<int>& demand, int& standard_dim,
	const float& eps_for_early_termination,
	const int& gcd,
	const int& no_teeth, const float& eps_for_int, const float& eps_for_violation, double& time, GRBEnv env
) {
	std::vector<std::vector<GRBVar>> tau;
	std::vector<std::vector<GRBVar>> mu;

	for (int t = 2; t <= no_teeth; t++) {
		std::vector<std::set<int>> teeth;

		// Setup separation model
		GRBModel* model = setup_sc_model(tau, mu, solution, shrunk_mapping, edges_subset, customers, vertices, capacity, demand, gcd, t, standard_dim, env);
		if (eps_for_early_termination > 0) model->getEnv().set(GRB_DoubleParam_BestObjStop, eps_for_early_termination);

		auto start = std::chrono::high_resolution_clock::now();

		model->optimize();

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
		time += duration.count();

		if (model->get(GRB_DoubleAttr_ObjVal) > max_violation) max_violation = model->get(GRB_DoubleAttr_ObjVal);

		for (int sol = 0; sol < model->get(GRB_IntAttr_SolCount); ++sol) {
			model->set(GRB_IntParam_SolutionNumber, sol);
			if (model->get(GRB_DoubleAttr_PoolObjVal) >= eps_for_violation) {
				teeth.resize(t + 1);
				for (int i : customers) {
					if (tau[i][0].get(GRB_DoubleAttr_Xn) >= 1 - eps_for_int) {
						teeth[0].insert(i);
					}
				}
				for (int k = 1; k <= t; k++) {
					for (int i : vertices) {
						if (tau[i][k].get(GRB_DoubleAttr_Xn) >= 1 - eps_for_int) {
							teeth[k].insert(i);
						}
					}
				}

				double rhs = 1;
				for (int k = 1; k <= t; k++) {
					for (int j = 0; j <= 2; j++) {
						rhs += mu[j][k].get(GRB_DoubleAttr_Xn);
					}
			
				}
				new_cuts.emplace_back(SC_Inequality(teeth, int(rhs), model->get(GRB_DoubleAttr_PoolObjVal) < 0.1, true));
				teeth.clear();
			}
		}
		if (max_violation >= 0.1) break;
	}
}
