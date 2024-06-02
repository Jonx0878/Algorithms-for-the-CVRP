#pragma once
#include <iostream>

#include "cvrpsep/cnstrmgr.h"
#include "cvrpsep/capsep.h"
#include "cvrpsep/compress.h"

#include "utils.h"

struct RC_Inequality {
	std::set<int> set;
	int rhs;
	bool check_violation;
	bool add_to_cnstrmgr;
};


void shrink_rcc_graph(
	ReachPtr SuperNodesRPtr,
	std::vector<int>& shrunk_demand,
	std::set<int>& shrunk_vertices,
	std::set<int>& shrunk_customers,
	std::vector<std::tuple<int, int, double>>& shrunk_solution,
	const std::vector<std::tuple<int, int, double>>& solution,
	const int& dimension,
	std::vector<int>& demand,
	const int& capacity,
	double& time,
	CnstrMgrPointer MyOldCutsCMP
) {
	auto start = std::chrono::high_resolution_clock::now();

	int i, j, NoOfEdges, MaxNoOfCuts, NoOfCustomers, NoOfV1Cuts, ShrunkGraphCustNodes;
	int* Demand, * EdgeTail, * EdgeHead, *SuperNodeSize, *NodeList;
	double CAP;
	double* EdgeX, * SuperDemand, *XInSuperNode;
	double** XMatrix, ** SMatrix;
	ReachPtr V1CutsPtr, SupportPtr, SAdjRPtr;

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
	ReachInitMem(&SAdjRPtr, NoOfCustomers + 1);

	SuperDemand = MemGetDV(NoOfCustomers + 1);
	SuperNodeSize = MemGetIV(NoOfCustomers + 1);
	NodeList = MemGetIV(NoOfCustomers + 1);
	XInSuperNode = MemGetDV(NoOfCustomers + 1);

	SMatrix = MemGetDM(NoOfCustomers + 2, NoOfCustomers + 2);
	XMatrix = MemGetDM(NoOfCustomers + 2, NoOfCustomers + 2);
	for (i = 1; i <= NoOfCustomers + 1; i++)
		for (j = 1; j <= NoOfCustomers + 1; j++)
			XMatrix[i][j] = 0.0;


	for (i = 1; i <= NoOfEdges; i++)
	{
		ReachAddForwArc(SupportPtr, EdgeTail[i], EdgeHead[i]);
		ReachAddForwArc(SupportPtr, EdgeHead[i], EdgeTail[i]);

		XMatrix[EdgeTail[i]][EdgeHead[i]] = EdgeX[i];
		XMatrix[EdgeHead[i]][EdgeTail[i]] = EdgeX[i];
	}

	V1CutsPtr = NULL;
	CAPSEP_GetOneVehicleCapCuts(MyOldCutsCMP,
		&V1CutsPtr,
		&NoOfV1Cuts);

	COMPRESS_ShrinkGraph(SupportPtr,
		NoOfCustomers,
		XMatrix,
		SMatrix,
		NoOfV1Cuts,
		V1CutsPtr,
		SAdjRPtr,
		SuperNodesRPtr,
		&ShrunkGraphCustNodes);

	ReachFreeMem(&V1CutsPtr);


	//Retrieve Shrunk Graph
	shrunk_demand.resize(ShrunkGraphCustNodes + 1);

	shrunk_vertices.insert(0);
	for (i = 1; i <= ShrunkGraphCustNodes; i++)
	{
		shrunk_vertices.insert(i);
		shrunk_customers.insert(i);
		shrunk_demand[i] = 0;
		for (j = 1; j <= SuperNodesRPtr->LP[i].CFN; j++)
		{
			k = SuperNodesRPtr->LP[i].FAL[j];
			shrunk_demand[i] += Demand[k];
		}
	}

	for (int i = 1; i <= ShrunkGraphCustNodes; i++) {
		for (j = 1; j <= SAdjRPtr->LP[i].CFN; j++)
		{
			k = SAdjRPtr->LP[i].FAL[j];
			if (k > i) {
				if (k == ShrunkGraphCustNodes + 1) {
					shrunk_solution.push_back({ 0, i, SMatrix[i][k] });
				}
				else {
					shrunk_solution.push_back({ i, k, SMatrix[i][k] });
				}
			}
		}
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
	time += duration.count();
}


void convert_cuts_rcc(std::vector<RC_Inequality>& cuts_found, ReachPtr graph_mapping, double& time) {
	auto start = std::chrono::high_resolution_clock::now();

	std::set<int> set;
	std::set<int> new_set;
	int rhs;
	bool check_viol;
	bool from_exact;

	for (int i = 0; i < cuts_found.size(); i++) {
		set = cuts_found[i].set;
		rhs = cuts_found[i].rhs;
		check_viol = cuts_found[i].check_violation;
		from_exact = cuts_found[i].add_to_cnstrmgr;

		if (from_exact) {
			for (int j : set) {
				for (int l = 1; l <= graph_mapping->LP[j].CFN; l++)
				{
					int k = graph_mapping->LP[j].FAL[l];
					new_set.insert(k);
				}
			}
			cuts_found[i] = RC_Inequality( new_set, rhs + new_set.size() - set.size(), check_viol, from_exact );
			new_set.clear();
		}
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
	time += duration.count();
}


void determine_rcc_type1(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const std::set<int>& cut_set,
	const int& cut_size,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<int>& vertices,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset
) {
	for (const auto& edge : delta_edges(edges, vertices, cut_set)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += x[i][j];
		else edges_not_subset.insert({ i, j, 1.0 });
	}
	sense = GRB_GREATER_EQUAL;
	rhs = -2 * (rhs - cut_size);
}


void determine_rcc_type2(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const std::set<int>& cut_set,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset
) {
	for (const auto& edge : edges_in(edges, cut_set)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += x[i][j];
		else edges_not_subset.insert({ i, j, 1.0 });
	}
	sense = GRB_LESS_EQUAL;
}


void determine_rcc_type3(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const std::set<int>& cut_set,
	const int cut_size,
	const std::set<int> set_bar,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<int>& customers,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset
) {
	for (const auto& edge : edges_in(edges, set_bar)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += 2*x[i][j];
		else edges_not_subset.insert({ i, j, 2.0 });
	}
	for (const auto& j : set_bar) {
		if (edges_subset.count({ 0, j }) == 1) lhs += x[0][j];
		else edges_not_subset.insert({ 0, j, 1.0 });
	}
	for (const auto& j : cut_set) {
		if (edges_subset.count({ 0, j }) == 1) lhs -= x[0][j];
		else edges_not_subset.insert({ 0, j, -1.0 });
	}
	sense = GRB_LESS_EQUAL;
	rhs += -cut_size + set_bar.size();
	rhs *= 2;
}


void determine_rcc(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const std::set<int>& cut_set,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<int>& vertices,
	const std::set<int>& customers,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset,
	const float &type2_max,
	const float &type3_min,
	const float &type3_max,
	const int& cut_implementation
) {
	int cut_size = static_cast<int>(cut_set.size());
	int cut_type;
	std::set<int> set_bar;

	if (cut_implementation == 1) {
		// Implementation by fewest variables - Only type 2 and 3
		std::set_difference(customers.begin(), customers.end(), cut_set.begin(), cut_set.end(), std::inserter(set_bar, set_bar.begin()));
		size_t no_vars_type2 = edges_in(edges_subset, cut_set).size();
		size_t no_vars_type3 = edges_in(edges_subset, set_bar).size() + edges_from_to(edges_subset, std::set<int>({ 0 }), customers).size();
		if (no_vars_type2 <= no_vars_type3) cut_type = 2;
		else cut_type = 3;
		}
	else if (cut_implementation == 2) {
		// Implementation by fewest variables
		std::set_difference(customers.begin(), customers.end(), cut_set.begin(), cut_set.end(), std::inserter(set_bar, set_bar.begin()));
		size_t no_vars_type1 = delta_edges(edges_subset, vertices, cut_set).size();
		size_t no_vars_type2 = edges_in(edges_subset, cut_set).size();
		size_t no_vars_type3 = edges_in(edges_subset, set_bar).size() + edges_from_to(edges_subset, std::set<int>({ 0 }), customers).size();
		if (no_vars_type2 <= no_vars_type3 && no_vars_type2 <= no_vars_type1) cut_type = 2;
		else if (no_vars_type3 <= no_vars_type1) cut_type = 3;
		else cut_type = 1;
		}
	else {
		// Implementation by cut size
		if (cut_size <= type2_max) cut_type = 2;
		else if (type3_min <= cut_size && cut_size <= type3_max) {
			std::set_difference(customers.begin(), customers.end(), cut_set.begin(), cut_set.end(), std::inserter(set_bar, set_bar.begin()));
			cut_type = 3;
		}
		else cut_type = 1;
	}


	if (cut_type == 2) determine_rcc_type2(lhs, sense, rhs, cut_set, x, edges, edges_subset, edges_not_subset);
	else if (cut_type == 3) determine_rcc_type3(lhs, sense, rhs, cut_set, cut_size, set_bar, x, customers, edges, edges_subset, edges_not_subset);
	else determine_rcc_type1(lhs, sense, rhs, cut_set, cut_size, x, vertices, edges, edges_subset, edges_not_subset);
}


void add_rc_cuts(
	const std::vector<RC_Inequality>& cuts,
	GRBModel* model,
	const std::vector<std::vector<GRBVar>>& x,
	const int& dimension,
	int& total_rcc,
	int& current_rcc,
	const std::set<int>& vertices,
	const std::set<int>& customers,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::vector<std::vector<GRBColumn>>& var_columns,
	const float& type2_max,
	const float& type3_min,
	const float& type3_max,
	const int& implementation_type,
	const float& eps_for_viol,
	double& time,
	CnstrMgrPointer MyCutsCMP
) {
	for (const auto& cut : cuts) {
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
			type2_max, type3_min, type3_max, implementation_type);

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
		time += duration.count();

		if (check_violation) {
			if ((sense == *("<") && lhs.getValue() >= rhs + eps_for_viol) || (sense == *(">") && lhs.getValue() <= rhs - eps_for_viol))
				goto addcut;
		}
		else {
		addcut:
			total_rcc++;
			current_rcc++;
			std::string name = "RCC" + itos(total_rcc);
			GRBConstr constr = model->addConstr(lhs, sense, rhs, name);
			for (const auto& var : edges_not_subset) {
				int i = std::get<0>(var);
				int j = std::get<1>(var);
				double coeff = std::get<2>(var);
				var_columns[i][j].addTerm(coeff, constr);
			}

			if (add_to_cnstrmgr) {
				std::vector<int> set_list(set.size() + 1);
				int index = 1;
				for (int i : set) {
					set_list[index] = i;
					index++;
				}
				CMGR_AddCnstr(MyCutsCMP, CMGR_CT_CAP, 0, int(set.size()), &set_list[0], rhs);
			}
		}
	}
}


void separate_rcc_heuristically(
	std::vector<RC_Inequality>& new_cuts,
	double& max_violation,
	const std::vector<std::tuple<int, int, double>>& solution,
	const int& dimension,
	std::vector<int>& demand,
	const int& capacity,
	const double& EpsForIntegrality,
	const double& EpsForViolation,
	CnstrMgrPointer MyOldCutsCMP,
	CnstrMgrPointer MyCutsCMP
) {
	char IntegerAndFeasible;
	int NoOfCustomers = dimension - 1;
	std::vector<double> demand_double(demand.begin(), demand.end());
	const double* Demand = &demand_double[0];
	double CAP = static_cast<double>(capacity);
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

	CAPSEP_SeparateCapCuts(NoOfCustomers,
		Demand,
		CAP,
		NoOfEdges,
		EdgeTail,
		EdgeHead,
		EdgeX,
		MyOldCutsCMP,
		MaxNoOfCuts,
		EpsForIntegrality,
		EpsForViolation,
		&IntegerAndFeasible,
		&max_violation,
		MyCutsCMP);

	// Deallocate memory
	delete[] EdgeTail;
	delete[] EdgeHead;
	delete[] EdgeX;

	int ListSize;
	double RHS;
	int* List;

	List = new int[dimension];
	for (int i = 0; i < MyCutsCMP->Size; i++)
	{
		if (MyCutsCMP->CPL[i]->CType == CMGR_CT_CAP) {
			std::set<int> S;
			ListSize = 0;
			for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; j++)
			{
				List[++ListSize] = MyCutsCMP->CPL[i]->IntList[j];
			}
			/* Now List contains the customer numbers defining the cut. */
			/* The right-hand side of the cut, */
			/* in the form x(S:S) <= |S| - k(S), is RHS. */
			RHS = MyCutsCMP->CPL[i]->RHS;
			for (int j = 1; j <= ListSize; j++) {
				S.insert(List[j]);
			}
			new_cuts.emplace_back(RC_Inequality(S, RHS, false, false));
		}
	}
}


GRBModel* setup_rcc_model(
	GRBVar& lambda,
	std::vector<GRBVar>& xi,
	const std::vector<std::tuple<int, int, double>>& solution,
	const std::set<int>& customers,
	const int& capacity,
	const std::vector<int>& demand,
	const int& gcd,
	GRBEnv env,
	ReachPtr supervertices
) {
	// Create Model
	GRBModel* model = new GRBModel(env);

	xi.resize(customers.size() + 1);
	std::vector<std::vector<GRBVar>> gamma(customers.size() + 2);
	for (int i = 0; i <= customers.size(); i++) {
		// Resize each inner vector
		gamma[i].resize(customers.size() + 1);
	}

	// Change parameters
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);
	//model->set(GRB_DoubleParam_Cutoff, 0.0);
	//model->set(GRB_IntParam_PoolSearchMode, 2);
	//model->set(GRB_IntParam_PoolSolutions, 20);

	// Add Variables
	lambda = model->addVar(0, GRB_INFINITY, NULL, GRB_INTEGER, "lambda");
	for (int i : customers) {
		xi[i] = model->addVar(0, 1, NULL, GRB_BINARY, "xi[" + std::to_string(i) + "]");
	}

	for (auto& edge : solution) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (i > 0) {
			gamma[i][j] = model->addVar(0, 1, NULL, GRB_CONTINUOUS, "gamma[" + std::to_string(i) + "," + std::to_string(j) + "]");
		}
	}
	model->update();

	// Add Constraints
	GRBLinExpr expr;
	for (int i : customers) {
		expr += demand[i] * xi[i];
	}
	model->addConstr(capacity * (lambda - 1) + gcd, GRB_LESS_EQUAL, expr);

	expr.clear();

	if (supervertices == NULL) {
		for (int i : customers) {
			expr += xi[i];
		}
	}
	else {
		for (int i : customers) {
			expr += supervertices->LP[i].CFN * xi[i];
		}
	}
	model->addConstr(expr, GRB_GREATER_EQUAL, 2);
	expr.clear();

	for (auto& edge : solution) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (i > 0) {
			model->addConstr(gamma[i][j], GRB_LESS_EQUAL, xi[i]);
			model->addConstr(gamma[i][j], GRB_LESS_EQUAL, xi[j]);
		}
	}
	// Set Objective
	GRBLinExpr obj = 2.0 * lambda;
	for (const auto& edge : solution) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		double value = std::get<2>(edge);
		if (i == 0) {
			obj -= value * xi[j];
		}
		else {
			obj -= value * (xi[i] + xi[j] - 2.0 * gamma[i][j]);
		}
	}
	model->setObjective(obj, GRB_MAXIMIZE);

	return model;
}


void separate_rcc_exactly(
	std::vector<RC_Inequality>& new_cuts,
	double& max_violation,
	const std::vector<std::tuple<int, int, double>>& solution,
	const std::set<int>& customers,
	const int& capacity,
	const std::vector<int>& demand,
	const float& eps_for_early_termination,
	const int& gcd,
	const float& eps_for_int,
	const float& eps_for_violation,
	double &time,
	GRBEnv env,
	ReachPtr shrunk_mapping
) {
	GRBVar lambda;
	std::vector<GRBVar> xi;

	// Setup separation model
	GRBModel* model = setup_rcc_model(lambda, xi, solution, customers, capacity, demand, gcd, env, shrunk_mapping);
	if (eps_for_early_termination > 0) model->getEnv().set(GRB_DoubleParam_BestObjStop, eps_for_early_termination);

	auto start = std::chrono::high_resolution_clock::now();

	model->optimize();

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
	time += duration.count();

	max_violation = std::fmax(max_violation, model->get(GRB_DoubleAttr_ObjVal));

	for (int sol = 0; sol < model->get(GRB_IntAttr_SolCount); ++sol) {
		model->set(GRB_IntParam_SolutionNumber, sol);
		if (model->get(GRB_DoubleAttr_PoolObjVal) >= eps_for_violation) {
			std::set<int> S;
			for (int i : customers) {
				if (xi[i].get(GRB_DoubleAttr_Xn) >= 1 - eps_for_int) {
					S.insert(i);
				}
			}
			new_cuts.emplace_back(RC_Inequality(S, S.size() - lambda.get(GRB_DoubleAttr_Xn), model->get(GRB_DoubleAttr_PoolObjVal) < 0.1, true));
		}
	}
	delete model;
}