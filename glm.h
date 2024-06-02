#pragma once
#include "cvrpsep/glmsep.h"

struct GLM_Inequality {
	std::set<int> nucleus;
	bool check_violation = true;
};


// Determines lhs, sense, and rhs of an GLM inequality of type 1
void determine_glm_type1(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const GLM_Inequality& glm,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset,
	const std::set<int>& customers,
	const std::vector<int>& demand,
	const int& capacity
) {
	for (const auto& edge : edges_in(edges, glm.nucleus)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += capacity * x[i][j];
		else edges_not_subset.insert({ i, j, capacity });
	}
	
	for (const int& k : customers) {
		if (glm.nucleus.count(k) == 0) {
			std::set<int> cust_set;
			cust_set.insert(k);
			for (const auto& edge : edges_from_to(edges, glm.nucleus, cust_set)) {
				int i = std::get<0>(edge);
				int j = std::get<1>(edge);
				if (edges_subset.count({ i, j }) == 1) lhs += demand[k] * x[i][j];
				else edges_not_subset.insert({ i, j, demand[k]});
			}
		}
	}
	sense = GRB_LESS_EQUAL;
	rhs = capacity * glm.nucleus.size();
	for (const int& i: glm.nucleus) {
		rhs -= demand[i];
	}
}


// Determines lhs, sense, and rhs of an GLM inequality of type 2
void determine_glm_type2(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const GLM_Inequality& glm,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset,
	const std::set<int>& vertices,
	const std::set<int>& customers,
	const std::vector<int>& demand,
	const int& capacity
) {
	std::map<std::tuple<int, int>, double> edge_count_map;

	for (const auto& edge : delta_edges(edges, vertices, glm.nucleus)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += capacity * x[i][j];
		else edge_count_map[edge] += capacity;
	}

	for (const int& k : customers) {
		if (glm.nucleus.count(k) == 0) {
			std::set<int> cust_set;
			cust_set.insert(k);
			for (const auto& edge : edges_from_to(edges, glm.nucleus, cust_set)) {
				int i = std::get<0>(edge);
				int j = std::get<1>(edge);
				if (edges_subset.count({ i, j }) == 1) lhs -= 2* demand[k] * x[i][j];
				else edge_count_map[{i, j}] -= 2 * demand[k];
			}
		}
	}

	for (const auto& [edge, coeff] : edge_count_map) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		edges_not_subset.insert({ i, j, coeff });
	}
	sense = GRB_GREATER_EQUAL;
	rhs = 0;
	for (const int& i : glm.nucleus) {
		rhs += 2 * demand[i];
	}
}


// Determines lhs, sense, and rhs of an GLM inequality
void determine_glm(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const GLM_Inequality& glm,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<int>& vertices,
	const std::set<int>& customers,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset,
	const std::vector<int>& demand,
	const int& capacity,
	const float& type1_max
) {
	if (glm.nucleus.size() <= type1_max) {
		determine_glm_type1(lhs, sense, rhs, glm, x, edges, edges_subset, edges_not_subset, customers, demand, capacity);
	}
	else determine_glm_type2(lhs, sense, rhs, glm, x, edges, edges_subset, edges_not_subset, vertices, customers, demand, capacity);
}


// Adds GLM inequalities to the model
void add_glm(
	int& no_new_cuts,
	const std::vector<GLM_Inequality>& cuts,
	GRBModel* model,
	const std::vector<std::vector<GRBVar>>& x,
	const int& dimension,
	int& total_glm,
	int& current_glm,
	const std::set<int>& vertices,
	const std::set<int>& customers,
	const std::vector<int>& demand,
	const int& capacity,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::vector<std::vector<GRBColumn>>& var_columns,
	const float& type1_max,
	const float& eps_for_viol,
	double& time
) {
	GRBLinExpr lhs;
	char sense = GRB_LESS_EQUAL;
	int rhs = 0;
	for (const auto& glm : cuts) {
		lhs.clear();
		rhs = 0;
		bool check_violation = glm.check_violation;
		std::set<std::tuple<int, int, double>> edges_not_subset;

		auto start = std::chrono::high_resolution_clock::now();

		determine_glm(lhs, sense, rhs, glm, x, vertices, customers,
			edges, edges_subset, edges_not_subset, demand, capacity, type1_max);
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
		time += duration.count();

		if (check_violation) {
			if ((sense == *("<") && lhs.getValue() >= rhs + eps_for_viol) || (sense == *(">") && lhs.getValue() <= rhs - eps_for_viol))
				goto addcut;
		}
		else {
		addcut:
			total_glm++;
			current_glm++;
			no_new_cuts++;
			std::string name = "GLM" + itos(total_glm);
			GRBConstr constr = model->addConstr(lhs, sense, rhs, name);
			for (const auto& var : edges_not_subset) {
				int i = std::get<0>(var);
				int j = std::get<1>(var);
				double coeff = std::get<2>(var);
				var_columns[i][j].addTerm(coeff, constr);
			}
		}
	}
}


// Separates GLM Inequalities using CVRPSEP
void separate_glm_CVRPSEP(
	double& max_violation,
	std::vector<GLM_Inequality>& new_cuts,
	const std::vector<std::tuple<int, int, double>>& solution,
	const int& dimension,
	std::vector<int>& demand,
	const int& capacity
	)
{
	int NoOfCustomers = dimension - 1;
	int NoOfEdges, MaxNoOfCuts, NListSize;
	int* EdgeTail, * EdgeHead, * NList;
	double* EdgeX;
	NList = new int[dimension];

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

	GLMSEP_SeparateGLM(
		NoOfCustomers,
		&demand[0],
		capacity,
		NoOfEdges,
		EdgeTail,
		EdgeHead,
		EdgeX,
		NList,
		&NListSize,
		&max_violation
	);

	//Retrieve Cut
	if (NListSize > 0) {
		std::set<int> nucleus;
		for (int i = 1; i <= NListSize; i++) {
			nucleus.insert(NList[i]);
		}
		new_cuts.emplace_back(GLM_Inequality(nucleus, true));
	}
}
