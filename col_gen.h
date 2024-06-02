#pragma once
#include <chrono>
#include <set>
#include <vector>
#include <unordered_set>


#include "cvrp.h"
#include "gurobi_c++.h"
#include "utils.h"
struct RCpair {
	std::tuple<int, int> edge;
	int RC;
};


// Calculates the reduced cost of a variable not in the model
double calculate_reduced_cost(const double& cost, GRBColumn& column) {
	double reduced_cost = cost;

	for (int k = 0; k < column.size(); k++) {
		GRBConstr constr = column.getConstr(k);
		double coeff = column.getCoeff(k);
		double pi = constr.get(GRB_DoubleAttr_Pi);
		reduced_cost -= coeff * pi;
	}
	return reduced_cost;
}


// Determines which variables to add to the model
std::set<std::tuple<int, int>> check_reduced_costs(
	GRBModel* model,
	const int type,
	const int no_edges,
	const float epsilon,
	std::set<std::tuple<int, int>>& edges_not_subset,
	std::vector<std::vector<int>>& cost,
	std::vector<std::vector<GRBColumn>>& var_columns
) {

	std::set<std::tuple<int, int>> new_vars;

	if (type == 1) {
		for (const auto& edge : edges_not_subset) {
			int i = std::get<0>(edge);
			int j = std::get<1>(edge);
			if (calculate_reduced_cost(cost[i][j], var_columns[i][j]) <= -epsilon) {
				new_vars.insert(edge);
			}
		}
	}
	else if (type == 2) {
		int no_new_edges = 0;
		for (const auto& edge : edges_not_subset) {
			int i = std::get<0>(edge);
			int j = std::get<1>(edge);
			if (calculate_reduced_cost(cost[i][j], var_columns[i][j]) <= -epsilon) {
				new_vars.insert(edge);
				no_new_edges++;
				if (no_new_edges == no_edges) break;
			}
		}
	}
	else {
		std::vector<RCpair> reduced_costs;
		for (const auto& edge : edges_not_subset) {
			int i = std::get<0>(edge);
			int j = std::get<1>(edge);
			double reduced_cost = calculate_reduced_cost(cost[i][j], var_columns[i][j]);
			if (reduced_cost <= epsilon) {
				reduced_costs.push_back(RCpair(edge, reduced_cost));
			}
		}

		std::sort(reduced_costs.begin(), reduced_costs.end(), [](const RCpair& a, const RCpair& b) {
			return a.RC < b.RC;
			});

		for (int k = 0; k < std::fmin(no_edges, reduced_costs.size()); k++) {
			new_vars.insert(reduced_costs[k].edge);
		}
	}
	return new_vars;
}


// Adds variables to the given model
bool add_variables(
	const int type,
	const int no_edges,
	const float epsilon,
	std::set<std::tuple<int, int>>& edges_not_subset,
	std::vector<std::vector<GRBColumn>>& var_columns,
	GRBModel* model,
	std::vector<std::vector<GRBVar>>& x,
	std::set<std::tuple<int, int>>& edges_subset,
	std::vector<std::vector<int>>& cost,
	double& time
) {

	auto start = std::chrono::high_resolution_clock::now();

	double ub;
	double obj;
	std::string name;
	int current_no_of_edges = edges_subset.size();
	std::set<std::tuple<int, int>> new_vars;

	while (true) {
		if (type == 0) {
			new_vars = check_reduced_costs(model, 1, no_edges, epsilon, edges_not_subset, cost, var_columns);
		}
		else {
			new_vars = check_reduced_costs(model, type, no_edges, epsilon, edges_not_subset, cost, var_columns);
		}

		for (const auto& edge : new_vars) {
			int i = std::get<0>(edge);
			int j = std::get<1>(edge);
			name = "x[" + itos(i) + "][" + itos(j) + "]";
			ub = (i == 0) ? 2.0 : 1.0; // UB = 2 if i = 0 else 1
			obj = (i < j) ? cost[i][j] : cost[j][i];

			x[i][j] = model->addVar(0.0, ub, obj, GRB_CONTINUOUS, var_columns[i][j], name);
			edges_not_subset.erase(edge);
		}
		edges_subset.insert(new_vars.begin(), new_vars.end());

		if (type != 0 || new_vars.size() == 0) break;
		else model->optimize();
		
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
	time += duration.count();
	return (edges_subset.size() > current_no_of_edges);
}


// Removes constraints from that model with slack larger that given epsilon
void remove_constraints_by_slack(
	const float slack_eps,
	int dimension,
	GRBModel* model,
	std::set<std::tuple<int, int>>& edges_not_subset,
	std::vector<std::vector<GRBColumn>>& var_columns,
	int& current_rcc,
	int & current_combs,
	int & current_mstars,
	int& current_glm,
	double& time
) {

	auto start = std::chrono::high_resolution_clock::now();

	GRBConstr* constrs = model->getConstrs();
	std::vector<GRBConstr> constrs_to_remove;
	std::unordered_set<std::string> constr_names_to_remove;
	for (int k = dimension - 1; k < model->get(GRB_IntAttr_NumConstrs); k++) {
		try {
			GRBConstr constr = constrs[k];
			double slack = constr.get(GRB_DoubleAttr_Slack);
			//td::cout << "Constraint index: " << k << ", Name: " << constr.get(GRB_StringAttr_ConstrName) << ", Slack: " << slack << std::endl;

			if (abs(slack) >= slack_eps) {
				constrs_to_remove.push_back(constr);
				constr_names_to_remove.insert(constr.get(GRB_StringAttr_ConstrName));
			}
		}
		catch (GRBException& e) {
			std::cerr << "Error accessing constraint: " << e.getMessage() << std::endl;
			// Handle the exception or log the error as needed
		}
	}

	for (const auto& edge : edges_not_subset) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		GRBColumn column = var_columns[i][j];
		std::vector<GRBConstr> constrs_in_column;
		std::vector<double> coeffs_in_column;
		for (int k = 0; k < column.size(); k++) {
			GRBConstr constr = column.getConstr(k);
			if (constr_names_to_remove.find(constr.get(GRB_StringAttr_ConstrName)) == constr_names_to_remove.end()) {
				constrs_in_column.push_back(constr);
				coeffs_in_column.push_back(column.getCoeff(k));
			}
		}
		GRBColumn new_column;
		new_column.addTerms(coeffs_in_column.data(), constrs_in_column.data(), coeffs_in_column.size());
		var_columns[i][j] = new_column;
	}

	for (GRBConstr constr : constrs_to_remove) {
		std::string cname = constr.get(GRB_StringAttr_ConstrName);
		if (cname.starts_with("RCC")) current_rcc--;
		else if (cname.starts_with("SC")) current_combs--;
		else if (cname.starts_with("MSTAR")) current_mstars--;
		else if (cname.starts_with("GLM")) current_glm--;
		model->remove(constr);
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
	time += duration.count();
}


// Eliminates variables based on reduced costs
void eliminate_variables(
	GRBModel* model,
	std::vector<std::vector<GRBVar>>& x,
	std::set<std::tuple<int, int>>& edges,
	std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int>>& edges_not_subset,
	std::vector<std::vector<int>>& cost,
	std::vector<std::vector<GRBColumn>>& var_columns,
	int upper_bound,
	float eps_for_elim,
	float eps_for_integrality,
	int& variables_eliminated,
	double& time
) {

	double obj = model->get(GRB_DoubleAttr_ObjVal);
	double gap = upper_bound - obj;

	for (const auto& edge : edges_not_subset) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		double reduced_cost = calculate_reduced_cost(cost[i][j], var_columns[i][j]);
		if (reduced_cost >= gap + eps_for_elim) {
			edges.erase(edge);
			edges_not_subset.erase(edge);
			variables_eliminated++;
		}
	}

	for (const auto& edge : edges_subset) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		double xval = x[i][j].get(GRB_DoubleAttr_X);
		if (xval <= eps_for_integrality) {
			double reduced_cost = x[i][j].get(GRB_DoubleAttr_RC);
			if (reduced_cost >= gap + eps_for_elim) {
				model->remove(x[i][j]);
				edges.erase(edge);
				edges_subset.erase(edge);
				variables_eliminated++;
			}
		}
		else if (xval >= 1 - eps_for_integrality) {
			double reduced_cost = x[i][j].get(GRB_DoubleAttr_RC);
			if (-reduced_cost >= gap + eps_for_elim) {
				double lb = (i > 0) ? 1.0 : 2.0;
				x[i][j].set(GRB_DoubleAttr_LB, lb);
				variables_eliminated++;
			}
		}
	}
}


// Adds all varaibles to the model
void add_all_variables(
	std::set<std::tuple<int, int>>& edges_not_subset,
	std::vector<std::vector<GRBColumn>>& var_columns,
	GRBModel* model,
	std::vector<std::vector<GRBVar>>& x,
	std::set<std::tuple<int, int>>& edges_subset,
	std::vector<std::vector<int>>& cost,
	double& time) {

	auto start = std::chrono::high_resolution_clock::now();

	double ub;
	double obj;
	std::string name;
	std::set<std::tuple<int, int>> new_vars;

	for (const auto& edge : edges_not_subset) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		name = "x[" + itos(i) + "][" + itos(j) + "]";
		ub = (i == 0) ? 2.0 : 1.0; // UB = 2 if i = 0 else 1
		obj = (i < j) ? cost[i][j] : cost[j][i];
		x[i][j] = model->addVar(0.0, ub, obj, GRB_CONTINUOUS, var_columns[i][j], name);
	}
}


// Adds all varaibles with negative reduced costs to the model. This function only works when flow varaiables are present
bool add_variables_with_flow(
	const float epsilon,
	std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int>>& edges_not_subset,
	std::set<std::tuple<int, int>>& asym_edges,
	std::vector<std::vector<GRBVar>>& x,
	std::vector<std::vector<std::vector<GRBVar>>>& f,
	std::set<int>& customers,
	double& time
) {
	auto start = std::chrono::high_resolution_clock::now();

	double ub;
	int new_vars = 0;
	for (const auto& edge : edges_not_subset) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (x[i][j].get(GRB_DoubleAttr_RC) <= -epsilon) {
			ub = (i == 0) ? 2.0 : 1.0; // UB = 2 if i = 0 else 1

			x[i][j].set(GRB_DoubleAttr_UB, ub);
			edges_not_subset.erase(edge);
			edges_subset.insert(edge);
			new_vars++;
		}
	}

	for (const auto& edge : asym_edges) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		for (int k : customers) {
			if (f[i][j][k].get(GRB_DoubleAttr_UB) == 0.0 && f[i][j][k].get(GRB_DoubleAttr_RC) <= -epsilon) {
				double ub = (i == 0) ? 2.0 : 1.0;
				f[i][j][k].set(GRB_DoubleAttr_UB, ub);
				new_vars++;
			}
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
	time += duration.count();
	return new_vars > 0;
}
