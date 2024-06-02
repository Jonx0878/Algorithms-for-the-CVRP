#pragma once
#include "cvrpsep/mstarsep.h"

struct MSTAR_Inequality {
	std::set<int> nucleus;
	std::set<int> connectors;
	std::set<int> satellites;
	int A = 1; 
	int B = 1;
	int L = 1;
	bool check_violation = true;
	bool from_exact = false;
};

// Convert sets obtained from shrunk graph to base graph for HPM inequalities
void convert_cuts_mstar(std::vector<MSTAR_Inequality>& cuts_found, ReachPtr graph_mapping) {
	for (int i = 0; i < cuts_found.size(); i++) {
		std::set<int> new_nucleus, new_connectors, new_satellites;
		MSTAR_Inequality mstar = cuts_found[i];
		if (mstar.from_exact) {
			for (int j : mstar.nucleus) {
				for (int l = 1; l <= graph_mapping->LP[j].CFN; l++)
				{
					int k = graph_mapping->LP[j].FAL[l];
					new_nucleus.insert(k);

				}
			}
			for (int j : mstar.connectors) {
				for (int l = 1; l <= graph_mapping->LP[j].CFN; l++)
				{
					int k = graph_mapping->LP[j].FAL[l];
					new_connectors.insert(k);
				}
			}

			for (int j : mstar.satellites) {
				for (int l = 1; l <= graph_mapping->LP[j].CFN; l++)
				{
					int k = graph_mapping->LP[j].FAL[l];
					new_satellites.insert(k);
				}
			}

			int A_mu = 2 * mstar.B * mstar.nucleus.size() - mstar.L;
			int new_L = 2 * mstar.B * new_nucleus.size() - A_mu;

			cuts_found[i] = MSTAR_Inequality(new_nucleus, new_connectors, new_satellites, mstar.A, mstar.B, new_L, mstar.check_violation, true);
		}
	}
}


// Determines lhs, sense, and rhs of an HPM inequality of type 1
void determine_mstar_type1(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const MSTAR_Inequality& mstar,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset
) {
	for (const auto& edge : edges_in(edges, mstar.nucleus)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += 2 * mstar.B * x[i][j];
		else edges_not_subset.insert({ i, j, 2.0 * mstar.B });
	}

	for (const auto& edge : edges_from_to(edges, mstar.connectors, mstar.satellites)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += mstar.A * x[i][j];
		else edges_not_subset.insert({ i, j, 1.0*mstar.A });
	}

	sense = GRB_LESS_EQUAL;
	rhs = 2.0f * mstar.B * int(mstar.nucleus.size()) - mstar.L;
}


// Determines lhs, sense, and rhs of an HPM inequality of type 2
void determine_mstar_type2(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const MSTAR_Inequality& mstar,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<int>& vertices,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset
) {

	std::map<std::tuple<int, int>, double> edge_count_map;

	for (const auto& edge : delta_edges(edges, vertices, mstar.nucleus)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += mstar.B * x[i][j];
		else edge_count_map[edge] += mstar.B;
	}

	for (const auto& edge : edges_from_to(edges, mstar.connectors, mstar.satellites)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs -= mstar.A * x[i][j];
		else edge_count_map[edge] -= mstar.A;
	}
	
	for (const auto& [edge, coeff] : edge_count_map) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		edges_not_subset.insert({ i, j, coeff });
	}
	sense = GRB_GREATER_EQUAL;
	rhs = mstar.L;
}


// Determines lhs, sense, and rhs of an HPM inequality of type 3
void determine_mstar_type3(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const MSTAR_Inequality& mstar,
	const std::set<int>& nucleus_bar,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset
) {
	for (const auto& edge : edges_in(edges, nucleus_bar)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += 2*mstar.B*x[i][j];
		else edges_not_subset.insert({ i, j, mstar.B });
	}
	for (const auto& j : nucleus_bar) {
		if (edges_subset.count({ 0, j }) == 1) lhs += mstar.B * x[0][j];
		else edges_not_subset.insert({ 0, j, mstar.B });
	}
	for (const auto& j : mstar.nucleus) {
		if (edges_subset.count({ 0, j }) == 1) lhs -= mstar.B * x[0][j];
		else edges_not_subset.insert({ 0, j, -mstar.B });
	}

	for (const auto& edge : edges_from_to(edges, mstar.connectors, mstar.satellites)) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (edges_subset.count({ i, j }) == 1) lhs += mstar.A * x[i][j];
		else edges_not_subset.insert({ i, j, mstar.A });
	}

	sense = GRB_LESS_EQUAL;
	rhs = 2.0f*mstar.B*int(nucleus_bar.size()) - mstar.L;
}


// Determines lhs, sense, and rhs of an HPM inequality
void determine_mstar(
	GRBLinExpr& lhs, char& sense, int& rhs,
	const MSTAR_Inequality& mstar,
	const std::vector<std::vector<GRBVar>>& x,
	const std::set<int>& vertices,
	const std::set<int>& customers,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::set<std::tuple<int, int, double>>& edges_not_subset,
	const float& type2_max,
	const float& type3_min,
	const float& type3_max
) {
	int nucleus_size = static_cast<int>(mstar.nucleus.size());
	std::set<int> nucleus_bar;

	if (nucleus_size <= type2_max) {
		determine_mstar_type2(lhs, sense, rhs, mstar, x, vertices, edges, edges_subset, edges_not_subset);
	}
	else if (type3_min <= nucleus_size && nucleus_size <= type3_max) {
		std::set_difference(customers.begin(), customers.end(), mstar.nucleus.begin(), mstar.nucleus.end(), std::inserter(nucleus_bar, nucleus_bar.begin()));
		determine_mstar_type3(lhs, sense, rhs, mstar, nucleus_bar, x, edges, edges_subset, edges_not_subset);
	}
	else determine_mstar_type1(lhs, sense, rhs, mstar, x, edges, edges_subset, edges_not_subset);
}


// Adds HPM inequalities to the model
void add_mstars(
	int& no_new_cuts,
	const std::vector<MSTAR_Inequality>& cuts,
	GRBModel* model,
	const std::vector<std::vector<GRBVar>>& x,
	const int& dimension,
	int& total_mstars,
	int& current_mstars,
	const std::set<int>& vertices,
	const std::set<int>& customers,
	const std::set<std::tuple<int, int>>& edges,
	const std::set<std::tuple<int, int>>& edges_subset,
	std::vector<std::vector<GRBColumn>>& var_columns,
	const float& type2_max,
	const float& type3_min,
	const float& type3_max,
	const float& eps_for_viol,
	double& time,
	CnstrMgrPointer MyCutsCMP
) {
	GRBLinExpr lhs;
	char sense = GRB_LESS_EQUAL;
	int rhs = 0;
	for (const auto& mstar : cuts) {
		lhs.clear();
		bool check_violation = mstar.check_violation;
		bool add_to_cnstrmgr = mstar.from_exact;
		std::set<std::tuple<int, int, double>> edges_not_subset;

		auto start = std::chrono::high_resolution_clock::now();

		determine_mstar(lhs, sense, rhs, mstar, x, vertices, customers,
			edges, edges_subset, edges_not_subset,
			type2_max, type3_min, type3_max);
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
		time += duration.count();

		if (check_violation) {
			if ((sense == *("<") && lhs.getValue() >= rhs + eps_for_viol) || (sense == *(">") && lhs.getValue() <= rhs - eps_for_viol))
				goto addcut;
		}
		else {
		addcut:
			total_mstars++;
			current_mstars++;
			no_new_cuts++;
			std::string name = "MSTAR" + itos(total_mstars);
			GRBConstr constr = model->addConstr(lhs, sense, rhs, name);
			for (const auto& var : edges_not_subset) {
				int i = std::get<0>(var);
				int j = std::get<1>(var);
				double coeff = std::get<2>(var);
				var_columns[i][j].addTerm(coeff, constr);
			}


			if (add_to_cnstrmgr) {
				std::vector<int> nucleus_list(mstar.nucleus.size() + 1);
				std::vector<int> satellites_list(mstar.satellites.size() + 1);
				std::vector<int> connectors_list(mstar.connectors.size() + 1);

				int index = 1;
				for (int i : mstar.nucleus) {
					nucleus_list[index] = i;
					index++;
				}
				index = 1;
				for (int i : mstar.satellites) {
					satellites_list[index] = i;
					index++;
				}
				index = 1;
				for (int i : mstar.connectors) {
					connectors_list[index] = i;
					index++;
				}

				CMGR_AddPartialMStar(
					MyCutsCMP,
					CMGR_CT_MSTAR, 0,
					nucleus_list.size()-1, &nucleus_list[0],
					satellites_list.size()-1, &satellites_list[0],
					connectors_list.size()-1, &connectors_list[0],
					mstar.A, mstar.B, mstar.L
				);
			}
		}
	}
}


// Separates HPM inequalities using CVRPSEP
void separate_mstar_heuristically(
	double& max_violation,
	std::vector<MSTAR_Inequality>& new_cuts,
	const std::vector<std::tuple<int, int, double>>& solution,
	const int& dimension,
	std::vector<int>& demand,
	const int& capacity,
	CnstrMgrPointer MyOldCutsCMP,
	CnstrMgrPointer MyCutsCMP
)
{
	int NoOfCustomers = dimension - 1;
	std::vector<double> demand_double(demand.begin(), demand.end());
	const double* Demand = &demand_double[0];
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

	MSTARSEP_SeparateMultiStarCuts(
		NoOfCustomers,
		Demand,
		capacity,
		NoOfEdges,
		EdgeTail,
		EdgeHead,
		EdgeX,
		MyOldCutsCMP,
		MaxNoOfCuts,
		&max_violation,
		MyCutsCMP);

	//Retrieve Cuts
	int A, B, L;
	std::set<int> nucleus, connectors, satellites;
	int* NList, * TList, * CList;
	/* Allocate memory for n+1 integers in each of */
	/* the vectors NList, TList, CList */

	NList = new int [dimension];
	TList = new int[dimension];
	CList = new int[dimension];

	for (int i = 0; i < MyCutsCMP->Size; i++){
		if (MyCutsCMP->CPL[i]->CType == CMGR_CT_MSTAR)
		{
			int Nsize = 0;
			/* Nucleus: */
			for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; j++) {
				NList[j] = MyCutsCMP->CPL[i]->IntList[j];
				Nsize++;
				nucleus.insert(MyCutsCMP->CPL[i]->IntList[j]);
			}
			/* Satellites: */
			for (int j = 1; j <= MyCutsCMP->CPL[i]->ExtListSize; j++) {
					TList[j] = MyCutsCMP->CPL[i]->ExtList[j];
					satellites.insert(MyCutsCMP->CPL[i]->ExtList[j]);
			}
			/* Connectors: */
			for (int j = 1; j <= MyCutsCMP->CPL[i]->CListSize; j++) {
					CList[j] = MyCutsCMP->CPL[i]->CList[j];
					connectors.insert(MyCutsCMP->CPL[i]->CList[j]);
			}
			/* Coefficients of the cut: */
			A = MyCutsCMP->CPL[i]->A;
			B = MyCutsCMP->CPL[i]->B;
			L = MyCutsCMP->CPL[i]->L;
			/* Lambda=L/B, Sigma=A/B */
			/* Add the cut to the LP */

			float sigma = float(A) / B;
			float rho = float(L) / B;
			float lambda = 1 / sigma;
			float mu = lambda * (Nsize - rho);

			new_cuts.emplace_back(MSTAR_Inequality(nucleus, connectors, satellites, A, B, L, (max_violation < 0.05f), false));
			nucleus.clear();
			connectors.clear();
			satellites.clear();
		}
	}
}


// Sets up the model used for separating HPM inequalities "exactly"
GRBModel* setup_mstar_model(
	std::vector<GRBVar>& n,
	std::vector<GRBVar>& c,
	std::vector<GRBVar>& s_true,
	GRBVar& lambda,
	GRBVar& mu,
	const std::vector<std::tuple<int, int, double>>& solution,
	const std::set<int>& customers,
	const std::set<int>& vertices,
	const int& capacity,
	const std::vector<int>& demand,
	const int& gcd,
	const int& standard_dim,
	GRBEnv env,
	ReachPtr supervertices
) {
	// Create Model
	GRBModel* model = new GRBModel(env);
	std::string name;

	// Calculate constants
	int no_of_locations = standard_dim;
	int max_alpha = no_of_locations - 2;

	std::vector<int> locations_smaleq(vertices.size() + 1);
	for (int i : customers) {
		locations_smaleq[i] = 0;
		for (int j : customers) {
			if (demand[i] >= demand[j]) locations_smaleq[i]++;
		}
	}

	double max_vehicles = 0;
	for (int i : customers) {
		max_vehicles += demand[i];
	}

	int max_demand = *std::max_element(demand.begin(), demand.end());

	std::vector<int> Log_UB(max_alpha+1);
	double LB;
	for (int alpha = 0; alpha <= max_alpha; alpha++) {
		Log_UB[alpha] = ceil(log2(no_of_locations - 2 + alpha));
	}

	// Setup variables
	n.resize(vertices.size() + 1);
	c.resize(vertices.size() + 1);
	std::vector<std::vector<GRBVar>> s(no_of_locations+1);
	GRBVar nu_CS;
	std::vector<GRBVar> nu(no_of_locations+1);
	std::vector<GRBVar> alpha_CS(no_of_locations+1);
	std::vector<GRBVar> alpha_S(no_of_locations+1);

	GRBVar UB_CS;
	GRBVar UB_CS1;
	GRBVar UB_CS2;
	GRBVar UB_CS3;
	std::vector<GRBVar> UB_alpha(no_of_locations+1);

	std::vector<GRBVar> lam_n(no_of_locations+1);
	std::vector<GRBVar> lam_alpha_CS(no_of_locations+1);

	s_true.resize(vertices.size()+1);

	std::vector<std::vector<GRBVar>> UBALPHA_LIN(no_of_locations+1);
	std::vector<std::vector<GRBVar>> LAMUBALP(no_of_locations+1);

	std::vector<std::vector<GRBVar>> ALPHASS(no_of_locations+1);


	std::vector<GRBVar> UB_alpha_y(no_of_locations+1);
	
	for (int i : vertices) {
		s[i].resize(locations_smaleq[i] + 1);
		ALPHASS[i].resize(max_alpha + 1);
	}

	for (int alpha = 0; alpha <= max_alpha; alpha++){
		UBALPHA_LIN[alpha].resize(Log_UB[alpha] + 1);
		LAMUBALP[alpha].resize(Log_UB[alpha] + 1);
	}

	// Change parameters
	//model->getEnv().set(GRB_IntParam_OutputFlag, 0);
	//model->set(GRB_DoubleParam_Cutoff, 0.0);
	//model->set(GRB_DoubleParam_BestObjStop, 0.2);
	double M = max_alpha;

	// Add Variables
	lambda = model->addVar(1.0, M, NULL, GRB_CONTINUOUS, "lambda");
	mu = model->addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS, "mu");

	nu_CS = model->addVar(1, max_vehicles, NULL, GRB_INTEGER, "nu_CS");
	for (int alpha = 0; alpha <= max_alpha; alpha++) {
		nu[alpha] = model->addVar(1, max_vehicles, NULL, GRB_INTEGER, "nu[" + std::to_string(alpha) + "]");
	}

	for (int i : customers) {
		n[i] = model->addVar(0, 1, NULL, GRB_BINARY, "n[" + std::to_string(i) + "]");
		c[i] = model->addVar(0, 1, NULL, GRB_BINARY, "c[" + std::to_string(i) + "]");
		s_true[i] = model->addVar(0, 1, NULL, GRB_BINARY, "s_true[" + std::to_string(i) + "]");
		for (int k = 1; k <= locations_smaleq[i]; k++) {
			s[i][k] = model->addVar(0, 1, NULL, GRB_BINARY, "s[" + std::to_string(i) + "," + std::to_string(k) + "]");
		}
		lam_n[i] = model->addVar(0, M, NULL, GRB_CONTINUOUS, "lam_n[" + std::to_string(i) + "]");

		for (int alpha = 2; alpha <= max_alpha; alpha++) {
			ALPHASS[i][alpha] = model->addVar(0, 1, NULL, GRB_CONTINUOUS, "ALPHA_SS[" + std::to_string(i) + "," + std::to_string(alpha) + "]");
		}
	}

	for (int alpha = 2; alpha <= max_alpha; alpha++) {
		alpha_CS[alpha] = model->addVar(0, 1, NULL, GRB_BINARY, "alpha_CS[" + std::to_string(alpha) + "]");
		alpha_S[alpha] = model->addVar(0, 1, NULL, GRB_BINARY, "alpha_S[" + std::to_string(alpha) + "]");

		lam_alpha_CS[alpha] = model->addVar(0, M, NULL, GRB_CONTINUOUS, "lam_alpha_CS[" + std::to_string(alpha) + "]");
	}

	UB_CS = model->addVar(1, max_alpha, NULL, GRB_CONTINUOUS, "UB_CS");
	UB_CS1 = model->addVar(2, 2 * (no_of_locations-2), NULL, GRB_CONTINUOUS, "UB_CS1");
	UB_CS2 = model->addVar(2, 2 * (no_of_locations-2), NULL, GRB_CONTINUOUS, "UB_CS2");
	UB_CS3 = model->addVar(1, no_of_locations - 2, NULL, GRB_CONTINUOUS, "UB_CS3");

	for (int alpha = 0; alpha <= max_alpha; alpha++) {
		alpha % 2 == 0 ? LB = 1 - alpha / 2.0f : 1 - (alpha + 1) / 2.0f;
		UB_alpha[alpha] = model->addVar(-alpha, no_of_locations - 2, NULL, GRB_CONTINUOUS, "UB_alpha[" + std::to_string(alpha) + "]");
		if (alpha >= 2) {
			UB_alpha_y[alpha] = model->addVar(0, 1, NULL, GRB_BINARY, "UB_alpha_y[" + std::to_string(alpha) + "]");
		}

		for (int i = 0; i <= Log_UB[alpha]; i++) {
			UBALPHA_LIN[alpha][i] = model->addVar(0, 1, NULL, GRB_BINARY, "UBALPHA_IN[" + std::to_string(alpha) + "," + std::to_string(i) + "]");
			LAMUBALP[alpha][i] = model->addVar(0, M, NULL, GRB_CONTINUOUS, "LAMUBALP[" + std::to_string(alpha) + "," + std::to_string(i) + "]");
		}
	}
	UB_alpha[0].set(GRB_DoubleAttr_LB, 1.0);
	model->update();

	// Add Constraints
	GRBLinExpr Lexpr;
	GRBLinExpr Lexpr2;
	GRBLinExpr Lexpr3;
	GRBLinExpr Lexpr4;


	//Set Constraints
	if (supervertices == NULL) {
		for (int i : customers) {
			Lexpr += n[i];
		}
	}
	else {
		for (int i : customers) {
			Lexpr += supervertices->LP[i].CFN * n[i];
		}
	}
	for (int i : customers) {
		Lexpr2 += c[i];
		model->addConstr(c[i] <= n[i], "c_bound" + std::to_string(i));
	}
	model->addConstr(Lexpr >= 1);
	model->addConstr(Lexpr <= no_of_locations - 2);

	model->addConstr(Lexpr2 >= 1);
	Lexpr.clear();
	Lexpr2.clear();

	for (int i : customers) {
		Lexpr += s_true[i];
		for (int k = 1; k <= locations_smaleq[i]; k++) {
			Lexpr2 += s[i][k];
		}
		model->addConstr(s_true[i] == Lexpr2);
		model->addConstr(s_true[i] <= 1 - n[i], "s_bound" + std::to_string(i));
		Lexpr2.clear();
	}
	model->addConstr(Lexpr >= 1);
	Lexpr.clear();

	for (int k = 1; k < no_of_locations - 1; k++) {
		for (int i : customers) {
			if (k <= locations_smaleq[i]) {
				Lexpr += demand[i] * s[i][k];
				Lexpr3 += s[i][k];
			}
			if (k + 1 <= locations_smaleq[i]) {
				Lexpr2 += demand[i] * s[i][k + 1];
				Lexpr4 += s[i][k + 1];
			}
		}
		model->addConstr(Lexpr - Lexpr3 * max_demand <= Lexpr2 - Lexpr4 * max_demand);
		model->addConstr(Lexpr3 <= 1);
		Lexpr.clear();
		Lexpr2.clear();
		Lexpr3.clear();
		Lexpr4.clear();
	}

	// Vehicle Constraints
	for (int i : customers) {
		Lexpr += demand[i] * c[i];
		for (int k = 1; k <= locations_smaleq[i]; k++) {
			Lexpr += demand[i] * s[i][k];
		}
	}
	model->addConstr(capacity * (nu_CS - 1) + gcd, GRB_LESS_EQUAL, Lexpr, "nu_CS_UB");
	model->addConstr(capacity * (nu_CS), GRB_GREATER_EQUAL, Lexpr, "nu_CS_LB");
	Lexpr.clear();

	for (int i : customers) {
		Lexpr += demand[i] * n[i];
	}
	model->addConstr(capacity * (nu[0] - 1) + gcd, GRB_LESS_EQUAL, Lexpr, "nu_" + std::to_string(0) + "_UB");
	model->addConstr(capacity * (nu[0]), GRB_GREATER_EQUAL, Lexpr, "nu_" + std::to_string(0) + "_LB");
	Lexpr.clear();

	for (int alpha = 1; alpha <= max_alpha; alpha++) {
		for (int i : customers) {
			Lexpr += demand[i] * n[i];
			for (int k = 1; k <= std::min(alpha, locations_smaleq[i]); k++) {
				Lexpr += demand[i] * s[i][k];
			}
		}
		model->addConstr(capacity * (nu[alpha] - 1) + gcd, GRB_LESS_EQUAL, Lexpr, "nu_" + std::to_string(alpha) + "_UB");
		model->addConstr(capacity * (nu[alpha]), GRB_GREATER_EQUAL, Lexpr, "nu_" + std::to_string(alpha) + "_LB");
		model->addConstr(nu[alpha] >= nu[alpha - 1]);
		Lexpr.clear();
	}

	// UB_CS Constraints
	if (supervertices == NULL) {
		for (int i : customers) {
			Lexpr += c[i];
			Lexpr2 += s_true[i];
		}
	}
	else {
		for (int i : customers) {
			Lexpr += supervertices->LP[i].CFN * c[i];
			Lexpr2 += supervertices->LP[i].CFN * s_true[i];
		}
	}

	model->addConstr(UB_CS1 == 2 * Lexpr);
	model->addConstr(UB_CS2 == 2 * Lexpr2);
	model->addConstr(UB_CS3 == Lexpr + Lexpr2 - nu_CS);

	GRBVar* vars = new GRBVar[3];
	vars[0] = UB_CS1;
	vars[1] = UB_CS2;
	vars[2] = UB_CS3;
	model->addGenConstrMin(UB_CS, vars, 3);

	Lexpr.clear();
	Lexpr2.clear();

	// alpha_CS constraints
	for (int alpha = 2; alpha <= max_alpha; alpha++) {
		model->addQConstr(alpha_CS[alpha] * UB_CS >= UB_CS - alpha + 1, "alpha_CS_UB_" + std::to_string(alpha));
		model->addConstr(alpha * alpha_CS[alpha] <= UB_CS, "alpha_CS_LB_" + std::to_string(alpha));
	}

	// alpha_S constraints
	if (supervertices == NULL) {
		for (int i : customers) {
			Lexpr2 += s_true[i];
		}
	}
	else {
		for (int i : customers) {
			Lexpr2 += supervertices->LP[i].CFN * s_true[i];
		}
	}
	for (int alpha = 2; alpha <= max_alpha; alpha++) {
		for (int i : customers) {
			Lexpr2 += ALPHASS[i][alpha];
		}
		model->addQConstr(alpha_S[alpha] * Lexpr >= Lexpr - alpha + 1, "alpha_S_UB_" + std::to_string(alpha));
		model->addConstr(alpha * alpha_S[alpha] <= Lexpr, "alpha_S_LB_" + std::to_string(alpha));
		Lexpr2.clear();
	}
	Lexpr.clear();


	// UB_alpha constraints
	if (supervertices == NULL) {
		for (int i : customers) {
			Lexpr += n[i];
			Lexpr2 += s_true[i];
		}
	}
	else {
		for (int i : customers) {
			Lexpr += supervertices->LP[i].CFN * n[i];
			Lexpr2 += supervertices->LP[i].CFN * s_true[i];
		}
	}

	model->addConstr(UB_alpha[0] == Lexpr - nu[0], "UB_alpha_" + std::to_string(0));
	model->addConstr(UB_alpha[1] == Lexpr - nu[1], "UB_alpha_" + std::to_string(1));
	for (int alpha = 2; alpha <= max_alpha; alpha++) {
		for (int i : customers) {
			Lexpr3 += ALPHASS[i][alpha];
		}

		double alpha_round_up = (alpha % 2 == 0) ? alpha / 2.0f : (alpha + 1) / 2.0f;

		model->addConstr(UB_alpha[alpha] <= Lexpr - alpha_round_up);
		model->addConstr(UB_alpha[alpha] <= Lexpr + Lexpr2 - Lexpr3 - nu[alpha] - (1 - alpha_S[alpha]) * alpha);
		model->addConstr(UB_alpha[alpha] >= Lexpr - alpha_round_up - (no_of_locations - 2 - alpha_round_up + alpha)*(1 - UB_alpha_y[alpha]));
		model->addConstr(UB_alpha[alpha] >= Lexpr + Lexpr2 - Lexpr3 - nu[alpha] - (1 - alpha_S[alpha]) * alpha - (no_of_locations - 3 + alpha) * UB_alpha_y[alpha]);

		Lexpr3.clear();
	}
	Lexpr.clear();
	Lexpr2.clear();

	model->addConstr(lambda <= UB_CS);

	// main constraints
	for (int i = 0; i <= Log_UB[0]; i++) {
		Lexpr += pow(2, i) * LAMUBALP[0][i];
	}
	model->addConstr(Lexpr <= mu, "main_0");
	Lexpr.clear();
	for (int i = 0; i <= Log_UB[1]; i++) {
		Lexpr += pow(2, i) * LAMUBALP[1][i];
	}
	model->addConstr(Lexpr - lambda + 1 <= mu, "main_1");
	Lexpr.clear();

	for (int alpha = 2; alpha <= max_alpha; alpha++) {
		for (int i = 0; i <= Log_UB[alpha]; i++) {
			Lexpr += pow(2, i) * LAMUBALP[alpha][i];
		}
		model->addConstr(Lexpr - lam_alpha_CS[alpha] * alpha + alpha_CS[alpha] * alpha <= mu, "main_" + std::to_string(alpha));
		Lexpr.clear();
	}

	//Lam * n constraints
	for (int i : customers) {
		model->addConstr(lam_n[i] <= n[i] * M);
		//model->addConstr(lam_n[i] >= lambda - (1 - n[i]) * M);
		model->addConstr(lam_n[i] <= lambda);
	}

	// Lam * alpha_CS
	for (int alpha = 2; alpha <= max_alpha; alpha++) {
		model->addConstr(lam_alpha_CS[alpha] <= alpha_CS[alpha] * M);
		model->addConstr(lam_alpha_CS[alpha] >= lambda - (1 - alpha_CS[alpha]) * M);
		model->addConstr(lam_alpha_CS[alpha] <= lambda);
	}

	// LAMUBALP
	for (int alpha = 0; alpha <= max_alpha; alpha++) {
		for (int i = 0; i <= Log_UB[alpha]; i++) {
			Lexpr += pow(2, i) * UBALPHA_LIN[alpha][i];
			if (alpha < 2) {
				model->addConstr(LAMUBALP[alpha][i] <= UBALPHA_LIN[alpha][i] * max_alpha);
				model->addConstr(LAMUBALP[alpha][i] >= lambda - (1 - UBALPHA_LIN[alpha][i]) * max_alpha);
				model->addConstr(LAMUBALP[alpha][i] <= lambda);
			}
			else {
				model->addConstr(LAMUBALP[alpha][i] <= UBALPHA_LIN[alpha][i] * max_alpha);
				model->addConstr(LAMUBALP[alpha][i] >= lam_alpha_CS[alpha] - (1 - UBALPHA_LIN[alpha][i]) * max_alpha);
				model->addConstr(LAMUBALP[alpha][i] <= lam_alpha_CS[alpha]);
			}
		}
		model->addConstr(UB_alpha[alpha] == Lexpr - alpha);
		Lexpr.clear();
	}

	// ALPHASS
	for (int i : customers) {
		for (int alpha = 2; alpha <= max_alpha; alpha++) {
			model->addConstr(ALPHASS[i][alpha] >= alpha_S[alpha] + s_true[i] - 1);
			model->addConstr(ALPHASS[i][alpha] <= alpha_S[alpha]);
			model->addConstr(ALPHASS[i][alpha] <= s_true[i]);
		}
	}

	// Set Objective
	GRBQuadExpr obj = -mu;

	for (const auto& edge : solution) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		double value = std::get<2>(edge);

		if (i > 0) {
			obj += value * (lam_n[i] * n[j] + c[i] * s_true[j] + c[j] * s_true[i]);
			Lexpr.clear();
			Lexpr2.clear();
		}
	}
	model->setObjective(obj, GRB_MAXIMIZE);

	return model;
}


// Separates HPM inequalities "exactly"
void separate_mstar_exactly(
	double& max_violation,
	std::vector<MSTAR_Inequality>& new_cuts,
	ReachPtr shrunk_mapping,
	const std::vector<std::tuple<int, int, double>>& solution,
	const std::set<int>& vertices,
	const std::set<int>& customers,
	const int& capacity,
	const std::vector<int>& demand,
	int& standard_dim,
	const float& eps_for_early_termination,
	const int& gcd,
	const float& eps_for_int, const float& eps_for_violation, double& time, GRBEnv env
) {
	std::vector<GRBVar> n;
	std::vector<GRBVar> c;
	std::vector<GRBVar> s;
	GRBVar lambda;
	GRBVar mu;
	int A, B, L;
	double lam_value;
	double mu_value;

	// Setup separation model
	GRBModel* model = setup_mstar_model(n, c, s, lambda, mu, solution, customers, vertices, capacity, demand, gcd, standard_dim, env, shrunk_mapping);
	if (eps_for_early_termination > 0) model->getEnv().set(GRB_DoubleParam_BestObjStop, eps_for_early_termination);

	auto start = std::chrono::high_resolution_clock::now();

	model->optimize();

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
	time += duration.count();

	if (model->get(GRB_DoubleAttr_ObjVal) > max_violation) max_violation = model->get(GRB_DoubleAttr_ObjVal);

	for (int sol = 0; sol < model->get(GRB_IntAttr_SolCount); ++sol) {
		std::set<int> nucleus, connectors, satellites;
		model->set(GRB_IntParam_SolutionNumber, sol);
		if (model->get(GRB_DoubleAttr_PoolObjVal) >= eps_for_violation) {
			std::set<int> S;
			for (int i : customers) {
				if (n[i].get(GRB_DoubleAttr_Xn) >= 1 - eps_for_int) {
					nucleus.insert(i);
				}
				if (c[i].get(GRB_DoubleAttr_Xn) >= 1 - eps_for_int) {
					connectors.insert(i);
				}
				if (s[i].get(GRB_DoubleAttr_Xn) >= 1 - eps_for_int) {
					satellites.insert(i);
				}
			}
			A = 1;
			lam_value = lambda.get(GRB_DoubleAttr_Xn);
			mu_value = mu.get(GRB_DoubleAttr_Xn);

			while (abs(remainder(lam_value, 2)) >= 0.001f || abs(remainder(mu_value, 1)) >= 0.001f) {
				A++;
				lam_value = A * lambda.get(GRB_DoubleAttr_Xn);
				mu_value = A * mu.get(GRB_DoubleAttr_Xn);
			}
			B = round(lam_value / 2.0f);
			L = 2.0f * B * nucleus.size() - mu_value;

			new_cuts.emplace_back(MSTAR_Inequality(nucleus, connectors, satellites, A, B, L, (model->get(GRB_DoubleAttr_PoolObjVal) < 0.1), true));
		}
	}
}
