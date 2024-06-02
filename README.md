# Algorithms-for-the-CVRP
Code used for the separation of inequalities and solving large-scale CVRP instances.

This code requires Gurobi and cURL. For configuring Gurobi in Visual Studio see https://support.gurobi.com/hc/en-us/articles/360013194392-How-do-I-configure-a-new-Gurobi-C-project-with-Microsoft-Visual-Studio.
As I do not know how to ensure these are included in the repository, these must be fixed manually. Moreover, the corresponding #include statements may need to be changed.

Instances are loaded from "/Instance Data/Instances/" whilst the best known solutions are loaded from "/Instance Data/Solutions/".

LP models are saved in "/Results/LP_models/" and optimal solutions found are saved in "/Results/Optimal Solutions/".
