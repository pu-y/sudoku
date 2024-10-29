#include "gurobi_c++.h"
#include <iostream>
#include <fstream>
using namespace std;

int sd, n;
char* fileIn, * fileOut;

int main(int argc, char* argv[])
{
    if (argc == 3) {
        fileIn = argv[1];
        fileOut = argv[2];
    }
    else {
        cout << "Error!" << endl;
        return 0;
    }
    ifstream FIC;
    FIC.open(fileIn);
    FILE* fp;
    fp = fopen(fileOut, "w+");
    try {
        GRBEnv env = GRBEnv();
        env.set(GRB_IntParam_Threads, 1);
        GRBModel model = GRBModel(env);
        model.set(GRB_DoubleParam_TimeLimit, 3600.0);

        // Create 3-D array of model variables
        GRBVar vars[81][81][81];

        int i, j, v;
        FIC >> sd;
        n = sd * sd;

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (v = 0; v < n; v++) {
                    vars[i][j][v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                }
            }
        }

        // Add constraints

        // Each cell must take one value

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                GRBLinExpr expr = 0;
                for (v = 0; v < n; v++)
                    expr += vars[i][j][v];
                model.addConstr(expr, GRB_EQUAL, 1.0);
            }
        }

        // Each value appears once per row

        for (i = 0; i < n; i++) {
            for (v = 0; v < n; v++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < n; j++)
                    expr += vars[i][j][v];
                model.addConstr(expr == 1.0);
            }
        }

        // Each value appears once per column

        for (j = 0; j < n; j++) {
            for (v = 0; v < n; v++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++)
                    expr += vars[i][j][v];
                model.addConstr(expr == 1.0);
            }
        }

        // Each value appears once per sub-grid

        for (v = 0; v < n; v++) {
            for (int i0 = 0; i0 < sd; i0++) {
                for (int j0 = 0; j0 < sd; j0++) {
                    GRBLinExpr expr = 0;
                    for (int i1 = 0; i1 < sd; i1++) {
                        for (int j1 = 0; j1 < sd; j1++) {
                            expr += vars[i0 * sd + i1][j0 * sd + j1][v];
                        }
                    }
                    model.addConstr(expr == 1.0);
                }
            }
        }

        // Fix variables associated with pre-specified cells

        int input;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                FIC >> input;
                int val = input - 1; // 0-based
                if (val >= 0)
                    vars[i][j][val].set(GRB_DoubleAttr_LB, 1.0);
            }
        }

        // Optimize model

        model.optimize();
        double solveTime = model.get(GRB_DoubleAttr_Runtime);

        fprintf(fp, "Time: %f\n", solveTime);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (v = 0; v < n; v++) {
                    if (vars[i][j][v].get(GRB_DoubleAttr_X) > 0.5)
                        fprintf(fp, "%d ", v + 1);
                }
            }
            fprintf(fp, "\n");
        }
    }
    catch (GRBException e) {
        fprintf(fp, "Error code = %d\n", e.getErrorCode());
        fprintf(fp, "%s\n", e.getMessage().c_str());
    }
    catch (...) {
        fprintf(fp, "Error during optimization\n");
    }

    return 0;
}