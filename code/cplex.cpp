#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <fstream>
#include <ilcplex/ilocplex.h>
using namespace std;
const int MAXN = 150;
int order, n, a[MAXN][MAXN];
char* fileIn, *fileOut;
int main(int argc, char** argv) {
    if (argc == 3){ 
        fileIn = argv[1];
        fileOut = argv[2];
    }
    ifstream FIC;
    FIC.open(fileIn);
    FILE* fp;
    fp = fopen(fileOut, "w+");
    FIC >> order;
    n = order * order;
    for (int i = 0;i < n; ++i)
        for (int j = 0; j < n; ++j)
            FIC >> a[i][j];
    IloEnv env;
    IloModel model(env);
    IloArray<IloArray<IloBoolVarArray>> x(env, n);
    for (int i = 0; i < n; ++i) {
        x[i] = IloArray<IloBoolVarArray>(env, n);
        for (int j = 0; j < n; ++j) {
            x[i][j] = IloBoolVarArray(env, n);
            for (int k = 0; k < n; ++k) {
                x[i][j][k] = IloBoolVar(env);
                model.add(x[i][j][k]);
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            IloExpr expr(env);
            for (int k = 0; k < n; ++k) {
                expr += x[i][j][k];
            }
            model.add(expr == 1);
            expr.end();
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            IloExpr exprRow(env);
            IloExpr exprCol(env);
            for (int j = 0; j < n; ++j) {
                exprRow += x[i][j][k];
                exprCol += x[j][i][k];
            }
            model.add(exprRow == 1);
            model.add(exprCol == 1);
            exprRow.end();
            exprCol.end();
        }
    }
    for (int i = 0; i < order; ++i) {
        for (int j = 0; j < order; ++j) {
            for (int k = 0; k < n; ++k) {
                IloExpr exprBox(env);
                for (int d1 = 0; d1 < order; ++d1)
                    for (int d2 = 0; d2 < order; ++d2) {
                        int ni = i * order + d1, nj = j * order + d2;
                        exprBox += x[ni][nj][k];
                    }
                model.add(exprBox == 1);
                exprBox.end();
            }
        }
    }
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (a[i][j] != 0)
                model.add(x[i][j][a[i][j] - 1] == 1);
    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, 3600);
    cplex.setParam(IloCplex::Param::Threads, 1);
    double startTime = cplex.getCplexTime();
    cplex.solve();
    double endTime = cplex.getCplexTime();
    if (cplex.getStatus() == IloAlgorithm::Optimal) {
        fprintf(fp, "Success! Time: %f\n", endTime - startTime);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    if (cplex.getValue(x[i][j][k]) > 0.5) {
                        //std::cout << k + 1 << " ";
                        fprintf(fp, "%d ", k + 1);
                    }
                }
            }
            fprintf(fp, "\n");
        }
    }
    else {
        fprintf(fp, "No solution found.\n");
    }
    env.end();
    return 0;
}