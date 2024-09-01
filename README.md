## Project Overview
This project implements a Iterated Local Search metaheuristic so solve a variation of the TSP. In more depth, we are given a list of VIP customers that we must visit before others to ensure we have enough products for them. The classic TSP does not accommodate these precedence constraints, as such a new heuristic solver is developed the CTSP.

This project is my submission for the course 42137 Optimization using Metaheuristics, taught by Thomas Jacob Riis Stidsen and Dario Pacino at DTU.

### Key aspects of this algorithm include:
-   **Neighborhood Operator**: The 3-opt move is employed as the neighborhood operator depending on the file you choose.
-   **Termination Criteria**: The algorithm uses a time limit as the primary termination criterion.
-   **Step Criterion**: Adopts the "Best Improvement" strategy.
-   **Cost**: The algorithm operates based on a cost matrix, which is specified through a filename.
-   **Initialization**: Implements a nearest neighbor construction heuristic for generating an initial solution.
-   **Perturbation**: Incorporates the Double Bridge Perturbation method.
