# <42137 Optimization using Metaheuristics -- Assignment 01>
# ILS for the Clever Traveling Salesperson Problem
#
# Key aspects of this algorithm include:
#   Neighborhood Operator: The 3-opt move is employed as the neighborhood operator.
#   Termination Criteria: The algorithm uses a time limit as the primary termination criterion.
#   Step Criterion: Adopts the "Best Improvement" strategy is utilized.
#   Cost: The algorithm operates based on a cost matrix, which is specified through a filename.
#   Initialization: Implements a nearest neighbor construction heuristic for generating an initial solution.
#   Perturbation: Incorporates the  Double Bridge Perturbation method.
#**************************************************************************************************************************************************
using Random

# Read an instance of the Clever Traveling Salesperson Problem --> Input: filename = path + filename of the instance
function read_instance(filename)
    # Opens the file
    f = open(filename)
    # Reads the name of the instance
    name = split(readline(f))[2]
    # reads the upper bound value
    upper_bound = parse(Int64,split(readline(f))[2])
    readline(f)#type
    readline(f)#comment
    # reads the dimentions of the problem
    dimention = parse(Int64,split(readline(f))[2]);
    readline(f)#Edge1
    readline(f)#Edge2
    readline(f)#Edge3
    readline(f)#Dimention 2
    # Initialises the cost matrix
    cost = zeros(Int64,dimention,dimention)
    # Reads the cost matrix
    for i in 1:dimention
        data = parse.(Int64,split(readline(f)))
        cost[i,:]=data
    end
    # Closes the file
    close(f)
    # Returns the input data
    return name, upper_bound, dimention, cost
end
#**************************************************************************************************************************************************


#**************************************************************************************************************************************************
# CSPSolver Struct Definition
# This struct defines the main components required for solving the Clever Traveling Salesperson Problem (CTSP) using Iterated Local Search (ILS).
# It encapsulates the cost matrix, problem dimensions, precedence constraints, timing information, and identifiers for the starting and ending customers.
struct CSPSolver
    costs::Array{Int64, 2}
    dimension::Int32
    precedence::Dict{Int32, Set{Int32}}
    startTime::UInt64
    timelimit::UInt32
    startingCustomer::Int32 
    endingCustomer::Int32 
end

# CSPSolution Struct Definition
# This mutable struct represents a potential solution to the CTSP.
# It holds the sequence in which customers are visited and the objective value of this sequence.
mutable struct CSPSolution
    sequence::Array{Int32,1}
    objective::Int64
    # Constructs a new CSPSolution with a sequence of given dimension filled with zeros and an initial objective value of 0.
    CSPSolution(dimension::Int32) = new(zeros(Int32, dimension), 0)
    # Constructs a new CSPSolution with a specified sequence and its associated objective value.
    CSPSolution(sequence::Array{Int32,1}, objective::Int64) = new(sequence, objective)
end
#**************************************************************************************************************************************************


#**************************************************************************************************************************************************
# This function constructs a dictionary that maps each customer to a set of customers that must be visited before it, based on the presence of -1 in the cost matrix.
function build_precedence_dictionary(solver::CSPSolver)
    precedence = Dict{Int32, Set{Int32}}()
    for i in 1:solver.dimension
        # Find indices with -1 and create a set for those indices
        restricted_indices = findall(x -> x == -1, solver.costs[i, :])
        # Only add to the dictionary if there are restrictions
        if !isempty(restricted_indices)
            precedence[i] = Set(restricted_indices)
        end
    end
    return precedence
end

# This function determine the starting for the tour,based on the special conditions in the cost matrix. 
function identify_starting_customer(costs::Array{Int64, 2}, dimension::Int32)
    for j in 1:dimension
        # all column entries are -1 (nxn customers) but one and that value is 0 
        if count(x -> x == -1, costs[:, j]) == dimension - 1 && costs[j, j] == 0
            return j
        end
    end
    return 1 # If no specific starting pattern is found, return a default value.
end

# This function determine the ending for the tour,based on the special conditions in the cost matrix. 
function identify_ending_customer(costs::Array{Int64, 2}, dimension::Int32)
    for i in 1:dimension
        # all row entries are -1 (nxn customers) but one and that value is 0
        if count(x -> x == -1, costs[i, :]) == dimension - 1 && costs[i, i] == 0
            return i
        end
    end
    return dimension # If no specific ending pattern is found, return the last customer as default.
end

# This function implements a contruction heuristic (nearest neighbor heuristic), modified to respect the precedence constraints defined in your problem. 
# By ensuring that only accessible and allowed customers are considered for the next step in the path, this function generates a feasible initial solution.
function nearest_neighbor_with_custom_restrictions(costs, dimension, starting_customer, ending_customer, precedence)
    n = dimension  # Number of nodes
    visited = Bool[false for _ in 1:n]  # Track visited nodes
    path = [starting_customer]  # Starting from starting_customer
    visited[starting_customer] = true  # Mark the starting node as visited

    function can_visit(node, visited, precedence)
        # Check if all precedence for the node are satisfied
        for req in get(precedence, node, Set())
            if !visited[req]
                return false
            end    
        end    
        return true
    end    
    
    while sum(visited) < n
        last_node = path[end]
        nearest_node = 0
        nearest_distance = Inf  # A large number
        for j in 1:n
            # Check if the node is unvisited, reachable, and satisfies precedence
            if !visited[j] && costs[last_node, j] != -1 && costs[last_node, j] < nearest_distance && can_visit(j, visited, precedence)
                nearest_distance = costs[last_node, j]
                nearest_node = j
            end    
        end    
        if nearest_node == 0 # Impossible to proceed without violating constraints.
            println("No further nodes can be visited due to restrictions.")
            break
        else
            push!(path, nearest_node)  # Add the nearest node to the path
            visited[nearest_node] = true  # Mark the nearest node as visited
        end    
    end    
    # Add the ending customer if it's not already in the sequence and if specified
    if ending_customer !== nothing && !(ending_customer in path)
        push!(path, ending_customer)
    end    
    # Calculate the total cost of the sequence
    total_cost = 0
    for i in 1:(length(path) - 1)
        total_cost += costs[path[i], path[i + 1]]
    end

    # Returns a tuple (path, total_cost) where `path` is the sequence of visited customers and `total_cost` is the path's total cost.
    return path, total_cost
end
#**************************************************************************************************************************************************


#**************************************************************************************************************************************************
# This function returns the elapsed time in seconds, given the solve which stores the start time.
function elapsed_time(m::CSPSolver)
    return round((time_ns()-m.startTime)/1e9,digits=3)
end

# Computes the total cost of a given solution based on the cost matrix
function calculate_objective(sol::CSPSolution, costs::Array{Int64, 2})::Int64
    total_cost = 0
    for i in 1:(length(sol.sequence) - 1)
        total_cost += costs[sol.sequence[i], sol.sequence[i+1]]
    end
    return total_cost
end

# Computes the change in objective for a 3-opt move
function calculate_3opt_Objective_Delta(m::CSPSolver, sol::CSPSolution, i::Int, j::Int, k::Int)::Int64
    before_move_cost = calculate_objective(sol, m.costs)
    perform_3opt_move!(sol, i, j, k)  # Apply the 3-opt move temporarily
    
    after_move_cost = calculate_objective(sol, m.costs)
    perform_3opt_move!(sol, i, j, k)  # Revert the 3-opt move
    
    return after_move_cost - before_move_cost
end
#**************************************************************************************************************************************************


#**************************************************************************************************************************************************
# This function checks if a given solution's path is valid considering precedence constraints.
# Inputs:
#   - sol: A CSPSolution instance containing a sequence of nodes.
#   - precedence: Dict mapping nodes to their respective sets of nodes that must precede them.
# Output:
#   - Boolean value indicating whether the path is valid (true) or not (false).
function is_path_valid(sol::CSPSolution, precedence::Dict{Int32, Set{Int32}})::Bool
    visited = Set{Int32}()  # Tracks nodes visited in the sequence.
    
    # Iterate through each node in the sequence.
    for (index, node) in enumerate(sol.sequence)
        # Retrieve predecessors for the current node and check if all have been visited.
        predecessors = get(precedence, node, Set())
        if !issubset(predecessors, visited)
            return false # Return false if any predecessor hasn't been visited.
        end

        # Additional check for successors listed in the precedence dict to ensure no forward references.
        for successor in keys(precedence)
            if node in get(precedence, successor, Set()) && successor in Set(sol.sequence[1:index-1])
                return false
            end
        end
        
        # Mark the current node as visited
        push!(visited, node)
    end
    
    return true # Return true if all checks pass, indicating a valid path.
end

# Generates all 3-opt reconnections of a given sequence that respect precedence constraints.
# Inputs:
#   - sequence: Current sequence of nodes.
#   - i, j, k: Indices defining the boundaries of three segments to be reconnected.
#   - precedence: Dict mapping nodes to their respective sets of nodes that must precede them.
# Output:
#   - valid_reconnections: List of valid sequences generated by applying 3-opt moves to the original sequence.
function generate_3opt_reconnections(sequence, i, j, k, precedence)
    valid_reconnections = []

    # Define segments A, B, C, and D based on indices i, j, k.    
    A = sequence[1:i]
    B = sequence[i+1:j]
    C = sequence[j+1:k]
    D = sequence[k+1:end]

    # Define all possible reconnections using 3-opt moves.
    possibilities = [
        vcat(A, reverse(B), reverse(C), D),  # Reverses segments B and C.
        vcat(A, C, reverse(B), D),           # Places segment C before B and reverses B.
        vcat(A, reverse(C), B, D),           # Reverses segment C and places it before B.
        vcat(reverse(A), reverse(B), C, D),  # Reverses segments A and B, maintaining C and D.
        vcat(reverse(A), B, reverse(C), D),  # Reverses segments A and C, maintaining B and D.
        vcat(reverse(A), C, B, D),           # Reverses segment A and swaps B and C.
        vcat(A, B, C, D),                    # Original sequence, included for completeness.
    ]

    # Check each new sequence for validity considering the precedence constraints.
    for new_sequence in possibilities
        if is_path_valid(CSPSolution(new_sequence, 0), precedence)
            push!(valid_reconnections, new_sequence)
        end
    end
    
    return valid_reconnections
end
#**************************************************************************************************************************************************


#**************************************************************************************************************************************************
# This function implements a "Double-Bridge Perturbation" to the given solution's sequence.
# A Double-Bridge Perturbation breaks the path into four segments and then reconnects them in a new order.
# This method is beneficial for escaping local optima by significantly altering the solution's structure.

function doubleBridgePerturbation(sol::CSPSolution)
    path = sol.sequence
    # Ensure there are enough elements for double bridge perturbation
    if length(path) < 5
        println("Path too short for double bridge perturbation.")
        return
    end

    first = path[1] #fix start 
    last = path[end] # fix end
    reduced_path = path[2:end-1]

    # Calculate the maximum segment size K
    K = length(reduced_path) ÷ 4
    if K == 0
        println("Segment too short for meaningful perturbation.")
        return
    end
    # Generate segment offsets A, B, C with the given constraints
    A = rand(1:K)
    B = A + rand(1:K)
    C = B + rand(1:K)

    # Create the new path by rearranging segments
    newPath = vcat(first, reduced_path[1:A], reduced_path[C+1:end], reduced_path[B+1:C], reduced_path[A+1:B], last)
    sol.sequence = newPath # Update the solution's sequence
end

# This function applies the doubleBridgePerturbation to the solution's sequence and checks for precedence constraint violations.
function perturb!(sol::CSPSolution, m::CSPSolver)
    original_sequence = copy(sol.sequence)  # Keep the original sequence for potential rollback
    # Apply the double bridge perturbation on the solution's sequence directly
    doubleBridgePerturbation(sol)
    # Check if the new sequence respects the precedence constraints
    if !is_path_valid(sol, m.precedence)
        # If the perturbed solution violates precedence constraints, revert to the original sequence
        sol.sequence = original_sequence
        # println("Perturbation reverted due to precedence violation.")
    else
        # Recalculate the solution's objective if the perturbation is valid
        sol.objective = calculate_objective(sol, m.costs)
    end
end
#**************************************************************************************************************************************************


#**************************************************************************************************************************************************
# EnhancedLocalSearch! performs an iterative improvement process on the given solution using 3-opt moves (Neighborhood Operator)
#   It aims to find a better solution by systematically exploring possible reconnections that respect precedence constraints.
#   The search continues until no further improvement can be found or the time limit is reached.
# Inputs:
#   - m: An instance of CSPSolver containing problem details and algorithm settings.
#   - sol: A CSPSolution object representing the current solution to be improved.
# The function updates `sol` in place if a better solution is found. It also ensures that the search does not exceed the specified time limit.
function EnhancedLocalSearch!(m::CSPSolver, sol::CSPSolution)
    improvement_found = true
    while improvement_found
        # Check for time limit at the start of each iteration to ensure the search does not overrun.
        if elapsed_time(m) >= m.timelimit
            println("Time limit reached during enhanced local search. Terminating early.")
            break  # Terminate the search if the time limit is reached.
        end
        improvement_found = false # Reset the flag for detecting improvements in this iteration.
        best_improvement = 0 # Track the best improvement found in this iteration.
        best_sequence = copy(sol.sequence) # Copy the current solution sequence for modification.
        
        # Iterate over all possible triplets of edges in the sequence to apply 3-opt moves.
        n = length(sol.sequence)
        for i in 2:n-4
            for j in i+2:n-2
                for k in j+2:n
                    # Perform a time check within the deepest loop to maintain responsiveness to the time limit.
                    if elapsed_time(m) >= m.timelimit
                        println("Time Limit Reached.")
                        sol.sequence = best_sequence  # Apply the best improvement found before terminating.
                        sol.objective = calculate_objective(CSPSolution(best_sequence, 0), m.costs)
                        return  # Exit the function forcefully due to time constraint.
                    end

                    # Generate all valid 3-opt reconnections for the current triplet and evaluate their objective values.
                    reconnections = generate_3opt_reconnections(sol.sequence, i, j, k, m.precedence)
                    for reconnection in reconnections
                        new_objective = calculate_objective(CSPSolution(reconnection, 0), m.costs)
                        # if new objective > objective --> the difference is + --> indicating a reconnection offers a better (lower cost) solution.
                        improvement = sol.objective - new_objective 
                        # Checks if the current reconnection offers a more significant reduction in the objective value than any previously considered reconnections in the same iteration
                        if improvement > best_improvement
                            best_improvement = improvement
                            best_sequence = reconnection
                            improvement_found = true # Flag that an improvement was found in this iteration.
                        end
                    end

                end
            end
        end

        # After exploring all reconnections, apply the best one found to the solution.
        if improvement_found
            sol.sequence = best_sequence # Update the solution with the improved sequence.
            sol.objective = calculate_objective(CSPSolution(best_sequence, 0), m.costs) # Recalculate the solution's objective.
            println("3-opt improvement found with objective: ", sol.objective)
        end
    end
end

# CTSP_ILS implements the Iterated Local Search (ILS) algorithm for the Clever Traveling Salesperson Problem (CTSP).
function CTSP_ILS(m::CSPSolver, initialSol::CSPSolution)
    s = deepcopy(initialSol)  # Create a deep copy of the initial solution to avoid altering the original.
    s.objective = calculate_objective(s, m.costs)  # Recalculate the initial objective for accuracy.
    println("Starting Enhanced Local Search with initial objective: ", s.objective)

    # First phase of ILS: Apply an enhanced local search (3-opt) to improve the initial solution.
    EnhancedLocalSearch!(m, s)
    best_cost = s.objective # Initialize the best known cost as the objective of the improved initial solution.
    
    # Main ILS loop: Continues until the time limit is reached.
    while elapsed_time(m) < m.timelimit
        ss = deepcopy(s)  # Create a deep copy of the current best solution for perturbation.
        perturb!(ss, m)   # Apply a perturbation to the copied solution to explore new parts of the solution space.

        # Check the elapsed time to ensure the algorithm does not exceed the time limit.
        if elapsed_time(m) >= m.timelimit
            println("Approached timie limit inside of ILS function.")
            break  # Exit the loop if the time limit has been reached.
        end

        # Second phase of ILS: Apply the enhanced local search to the perturbed solution.
        EnhancedLocalSearch!(m, ss) 
        # If the perturbed and locally searched solution is better, update the current best solution.
        if ss.objective < best_cost  # Check for improvement -- step funciton 
            s = ss # Update the best solution.
            best_cost = ss.objective # Update the best known cost.
            println("ILS improvement found with objective: ", ss.objective)
        end
        # The loop will continue to alternate between perturbation and local search, seeking further improvements.
    end
    
    # Check the time again before starting the next iteration

    println("Final best objective: ", best_cost)  # Announce the best objective value found within the time limit.
    return s # Return the best solution found.
end
#**************************************************************************************************************************************************


#**************************************************************************************************************************************************
# This function reduces the indices of the optimal path by 1 in accordance to the Solution Checker
function reduce_by_one(path::Vector{<:Integer})
    new_path = [i - 1 for i in path]
    return new_path
end
#**************************************************************************************************************************************************


#**************************************************************************************************************************************************
# Entry point of the CTSP solver program.
# The program requires three command-line arguments to run: 
#    -- the path to the instance file
#    -- the path where the solution should be saved 
#    -- the time limit for the algorithm in seconds.

function main()
    # Verify the correct number of command-line arguments.
    if length(ARGS) != 3
        println("Usage: julia script.jl <instance_file> <solution_file> <time_limit>")
        return
    end
    # Extract command-line arguments.
    instance_file = ARGS[1]
    solution_file = ARGS[2]
    time_limit = parse(Int, ARGS[3])
    # Initialize the random number generator for consistent, reproducible results.
    Random.seed!(1234)
    # Load problem instance data.
    name, upper_bound, dimension, costs = read_instance(instance_file)

    # Determine start and end points for the tour, if specified within the problem instance.
    starting_customer = identify_starting_customer(costs, Int32(dimension))
    ending_customer = identify_ending_customer(costs, Int32(dimension))

    # Initialize the solver with problem data and configuration.
    solver = CSPSolver(costs, Int32(dimension), Dict{Int32, Set{Int32}}(), time_ns(), UInt32(time_limit), Int32(starting_customer), Int32(ending_customer))

    # Establish precedence constraints based on the problem instance.
    precedence = build_precedence_dictionary(solver)

    # Update the solver with the precedence 
    solver = CSPSolver(costs, Int32(dimension), precedence, time_ns(), UInt32(time_limit), Int32(starting_customer), Int32(ending_customer))

    # Generate an initial solution that respects any precedence constraints.
    initialSequence, initialObjective = nearest_neighbor_with_custom_restrictions(costs, Int32(dimension), Int32(starting_customer), Int32(ending_customer), precedence)
    println("Initial path generated: ", initialSequence)
    println("Initial cost : ", initialObjective)
    initialSolution = CSPSolution(initialSequence, initialObjective)

    # Execute the Iterated Local Search algorithm to improve upon the initial solution.
    finalSolution = CTSP_ILS(solver, initialSolution)

    # Adjust the final solution to match the expected output format.
    modified_final_sequence = reduce_by_one(finalSolution.sequence)

    # Save the final solution sequence to the specified output file.
    open(solution_file, "w") do file
        write(file, join(modified_final_sequence, " "))
    end
end
#**************************************************************************************************************************************************


#**************************************************************************************************************************************************
# Call main function if this script is executed as a program
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
#**************************************************************************************************************************************************
