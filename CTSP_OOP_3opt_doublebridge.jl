using Random
using Dates  # Ensure the Dates module is available for time calculations

# Read an instance of the Clever Traveling Salesperson Problem
# Input: filename = path + filename of the instance
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

##########################################################################################################

struct CSPSolver
    costs::Array{Int64, 2}
    dimension::Int32
    precedence::Dict{Int32, Set{Int32}}
    startTime::UInt64
    timelimit::UInt32
    startingCustomer::Int32 
    endingCustomer::Int32 
end

# struct representing a solution
mutable struct CSPSolution
    sequence::Array{Int32,1}
    objective::Int64

    # Original constructor
    CSPSolution(dimension::Int32) = new(zeros(Int32, dimension), 0)

    # Additional constructor
    CSPSolution(sequence::Array{Int32,1}, objective::Int64) = new(sequence, objective)
end

##########################################################################################################

# The function iterates over each row of the costs matrix and identifies indices with -1. These indices are added to a set, and the set is then associated with the row index in the precedence dictionary. Only non-empty sets are added to the dictionary, which means only customers with restrictions will have an entry.
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

# allazo logiki me iakovo
function identify_starting_customer(costs::Array{Int64, 2}, dimension::Int32)
    for j in 1:dimension
        # all column entries are -1 (nxn customers) but one and that value is 0 -- no travel needed
        if count(x -> x == -1, costs[:, j]) == dimension - 1 && costs[j, j] == 0
            return j
        end
    end
    return 1 # If no specific starting pattern is found, return a default value.
end

function identify_ending_customer(costs::Array{Int64, 2}, dimension::Int32)
    for i in 1:dimension
        # all row entries are -1 (nxn customers) but one and that value is 0
        if count(x -> x == -1, costs[i, :]) == dimension - 1 && costs[i, i] == 0
            return i
        end
    end
    return dimension # If no specific ending pattern is found, return the last customer as default.
end

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
        
        if nearest_node == 0
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

    return path, total_cost
end
##########################################################################################################
# Part B: Progressing the metaheuristic

# Given the solver, which stores the start time, this function returns the elapsed time in seconds.
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

function generate_3opt_reconnections(sequence, i, j, k, precedence)
    # This function will generate all valid reconnections considering the precedence constraints
    valid_reconnections = []

    # The segments are defined as follows: A = sequence[1:i], B = sequence[i+1:j], C = sequence[j+1:k]
    A = sequence[1:i]
    B = sequence[i+1:j]
    C = sequence[j+1:k]
    D = sequence[k+1:end]

    # Possible reconnections
    possibilities = [
        vcat(A, reverse(B), reverse(C), D),  # 1
        vcat(A, C, reverse(B), D),          # 2
        vcat(A, reverse(C), B, D),          # 3
        vcat(reverse(A), reverse(B), C, D), # 4
        vcat(reverse(A), B, reverse(C), D), # 5
        vcat(reverse(A), C, B, D),          # 6
        vcat(A, B, C, D),                   # 7 - This is the original sequence
    ]

    # Check if the reconnections are valid considering precedence
    for new_sequence in possibilities
        if is_path_valid(CSPSolution(new_sequence, 0), precedence)
            push!(valid_reconnections, new_sequence)
        end
    end
    
    return valid_reconnections
end


function is_path_valid(sol::CSPSolution, precedence::Dict{Int32, Set{Int32}})::Bool
    visited = Set{Int32}()  # Keep track of visited nodes

    # Iterate through the path
    for (index, node) in enumerate(sol.sequence)
        # Check if all predecessors of the current node have been visited
        predecessors = get(precedence, node, Set())
        if !issubset(predecessors, visited)
            return false
        end

        # Check if the current node is a successor of any node that comes after it
        for successor in keys(precedence)
            if node in get(precedence, successor, Set()) && successor in Set(sol.sequence[1:index-1])
                return false
            end
        end

        # Mark the current node as visited
        push!(visited, node)
    end

    return true
end
##########################################################################################################
# Part C: Performing ILS

# The goal here is to systematically explore the solution space using 3-opt moves that respect precedence constraints and to apply improvements as they are discovered
function EnhancedLocalSearch!(m::CSPSolver, sol::CSPSolution)
    improvement_found = true
    while improvement_found
        if elapsed_time(m) >= m.timelimit
            println("Time limit reached during enhanced local search. Terminating early.")
            break  # Exit the while loop if time limit is reached
        end
        improvement_found = false
        best_improvement = 0
        best_sequence = copy(sol.sequence)
        
        # Iterate over all triplets of edges to consider for 3-opt moves
        n = length(sol.sequence)
        for i in 2:n-4
            for j in i+2:n-2
                for k in j+2:n
                    # Time check inside the most nested loop, to ensure responsiveness
                    if elapsed_time(m) >= m.timelimit
                        println("Time Limit Reached.")
                        sol.sequence = best_sequence  # Apply any pending best sequence change
                        sol.objective = calculate_objective(CSPSolution(best_sequence, 0), m.costs)
                        return  # Return early from the function
                    end
                    # Generate and evaluate 3-opt reconnections
                    reconnections = generate_3opt_reconnections(sol.sequence, i, j, k, m.precedence)
                    for reconnection in reconnections
                        new_objective = calculate_objective(CSPSolution(reconnection, 0), m.costs)
                        improvement = sol.objective - new_objective
                        if improvement > best_improvement
                            best_improvement = improvement
                            best_sequence = reconnection
                            improvement_found = true
                        end
                    end
                end
            end
        end

        # Apply the best 3-opt improvement found
        if improvement_found
            sol.sequence = best_sequence
            sol.objective = calculate_objective(CSPSolution(best_sequence, 0), m.costs)
            println("3-opt improvement found with objective: ", sol.objective)
        end
    end
end
# This approach ensures that each iteration of the local search makes the best possible improvement using 3-opt moves while respecting precedence constraints. The search continues until no further improvements can be found.

function doubleBridgePerturbation(sol::CSPSolution)
    path = sol.sequence
    # Ensure there are enough elements for double bridge perturbation
    if length(path) < 5
        println("Path too short for double bridge perturbation.")
        return
    end

    first = path[1]
    last = path[end]
    path_2 = path[2:end-1]

    # Calculate the maximum segment size K
    K = length(path_2) ÷ 4
    if K == 0
        println("Segment too short for meaningful perturbation.")
        return
    end

    # Generate segment offsets A, B, C with the given constraints
    A = rand(1:K)
    B = A + rand(1:K)
    C = B + rand(1:K)

    # Create the new path by rearranging segments
    newPath = vcat(first, path_2[1:A], path_2[C+1:end], path_2[B+1:C], path_2[A+1:B], last)

    # Update the solution's sequence
    sol.sequence = newPath
    # println("Double bridge perturbation applied.")
end


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


function CTSP_ILS(m::CSPSolver, initialSol::CSPSolution)
    s = deepcopy(initialSol)  # Assuming initialSol is generated considering precedence constraints
    s.objective = calculate_objective(s, m.costs)  # Calculate initial objective
    
    println("Starting Enhanced Local Search with initial objective: ", s.objective)

    # Apply enhanced local search (3-opt) to the initial solution
    EnhancedLocalSearch!(m, s)
    
    best_cost = s.objective
    
    while elapsed_time(m) < m.timelimit
        ss = deepcopy(s)  # Deep copy the solution for perturbation
        perturb!(ss, m)  # Apply perturbation

        # Before proceeding with another round of local search, check if time limit has been reached
        if elapsed_time(m) >= m.timelimit
            println("Approached timie limit inside of ILS function")
            break
        end

        EnhancedLocalSearch!(m, ss)  # Apply enhanced local search (3-opt) to the perturbed solution
        
        if ss.objective < best_cost  # Check for improvement -- step funciton 
            s = ss
            best_cost = ss.objective
            println("Improvement found with objective: ", ss.objective)
        end
    end
    
    # Check the time again before starting the next iteration

    println("Final best objective: ", best_cost)  # Confirm the final best objective
    return s
end


##########################################################################################################
# Part D: Output

function reduce_by_one(path::Vector{<:Integer})
    new_path = [i - 1 for i in path]
    return new_path
end


function main()
    # Check if the correct number of arguments are passed
    if length(ARGS) != 3
        println("Usage: julia script.jl <instance_file> <solution_file> <time_limit>")
        return
    end

    # Parse command-line arguments
    instance_file = ARGS[1]
    solution_file = ARGS[2]
    time_limit = parse(Int, ARGS[3])

    # Seed the random number generator for reproducibility
    Random.seed!(1234)

    # Read the instance data
    name, upper_bound, dimension, costs = read_instance(instance_file)

    # Identify starting and ending customers based on the cost matrix
    starting_customer = identify_starting_customer(costs, Int32(dimension))
    ending_customer = identify_ending_customer(costs, Int32(dimension))


    # Now that starting_customer and ending_customer are defined, initialize the solver 
    solver = CSPSolver(costs, Int32(dimension), Dict{Int32, Set{Int32}}(), time_ns(), UInt32(time_limit), Int32(starting_customer), Int32(ending_customer))

    # Build precedence dictionary from the cost matrix
    precedence = build_precedence_dictionary(solver)

    # Update the solver with the precedence 
    solver = CSPSolver(costs, Int32(dimension), precedence, time_ns(), UInt32(time_limit), Int32(starting_customer), Int32(ending_customer))

    # Generate an initial solution respecting precedence constraints
    initialSequence, initialObjective = nearest_neighbor_with_custom_restrictions(costs, Int32(dimension), Int32(starting_customer), Int32(ending_customer), precedence)
    println("Initial path generated: ", initialSequence)
    println("Initial cost : ", initialObjective)
    initialSolution = CSPSolution(initialSequence, initialObjective)

    # Perform Iterated Local Search and get the final solution
    finalSolution = CTSP_ILS(solver, initialSolution)

    # Modify the final solution sequence to start customer numbers from 0
    modified_final_sequence = reduce_by_one(finalSolution.sequence)

    # Write the modified solution sequence to the solution file
    open(solution_file, "w") do file
        write(file, join(modified_final_sequence, " "))
    end
end

# Call main function if this script is executed as a program
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

