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
    CSPSolution(sequence::Array{Int64,1}, objective::Int64) = new(sequence, objective)
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

# Computes the total cost of a given solution based on the cost matrix
function calculate_objective(sol::CSPSolution, costs::Array{Int64, 2})::Int64
    total_cost = 0
    for i in 1:(length(sol.sequence) - 1)
        total_cost += costs[sol.sequence[i], sol.sequence[i+1]]
    end
    return total_cost
end

function calculate_Objective_Delta(m::CSPSolver, sol::CSPSolution, i::Int, j::Int)::Int64
    before_swap_cost = calculate_objective(sol, m.costs)
    swap!(m, sol, i, j)  # Apply the swap temporarily

    after_swap_cost = calculate_objective(sol, m.costs)
    swap!(m, sol, i, j)  # Revert the swap
    # println("Swap between $i and $j, before: $before_swap_cost, after: $after_swap_cost, delta: ", after_swap_cost - before_swap_cost)
    return after_swap_cost - before_swap_cost
end

# Given the solver, which stores the start time, this function returns the elapsed time in seconds.
function elapsed_time(m::CSPSolver)
    return round((time_ns()-m.startTime)/1e9,digits=3)
end

# swaps two elements in the solution's sequence
function swap!(m::CSPSolver, sol::CSPSolution, i::Int, j::Int)
    sol.sequence[i], sol.sequence[j] = sol.sequence[j], sol.sequence[i]
end

# Evaluates if a swap between two positions in the sequence respects the problem's precedence constraints.
# t creates a temporary copy of the solution sequence and reverses the subsequence between indices i and j. Then, for each node within this swapped subsequence, it checks whether all required precedents (as specified in the m.precedence dictionary) appear before the node in the temporary sequence. If any required precedent is missing, the function returns false, indicating the swap would violate precedence constraints.
function is_valid_swap(m::CSPSolver, sol::CSPSolution, i::Int, j::Int)::Bool
    # Create a temporary sequence with the proposed swap
    temp_sequence = copy(sol.sequence)
    temp_sequence[i:j] = reverse(temp_sequence[i:j])

    # Check the entire sequence for precedence constraint violations
    for idx in 1:length(temp_sequence)
        node = temp_sequence[idx]
        for req in get(m.precedence, node, Set())
            if !(req in temp_sequence[1:idx-1])
                return false  # Precedence violation found
            end
        end
    end
    return true  # No violations found
end


##########################################################################################################
# Part C: Performing ILS

# Implements a local search optimization that iteratively looks for the best swap operation that improves the solution while respecting precedence constraints
function LocalSearch!(m::CSPSolver, sol::CSPSolution)
    improved = true
    while improved
        improved = false
        best_delta = 0
        best_i, best_j = 0, 0
        for i in 2:(length(sol.sequence) - 1)
            for j in (i + 1):length(sol.sequence)
                if is_valid_swap(m, sol, i, j)
                    delta = calculate_Objective_Delta(m, sol, i, j)
                    # println("Delta for swap between $i and $j: $delta")  # Added print
                    if delta < best_delta
                        best_delta = delta
                        best_i = i
                        best_j = j
                        improved = true
                    end
                end
            end
        end
        # calculates the best swap
        if improved
            # Perform the best swap found and update the objective accordingly
            swap!(m, sol, best_i, best_j)
            sol.objective += best_delta
            # println("Performed swap between $best_i and $best_j, new objective: ", sol.objective)
        end
    end
end


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
    if !is_path_valid(sol, m)
        # If the perturbed solution violates precedence constraints, revert to the original sequence
        sol.sequence = original_sequence
        # println("Perturbation reverted due to precedence violation.")
    else
        # Recalculate the solution's objective if the perturbation is valid
        sol.objective = calculate_objective(sol, m.costs)

    end
end

function is_path_valid(sol::CSPSolution, m::CSPSolver)::Bool
    visited = Set{Int32}()

    for node in sol.sequence
        if haskey(m.precedence, node)
            for prereq in m.precedence[node]
                if !(prereq in visited)
                    return false
                end
            end
        end
        push!(visited, node)
    end

    return true
end



function CTSP_ILS(m::CSPSolver, initialSol::CSPSolution)
    start_time = Dates.now()
    s = deepcopy(initialSol)  # Assuming initialSol is generated considering precedence constraints
    s.objective = calculate_objective(s, m.costs)  # Calculate initial objective
    
    println("Starting Local Search with initial objective: ", s.objective) # Add this line to confirm entry

    # Apply local search to the initial solution
    LocalSearch!(m, s)
    
    noImprovementIterations = 0
    best_cost = s.objective
    
    while (Dates.now() - start_time) < Dates.Second(m.timelimit)
        ss = deepcopy(s)  # Deep copy the solution for perturbation
        perturb!(ss, m)  # Apply perturbation
        LocalSearch!(m, ss)  # Apply local search to the perturbed solution
        
        if ss.objective < best_cost  # Check for improvement Best swap selection
            s = ss
            best_cost = ss.objective
            noImprovementIterations = 0
            println("Improvement found with objective: ", ss.objective)
        else
            noImprovementIterations += 1
        end
    end
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
    solver = CSPSolver(costs, Int32(dimension), Dict{Int32, Set{Int32}}(), UInt64(0), UInt32(time_limit), Int32(starting_customer), Int32(ending_customer))


    # Build precedence dictionary from the cost matrix
    precedence = build_precedence_dictionary(solver)

    # Initialize the solver with the problem data and the time limit
    # This line is redundant since the solver is already initialized above with the same parameters, so it can be removed or commented out.
    solver = CSPSolver(costs, Int32(dimension), precedence, UInt64(0), UInt32(time_limit), Int32(starting_customer), Int32(ending_customer))

    # Generate an initial solution respecting precedence constraints
    initialSequence, initialObjective = nearest_neighbor_with_custom_restrictions(costs, dimension, starting_customer, ending_customer, precedence)
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

# Call main function if this script is executed as a program
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


# D:
# @ cd D:\Aλέξανδρος\Code_Warehouse\Julia\Metaheuristics\Project_1

#  julia s222717.jl instances/ESC07.sop sols/solution_ESC07.sol 60
# CTSPChecker.exe Instances/ESC07.sop sols/solution_ESC07.sol
# Objective = 2600 -- gap 22.3% --> 3OPT = 2150  -- gap = 1.1%

# julia s222717.jl instances/br17.10.sop sols/solution_br17.10.sol 60
# CTSPChecker.exe Instances/br17.10.sop sols/solution_br17.10.sol
# obj = 63 -- 32% ------ 3OPT = 58 -- gap 5.4%

# julia s222717.jl instances/br17.12.sop sols/solution_br17.12.sol 60
# CTSPChecker.exe Instances/br17.12.sop sols/solution_br17.12.sol
# Objective = 63 -- gap 129% two opt --> 3opt = 58 -- gap 5.4%

# julia s222717.jl instances/ESC12.sop sols/solution_ESC12.sol 60
# CTSPChecker.exe Instances/ESC12.sop sols/solution_ESC12.sol
# Objective =  1752 -- gap 0%  --- 3opt= 1793 -- gap  7.0%



# julia s222717.jl instances/ESC63.sop sols/solution_ESC63.sol 60
# CTSPChecker.exe Instances/ESC63.sop sols/solution_ESC63.sol
# Objective =  72 -- gap 33,8% ------ 3OPT =  75  -- gap <22.5%

# julia s222717.jl instances/ESC78.sop sols/solution_ESC78.sol 60
# CTSPChecker.exe Instances/ESC78.sop sols/solution_ESC78.sol
# Objective = 22600 -- gap 43.8% -------- 3opt = 22300 -- gap = <23.9%  

# julia s222717.jl instances/ft53.1.sop sols/solution_ft53.1.sol 60
# CTSPChecker.exe Instances/ft53.1.sop sols/solution_ft53.1.sol
# Objective =  9609 -- gap 53.3% ---- 3opt = 9883  -- gap 31.2%



# julia s222717.jl instances/ft53.2.sop sols/solution_ft53.2.sol 60
# CTSPChecker.exe Instances/ft53.2.sop sols/solution_ft53.2.sol
# Objective =  12298 -- gap 71.3% -- 3opt = 11994 -- gap <56.3

# julia s222717.jl instances/ESC25.sop sols/solution_ESC25.sol 60
# CTSPChecker.exe Instances/ESC25.sop sols/solution_ESC25.sol
# Objective =  3075 -- gap 138% -- 2962 WITH 3opt -- gap 76.2%

# julia s222717.jl instances/ESC47.sop sols/solution_ESC47.sol 60
# CTSPChecker.exe Instances/ESC47.sop sols/solution_ESC47.sol
# Objective = 3489  -- gap 256%% ----- 3opt = 2993 -- gap 193%
