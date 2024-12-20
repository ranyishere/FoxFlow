
using Random

# Function to calculate the propensity for each reaction
function calculate_propensities(rates, state, stoich)
    # res = [stoich[r] for (r, rate) in enumerate(rates)]
    # println("res: ", res)
    println("state: ", state)
    println("rate: ", rates)
    println("stoich: ", stoich)
    res = [rate * state .^ stoich[r] for (r, rate) in enumerate(rates)]
    println(res)
    res
end

# Function to execute the Gillespie algorithm
function gillespie_algorithm(stoich, rates, initial_state, max_time)

    state = copy(initial_state)
    println("rates: ", rates)

    println("state: ", state)

    time = 0.0
    trajectory = [(time, copy(state))]  # Record time and state

    while time < max_time
        propensities = calculate_propensities(rates, state, stoich)
        total_rate = sum(propensities)

        if total_rate == 0
            break  # No more reactions can occur
        end

        # Time to the next reaction
        tau = -log(rand()) / total_rate

        # Determine which reaction occurs
        println("total_rate: ", total_rate)
        println("cumsum(propensities): ", cumsum(propensities))

	println("tau: ", tau)

        # r = searchsortedfirst(cumsum(propensities) .>= rand() * total_rate)
	total_prop = cumsum(propensities)

	tmp_rate = rand() * total_rate

	println("tmp_rate: ", tmp_rate)

	println("total_prop: ", total_prop)

        r = cumsum(propensities) .>= rand() * total_rate

        # Update the state
        state += stoich[r]

        time += tau

        # Record the time and state
        push!(trajectory, (time, copy(state)))
    end
    return trajectory
end

# Example usage
# 2A -> A (Input is 2, net is -1)
# A -> null (Input is 1, net is -1)
# NEED Input stochiometry and Net stochiometry
# Cha
stoich = [  # Stoichiometry matrix
    [-1, 1];   # Reaction 1: A -> B
    [1, -1]    # Reaction 2: B -> A
]

rates = [1.0, 0.5]  # Reaction rates for each reaction
initial_state = [100.0, 0.1]  # Initial state of species A and B
max_time = 10.0  # Maximum simulation time

# Run the Gillespie algorithm
trajectory = gillespie_algorithm(stoich, rates, initial_state, max_time)

# Output the result
for (t, state) in trajectory
    println("Time: $t, State: $state")
end
