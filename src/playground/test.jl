
include("../parser.jl")
include("../tokens.jl")
include("../semantics.jl")
# include("../simulate.jl")
# include("../simulator.jl")
include("../ir_generation.jl")

function main()

    source_code_1 = """
        grammar discretetime symbol1 (Node(x) -> Node(y)) {
            Node(x) -> Node(y) with f(a) where a < 10;
            Node(a) -> Node(b) with g(b) where b < 10;
    };
    """

    source_code_2 = """
    grammar discretetime PredatorPrey () {

        initial_conditions {
            rabbit = 10;
            fox = 5;
        };

        rabbit -> null with alpha where alpha = 0.3;
        fox -> null with beta where beta = 0.3;
        rabbit -> rabbit with gamma where gamma = 1.0;
    };
    """

    # TODO: Add observables.
    #       Infer graph grammar from data.
    #       What's the coarse graining of particles.
    #       Related to variables.
    #       Coarse graining of dynamics is different.
    """
    grammar discretetime Analytic () {

        initial_conditions {
            rabbit = 10;
            fox = 5;
        };

        observables {
            particle_count;
            moments;
        }

        rabbit -> fox with alpha;
        fox -> rabbit with beta;
    };
    """

    # sys.exec(./foxflow -m ".....................................")
    # ./foxflow -r "current_session.flow"


    # // Perhaps allow user to sample from probability distribution
    # // NOTE allow f to be defined outside the grammar
    source_code_3 = """
    grammar discretetime PredatorPrey () {

        types {
            int rabbit ;
            int fox;
            float alpha;
            float beta;
            float gamma;
        };

        initial_conditions {
            rabbit = 10;
            fox = 5;
            f(alpha) = 1;
        };

        rabbit -> null with f(alpha) where alpha < 0.3;
        fox -> null with f(alpha) where beta < 0.3;
        rabbit -> rabbit with f(alpha) where gamma < 1.0;

    };
    """

    # fox -> fox with f(fox, rabbit) where rabbit < 1.0;

    println("source_code_2: ", source_code_2)

    # tokens = tokenize(source_code_1)
    tokens = tokenize(source_code_2)

    # println("tokens: ", tokens)
    # exit(0)

    symbol_table = Dict()
    symbol_table["null"] = nothing
    grammars = parse!(tokens, symbol_table)
    generate_im(grammars, symbol_table)

    # println("done")
    # semantics!(grammars, symbol_table)
    # println("done")

end

main()
