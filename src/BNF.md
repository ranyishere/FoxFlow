# FoxFlow BNF Grammar

This document describes the grammar of the FoxFlow DSL as implemented by the
parser (`parser/main_parser_op.jl`) and AST (`parser/ast_nodes.jl`).

---

## Terminals & Lexical Elements

```
<symbol-name>       ::= <letter> { <letter> | <digit> | "_" }
<letter>            ::= [A-Z] | [a-z]
<digit>             ::= [0-9]
<integer>           ::= <digit> { <digit> }
<float>             ::= <digit> { <digit> } "." <digit> { <digit> }
<number>            ::= <integer> | <float>
<string-literal>    ::= '"' { <any-char> } '"'
<literal>           ::= <integer> | <float> | <symbol-name>
```

---

## Top-Level File Structure

A FoxFlow file consists of one or more sections in any order:

```
<file> ::= { <section> }

<section> ::= <type-section>
            | <parameters-section>
            | <functions-section>
            | <rule-section>
            | <simulations-section>
            | <observables-section>
```

---

## Symbol Parameters

Symbol parameters appear after `<<` and before `>>`, used to annotate types,
rule nodes, and rule graph patterns with typed attribute bindings.

```
<symbol-parameters> ::= "<<" <symbol-parameter-list> ">>"
                       | ε

<symbol-parameter-list> ::= <symbol-parameter> { "," <symbol-parameter> }

<symbol-parameter> ::= <named-parameter>
                     | <expression>

<named-parameter> ::= <symbol-name> ":" <type-signature>
                     | <symbol-name> ":" <vertex-attribute-access>

<vertex-attribute-access> ::= <symbol-name> "::" <symbol-name>
```

A `<vertex-attribute-access>` binds a user-chosen name to a specific
attribute of an already-declared graph node, using `<node-name> "::"
<attribute-name>`. This form is only meaningful in rule parameter blocks,
where the referenced node has been introduced in the rule graph pattern
(e.g. `rho1 : c1::rho` binds `rho1` to the `rho` attribute of node `c1`).

---

## Type Signatures

```
<type-signature> ::= <type-name> <symbol-parameters>
                   | <type-name> "->" <type-signature>
                   | <type-name>

<type-name> ::= <symbol-name> | <built-in-type>

<built-in-type> ::= "Float" | "Integer" | "FixedList" | "TFTensor" | "ODE" | "Type"
```

Two list-like types exist:

- `FixedList<<N, Float>>` is a fixed-size, stack-allocated 1-D array of `N`
  doubles (`TFArray<N>`, defined in the generated `tf_array.h`). It supports
  element-wise arithmetic with scalars and arrays, `exp/log/sqrt/abs/sin/cos/tanh`,
  `dot`, `sum/mean/norm/min/max`, `ones_matrix(N)` / `zeros_matrix(N)` with a
  literal `N`, and array literals. It is much cheaper than a tensor for small
  vectors. It must be 1-D and `Float`; ODE solving, `einsum`, `mat_inv`, nested
  literals and dynamic shapes are not available.
- `TFTensor<<d1, ..., Float|Integer>>` is a `torch::Tensor` with the full set of
  tensor operations. Conversion between the two is explicit only:
  `to_tensor(array)` and `to_array(tensor, N)`.

A local declared as `x : FixedList := ...` (no size) takes the type of its
value; dimension-less `FixedList` is not allowed in function signatures.

---

## Types Section

Declares the named types available in a namespace. Each type has a name,
optional generic parameters, and a type signature.

```
<type-section> ::= "types" <symbol-name> "{" <type-declaration-list> "}"

<type-declaration-list> ::= { <type-declaration> "\n" }

<type-declaration> ::= <symbol-name> <symbol-parameters> ":" <type-signature>
```

### Examples

```
types Microtubule {
    Intermediate << Position : FixedList<<3, Float>>,
                    Direction : FixedList<<3, Float>> >> : Type
    Positive     << Position : FixedList<<3, Float>>,
                    Direction : FixedList<<3, Float>> >> : Type
    Nucleator    << Position : FixedList<<3, Float>>,
                    Count : TFTensor<<2, Integer>> >> : Type
}
```

---

## Parameters Section

Declares simulation parameters as typed, named values.

```
<parameters-section> ::= "parameters" <symbol-name> "{" <parameter-list> "}"

<parameter-list> ::= { <parameter> "\n" }

<parameter> ::= <symbol-name> <symbol-parameters> ":" <type-signature> ":=" <parameter-value>

<parameter-value> ::= "{" <type-declaration-list> "}"
                    | <from-file>
                    | <expression>

<from-file> ::= "from_file" "(" <string-literal> ")"
```

A `TFTensor` parameter with `from_file` loads a raw binary file (saved via
`numpy().tofile()`) into a `torch::Tensor` at startup using memory-mapping:

```
images : TFTensor<<60000, 1, 28, 28, Float>>   := from_file("mnist_images.bin")
labels : TFTensor<<60000, Integer>>             := from_file("mnist_labels.bin")
```

The last type parameter of the `TFTensor` (`Float` or `Integer`) determines
the element dtype (`torch::kFloat64` / `torch::kInt64`). All preceding
integer parameters are the tensor dimensions.

### Examples

```
parameters Microtubule {
    creation_rate : Float := 0.5
    offset : Float := 10.0
    boundary_pts : Integer := 10
}

parameters NeuralNetwork {
    images : TFTensor<<60000, 1, 28, 28, Float>> := from_file("mnist_images.bin")
    labels : TFTensor<<60000, Integer>>           := from_file("mnist_labels.bin")
}
```

---

## Expressions
Operator precedence (lowest to highest):

1. Assignment: `=`, `:=`
2. Logical OR: `||`
3. Logical AND: `&&`
4. Equality: `==`, `!=`
5. Relational: `<`, `<=`, `>`, `>=`
6. Additive: `+`, `-`
7. Multiplicative: `*`, `/`
8. Exponential: `^` (right-associative)
9. Unary: `-`, `+`, `!`, `~`
10. Postfix: index access `[...]`, function call `(...)`

```
<expression> ::= <assignment>

<assignment> ::= <logical-or> [ ("=" | ":=") <assignment> ]

<logical-or> ::= <logical-and> { "||" <logical-and> }

<logical-and> ::= <equality> { "&&" <equality> }

<equality> ::= <relational> { ("==" | "!=") <relational> }

<relational> ::= <additive> { ("<" | "<=" | ">" | ">=") <additive> }

<additive> ::= <multiplicative> { ("+" | "-") <multiplicative> }

<multiplicative> ::= <exponential> { ("*" | "/") <exponential> }

<exponential> ::= <factor> { "^" <exponential> }

<factor> ::= <integer>
           | <float>
           | <identifier-or-call>
           | <unary-expression>
           | <grouped-expression>
           | <array-literal>

<identifier-or-call> ::= [ <symbol-name> "::" ] <symbol-name>
                          [ "(" <arg-list> ")" ]
                          { "[" <index-list> "]" }

<unary-expression> ::= "-" <factor>
                     | "+" <factor>
                     | "!" <factor>
                     | "~" <factor>

<grouped-expression> ::= "(" <expression> ")"

<array-literal> ::= "[" "]"
                  | "[" <expression> { "," <expression> } "]"

<arg-list> ::= <expression> { "," <expression> }

<index-list> ::= <index-element> { "," <index-element> }

<index-element> ::= <expression>
                  | <slice>

<slice> ::= [ <expression> ] ":" [ <expression> ] [ ":" <expression> ]
```

---

## Rules Section

```
<rule-section> ::= "rules" <symbol-name> "{" <rule-list> "}"

<rule-list> ::= { <rule> }
```

### Rule Structure

A rule has a name, a LHS graph pattern, an optional LHS parameter block,
a RHS graph pattern, an optional RHS parameter block, and a modify clause.

```
<rule> ::= <symbol-name> ":="
           <rule-graph-pattern> [ <symbol-parameters> ]
           "->"
           <rule-graph-pattern> [ <symbol-parameters> ]
           <modify-clause>
```

### Rule Graph Patterns (LHS and RHS)

The LHS and RHS of a rule are **graph patterns** — they describe the
structure of the sub-graph to match (LHS) or produce (RHS). Each element is
either a standalone typed node or an undirected edge connecting two typed
nodes. Parentheses around each node are syntactic delimiters consumed by the
parser.

**Important:** These are *not* assignments. There are no curly braces or
`:=` here — just `(name : Type)` declarations, optionally connected with
`--` edges.

```
<rule-graph-pattern> ::= <graph-element> { <graph-element> }

<graph-element> ::= <type-instance>
                  | <type-instance> "--" <type-instance>

```

Edges can chain: when the parser sees `(a : T) -- (b : T)` followed by
another `-- (c : T)`, the rightmost node of the previous edge (`b`) becomes
the left node of the new edge, producing two `UndirectedTypeEdgeNode`s
sharing node `b`.

### Rule Parameter Blocks

The `<< ... >>` blocks that follow the LHS and RHS graph patterns bind
user-chosen names to the attributes of each graph node, grouped by node in
the order the nodes appear in the pattern:

```
<< (attr1 : Type1, attr2 : Type2),   -- bindings for graph node 1
   (attr3 : Type3) >>                 -- bindings for graph node 2
```

Each named parameter here is a **binding** that gives the user a name to
refer to that node's attribute in the modify clause (where body or solving
body).

A binding may take one of two forms:

- **Typed binding** — `<name> ":" <type-signature>` declares a fresh name
  with an explicit attribute type (e.g. `im_pos : FixedList<<3, Float>>`).
- **Vertex attribute access** — `<name> ":" <node-name> "::" <attribute-name>`
  binds a name to an existing attribute of a node introduced in the rule
  graph pattern (e.g. `rho1 : c1::rho` binds `rho1` to the `rho` attribute of
  node `c1`). This lets a rule reference attributes from specific vertices by
  the node names used in the graph pattern.

### Examples

```
# One LHS node, two RHS nodes (no edges)
start_to_node := (start : StartType) << (start_pos : start::Position) >>
    -> (p1 : Nucleator) (b0 : CellBoundary)
       << (nuc_pos : p1::Position, nuc_count : p1::Count),
          (b_pos : b0::Position, b_unit : b0::Direction) >>
    with (heaviside(10, 1)) where { ... }

# Vertex attribute access — bind names to attributes of nodes c1 and c2
equalize_density := (c1 : Fluid) -- (c2 : Fluid)
    << (rho1 : c1::rho, rho2 : c2::rho) >>
    ->
    (c1 : Fluid) -- (c2 : Fluid)
    << (rho1 : c1::rho, rho2 : c2::rho) >>
    with (1.0) where {
        rho1 = (rho1 + rho2) / 2.0
    }

# Edge on both sides
growing_rule := (im : Intermediate) -- (pos : Positive)
    << (im_pos : im::Position), (p_pos : pos::Position) >>
    -> (new_im : Intermediate) -- (new_pos : Positive)
       << (im_pos : new_im::Position), (dpos : new_pos::Position) >>
    solving (...) { ... }

# Multiple disconnected components (two edges)
boundary_catastrophe := (im0 : Intermediate) -- (pos : Positive)
    (b0 : CellBoundary) -- (b1 : CellBoundary)
    << (im_pos : im0::Position), (p_pos : pos::Position),
       (b0_pos : b0::Position), (b1_pos : b1::Position) >>
    -> (im0 : Intermediate) -- (ret : Retraction)
       (b0 : CellBoundary) -- (b1 : CellBoundary)
       << ... >>
    with (...) where { ... }
```

---

## Modify Clauses

Every rule has exactly one modify clause: a `with` clause (stochastic /
propensity-driven rules), a `solving` clause (ODE-based rules) or an `every`
clause (clocked rules). A `with` clause can end in a `via` call instead of a
`where` block: a stochastic rule whose update runs another grammar as a nested
simulation.

```
<modify-clause> ::= <with-clause>
                  | <solve-clause>
                  | <every-clause>
```

### Every Clause (clocked rules)

A clocked rule is a deterministic, data-only rewrite that fires on a fixed
clock instead of with a propensity. It maps to `DGGML::EveryRule`.

```
<every-clause> ::= "every" "(" <expr> [ "," <expr> [ "," <expr> ] ] ")"
                   [ "additive" ] "where" "{" <where-body> "}"
```

The three arguments are `n_r` (fire every n_r-th tick, default 1), `p_r`
(phase, default 0) and `T_clock` (clock period; default 0, meaning the
simulation `DELTA`). They may only use numbers and parameters. `every` is a
contextual keyword, recognised only directly after a rule's RHS parameters.

The LHS and RHS must contain the same nodes (same alias and type, same order)
and the same edges. The first node is the principal node; only it is
committed. With `additive`, the engine adds `after - before` to the live
value of each bound attribute instead of overwriting it.

```
friction_relax := (w : Wall) << (fr : w::Friction) >> ->
    (w : Wall) << (frn : w::Friction) >>
    every (1, 0, 0.5) additive where {
        frn = fr * 0.9
    }
```

### Via Clause (sub-grammar call)

A `via` rule is a `with` rule (stochastic, fires with the given propensity) whose
update, instead of a `where` block, is a call to another grammar. On each firing it

1. **unpacks** the matched LHS into a fresh graph (the `unpack` pattern, filled in
   by the same shared-name copies / `where` assignments as any rule),
2. runs the grammar `g` on that graph for `t` time units (a nested, quiet,
   serial simulation in a small geoplex of its own), and
3. **packs** the final graph back: the `pack` pattern is matched in it and its
   bindings are written to the rule's RHS.

```
<with-via>     ::= "with" <propensity> "via" <identifier> "(" <expr> ")" "{" { <via-section> } "}"
<via-section>  ::= ("unpack" | "pack") ("=" | ":") <graph-pattern>
                   [ "<<" <parameters> ">>" ] [ "where" "{" <where-body> "}" ]
```

`<identifier>` names a `Grammar` declared in the simulations section; `t` may
only use numbers and parameters. Both sections are required. The LHS and RHS are
ordinary `with`-rule graphs: nodes missing from the RHS are destroyed, new RHS
nodes are created at the position of the first LHS node (unless the pack `where`
sets it) and nodes sharing an alias persist. Pack matching takes the first match
(nodes of equal type, pattern edges present, ascending key order); with no match
the RHS keeps its default values and a warning is printed once. Unpack/pack node
aliases must differ from the rule's own aliases. Not supported with
`domain_parallel`. Grammars with PDE rules are not supported as sub-grammars.

```
cell_signalling := (c : Cell) << (c1 : c::C1, c2 : c::C2) >> ->
    (c : Cell) << (c1_new : c::C1, c2_new : c::C2) >>
    with (4.0) via intra_cellular(0.5) {
        unpack = (pa : ProteinA) (pb : ProteinB) << (c1 : pa::Level, c2 : pb::Level) >>
        pack   = (qa : ProteinA) (qb : ProteinB) << (c1_new : qa::Level, c2_new : qb::Level) >>
    }
```

See `src/tests/via_cell`.

### With Clause

```
<with-clause> ::= "with" <propensity> "where" "{" <where-body> "}"

<propensity> ::= "(" <expression> ")"

<where-body> ::= { <where-entry> "\n" }

<where-entry> ::= <definition>
                | <type-instance-update>
```

#### Definition (intermediate variable)

Introduces a new local variable computed from an expression. Used in both
`where` bodies and `solving` bodies.

```
<definition> ::= <symbol-name> ":" <type-signature> ":=" <expression>
```

#### Type Instance Update (attribute assignment)

Assigns a value to a bound attribute name. The LHS target can be a bare name
(assigns the whole attribute) or an indexed name (assigns a single element of
a tensor attribute).

```
<type-instance-update> ::= <lhs-target> "=" <update-value>

<lhs-target> ::= <symbol-name>
               | <symbol-name> { "[" <index-list> "]" }

<update-value> ::= "{" <type-declaration-list> "}"
                 | <expression>
```

### Examples (where body)

```
where {
    mt_seg_len : Float := ~UniformDistribution(mt_min, mt_max)

    x_c : Float := nuc_pos[0] + ~UniformDistribution(-1*eps, eps)

    ret_pos = [x_l, y_l, z_l]
    im_pos = [x_c, y_c, z_c]

    nuc_count = [0, nuc_count[1]]

    np_unit = [ -p_unit[0], -p_unit[1], p_unit[2] ]
}
```

---

### Solve Clause

```
<solve-clause> ::= "solving" "(" <binding-variable-list> ")"
                   "{" <solve-body> "}"

<binding-variable-list> ::= <binding-variable> { "," <binding-variable> }

<binding-variable> ::= <binding-name> ":=" "D" "(" <ode-var-list> ")"

<binding-name> ::= <symbol-name>
                 | <symbol-name> "[" <expression> "]"

<ode-var-list> ::= <ode-var> { "," <ode-var> }

<ode-var> ::= <symbol-name>
            | <symbol-name> "[" <expression> "]"
```

The binding variables establish which derivatives are being solved. For
example, `dpos[0] := D(im_pos[0], t)` says "`dpos[0]` is the time
derivative of `im_pos[0]`".

```
<solve-body> ::= { <solve-entry> "\n" }

<solve-entry> ::= <ode-equation>
                | <definition>

<ode-equation> ::= <ode-name> ":" "ODE" "=" <expression>

<ode-name> ::= <symbol-name>
             | <symbol-name> "[" <index-list> "]"
```

A `<definition>` inside a solving body uses `:=` (define token) and creates
an intermediate local variable, while an `<ode-equation>` uses `=` (equal
token) and contributes to the right-hand side of the ODE system.

### Examples (solving)

```
solving (dpos[0] := D(im_pos[0], t), dpos[1] := D(im_pos[1], t)) {

    check : Float := HELP::distance(im_pos[0], im_pos[1], p_pos[0], p_pos[1])

    dpos[0] : ODE = 0.0615 * p_unit[0] * 4
    dpos[1] : ODE = 0.0615 * p_unit[1] * 4
}
```

### PDE Solve Clause (MFEM)

A solving body that contains an `Element` entry is a PDE rule. It works like an ODE
solving rule: one match of the LHS pattern is one finite element, so the mesh is the
set of matches in the graph and follows it when the graph changes. The PDE is written
in differential form; FoxFlow derives the MFEM weak form `M du/dt = L u + b` and
generates a `<rule>_Operator` class that the DGGML solver integrates together with the
grammar ODEs (implicitly when there is no source or Neumann term, explicitly otherwise).

```
<pde-entry> ::= <symbol-name> ":" <pde-kind> "=" <expression>

<pde-kind>  ::= "Element" | "PDE" | "Dirichlet" | "Neumann"
```

| Entry | Meaning |
|---|---|
| `c : Element = segment(a, b)` / `triangle(a, b, c)` / `quad(a, b, c, d)` / `tet(a, b, c, d)` / `hex(a, ..., h)` | element shape over LHS nodes (exactly one); the pattern edges must cover the element edges |
| `dua : PDE = <expr>` | right-hand side of `D(ua, t)` (exactly one; named after one binding variable) |
| `x : Dirichlet = boundary(<predicate>, g)` | hold `u = g` at every mesh vertex whose position satisfies the predicate |
| `x : Neumann = boundary(<predicate>, g)` | flux `du/dn = g` on boundary faces whose vertices satisfy the predicate |

Predicates and values may use `x`, `y`, `z`, `pi`, parameters, comparisons and `&&`/`||`.
Supported terms in the `PDE` expression: `c * laplacian(u)`, `div(c * grad(u))`,
`dot([bx, by], grad(u))`, `dx(u)`/`dy(u)`/`dz(u)`, `c * u`, and `u`-free terms (sources).

Every element node needs one binding variable on the same `Float` field (the unknown).
Nodes need no extra fields: the field lives on the graph nodes. The initial graph is
produced in `simulation.fflow` by a generator state:

```
initialState : State := grid(Vertex, u, <init expr of x, y, z>, [x0, x1, nx], [y0, y1, ny])
```

Each axis is `[a, b, n]` (n elements) or `[c]` (fixed); missing axes sit at the cell
centre. When MFEM is not found, the generated CMake fetches MFEM v4.7; pass
`-DMFEM_DIR=<mfem-build>` (or set the `MFEM_DIR` environment variable) to reuse a
local build. Examples: `src/tests/heat_rod_pde/`, `heat_plate_pde/`, `heat_coupled_pde/`.

```
heat := (a : Vertex) -- (b : Vertex) << (ua : a::u), (ub : b::u) >> ->
        (a : Vertex) -- (b : Vertex) << (dua : a::u), (dub : b::u) >>
    solving (dua := D(ua, t), dub := D(ub, t)) {
        cell  : Element   = segment(a, b)
        dua   : PDE       = kappa * laplacian(ua)
        left  : Dirichlet = boundary(x < rod_start + 0.000001, u_left)
        right : Dirichlet = boundary(x > rod_end - 0.000001, u_right)
    }
```

---

## Functions Section

```
<functions-section> ::= "functions" <symbol-name> "{" <function-list> "}"

<function-list> ::= { <function-declaration> "\n" }

<function-declaration> ::= <symbol-name> ":" <function-type-signature>
                           ":=" <function-body>

<function-type-signature> ::= "Function" "<<" "(" <function-arg-list> ")" ">>"
                              "->" <return-type>

<function-arg-list> ::= <function-arg> { "," <function-arg> }

<function-arg> ::= <symbol-name> ":" <arg-type>

<arg-type> ::= <type-name>
             | ("FixedList" | "TFTensor") "<<" <integer> "," <type-name> ">>"

<return-type> ::= <type-name>
                | ("FixedList" | "TFTensor") "<<" <integer> "," <type-name> ">>"
```

### Function Body

```
<function-body> ::= <regular-function-body>
                  | <model-load>

<regular-function-body> ::= "{" { <function-statement> "\n" }
                            "return" <expression> "}"

<function-statement> ::= <symbol-name> ":" <type-name> ":=" <expression>

<model-load> ::= "load" "(" <string-literal> ")"
```

### Examples

```
functions Microtubule {
    distance : Function << (x1 : Float, y1 : Float, x2 : Float, y2 : Float) >> -> Float := {
        dx : Float := x2 - x1
        dy : Float := y2 - y1
        return sqrt(dx * dx + dy * dy)
    }
}
```

---

## Simulations Section

```
<simulations-section> ::= "simulations" <symbol-name> "{"
                           <sim-declaration-list>
                          "}"

<sim-declaration-list> ::= { ( <sim-declaration> | <trace-statement> | <measure-statement> | <run-statement> ) "\n" }

<sim-declaration> ::= <symbol-name> ":" <sim-type> ":=" <sim-value>

<trace-statement> ::= "trace" "(" <symbol-name> { "," <string-literal> } ")"

<measure-statement> ::= "measure" "(" <symbol-name> { "," <symbol-name> } { "," <string-literal> } ")"

<sim-type> ::= "SimulationParameters"
             | "State"
             | "SimulationRules"
             | "SimulationTypes"
             | "SimulationObservables"
             | "SimulationFunctions"
             | "Grammar"
             | "Integer"
             | "Float"
             | "Simulation"

<sim-value> ::= <load-file>
              | <number>
              | <run-simulation>
              | <grammar-value>
              | <simulation-config>

<grammar-value> ::= "Grammar" "(" <symbol-name> { "," <symbol-name> } ")"

<simulation-config> ::= "Simulation" "(" [ <string-literal> "," ]
                        <symbol-name> ","   (* a Grammar *)
                        <symbol-name> ","   (* a State *)
                        <symbol-name> ","   (* SimulationParameters *)
                        <symbol-or-number>  (* simulation time *)
                        [ "," <integer-or-symbol> ]   (* checkpoint frequency *)
                        [ "," "dim" "=" ( "2" | "3" ) ]   (* cell complex dimension, default 3 *)
                        ")"

<run-statement> ::= "run" "(" <symbol-name> "," <integer>
                    { "," ( "domain_parallel" "=" ( "true" | "false" )
                          | "jobs" "=" <integer>
                          | "seed" "=" <integer>
                          | "name" "=" <string-literal>
                          | "listeners" "=" ( "true" | "false" )
                          | "log" "=" ( "true" | "false" ) ) }
                    ")"

<load-file> ::= "load" "(" <string-or-identifier> ")"

<string-or-identifier> ::= <string-literal> | <symbol-name>

<run-simulation> ::= "RunSimulation" "("
                     [ <string-literal> "," ]
                     <symbol-name> ","
                     <symbol-name> ","
                     <symbol-name> ","
                     <symbol-name> ","
                     <symbol-or-number>
                     [ "," <integer-or-symbol> ]
                     [ "," <symbol-name> ]
                     ")"

<symbol-or-number> ::= <symbol-name> | <number>
<integer-or-symbol> ::= <integer> | <symbol-name>
```

The optional leading `<string-literal>` sets the file name used for the
"latest" save-state bin file written by the checkpoint (defaults to
`simulation_state_latest.bin`).

The optional `<integer-or-symbol>` after the simulation time is the checkpoint
frequency: VTK/CSV output, observables and save-state files are written every N
steps (and always on the final step). It defaults to 1 (every step). The
observables argument, when present, always comes last; a single trailing
identifier bound to an integer is read as the checkpoint frequency.

The optional `dim = 2|3` keyword of `Simulation(...)` selects the expanded cell
complex the model runs on: `3` (default) uses `ExpandedComplex3D`, `2` uses
`ExpandedComplex2D`. Each `Simulation` picks its own dimension.

### Grammars, simulations and replicated runs

A `Grammar` bundles the pieces that define a model: its functions, rules, types
and (optionally) observables. A `Simulation(...)` then pairs a grammar with an
initial state, parameters, a time and a checkpoint frequency. `RunSimulation(...)`
stays as the older, flat form and behaves as before.

```
funcs   : SimulationFunctions   := load("functions.fflow")
rules   : SimulationRules       := load("rules_movement.fflow")
types   : SimulationTypes       := load("types.fflow")
obs     : SimulationObservables := load("observables.fflow")
g       : Grammar               := Grammar(funcs, rules, types, obs)
start   : State                 := load("my_results/boundary_latest.bin")
sim     : Simulation            := Simulation("movement_latest.bin", g, start, params, time, every)

run(sim, 8, domain_parallel = false)
```

`Grammar(...)` takes its parts in any order; each is a rules, types, observables
or functions declaration (a grammar needs rules and types). The functions file
defaults to `functions.fflow`; stages naming different functions files are an error.

`run(sim, N, ...)` runs `N` replicates of `sim`. A `Simulation(...)` declaration only
defines the simulation; it is executed only when a `run(...)` statement names it (a
simulation that is never run can still be measured from existing checkpoints).
`RunSimulation(...)` runs without a `run` statement, in-process and serially. Options:

- `domain_parallel = true|false`: compiles the simulation with the OpenMP geocell
  domain-decomposition runner (`true`) or the serial runner (`false`, the default).
  With `false` each replicate runs single-threaded and several replicates run at once.
- `jobs = K`: replicates running at once (default 1 when domain parallelism is on,
  otherwise the number of hardware threads; never more than `N`).
- `seed = S`: base seed; replicate `i` is seeded with `S + i` (default base 0).
- `name = "exp"`: names the experiment. Its runs go to `runs/<name>/run_<i>/` instead of
  `runs/<sim>/run_<i>/`, also for `N = 1` (which is then run in a child process instead of
  in `my_results/`). The name may contain letters, digits, `_`, `-` and `.`, and must be
  unique among the `run` statements. A `measure(...)` whose simulations all belong to the
  same named experiment writes its CSV to `runs/<name>/measure_<k>.csv`.
- `listeners = true|false`: with `false` the simulator writes no listener files
  (`PropensityListener.csv`, `RuleFireListener.csv`, `TimeProgressListener.csv`,
  `InvalidationValidationListener.csv`, `cell_listener.txt`). Default `true`. `trace(...)`
  output is independent of this option.
- `log = true|false`: with `false` the standard output of the run (progress and per-geocell
  timing) is discarded. For replicates `run.log` is then empty apart from errors, which
  are still written to it. Default `true`.

With `N > 1`, the program starts one child process per replicate (the RNG and
the OpenMP settings are process-global). Replicate `i` runs in
`runs/<sim>/run_<i>/` with its log in `run.log` and its own `my_results/`. Input
files named in the `.fflow` (`load("...")` states, torch models) are resolved
against the base directory first (`FOXFLOW_BASE_DIR`). A failing replicate makes
the program exit non-zero. `measure(sim, ...)` of a replicated simulation reads
every replicate and adds a `run` column (see below). Trace and listener files of a
replicated simulation are written inside the run directories.

### Rule tracing

`trace(<simulation>, "rule1", "rule2", ...)` records every fire of the named
rules in the simulation `<simulation>` (a `Simulation` declared in the same
section). `"*"` traces all rules. Each record holds the matched nodes before and
after the rewrite, their neighbours, edge changes, the propensity and the
`where` locals. Rule names are checked against the simulation's rules file at
generation time. Records go to `trace_<simulation>.jsonl` in the run directory;
when the run ends, `trace_<simulation>.html` is written next to it. Re-render with
filters, without rebuilding: `./run_example.sh trace <example> --rule <name> --step <n> --max <n>`.

### Measuring global observables over simulations

`measure(sim1, sim2, ..., "obs1", "obs2", ...)` re-evaluates the global observables
(see the `observables` section) on every checkpointed state of the listed
simulations. The simulation names come first; the optional trailing strings
restrict the measurement to those observables (default: all). Each `measure`
statement produces `my_results/measure_<k>.csv` (`k` counts the statements from 1)
with columns `simulation,run,step,time,<observable...>` (`run` is the replicate
index, 0 for a single run). Measuring several simulations
lists them in one statement; measure statements may overlap.

A measured simulation writes its checkpoints to `my_results/<sim>/state_<step>.bin`
(cleared when it starts) instead of the shared `simulation_state_<step>.bin`.
Observable names are checked at generation time. Stages that use the same rules
file share one generated rules header.

### Examples

```
simulations Microtubule {
    params : SimulationParameters := load("params.fflow")
    types  : SimulationTypes      := load("types.fflow")
    rules  : SimulationRules      := load("rules.fflow")
    time   : Float                := 100.0
    sim    : Simulation           := RunSimulation(types, params, rules, funcs, time)
}

simulations NeuralNetwork {
    params : SimulationParameters  := load("params.fflow")
    types  : SimulationTypes       := load("types.fflow")
    rules  : SimulationRules       := load("rules.fflow")
    obs    : SimulationObservables := load("observables.fflow")
    time   : Float                 := 100.0
    sim    : Simulation            := RunSimulation(types, params, rules, funcs, time, obs)
}

simulations Microtubule {
    # Optional first argument: name of the latest save-state bin file.
    sim : Simulation := RunSimulation("my_state_latest.bin", state, params, rules, types, time)

    # Optional checkpoint frequency after the time: output every 60 steps.
    every : Integer := 60
    sim2  : Simulation := RunSimulation("my_state_latest.bin", state, params, rules, types, time, every)

    # Trace every fire of the "crossover" rule in sim2.
    trace(sim2, "crossover")
}
```

---

## Observables Section

An `observables` section defines two kinds of measurements:

- an **`Observable`** is a functional of the whole state (a pattern plus a
  reduction over its matches). It is not evaluated during a run: it is only
  computed by `measure(...)` on the saved checkpoints, and nothing is generated
  for it when the program has no `measure` statement;
- an **`Attribute`** is a per-node value written as a ParaView column.

```
<observables-section> ::= "observables" <symbol-name> "{"
                           { <global-observable> | <observable-definition> }
                          "}"

<global-observable> ::= <symbol-name> ":" "Observable" ":=" <reduction>
                         <pattern> [ "within" <number> ]
                         [ "<<" <binding-list> ">>" ]
                         [ "return" "(" <expression> ")" ]

<reduction> ::= "sum" | "mean" | "count" | "min" | "max"

<observable-definition> ::= <symbol-name> ":" <obs-kind>
                             "<<" <obs-param-block> ">>"
                             "->" <return-type>
                             ":=" <obs-body>

<obs-kind> ::= "Attribute" | "Function"

<obs-param-block> ::= "(" <symbol-name> ":" <symbol-name> ")"
                      [ "<<" <function-arg-list> ">>" ]

<obs-body> ::= "{" { <function-statement> "\n" } "return" <expression> "}"
```

- The `(p : TypeName)` block is the **primary node argument** — the graph
  node instance being observed.
- The inner `<< alias : FieldType, ... >>` block **destructures** that node's
  attributes into named local variables inside the body.
- `Attribute` definitions are called automatically by the runtime for every
  matching node at each output step. (The former spelling
  `name : Observable << ... >>` is still accepted as an `Attribute`.)
- `Function` definitions are helper functions callable from `Attribute`
  bodies.

### Global observables

```
n_tips : Observable := count (p : Positive)

tip_segment_length : Observable := mean
    (a : Intermediate) -- (b : Positive)
    << (a_pos : a::Position), (b_pos : b::Position) >>
    return ( HELP::distance(a_pos[0], a_pos[1], b_pos[0], b_pos[1]) )

# pair observable: mean of cos^2 between the growth directions of all tip pairs closer than 0.1
local_alignment : Observable := mean
    (a : Positive) (b : Positive) within 0.1
    << (a_dir : a::Direction), (b_dir : b::Direction) >>
    return ( HELP::unit_dot3_abs(a_dir[0], a_dir[1], a_dir[2], b_dir[0], b_dir[1], b_dir[2]) *
             HELP::unit_dot3_abs(a_dir[0], a_dir[1], a_dir[2], b_dir[0], b_dir[1], b_dir[2]) )
```

- The pattern uses the rule-pattern syntax and must be connected, except for a
  **pair observable**: `(a : T) (b : U) within r` has exactly two nodes and no
  edges, and matches every ordered pair of distinct nodes whose spatial distance
  is at most `r` (found with a uniform grid, so the cost is linear in the node
  count). The optional `<< ... >>` binds node attributes to names, exactly as in
  a rule.
- `return ( expr )` is the value a matched configuration returns, evaluated once per
  match of the pattern (`with` is reserved for rule propensities and is rejected
  here). `sum`, `mean`, `min` and `max` combine those values (`mean`/`min`/`max` are `nan` when there
  is no match); `count` needs no expression and counts the matches (with an
  expression it counts the matches where it is nonzero).
- Matches are ordered: a pattern whose nodes have the same type, such as
  `(a : T) -- (b : T)`, matches each edge twice.
- Each observable is compiled to a hidden rule `obs_<name>` (header
  `observable_rules.h`); the model evaluates it at every checkpoint, so output
  follows the `RunSimulation` checkpoint frequency.

### Examples

```
observables NeuralNetwork {
    # Helper: run a forward pass and return the predicted class
    classify : Function
        << (node : InputLayer) << input_id : Integer, layer_id : Integer >> >>
        -> Integer := {
            img   : torch::Tensor := images[input_id]
            logits : torch::Tensor := forward(img)
            return argmax(logits)
        }

    # Attribute: emit (node_id, predicted_label, true_label) each step
    accuracy : Attribute
        << (node : InputLayer) << input_id : Integer >> >>
        -> Float := {
            pred  : Integer := classify(node)
            truth : Integer := labels[input_id]
            return indicator(pred == truth)
        }
}
```

---

## Built-in Functions

These are not declared in the grammar but are recognised by the code
generator:

### General

| Function | Description |
|----------|-------------|
| `indicator(pred)` | 1.0 if predicate is true, else 0.0 |
| `heaviside(x, n)` | Heaviside step function |
| `cos(x)`, `sin(x)`, `sqrt(x)` | Standard math |
| `HELP::distance(...)` | Namespaced helper (euclidean distance) |
| `HELP::minimum_distance_2d(...)` | Namespaced helper |
| `~UniformDistribution(lo, hi)` | Random sample (unary `~` prefix) |

### Tensor Operations

These built-ins operate on `TFTensor` values (`FixedList` supports the subset listed above) and are
available in expressions, rule bodies, and observable bodies:

| Function | Description |
|----------|-------------|
| `zeros_matrix(m, n)` | `m × n` zero tensor |
| `rand_matrix(m, n)` | `m × n` uniform-random tensor |
| `mat_dot(A, B)` | Matrix–matrix or matrix–vector product |
| `mat_add(A, B)` | Element-wise addition |
| `mat_mul(A, B)` | Element-wise (Hadamard) multiplication |
| `transpose(A)` | Matrix transpose |
| `einsum(expr, A, B)` | Einstein summation (passes through to `torch::einsum`) |
| `permute(A, d0, d1, ...)` | Permute tensor dimensions |
| `autodiff(loss, param)` | Compute gradient of `loss` w.r.t. `param` |
| `from_file(path)` | Load raw binary file into tensor (params only — see Parameters Section) |

Additional built-ins available in `functions` sections:

| Function | Description |
|----------|-------------|
| `exp(x)` | Element-wise exponential |
| `dot(a, b)` | Inner product of two 1-D tensors (returns a scalar) |
| `mat_inv(A)` | Matrix inverse |
| `einsum(spec, A, B, ...)` | Emitted as `torch::einsum(spec, {A, B, ...})` |

Any other function name inside a `functions` body is emitted as a raw call,
so functions of the same section (and earlier sections) can call each other,
provided they are defined first. Rules call them as `NAME::fn(args)`.

An empty subtype is allowed in a `types` section: `Fluid << >> : Site`.

---

## Namespace-Qualified Identifiers

Function calls and identifiers can be namespace-qualified with `::`:

```
<namespaced-call> ::= <symbol-name> "::" <symbol-name> "(" <arg-list> ")"
```
