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
```

---

## Type Signatures

```
<type-signature> ::= <type-name> <symbol-parameters>
                   | <type-name> "->" <type-signature>
                   | <type-name>

<type-name> ::= <symbol-name> | <built-in-type>

<built-in-type> ::= "Float" | "Integer" | "FixedList" | "ODE" | "Type"
```

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
                    Count : FixedList<<2, Integer>> >> : Type
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

A `FixedList` parameter with `from_file` loads a raw binary file (saved via
`numpy().tofile()`) into a `torch::Tensor` at startup using memory-mapping:

```
images : FixedList<<60000, 1, 28, 28, Float>>   := from_file("mnist_images.bin")
labels : FixedList<<60000, Integer>>             := from_file("mnist_labels.bin")
```

The last type parameter of the `FixedList` (`Float` or `Integer`) determines
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
    images : FixedList<<60000, 1, 28, 28, Float>> := from_file("mnist_images.bin")
    labels : FixedList<<60000, Integer>>           := from_file("mnist_labels.bin")
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

<type-instance> ::= "(" <symbol-name> ":" <type-signature> ")"
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

### Examples

```
# One LHS node, two RHS nodes (no edges)
start_to_node := (start : StartType) << (start_pos : FixedList<<3, Float>>) >>
    -> (p1 : Nucleator) (b0 : CellBoundary)
       << (nuc_pos : FixedList<<3, Float>>, nuc_count : FixedList<<2, Integer>>),
          (b_pos : FixedList<<3, Float>>, b_unit : FixedList<<3, Float>>) >>
    with (heaviside(10, 1)) where { ... }

# Edge on both sides
growing_rule := (im : Intermediate) -- (pos : Positive)
    << (im_pos : FixedList<<3, Float>>), (p_pos : FixedList<<3, Float>>) >>
    -> (new_im : Intermediate) -- (new_pos : Positive)
       << (im_pos : FixedList<<3, Float>>), (dpos : FixedList<<3, Float>>) >>
    solving (...) { ... }

# Multiple disconnected components (two edges)
boundary_catastrophe := (im0 : Intermediate) -- (pos : Positive)
    (b0 : CellBoundary) -- (b1 : CellBoundary)
    << (im_pos : FixedList<<3, Float>>), (p_pos : FixedList<<3, Float>>),
       (b0_pos : FixedList<<3, Float>>), (b1_pos : FixedList<<3, Float>>) >>
    -> (im0 : Intermediate) -- (ret : Retraction)
       (b0 : CellBoundary) -- (b1 : CellBoundary)
       << ... >>
    with (...) where { ... }
```

---

## Modify Clauses

Every rule has exactly one modify clause, either a `with` clause (for
stochastic / propensity-driven rules) or a `solving` clause (for ODE-based
rules).

```
<modify-clause> ::= <with-clause>
                  | <solve-clause>
```

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
             | "FixedList" "<<" <integer> "," <type-name> ">>"

<return-type> ::= <type-name>
                | "FixedList" "<<" <integer> "," <type-name> ">>"
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

<sim-declaration-list> ::= { <sim-declaration> "\n" }

<sim-declaration> ::= <symbol-name> ":" <sim-type> ":=" <sim-value>

<sim-type> ::= "SimulationParameters"
             | "State"
             | "SimulationRules"
             | "SimulationTypes"
             | "SimulationObservables"
             | "Integer"
             | "Float"
             | "Simulation"

<sim-value> ::= <load-file>
              | <number>
              | <run-simulation>

<load-file> ::= "load" "(" <string-or-identifier> ")"

<string-or-identifier> ::= <string-literal> | <symbol-name>

<run-simulation> ::= "RunSimulation" "("
                     <symbol-name> ","
                     <symbol-name> ","
                     <symbol-name> ","
                     <symbol-name> ","
                     <symbol-or-number>
                     [ "," <symbol-name> ]
                     ")"

<symbol-or-number> ::= <symbol-name> | <number>
```

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
```

---

## Observables Section

An `observables` section defines per-node measurement functions that are
evaluated after each simulation step and written to output.

```
<observables-section> ::= "observables" <symbol-name> "{"
                           { <observable-definition> }
                          "}"

<observable-definition> ::= <symbol-name> ":" <obs-kind>
                             "<<" <obs-param-block> ">>"
                             "->" <return-type>
                             ":=" <obs-body>

<obs-kind> ::= "Observable" | "Function"

<obs-param-block> ::= "(" <symbol-name> ":" <symbol-name> ")"
                      [ "<<" <function-arg-list> ">>" ]

<obs-body> ::= "{" { <function-statement> "\n" } "return" <expression> "}"
```

- The `(p : TypeName)` block is the **primary node argument** — the graph
  node instance being observed.
- The inner `<< alias : FieldType, ... >>` block **destructures** that node's
  attributes into named local variables inside the body.
- `Observable` definitions are called automatically by the runtime for every
  matching node at each output step.
- `Function` definitions are helper functions callable from `Observable`
  bodies.

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

    # Observable: emit (node_id, predicted_label, true_label) each step
    accuracy : Observable
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

These built-ins operate on `FixedList` / `torch::Tensor` values and are
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

---

## Namespace-Qualified Identifiers

Function calls and identifiers can be namespace-qualified with `::`:

```
<namespaced-call> ::= <symbol-name> "::" <symbol-name> "(" <arg-list> ")"
```
