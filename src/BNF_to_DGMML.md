
## Symbol
<symbol-name> ::= [Aa-Zz] | [0-9] | _ | <symbol-name>

<symbol-parameters> ::= "<<" <symbol-parameter> ">> | ""

<symbol-parameter> ::= "(" <identifier-list> ")" | <identifier-list>
<identifier-list> ::= <symbol-name> | <symbol-name>, <identifier-list>

##  Types
<type-section> ::= "types" <symbol-name> "{" <type-declaration-list> "}"

<type-declaration-list> ::= <type-declaration> "\n" | <type-declaration> <type-declaration-list>

<type-declaration> ::= <symbol-name> <symbol-parameters> ":" <type-signature-list>
                    
<type-assignment> ::= <symbol-name> <symbol-parameters> ":" <type-signature-list> ":=" <literal>
                    | <symbol-name> <symbol-parameters> ":" <type-signature-list> ":=" "{" <type-declaration-list> "}"
                    | <symbol-name> ":" <type-signature-list> ":=" "{" <type-declaration-list> "}"
                    | <symbol-name> ":" <type-signature-list> ":=" <literal>
                    | <symbol-name> ":" <type-signature-list> ":=" <expression>


<type-signature-list> ::= <symbol-name> "->" <type-signature-list> 
                        | <symbol-name>
                        | <built-in-type>

<built-in-type> ::= "Float" | "Integer"

## Parameters
<parameters> ::= "parameters" <symbol-name> "{" <parameter-list> "}"
<parameter-list> ::= <parameter> "\n" | <parameter> <parameter-list>
<parameter> ::= <type-declaration> | <type-assignment>

## Functions
<functions> ::= "functions" <symbol-name> "{" <function-list> "}"
<function-list> ::= <function> <function-list> | <function> "\n"
<function> ::= <type-assignment>

## Expressions
<expression> ::= <term> | <expression> "+" <term> | <expression> "-" <term>
<term> ::= <factor> | <term> "*" <factor> | <term> "/" <factor>
<factor> ::= <number> | <variable> | "(" <expression> ")"
<number> ::= <digit> | <digit> <number>
<digit> ::= [0-9]
<variable> ::= <letter> | <letter> <variable>
<letter> ::= [Aa-Zz] 

## Rules
<rules> ::=  "rules" <symbol-name> "{" <rules-list> "}"
<rule-list> ::= <rule> <rules-list> | <rule>

<rule> ::= <parameterized-types-left>  <symbol-parameters> 
        "->" <parameterized-types-right> <symbol-parameters> "with" "(" <with-clause> ")" "where" "{" <where-clause> "}" "\n"

<solve-rule> ::= <parameterized-type-left> "->" <parameterized-type-right> "solving" <solve-clause> "where" "{" <where-clause> "}" "\n"

<parameterized-types-left> ::= <parameterized-types> "," | <parameterized-type>
<parameterized-types-right> ::= <parameterized-types> "," | <parameterized-type>
<parameterized-types> ::= <parameterized-types> "," | <parameterized-type>

<parameterized-type> ::= <symbol-name> ":" <type-signature-list> 

# TODO:  define predicate and expressions to be a list of them.
<with-clause> ::= <function-type> "where" "{" <predicate> "}" "\n" | <function-type> "\n"
<solve-clause> ::= <function-type> "solving" "{" <expression> "}" "\n" | <function-type> "\n"
<where-clause> ::= <expressions>

## Observables
<observables> ::=  "observables" <symbol-name> "{" <observables-list> "}"
<observable-list> ::= <observable> <observable-list> | <observable>
<observable> ::= <type-declaration>

## State
<states> ::=  "states" <symbol-name> "{" <states-list> "}"
<states-list> ::= <state> <states-list> | <state>
<state> ::= <type-declaration>

## Grammar
<grammar> ::= "grammar" <symbol-name> ":" <time-type> 
          "<<" <parameterized-types> ">>" "{" <rules> "}"
<time-type> ::= "DiscreteTime)" | "ContinuousTime"

## Simulations
<simulations> ::=  "simulations" <symbol-name> "{" <simulations-list> "}"
<simulations-list> ::= <simulation> <simulation-list> | <simulation>
<simulation> ::= "ApproximateSimulation" "<<" <parameterized-types> ">>" "\n"


<rules> ::= <rules> | <rule> | <solve-rule>
<grammar-signature> ::= <parameterized-types-left> "->" <parameterized-types-right>

 <predicate> ::= <symbol-name> "!=" <symbol-name> 
                | <symbol-name> "<" <symbol-name> 
                | <symbol-name> ">" <symbol-name> 
                | <symbol-name> "==" <symbol-name>

<function-type> ::= <symbol-name> "(" <symbol-name> ")"
<literal> ::= [0-9] | [Aa-Zz] | <literal>

