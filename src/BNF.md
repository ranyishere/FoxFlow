## Symbol
<symbol-name> ::= [Aa-Zz] | [0-9] | _ | <symbol-name>
<symbol-parameters> ::= "<<" <symbol-parameter> ">> | ""
<symbol-parameter> ::= "(" <identifier-list> ")" | <identifier-list> | <literal>
<identifier-list> ::= <symbol-name> | <symbol-name>, <identifier-list>

##  Types
<type-section> ::= "types" <symbol-name> "{" <type-declaration-list> "}"

<type-declaration-list> ::= <type-declaration> "\n" | <type-declaration> <type-declaration-list>

<type-declaration> ::= <symbol-name> <symbol-parameters> ":" <type-signature-list>

<type-assignment-list> ::= <type-assignment> | <type-assignment> "\n" <type-assignment-list>
                    
<type-assignment> ::= <symbol-name> <symbol-parameters> ":" <type-signature-list> ":=" <literal>
                    | <symbol-name> <symbol-parameters> ":" <type-signature-list> ":=" "{" <type-declaration-list> "}"
                    | <symbol-name> ":" <type-signature-list> ":=" "{" <type-declaration-list> "}"
                    | <symbol-name> ":" <type-signature-list> ":=" <literal>
                    | <symbol-name> ":" <type-signature-list> ":=" <expression>

<type-update-list> ::= <type-update> | <type-update> "\n" <type-update-list>
<type-update> ::= <symbol-name> <symbol-parameters> ":" <type-signature-list> "=" <literal>
                    | <symbol-name> <symbol-parameters> ":" <type-signature-list> "=" "{" <type-declaration-list> "}"
                    | <symbol-name> ":" <type-signature-list> "=" "{" <type-declaration-list> "}"
                    | <symbol-name> ":" <type-signature-list> "=" <literal>
                    | <symbol-name> ":" <type-signature-list> "=" <expression>


<type-signature-list> ::= <symbol-name> "->" <type-signature-list> 
                        | <symbol-name>
                        | <built-in-type>

<built-in-type> ::= "Float" | "Integer"

## Parameters
<parameters-section> ::= "parameters" <symbol-name> "{" <parameter-list> "}"
<parameter-list> ::= <parameter> "\n" | <parameter> <parameter-list> | nothing
<parameter> ::= <type-declaration> | <type-assignment>

## Functions
<functions> ::= "functions" <symbol-name> "{" <function-list> "}"
<function-list> ::= <function> <function-list> | <function> "\n"
<function> ::= <type-assignment>

## Expressions
<expression> ::= <term>
                | <expression> "+" <term> 
                | <expression> "-" <term>
<term> ::= <factor> 
            | <term> "\**" <factor> 
            | <term> "/" <factor>

<factor> ::= <number> 
            | <variable> 
            | "(" <expression> ")"

<number> ::= <digit> | <digit> <number>
<digit> ::= [0-9]
<variable> ::= <letter> | <letter> <variable>
<letter> ::= [Aa-Zz] 

## Rules
<rule-section> ::=  "rules" <symbol-name> "{" <rules-list> "}"
<rule-list> ::= <rule> <rules-list> | <rule>

<rule> ::= <symbol-name> ":=" <parameterized-types-left>  <symbol-parameters> 
        "->" <parameterized-types-right> <symbol-parameters> "with" "(" <with-clause> ")" "where" "{" <where-clause> "}" "\n"

<solve-rule> ::= <parameterized-type-left> "->" <parameterized-type-right> "solving" <solve-clause> "where" "{" <where-clause> "}" "\n"

<parameterized-types-left> ::= <parameterized-types> "," | <parameterized-type>
<parameterized-types-right> ::= <parameterized-types> "," | <parameterized-type>
<parameterized-types> ::= <parameterized-types> "," | <parameterized-type>

<parameterized-type> ::= <symbol-name> ":" <type-signature-list> 

# TODO:  define predicate and expressions to be a list of them.
<with-clause> ::= "indicator" "(" <predicate> ")" | <function-type> "\n"
<solve-clause> ::= <function-type> "solving" "{" <expression> "}" "\n" | <function-type> "\n"

<where-clause> ::= <expressions> | <type-assignment-list>

<expressions> ::= <expression> "\n" <expressions> 
                | <expression>

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
<time-type> ::= "DiscreteTime" | "ContinuousTime"


## Functions
<function-parameterized-types> ::= <function-parameterized-types> | ""
<function-parameterized-type> ::= <symbol-name> ":" <type-signature-list>

# Function Section
<functions-section> ::= "functions" <symbol-name> "{" <function-list> "}"
<function-list> ::= <function-declaration> <function-list> | <function-declaration>
<function-type-signature> ::= "<<" <parameterized-types> ">>" <symbol-name> "->"
                        | <symbol-name>

<function-declaration> ::= <symbol-name>  ":" "Function" <function-type-signature> ":=" "{" <function-expression-list> "return" <expression> "}" "\n"

<function-expression-list> ::= <expressions> | <type-assignment-list>
<function-expression> ::= "return" <expression>

## Simulations
<simulations> ::=  "simulations" <symbol-name> "{" <simulations-list> "}"
<simulations-list> ::= <simulation> <simulation-list> | <simulation>
<simulation> ::= "ApproximateSimulation" "<<" <parameterized-types> ">>" "\n"

<rules> ::= <rules> | <rule> | <solve-rule>
<grammar-signature> ::= <parameterized-types-left> "->" <parameterized-types-right>

<predicate> ::= <expression> "!=" <expression> 
            | <expression> "<" <expression> 
            | <expression> ">" <expression> 
            | <expression> "==" <expression>

<function-args> ::= <symbol-name> | <symbol-name> "," <function-args>
<function-type> ::= <symbol-name> "(" <function-args> ")"
<literal> ::= [0-9] | [Aa-Zz] | <literal>
