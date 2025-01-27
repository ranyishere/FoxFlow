 
## Symbol
<symbol-name> ::= [Aa-Zz] | [0-9] | _ | <symbol-name>
 <symbol-parameters> ::= "<<" <symbol-parameter> ">>
 <symbol-parameter> ::= "(" <identifier-list> ")"
 <identifier-list> ::= <symbol-name> | <symbol-name>, <identifier-list>

##  Types
<type-section> ::= "types" <symbol-name> "{" <type-declaration-list> "}"

<type-declaration-list> ::= <type-declaration> | <type-declaration> ; <type-declaration-list>

<type-declaration> ::= <symbol-name> <symbol-parameters> "::" <type-signature-list>;
                    | <symbol-name> <symbol-parameters> "::" <type-signature-list> "=" <literal>
                    | <symbol-name> <symbol-parameters> "::" <type-signature-list> "=" "{" <type-declaration-list> "}"
                    | <symbol-name> "::" <type-signature-list> "=" "{" <type-declaration-list> "}"
                    | <symbol-name> "::" <type-signature-list> "=" <literal>

<type-signature-list> ::= <symbol-name> "->" <type-signature-list> 
                        | <symbol-name>
                        | <built-in-type>

<built-in-type> ::= "Float" | "Integer"

## Parameters
<parameters> ::= "parameters" <symbol-name> "{" <parameter-list> "}"
<parameter-list> ::= <parameter> ";" | <parameter> <parameter-list>
<parameter> ::= <type-declaration> "=" <literal> | <type-declaration>

## Functions
<functions> ::= "functions" <symbol-name> "{" <function-list> "}"
<function-list> ::= <function> <function-list> | <function> ";"
<function> ::= <type-declaration> "=" <expression>

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

<rule> ::= <parameterized-types-left> "->" <parameterized-types-right> "with" <with-clause> ";"
<solve-rule> ::= <parameterized-type-left> "->" <parameterized-type-right> "solving" <solve-clause> ";"
<parameterized-types-left> ::= <parameterized-types> "," | <parameterized-type>
<parameterized-types-right> ::= <parameterized-types> "," | <parameterized-type>
<parameterized-types> ::= <parameterized-types> "," | <parameterized-type>

<parameterized-type> ::= <symbol-name> <symbol-parameters> 
                      | "{" <parameterized-types> "|" <predicate> "}"

<with-clause> ::= <function-type> "where" <predicate> ";" | <function-type> ";"
<solve-clause> ::= <function-type> "solving" <expression> ";" | <function-type> ";"

## Observables


## State

## Grammar
<grammar> ::= "grammar" <time-type> <symbol-name> "(" <grammar-signature>) "{" <initial-conditions> <rules> "}"
<time-type> ::= "(discretetime)" | "(continuoustime)"

## Meta-Grammar

## State Type
<initial-conditions> ::= "initial_conditions" "{" <initial-condition-list> "}"
<initial-condition-list> ::= <initial-condition> | <initial-condition> <initial-condition-list>
<initial-condition> ::= <symbol-name> "=" <literal> ";"

<rules> ::= <rules> | <rule> | <solve-rule>

<grammar-signature> ::= <parameterized-types-left> "->" <parameterized-types-right>


 
 <predicate> ::= <symbol-name> "!=" <symbol-name> 
                 | <symbol-name> "<" <symbol-name> 
                 | <symbol-name> ">" <symbol-name> 
                 | <symbol-name> "==" <symbol-name>

<function-type> ::= <symbol-name> "(" <symbol-name> ")"
<literal> ::= [0-9] | [Aa-Zz] | <literal>
