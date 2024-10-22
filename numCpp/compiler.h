#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cctype>
#include <stack>
#include <unordered_map>
#include <functional>

using namespace std;

enum TokenType {
  AND,
  OR,
  NOT,
  LPAREN,
  RPAREN,
  IDENTIFIER,
  PLUS,
  MINUS,
  MULTIPLY,
  DIVIDE,
  VARIABLE,
  EQUAL,
  NOT_EQUAL,
  LESS_THAN,
  GREATER_THAN,
  LESS_EQUAL,
  GREATER_EQUAL,
  END
};

struct Token {
  TokenType type;
  string value;
};

vector<Token> tokenize(const string& input) {
  vector<Token> tokens;
  size_t i = 0;
  while (i < input.size()) {
    char current = input[i];
    switch (current) {
      case ' ':
        ++i;
        break;
      case '&':
        tokens.push_back({AND, "&"});
        ++i;
        break;
      case '|':
        tokens.push_back({OR, "|"});
        ++i;
        break;
      case '!':
        if (i + 1 < input.size() && input[i + 1] == '=') {
          tokens.push_back({NOT_EQUAL, "!="});
          i += 2;
        } else {
          tokens.push_back({NOT, "!"});
          ++i;
        }
        break;
      case '(':
        tokens.push_back({LPAREN, "("});
        ++i;
        break;
      case ')':
        tokens.push_back({RPAREN, ")"});
        ++i;
        break;
      case '+':
        tokens.push_back({PLUS, "+"});
        ++i;
        break;
      case '-':
        tokens.push_back({MINUS, "-"});
        ++i;
        break;
      case '*':
        tokens.push_back({MULTIPLY, "*"});
        ++i;
        break;
      case '/':
        tokens.push_back({DIVIDE, "/"});
        ++i;
        break;
      case 'x':
        tokens.push_back({VARIABLE, "x"});
        ++i;
        break;
      case '=':
        if (i + 1 < input.size() && input[i + 1] == '=') {
          tokens.push_back({EQUAL, "=="});
          i += 2;
        } else {
          throw runtime_error("Unknown character");
        }
        break;
      case '<':
        if (i + 1 < input.size() && input[i + 1] == '=') {
          tokens.push_back({LESS_EQUAL, "<="});
          i += 2;
        } else {
          tokens.push_back({LESS_THAN, "<"});
          ++i;
        }
        break;
      case '>':
        if (i + 1 < input.size() && input[i + 1] == '=') {
          tokens.push_back({GREATER_EQUAL, ">="});
          i += 2;
        } else {
          tokens.push_back({GREATER_THAN, ">"});
          ++i;
        }
        break;
      default:
        if (isdigit(current) || current == '.') {
          string number;
          while (i < input.size() && (isdigit(input[i]) || input[i] == '.')) {
            number += input[i++];
          }
          tokens.push_back({IDENTIFIER, number});
        } else {
          throw runtime_error("Unknown character");
        }
        break;
    }
  }
  tokens.push_back({END, ""});
  return tokens;
}

struct ASTNode {
  TokenType type;
  string value;
  vector<ASTNode*> children;
};

class Parser {
  vector<Token> tokens;
  size_t pos;
  

  Token currentToken() { return tokens[pos]; }

  void consume(TokenType type) {
    if (currentToken().type == type) {
      ++pos;
    } else {
      ok=false; // throw runtime_error("Unexpected token");
    }
  }

  ASTNode* parseExpression() { return parseOr(); }

  ASTNode* parseOr() {
    ASTNode* node = parseAnd();
    while (currentToken().type == OR) {
      Token token = currentToken();
      consume(OR);
      ASTNode* newNode = new ASTNode{token.type, token.value};
      newNode->children.push_back(node);
      newNode->children.push_back(parseAnd());
      node = newNode;
    }
    return node;
  }

  ASTNode* parseAnd() {
    ASTNode* node = parseEquality();
    while (currentToken().type == AND) {
      Token token = currentToken();
      consume(AND);
      ASTNode* newNode = new ASTNode{token.type, token.value};
      newNode->children.push_back(node);
      newNode->children.push_back(parseEquality());
      node = newNode;
    }
    return node;
  }

  ASTNode* parseEquality() {
    ASTNode* node = parseRelational();
    while (currentToken().type == EQUAL || currentToken().type == NOT_EQUAL) {
      Token token = currentToken();
      consume(token.type);
      ASTNode* newNode = new ASTNode{token.type, token.value};
      newNode->children.push_back(node);
      newNode->children.push_back(parseRelational());
      node = newNode;
    }
    return node;
  }

  ASTNode* parseRelational() {
    ASTNode* node = parseAdditive();
    while (currentToken().type == LESS_THAN ||
           currentToken().type == GREATER_THAN ||
           currentToken().type == LESS_EQUAL ||
           currentToken().type == GREATER_EQUAL) {
      Token token = currentToken();
      consume(token.type);
      ASTNode* newNode = new ASTNode{token.type, token.value};
      newNode->children.push_back(node);
      newNode->children.push_back(parseAdditive());
      node = newNode;
    }
    return node;
  }

  ASTNode* parseAdditive() {
    ASTNode* node = parseMultiplicative();
    while (currentToken().type == PLUS || currentToken().type == MINUS) {
      Token token = currentToken();
      consume(token.type);
      ASTNode* newNode = new ASTNode{token.type, token.value};
      newNode->children.push_back(node);
      newNode->children.push_back(parseMultiplicative());
      node = newNode;
    }
    return node;
  }

  ASTNode* parseMultiplicative() {
    ASTNode* node = parseUnary();
    while (currentToken().type == MULTIPLY || currentToken().type == DIVIDE) {
      Token token = currentToken();
      consume(token.type);
      ASTNode* newNode = new ASTNode{token.type, token.value};
      newNode->children.push_back(node);
      newNode->children.push_back(parseUnary());
      node = newNode;
    }
    return node;
  }

  ASTNode* parseUnary() {
    if (currentToken().type == NOT) {
      Token token = currentToken();
      consume(NOT);
      ASTNode* node = new ASTNode{token.type, token.value};
      node->children.push_back(parseUnary());
      return node;
    } else if (currentToken().type == LPAREN) {
      consume(LPAREN);
      ASTNode* node = parseExpression();
      consume(RPAREN);
      return node;
    } else if (currentToken().type == IDENTIFIER ||
               currentToken().type == VARIABLE) {
      Token token = currentToken();
      consume(token.type);
      return new ASTNode{token.type, token.value};
    } else {
      ok=false; // throw runtime_error("Unexpected token");
      return nullptr;
    }
  }

 public:
  Parser(const vector<Token>& tokens) : tokens(tokens), pos(0) {}
  
  ASTNode* parse() { return parseExpression(); }
  bool ok=true;
};

void generateCode(ASTNode* node) {
  if (node->type == IDENTIFIER || node->type == VARIABLE) {
    cout << "PUSH " << node->value << endl;
  } else if (node->type == AND) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "AND" << endl;
  } else if (node->type == OR) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "OR" << endl;
  } else if (node->type == NOT) {
    generateCode(node->children[0]);
    cout << "NOT" << endl;
  } else if (node->type == PLUS) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "ADD" << endl;
  } else if (node->type == MINUS) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "SUB" << endl;
  } else if (node->type == MULTIPLY) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "MUL" << endl;
  } else if (node->type == DIVIDE) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "DIV" << endl;
  } else if (node->type == EQUAL) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "EQ" << endl;
  } else if (node->type == NOT_EQUAL) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "NEQ" << endl;
  } else if (node->type == LESS_THAN) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "LT" << endl;
  } else if (node->type == GREATER_THAN) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "GT" << endl;
  } else if (node->type == LESS_EQUAL) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "LE" << endl;
  } else if (node->type == GREATER_EQUAL) {
    generateCode(node->children[0]);
    generateCode(node->children[1]);
    cout << "GE" << endl;
  }
}

bool execute(ASTNode* node, double x) {
  stack<double> stack;

  function<void(ASTNode*)> eval = [&](ASTNode* node) {
    switch (node->type) {
      case IDENTIFIER:
        stack.push(stod(node->value));
        break;
      case VARIABLE:
        stack.push(x);
        break;
      case AND:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a && b);
        }
        break;
      case OR:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a || b);
        }
        break;
      case NOT:
        eval(node->children[0]);
        {
          double a = stack.top();
          stack.pop();
          stack.push(!a);
        }
        break;
      case PLUS:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a + b);
        }
        break;
      case MINUS:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a - b);
        }
        break;
      case MULTIPLY:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a * b);
        }
        break;
      case DIVIDE:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a / b);
        }
        break;
      case EQUAL:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a == b);
        }
        break;
      case NOT_EQUAL:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a != b);
        }
        break;
      case LESS_THAN:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a < b);
        }
        break;
      case GREATER_THAN:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a > b);
        }
        break;
      case LESS_EQUAL:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a <= b);
        }
        break;
      case GREATER_EQUAL:
        eval(node->children[0]);
        eval(node->children[1]);
        {
          double b = stack.top();
          stack.pop();
          double a = stack.top();
          stack.pop();
          stack.push(a >= b);
        }
        break;
      default:
        throw runtime_error("Unknown node type");
    }
  };

  eval(node);
  return stack.top();
}

class Evaluator {
 public:  // props.
  ASTNode* ast = nullptr;
  bool ok=true;

 public:
  Evaluator(string expr) {
    auto tokens = tokenize(expr);
    Parser parser(tokens);
    ast = parser.parse();
    
    ok=parser.ok;
  }
  ~Evaluator() { deleteAST(ast); }

  bool evaluate(double x) { return execute(ast, x); }

 private:
  void deleteAST(ASTNode* node) {
    if (node == nullptr) return;
    for (ASTNode* child : node->children) deleteAST(child);

    delete node;
  }
};


