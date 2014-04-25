/****************************************************************************
 * Copyright (C) 2014 by Brendan Duncan.                                    *
 *                                                                          *
 * This file is part of DartRay.                                            *
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License");          *
 * you may not use this file except in compliance with the License.         *
 * You may obtain a copy of the License at                                  *
 *                                                                          *
 * http://www.apache.org/licenses/LICENSE-2.0                               *
 *                                                                          *
 * Unless required by applicable law or agreed to in writing, software      *
 * distributed under the License is distributed on an "AS IS" BASIS,        *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 * See the License for the specific language governing permissions and      *
 * limitations under the License.                                           *
 *                                                                          *
 * This project is based on PBRT v2 ; see http://www.pbrt.org               *
 * pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.  *
 ****************************************************************************/
part of dartray;

/**
 * A lexer breaks up an incoming stream of characters into a sequence of
 * tokens.
 */
class PbrtLexer {
  static const int TOKEN_EOF = -1;
  static const int TOKEN_IDENTIFIER = -4;
  static const int TOKEN_NUMBER = -5;
  static const int TOKEN_DOUBLE_EQUAL = -6; // ==
  static const int TOKEN_NOT_EQUAL = -7; // !=
  static const int TOKEN_DOUBLE_OR = -8; // ||
  static const int TOKEN_DOUBLE_AND = -9; // &&
  static const int TOKEN_STRING = -10;

  static const int TOKEN_SPACE = 32; // ' '
  static const int TOKEN_TAB = 9; // \t
  static const int TOKEN_LINEFEED = 13; // \r
  static const int TOKEN_NEWLINE = 10; // \n
  static const int TOKEN_GREATER = 62; // >
  static const int TOKEN_LESS = 60; // <
  static const int TOKEN_LEFT_PAREN = 40; // (
  static const int TOKEN_RIGHT_PAREN = 41; // )
  static const int TOKEN_COMMA = 44; // ,
  static const int TOKEN_SEMICOLON = 59; // ;
  static const int TOKEN_PLUS = 43; // +
  static const int TOKEN_MINUS = 45; // -
  static const int TOKEN_MULTIPLY = 42; // *
  static const int TOKEN_DIVIDE = 47; // /
  static const int TOKEN_CARROT = 94; // ^
  static const int TOKEN_PERCENT = 37; // %
  static const int TOKEN_HASH = 35; // #
  static const int TOKEN_NOT = 33; // !
  static const int TOKEN_OR = 124; // |
  static const int TOKEN_AND = 38; // &
  static const int TOKEN_QUESTION = 63;// ?
  static const int TOKEN_COLON = 58; // :
  static const int TOKEN_DOT = 46; // .
  static const int TOKEN_DOUBLE_QUOTE = 34; // "
  static const int TOKEN_SINGLE_QUOTE = 39; // '
  static const int TOKEN_LEFT_BRACKET = 91; // [
  static const int TOKEN_RIGHT_BRACKET = 93; // ]
  static const int TOKEN_UNDERSCORE = 95; // _
  static const int TOKEN_EQUAL = 61;

  static const int TOKEN_0 = 48; // 0
  static const int TOKEN_9 = 57; // 9
  static const int TOKEN_a = 97; // a
  static const int TOKEN_e = 101; // e
  static const int TOKEN_z = 122; // z
  static const int TOKEN_A = 65; // A
  static const int TOKEN_E = 96; // E
  static const int TOKEN_Z = 90; // Z

  PbrtLexer(List<int> input, String path) {
    _inputStack.add(new _PbrtLexerInput(input, path));
  }

  void addInclude(List<int> input, String path) {
    _inputStack.add(new _PbrtLexerInput(input, path));
  }

  String get path => _inputStack.last.path;

  int get line => _inputStack.last.line;

  /**
   * Returns true if the parse position is at the end of the input.
   */
  bool isEof() => _curToken == TOKEN_EOF;

  /**
   * Returns the code of the last parsed token.
   */
  int get currentToken => _curToken;

  /**
   * Returns the string of the last parsed token.
   */
  String get currentTokenString => _curTokenStr;

  /**
   * If the parsed token was an identifier, this is the string name of the
   * identifier.
   */
  String get identifier => _identifierStr;

  /**
   * If the parsed token was a number, this is the number.
   */
  double get value => _numValue;

  /**
   * Returns the next parsed token.
   */
  int nextToken() {
    _identifierStr = '';
    _curTokenStr = '';

    // Ignore white-space (space, tabs, new-lines).
    while (_isWhitespace(_lastChar)) {
      if (_lastChar == TOKEN_NEWLINE) {
        _inputStack.last.line++;
      }
      _lastChar = _nextChar();
    }

    // Nothing more to do if we're at the end of the input.
    if (_lastChar == _EOF) {
      _curToken = TOKEN_EOF;
      return TOKEN_EOF;
    }

    // Double-quote string: "*"
    if (_lastChar == TOKEN_DOUBLE_QUOTE) {
      _lastChar = _nextChar();
      while (_lastChar != TOKEN_DOUBLE_QUOTE && _lastChar != _EOF) {
        _identifierStr += new String.fromCharCode(_lastChar);
        _lastChar = _nextChar();
      }

      _curTokenStr = _identifierStr;
      _curToken = TOKEN_STRING;

      // Eat the trailing quotation mark.
      _lastChar = _nextChar();

      return _curToken;
    }

    // Single-quote string: '*'
    if (_lastChar == TOKEN_SINGLE_QUOTE) {
      _lastChar = _nextChar();
      while (_lastChar != TOKEN_SINGLE_QUOTE && _lastChar != _EOF) {
        _identifierStr += new String.fromCharCode(_lastChar);
        _lastChar = _nextChar();
      }

      _curTokenStr = _identifierStr;
      _curToken = TOKEN_STRING;

      // Eat the trailing quotation mark.
      _lastChar = _nextChar();

      return _curToken;
    }

    // Identifier: [a-zA-Z][a-zA-Z0-9._]*
    if (_isAlpha(_lastChar)) {
      _identifierStr = new String.fromCharCode(_lastChar);
      _lastChar = _nextChar();
      while (_lastChar != _EOF && _isAlphaNum(_lastChar) ||
             _lastChar == TOKEN_DOT || _lastChar == TOKEN_UNDERSCORE) {
        _identifierStr += new String.fromCharCode(_lastChar);
        _lastChar = _nextChar();
      }

      _curTokenStr = _identifierStr;
      _curToken = TOKEN_IDENTIFIER;

      return _curToken;
    }

    if (_lastChar == TOKEN_MINUS &&
        !_isDigit(_peekNext()) && _peekNext() != TOKEN_DOT) {
      _curToken = TOKEN_MINUS;
      _curTokenStr = new String.fromCharCode(_lastChar);
      return _curToken;
    }

    // Number: [-](0-9)+[e[+-](0-9)+]
    if (_lastChar == TOKEN_MINUS || _isDigit(_lastChar) ||
        _lastChar == TOKEN_DOT) {
      String numStr = '';
      do {
        numStr += new String.fromCharCode(_lastChar);
        _lastChar = _nextChar();
        if (_lastChar == TOKEN_e || _lastChar == TOKEN_E) {
          numStr += new String.fromCharCode(_lastChar);
          _lastChar = _nextChar();
          if (_lastChar == TOKEN_PLUS || _lastChar == TOKEN_MINUS) {
            numStr += new String.fromCharCode(_lastChar);
            _lastChar = _nextChar();
          }
        }
      } while (_lastChar != _EOF && (_isDigit(_lastChar) ||
               _lastChar == TOKEN_DOT));

      _curTokenStr = numStr;
      _numValue = double.parse(numStr);
      _curToken = TOKEN_NUMBER;

      return _curToken;
    }

    // Skip # end-of-line comments
    if (_lastChar == TOKEN_HASH) {
      // Comment until end of line.
      do {
        _lastChar = _nextChar();
      } while (_lastChar != _EOF && _lastChar != TOKEN_NEWLINE &&
               _lastChar != TOKEN_LINEFEED);

      if (_lastChar != _EOF) {
        _curToken = nextToken();
        return _curToken;
      }
    }

    // Check for end of file.  Don't eat the EOF.
    if (_lastChar == _EOF) {
      _curToken = TOKEN_EOF;
      return _curToken;
    }

    // Otherwise, just return the character as its ascii value.
    int thisChar = _lastChar;
    _curTokenStr = new String.fromCharCode(thisChar);
    _lastChar = _nextChar();

    // Check for multi-charace operators.
    if (thisChar == TOKEN_EQUAL && _lastChar == TOKEN_EQUAL) {
      _lastChar = _nextChar();
      _curTokenStr = '==';
      _curToken = TOKEN_DOUBLE_EQUAL;
      return _curToken;
    } else if (thisChar == TOKEN_NOT && _lastChar == TOKEN_EQUAL) {
      _lastChar = _nextChar();
      _curTokenStr = '!=';
      _curToken = TOKEN_NOT_EQUAL;
      return _curToken;
    } else if (thisChar == TOKEN_OR && _lastChar == TOKEN_OR) {
      _lastChar = _nextChar();
      _curTokenStr = '||';
      _curToken = TOKEN_DOUBLE_OR;
      return _curToken;
    } else if (thisChar == TOKEN_AND && _lastChar == TOKEN_AND) {
      _lastChar = _nextChar();
      _curTokenStr = '&&';
      _curToken = TOKEN_DOUBLE_AND;
      return _curToken;
    }

    _curToken = thisChar;
    return _curToken;
  }

  /**
   * Get the next character from the input buffer.
   */
  int _nextChar() {
    if (_inputStack.last.isEOF) {
      if (_inputStack.length == 1) {
        return _EOF;
      }
      _inputStack.removeLast();
    }
    return _inputStack.last.nextChar();
  }

  int _peekNext() {
    int len = _inputStack.length;
    while (len > 0) {
      if (!_inputStack[len - 1].isEOF) {
        return _inputStack[len - 1].peekNext();
      }
    }
    return _EOF;
  }

  /**
   * Is the character in the range [a, z] or [A, Z]?
   */
  bool _isAlpha(int c) {
    return ((c >= TOKEN_a) && (c <= TOKEN_z)) ||
           ((c >= TOKEN_A) && (c <= TOKEN_Z));
  }

  /**
   * Is the character in the range [0, 9]
   */
  bool _isDigit(int c) {
    return (c >= TOKEN_0) && (c <= TOKEN_9);
  }

  /**
   * Is the character an letter or a number?
   */
  bool _isAlphaNum(int c) => _isAlpha(c) || _isDigit(c);

  /**
   * Is the character a space or tab?
   */
  bool _isWhitespace(int c) => (c == TOKEN_SPACE || c == TOKEN_TAB ||
                                c == TOKEN_NEWLINE || c == TOKEN_LINEFEED);

  static const int _EOF = 0;

  List<_PbrtLexerInput> _inputStack = [];
  String _identifierStr;
  double _numValue;
  int _lastChar = TOKEN_SPACE;
  int _curToken;
  String _curTokenStr;
}

class _PbrtLexerInput {
  List<int> input;
  String path;
  int line = 0;
  int position = 0;

  _PbrtLexerInput(this.input, this.path);

  bool get isEOF => position >= input.length;

  int nextChar() => input[position++];

  int peekNext() => input[position];
}
