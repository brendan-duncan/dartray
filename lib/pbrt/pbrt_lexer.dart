/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the "License");         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an "AS IS" BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of pbrt;

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
  static const int TOKEN_DOUBLE_BAR = -8; // ||
  static const int TOKEN_DOUBLE_AMPERSAND = -9; // &&
  static const int TOKEN_STRING = -10;

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
  static const int TOKEN_BANG = 33; // !
  static const int TOKEN_BAR = 124; // |
  static const int TOKEN_AMPERSAND = 38; // &
  static const int TOKEN_QUESTION = 63;// ?
  static const int TOKEN_COLON = 58; // :
  static const int TOKEN_DOT = 46; // .
  static const int TOKEN_DOUBLE_QUOTE = 34; // "
  static const int TOKEN_SINGLE_QUOTE = 39; // '
  static const int TOKEN_LEFT_BRACKET = 91;
  static const int TOKEN_RIGHT_BRACKET = 93;

  static const int TOKEN_0 = 48; // 0
  static const int TOKEN_9 = 57; // 9
  static const int TOKEN_a = 97; // a
  static const int TOKEN_z = 122; // z
  static const int TOKEN_A = 65; // A
  static const int TOKEN_Z = 90; // Z

  PbrtLexer(String input) {
    _inputStack.add(new _PbrtLexerInput(input));
  }

  void addInclude(String input) {
    _inputStack.add(new _PbrtLexerInput(input));
  }

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
      _lastChar = _nextChar();
    }

    // Nothing more to do if we're at the end of the input.
    if (_lastChar == _EOF) {
      _curToken = TOKEN_EOF;
      return TOKEN_EOF;
    }

    // Double-quote string: "*"
    if (_lastChar == '"') {
      _lastChar = _nextChar();
      while (_lastChar != '"' && _lastChar != _EOF) {
        _identifierStr += _lastChar;
        _lastChar = _nextChar();
      }

      _curTokenStr = _identifierStr;
      _curToken = TOKEN_STRING;

      // Eat the trailing quotation mark.
      _lastChar = _nextChar();

      return _curToken;
    }

    // Single-quote string: '*'
    if (_lastChar == "'") {
      _lastChar = _nextChar();
      while (_lastChar != "'" && _lastChar != _EOF) {
        _identifierStr += _lastChar;
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
      _identifierStr = _lastChar;
      while (_lastChar != _EOF &&
             (_isAlphaNum((_lastChar = _nextChar())) || _lastChar == '.' ||
              _lastChar == '_')) {
        _identifierStr = _identifierStr += _lastChar;
      }

      _curTokenStr = _identifierStr;
      _curToken = TOKEN_IDENTIFIER;

      return _curToken;
    }

    if (_lastChar == '-' && !_isDigit(_peekNext()) && _peekNext() != '.') {
      _curToken = TOKEN_MINUS;
      _curTokenStr = _lastChar;
      return _curToken;
    }

    // Number: [-](0-9)+[e[+-](0-9)+]
    if (_lastChar == '-' || _isDigit(_lastChar) || _lastChar == '.') {
      String numStr = '';
      do {
        numStr += _lastChar.toString();
        _lastChar = _nextChar();
        if (_lastChar == 'e') {
          numStr += _lastChar;
          _lastChar = _nextChar();
          if (_lastChar == '+' || _lastChar == '-') {
            numStr += _lastChar;
            _lastChar = _nextChar();
          }
        }
      } while (_lastChar != _EOF && (_isDigit(_lastChar) || _lastChar == '.'));

      _curTokenStr = numStr;
      _numValue = double.parse(numStr);
      _curToken = TOKEN_NUMBER;

      return _curToken;
    }

    // Skip # end-of-line comments
    if (_lastChar == '#') {
      // Comment until end of line.
      do {
        _lastChar = _nextChar();
      } while (_lastChar != _EOF && _lastChar != '\n' && _lastChar != '\r');

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
    String thisChar = _lastChar;
    _curTokenStr = thisChar;
    _lastChar = _nextChar();

    // Check for multi-charace operators.
    if (thisChar == '=' && _lastChar == '=') {
      _lastChar = _nextChar();
      _curTokenStr = '==';
      _curToken = TOKEN_DOUBLE_EQUAL;
      return _curToken;
    } else if (thisChar == '!' && _lastChar == '=') {
      _lastChar = _nextChar();
      _curTokenStr = '!=';
      _curToken = TOKEN_NOT_EQUAL;
      return _curToken;
    } else if (thisChar == '|' && _lastChar == '|') {
      _lastChar = _nextChar();
      _curTokenStr = '||';
      _curToken = TOKEN_DOUBLE_BAR;
      return _curToken;
    } else if (thisChar == '&' && _lastChar == '&') {
      _lastChar = _nextChar();
      _curTokenStr = '&&';
      _curToken = TOKEN_DOUBLE_AMPERSAND;
      return _curToken;
    }

    _curToken = thisChar.codeUnitAt(0);
    return _curToken;
  }

  /**
   * Get the next character from the input buffer.
   * TODO use the ord of the character instead of a substring.
   * _intput should be converted to a list of ord values...
   * Then, all character tests can be done as ints and not string
   * comparisons.
   */
  String _nextChar() {
    if (_inputStack.last.isEOF) {
      if (_inputStack.length == 1) {
        return _EOF;
      }
      _inputStack.removeLast();
    }
    return _inputStack.last.nextChar();
  }

  String _peekNext() {
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
  bool _isAlpha(String c) {
    if (c.isEmpty) {
      return false;
    }
    int cc = c.codeUnitAt(0);
    return ((cc >= TOKEN_a) && (cc <= TOKEN_z)) ||
           ((cc >= TOKEN_A) && (cc <= TOKEN_Z));
  }

  /**
   * Is the character in the range [0, 9]
   */
  bool _isDigit(String c) {
    if (c.isEmpty) {
      return false;
    }
    int cc = c.codeUnitAt(0);
    return (cc >= TOKEN_0) && (cc <= TOKEN_9);
  }

  /**
   * Is the character an letter or a number?
   */
  bool _isAlphaNum(String c) => _isAlpha(c) || _isDigit(c);

  /**
   * Is the character a space or tab?
   */
  bool _isWhitespace(String c) => (c == ' ' || c == '\t' ||
                                   c == '\n' || c == '\r');

  static const String _EOF = '';

  List<_PbrtLexerInput> _inputStack = [];
  String _identifierStr;
  double _numValue;
  String _lastChar = ' ';
  int _curToken;
  String _curTokenStr;
}

class _PbrtLexerInput {
  String input;
  int position = 0;

  _PbrtLexerInput(this.input);

  bool get isEOF => position >= input.length;

  String nextChar() => input[position++];

  String peekNext() => input[position];
}
