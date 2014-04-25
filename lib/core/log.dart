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
part of core;

const int LOG_INFO = 0;
const int LOG_WARNING = 1;
const int LOG_ERROR = 2;
const int LOG_SEVERE = 3;
const int LOG_DEBUG = 4;
const List<String> LOG_TYPES = const [
  'INFO',
  'WARNING',
  'ERROR',
  'SEVERE',
  'DEBUG'];

/**
 * A callback to call when a log message has been recieved.
 * [type] is one of the following: [LOG_INFO], [LOG_WARNING], [LOG_ERROR],
 * or [LOG_SEVERE].
 */
typedef LogCallback(int type, String msg);

void PrintLogger(int type, String msg) {
  print('${LOG_TYPES[type]}: $msg');
  if (type == LOG_SEVERE) {
    throw new Exception(msg);
  }
}

LogCallback Log = PrintLogger;

void SetLogger(LogCallback cb) {
  Log = cb;
}

void LogDebug(String msg) {
  Log(LOG_DEBUG, msg);
}

void LogInfo(String msg) {
  Log(LOG_INFO, msg);
}

void LogWarning(String msg) {
  Log(LOG_WARNING, msg);
}

void LogError(String msg) {
  Log(LOG_ERROR, msg);
}

void LogSevere(String msg) {
  Log(LOG_SEVERE, msg);
}
