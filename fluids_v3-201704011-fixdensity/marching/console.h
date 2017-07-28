#if defined(_MSC_VER)
#pragma once
#endif

#ifndef CONSOLE_H
#define CONSOLE_H

/************************************************************************/
/*					Console output routines                             */
/************************************************************************/

void Warning(const char *fmt,...);
void Error(const char *fmt,...);
void Debug(const char *fmt,...);
void Fatal(const char *fmt,...);

#endif // CONSOLE_H
