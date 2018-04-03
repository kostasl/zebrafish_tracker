#ifndef ERRORHANDLERS_H
#define ERRORHANDLERS_H

/// For Stack Trace Debugging //
#include <execinfo.h>
#include <cxxabi.h>
#include <signal.h>
#include <csetjmp>
#include <ucontext.h>
#include <unistd.h>

#include <iostream>



typedef struct _sig_ucontext {
 unsigned long     uc_flags;
 struct ucontext   *uc_link;
 stack_t           uc_stack;
 struct sigcontext uc_mcontext;
 sigset_t          uc_sigmask;
} sig_ucontext_t;



void on_sigabrt (int signum);
void crit_err_hdlr(int sig_num, siginfo_t * info, void * ucontext);
void handler(int sig);

#endif // ERRORHANDLERS_H

