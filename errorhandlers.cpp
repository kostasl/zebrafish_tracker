#include <larvatrack.h> //For CloseDAtaFile
#include <errorhandlers.h>
#include <GUI/mainwindow.h>


extern QString outfilename;
extern QFile outfishdatafile;
extern MainWindow* pwindow_main;




void on_sigabrt (int signum)
{
    void *array[10];
    size_t size;

    std::cerr << std::endl << "While Processing :"  << outfilename.toStdString() << " frame:" << pwindow_main->nFrame << std::endl;
    std::cerr << std::endl << ">>>> Simple SIG ABORT Handler Triggered <<<<<" << std::endl;
    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", signum);
    backtrace_symbols_fd(array, size, STDERR_FILENO);

    pwindow_main->LogEvent(QString("ERROR: ABORT SIGNAL RECEIVED AND HANDLED"));
    std::cerr << ">>>> Stop Execution <<<<<" << std::endl;
    // Skip  Code Block
    //std::cerr << ">>>> Skipping Execution <<<<<" << std::endl;

    closeDataFile(outfishdatafile);
    std::cerr << "Delete the output File" << std::endl;

    if (outfishdatafile.exists())
       outfishdatafile.remove();
    /*
     * The abort function causes abnormal program termination to occur, unless the signal
     * SIGABRT is being caught and the signal handler does not return. ...
     */
    //longjmp (env, 1);

}


/// Seg Fault Error Handler With Demangling - From stackoverflow//
void crit_err_hdlr(int sig_num, siginfo_t * info, void * ucontext)
{
    void *             array[50];
     void *             caller_address;
     char **            messages;
     int                size, i;
     sig_ucontext_t *   uc;

     uc = (sig_ucontext_t *)ucontext;

     std::cerr << "While Processing :"  << outfilename.toStdString() << " frame:" << pwindow_main->nFrame << std::endl;
     std::cerr << ">>>>  SIG SEG Handler with Demangling was Triggered <<<<<" << std::endl;

     closeDataFile(outfishdatafile);
     std::cerr << "Delete the output File" << std::endl;

     if (outfishdatafile.exists())
        outfishdatafile.remove();


     /* Get the address at the time the signal was raised */
    #if defined(__i386__) // gcc specific
     caller_address = (void *) uc->uc_mcontext.eip; // EIP: x86 specific
    #elif defined(__x86_64__) // gcc specific
     caller_address = (void *) uc->uc_mcontext.rip; // RIP: x86_64 specific
    #else
    #error Unsupported architecture. // TODO: Add support for other arch.
    #endif


    std::cerr << "signal " << sig_num;
    std::cerr << " (" << strsignal(sig_num) << "), address is ";
    std::cerr << info->si_addr << " from " << caller_address ;
    std::cerr << std::endl;



    size = backtrace(array, 50);
 /* overwrite sigaction with caller's address */
    array[1] = caller_address;

    messages = backtrace_symbols(array, size);

    // skip first stack frame (points here)
    for (int i = 1; i < size && messages != NULL; ++i)
    {
        char *mangled_name = 0, *offset_begin = 0, *offset_end = 0;

        // find parantheses and +address offset surrounding mangled name
        for (char *p = messages[i]; *p; ++p)
        {
            if (*p == '(')
            {
                mangled_name = p;
            }
            else if (*p == '+')
            {
                offset_begin = p;
            }
            else if (*p == ')')
            {
                offset_end = p;
                break;
            }
        }

        // if the line could be processed, attempt to demangle the symbol
        if (mangled_name && offset_begin && offset_end &&
            mangled_name < offset_begin)
        {
            *mangled_name++ = '\0';
            *offset_begin++ = '\0';
            *offset_end++ = '\0';

            int status;
            char * real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);

            // if demangling is successful, output the demangled function name
            if (status == 0)
            {
                std::cerr << "[bt]: (" << i << ") " << messages[i] << " : "
                          << real_name << "+" << offset_begin << offset_end
                          << std::endl;

            }
            // otherwise, output the mangled function name
            else
            {
                std::cerr << "[bt]: (" << i << ") " << messages[i] << " : "
                          << mangled_name << "+" << offset_begin << offset_end
                          << std::endl;
            }
            free(real_name);
        }
        // otherwise, print the whole line
        else
        {
            std::cerr << "[bt]: (" << i << ") " << messages[i] << std::endl;
        }
    }
    std::cerr << std::endl;

    free(messages);

    exit(EXIT_FAILURE);
}

void handler(int sig) {
  void *array[10];
  size_t size;

  std::cerr << ">>>> Simple SIG SEG Handler Triggered <<<<<" << std::endl;
  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);


  exit(1);
}
